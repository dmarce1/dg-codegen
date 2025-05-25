#pragma once
#include "Basis.hpp"
#include "ContainerArithmetic.hpp"
#include "Hdf5.hpp"
#include "Matrix.hpp"
#include "Quadrature.hpp"

#include <numeric>

template<typename Type>
inline constexpr Type minmod(Type const &a, Type const &b) {
	using namespace Math;
	Type const sgn = copysign(Type(0.5), a) + copysign(Type(0.5), b);
	Type const mag = min(abs(a), abs(b));
	return sgn * mag;
}

template<typename State, int cellsAcrossInterior, int basisOrder, typename RungeKutta>
class HyperGrid {

	static constexpr int dimensionCount = State::dimCount();
	static constexpr int ghostWidth = 2;
	static constexpr int fieldCount = State::fieldCount();
	static constexpr int rungeKuttaStageCount = RungeKutta::stageCount();
	static constexpr int basisSize = BasisIndexType<basisOrder, dimensionCount>::count();
	static constexpr int cellsAcrossExterior = cellsAcrossInterior + 2 * ghostWidth;
	static constexpr int exteriorVolume = ipow(cellsAcrossExterior, dimensionCount);
	static constexpr RungeKutta butcherTable { };
	static constexpr Range<int, dimensionCount> exteriorBox {
			repeat<dimensionCount>(-ghostWidth),
			repeat<dimensionCount>(cellsAcrossInterior + ghostWidth) };
	static constexpr Range<int, dimensionCount> interiorBox {
			repeat<dimensionCount>(0),
			repeat<dimensionCount>(cellsAcrossInterior) };
	static constexpr Quadrature<typename State::value_type, basisOrder, dimensionCount> volumeQuadrature { };
	static constexpr Basis<typename State::value_type, basisOrder, dimensionCount> orthogonalBasis { };

	using Type = State::value_type;
	using BasisIndex = BasisIndexType<basisOrder, dimensionCount>;
	using QuadratureType = Quadrature<Type, basisOrder, State::dimCount()>;
	using InteriorIndex = MultiIndex<exteriorBox, interiorBox>;

	std::vector<std::array<std::array<std::array<Type, exteriorVolume>, basisSize>, fieldCount>> stageDerivatives_;
	std::vector<std::array<std::array<Type, exteriorVolume>, basisSize>> nextState;
	std::vector<std::array<std::array<Type, exteriorVolume>, basisSize>> currentState;
	const Type cellWidth;
	const Type inverseCellWidth;

	static constexpr auto interiorIndexMap(int i) {
		static constexpr auto map = createMultiIndexMap<exteriorBox, interiorBox>();
		return map[i];
	}
	static constexpr int stride(int d) {
		return ipow(cellsAcrossExterior, dimensionCount - 1 - d);
	}

public:
	HyperGrid(Type const &xNint = Type(1)) :
			nextState(fieldCount), cellWidth(xNint / Type(cellsAcrossInterior)), inverseCellWidth(Type(cellsAcrossInterior) / xNint) {
	}
	void initialize(std::function<State(std::array<Type, dimensionCount> const&)> const &initialState) {
		Type const halfCellWidth = Type(0.5) * cellWidth;
		auto const massMatrix = orthogonalBasis.massMatrix();
		auto const inverseMassMatrix = matrixInverse(massMatrix);
		for (auto cellMultiIndex = InteriorIndex::begin(); cellMultiIndex != InteriorIndex::end(); cellMultiIndex++) {
			int const cellFlatIndex = cellMultiIndex;
			for (int basisIndex = 0; basisIndex < basisSize; basisIndex++) {
				for (int fieldIndex = 0; fieldIndex < fieldCount; fieldIndex++) {
					nextState[fieldIndex][basisIndex][cellFlatIndex] = Type(0);
				}
				for (int quadratureIndex = 0; quadratureIndex < volumeQuadrature.size(); quadratureIndex++) {
					auto const basis = orthogonalBasis(volumeQuadrature.point(quadratureIndex));
					auto const weight = volumeQuadrature.weight(quadratureIndex);
					auto const quadraturePoint = volumeQuadrature.point(quadratureIndex);
					std::array<Type, dimensionCount> position;
					for (int dimension = 0; dimension < dimensionCount; dimension++) {
						position[dimension] = (Type(2 * cellMultiIndex[dimension] + 1) + quadraturePoint[dimension]) * halfCellWidth;
					}
					auto const thisState = initialState(position);
					for (int fieldIndex = 0; fieldIndex < fieldCount; fieldIndex++) {
						nextState[fieldIndex][basisIndex][cellFlatIndex] += weight * basis[basisIndex] * thisState[fieldIndex];
					}
				}
				for (int fieldIndex = 0; fieldIndex < fieldCount; fieldIndex++) {
					nextState[fieldIndex][basisIndex][cellFlatIndex] *= inverseMassMatrix(basisIndex, basisIndex);
				}
			}
		}
	}
	void output(const char *filenameBase, int timeStepNumber, Type const &time) {
		std::string filename = std::string(filenameBase) + "." + std::to_string(timeStepNumber) + ".h5";
		writeHdf5<Type, dimensionCount, cellsAcrossInterior, basisOrder, ghostWidth>(filename, cellWidth, nextState, State::getFieldNames());
		writeList("X.visit", "!NBLOCKS 1\n", filename + ".xmf");
	}
	void enforceBoundaryConditions() {
		using Index = MultiIndex<exteriorBox>;
		for (auto ghostZoneIndex = Index::begin(); ghostZoneIndex != Index::end(); ghostZoneIndex++) {
			Index interiorIndex;
			bool isGhostZone = false;
			for (int dimensionIndex = 0; dimensionIndex < dimensionCount; dimensionIndex++) {
				if (ghostZoneIndex[dimensionIndex] < 0) {
					isGhostZone = true;
					interiorIndex[dimensionIndex] = ghostZoneIndex[dimensionIndex] + cellsAcrossInterior;
				} else if (ghostZoneIndex[dimensionIndex] >= cellsAcrossInterior) {
					isGhostZone = true;
					interiorIndex[dimensionIndex] = ghostZoneIndex[dimensionIndex] - cellsAcrossInterior;
				} else
					interiorIndex[dimensionIndex] = ghostZoneIndex[dimensionIndex];
			}
			if (!isGhostZone) {
				continue;
			}
			for (int fieldIndex = 0; fieldIndex < fieldCount; fieldIndex++) {
				for (int basisIndex = 0; basisIndex < basisSize; basisIndex++) {
					nextState[fieldIndex][basisIndex][ghostZoneIndex] = nextState[fieldIndex][basisIndex][interiorIndex];
				}
			}
		}
	}
	Type beginStep() {
		applyLimiter();
		using namespace Math;
		currentState = nextState;
		stageDerivatives_.resize(rungeKuttaStageCount);
		std::array<Type, dimensionCount> maximumEigenvalue;
		maximumEigenvalue.fill(Type(0));
		for (auto cellMultiIndex = InteriorIndex::begin(); cellMultiIndex != InteriorIndex::end(); cellMultiIndex++) {
			int const cellFlatIndex = cellMultiIndex;
			for (int quadratureIndex = 0; quadratureIndex < volumeQuadrature.size(); quadratureIndex++) {
				State thisState;
				auto const basis = orthogonalBasis(volumeQuadrature.point(quadratureIndex));
				for (int fieldIndex = 0; fieldIndex < fieldCount; fieldIndex++) {
					thisState[fieldIndex] = Type(0);
					for (int modeIndex = 0; modeIndex < basisSize; modeIndex++) {
						thisState[fieldIndex] += nextState[fieldIndex][modeIndex][cellFlatIndex] * basis[modeIndex];
					}
				}
				for (int dimension = 0; dimension < dimensionCount; dimension++) {
					auto const eigenvalues = thisState.eigenvalues(dimension);
					for (auto thisEigenvalue : eigenvalues) {
						maximumEigenvalue[dimension] = max(maximumEigenvalue[dimension], abs(thisEigenvalue));
					}
				}
			}
		}
		Type maximumEigenvalueSum = Type(0);
		for (int dimension = 0; dimension < dimensionCount; dimension++) {
			maximumEigenvalueSum += maximumEigenvalue[dimension];
		}
		Type const timeStepSize = (cellWidth * butcherTable.cfl()) / (Type(2 * basisOrder - 1) * maximumEigenvalueSum);
		return timeStepSize;
	}
	void subStep(Type const &timeStepSize, int stageIndex) {
		nextState = currentState;
		for (int fieldIndex = 0; fieldIndex < fieldCount; fieldIndex++) {
			for (int modeIndex = 0; modeIndex < basisSize; modeIndex++) {
				for (int thisStage = 0; thisStage < stageIndex; thisStage++) {
					for (auto cellMultiIndex = InteriorIndex::begin(); cellMultiIndex != InteriorIndex::end(); cellMultiIndex++) {
						int const cellFlatIndex = cellMultiIndex;
						nextState[fieldIndex][modeIndex][cellFlatIndex] +=
								butcherTable.a(stageIndex, thisStage) * stageDerivatives_[thisStage][fieldIndex][modeIndex][cellFlatIndex];
					}
				}
				stageDerivatives_[stageIndex][fieldIndex][modeIndex].fill(Type(0));
			}
		}
		enforceBoundaryConditions();
		applyLimiter();
		computeDudt(timeStepSize, stageDerivatives_[stageIndex], std::make_integer_sequence<int, dimensionCount> { });
	}
	void endStep() {
		for (int fieldIndex = 0; fieldIndex < fieldCount; fieldIndex++) {
			for (auto modeMultiIndex = BasisIndex::begin(); modeMultiIndex != BasisIndex::end(); modeMultiIndex++) {
				int const modeFlatIndex = modeMultiIndex;
				for (auto cellMultiIndex = InteriorIndex::begin(); cellMultiIndex != InteriorIndex::end(); cellMultiIndex++) {
					int const cellFlatIndex = cellMultiIndex;
					nextState[fieldIndex][modeFlatIndex][cellFlatIndex] = currentState[fieldIndex][modeFlatIndex][cellFlatIndex];
					for (int stageIndex = 0; stageIndex < rungeKuttaStageCount; stageIndex++) {
						nextState[fieldIndex][modeFlatIndex][cellFlatIndex] +=
								butcherTable.b(stageIndex) * stageDerivatives_[stageIndex][fieldIndex][modeFlatIndex][cellFlatIndex];
					}
				}
			}
		}
		stageDerivatives_ = { };
		currentState = { };
	}
	void applyLimiter() {
		constexpr auto limiterCoefficients = []() {
			std::array<Type, basisOrder> limiterCoefficients;
			for (int orderIndex = 0; orderIndex < basisOrder; orderIndex++) {
				limiterCoefficients[orderIndex] = Type(1) / Type(2 * orderIndex + 1);
			}
			return limiterCoefficients;
		}();

		constexpr Range<int, dimensionCount> limiterBox {
				repeat<dimensionCount>(-1),
				repeat<dimensionCount>(cellsAcrossInterior + 1) };

		using LimiterIndex = MultiIndex<exteriorBox, limiterBox>;

		for (int polynomialDegree = basisOrder - 1; polynomialDegree > 0; polynomialDegree--) {
			for (auto targetModeIndex = BasisIndex::begin(); targetModeIndex != BasisIndex::end(); targetModeIndex++) {
				for (int dimension = 0; dimension < dimensionCount; dimension++) {

					int totalDegree = 0;
					for (int dimension = 0; dimension < dimensionCount; dimension++) {
						totalDegree += targetModeIndex[dimension];
						if (totalDegree > polynomialDegree) {
							break;
						}
					}
					if ((totalDegree != polynomialDegree) || (targetModeIndex[dimension] == 0)) {
						continue;
					}

					auto lowerModeIndex = targetModeIndex;
					lowerModeIndex[dimension]--;
					auto const lowerModeFlatIndex = lowerModeIndex;
					auto const targetModeFlatIndex = targetModeIndex;
					auto const strideOffset = stride(dimension);

					for (auto cellMultiIndex = LimiterIndex::begin(); cellMultiIndex != LimiterIndex::end(); cellMultiIndex++) {
						int const cellFlatIndex = cellMultiIndex;
						State highOrderMode, referenceState, slopeDifference;

						for (int fieldIndex = 0; fieldIndex < fieldCount; fieldIndex++) {
							highOrderMode[fieldIndex] = nextState[fieldIndex][targetModeFlatIndex][cellFlatIndex];
							referenceState[fieldIndex] = nextState[fieldIndex][0][cellFlatIndex];
							auto const lowOrderMode = nextState[fieldIndex][lowerModeFlatIndex][cellFlatIndex];
							auto const differenceRight = nextState[fieldIndex][lowerModeFlatIndex][cellFlatIndex + strideOffset] - lowOrderMode;
							auto const differenceLeft = lowOrderMode - nextState[fieldIndex][lowerModeFlatIndex][cellFlatIndex - strideOffset];
							slopeDifference[fieldIndex] = minmod(differenceRight, differenceLeft);
						}

						auto const [_, rightEigenvectors] = referenceState.eigenSystem(dimension);
						auto const leftEigenvectors = matrixInverse(rightEigenvectors);
						auto const limiterCoefficient = limiterCoefficients[lowerModeIndex[dimension]];

						highOrderMode = leftEigenvectors * highOrderMode;
						slopeDifference = leftEigenvectors * slopeDifference;

						for (int fieldIndex = 0; fieldIndex < fieldCount; fieldIndex++) {
							highOrderMode[fieldIndex] = minmod(highOrderMode[fieldIndex], limiterCoefficient * slopeDifference[fieldIndex]);
						}

						highOrderMode = rightEigenvectors * highOrderMode;

						for (int fieldIndex = 0; fieldIndex < fieldCount; fieldIndex++) {
							nextState[fieldIndex][targetModeFlatIndex][cellFlatIndex] = highOrderMode[fieldIndex];
						}

					}
				}
			}
		}
	}

	template<int ... dimension>
	void computeDudt(Type timeStepSize, std::array<std::array<std::array<Type, exteriorVolume>, basisSize>, fieldCount> &stateDerivative,
			std::integer_sequence<int, dimension...>) {
		(computeDudtByDim<dimension>(timeStepSize, stateDerivative), ...);
	}
	template<int dimension>
	void computeDudtByDim(Type timeStepSize, std::array<std::array<std::array<Type, exteriorVolume>, basisSize>, fieldCount> &stateDerivative) {
		constexpr Range<int, dimensionCount> interiorBoxPlusOne {
				repeat<dimensionCount>(0),
				repeat<dimensionCount>(cellsAcrossInterior) + unit<dimensionCount>(dimension) };
		constexpr Quadrature<Type, basisOrder, dimensionCount - 1> surfaceQuadrature;
		using FluxIndexType = MultiIndex<exteriorBox, interiorBoxPlusOne>;
		Type const lambda = Type(2) * timeStepSize * inverseCellWidth;
		for (auto cellMultIndex = FluxIndexType::begin(); cellMultIndex != FluxIndexType::end(); cellMultIndex++) {
			int const rightInterfaceIndex = cellMultIndex;
			int const leftInterfaceIndex = rightInterfaceIndex - stride(dimension);
			for (int quadratureIndex = 0; quadratureIndex < surfaceQuadrature.size(); quadratureIndex++) {
				auto const leftPosition = insert<dimensionCount>(Type(+1), dimension, surfaceQuadrature.point(quadratureIndex));
				auto const rightPosition = insert<dimensionCount>(Type(-1), dimension, surfaceQuadrature.point(quadratureIndex));
				Type const quadratureWeight = surfaceQuadrature.weight(quadratureIndex);
				auto const basisLeft = orthogonalBasis(leftPosition);
				auto const basisRight = orthogonalBasis(rightPosition);
				State stateRight, stateLeft;
				for (int fieldIndex = 0; fieldIndex < fieldCount; fieldIndex++) {
					stateRight[fieldIndex] = Type(0);
					stateLeft[fieldIndex] = Type(0);
					for (int modeIndex = 0; modeIndex < basisSize; modeIndex++) {
						stateLeft[fieldIndex] += nextState[fieldIndex][modeIndex][leftInterfaceIndex] * basisLeft[modeIndex];
						stateRight[fieldIndex] += nextState[fieldIndex][modeIndex][rightInterfaceIndex] * basisRight[modeIndex];
					}
				}
				auto const riemannFlux = solveRiemannProblem(stateLeft, stateRight, dimension);
				if (cellMultIndex[dimension] > 0) {
					for (int fieldIndex = 0; fieldIndex < fieldCount; fieldIndex++) {
						for (int modeIndex = 0; modeIndex < basisSize; modeIndex++) {
							stateDerivative[fieldIndex][modeIndex][leftInterfaceIndex] -= quadratureWeight * lambda * riemannFlux[fieldIndex] * basisLeft[modeIndex];
						}
					}
				}
				if (cellMultIndex[dimension] < cellsAcrossInterior) {
					for (int fieldIndex = 0; fieldIndex < fieldCount; fieldIndex++) {
						for (int modeIndex = 0; modeIndex < basisSize; modeIndex++) {
							stateDerivative[fieldIndex][modeIndex][rightInterfaceIndex] += quadratureWeight * lambda * riemannFlux[fieldIndex] * basisRight[modeIndex];
						}
					}
				}
			}
		}
		for (auto cellMultiIndex = InteriorIndex::begin(); cellMultiIndex != InteriorIndex::end(); cellMultiIndex++) {
			int const cellFlatIndex = cellMultiIndex;
			for (int quadratureIndex = 0; quadratureIndex < volumeQuadrature.size(); quadratureIndex++) {
				auto const quadraturePoint = volumeQuadrature.point(quadratureIndex);
				auto const quadratureWeight = volumeQuadrature.weight(quadratureIndex);
				auto const basis = orthogonalBasis(quadraturePoint);
				auto const basisDerivative = orthogonalBasis.gradient(dimension, quadraturePoint);
				State state;
				for (int fieldIndex = 0; fieldIndex < fieldCount; fieldIndex++) {
					state[fieldIndex] = Type(0);
					for (int basisIndex = 0; basisIndex < basisSize; basisIndex++) {
						state[fieldIndex] += nextState[fieldIndex][basisIndex][cellFlatIndex] * basis[basisIndex];
					}
				}
				auto const flux = state.flux(dimension);
				for (int fieldIndex = 0; fieldIndex < fieldCount; fieldIndex++) {
					for (int basisIndex = 0; basisIndex < basisSize; basisIndex++) {
						stateDerivative[fieldIndex][basisIndex][cellFlatIndex] += quadratureWeight * lambda * flux[fieldIndex] * basisDerivative[basisIndex];
					}
				}
			}
		}
		for (auto cellMultiIndex = InteriorIndex::begin(); cellMultiIndex != InteriorIndex::end(); cellMultiIndex++) {
			int const cellFlatIndex = cellMultiIndex;
			for (int fieldIndex = 0; fieldIndex < fieldCount; fieldIndex++) {
				for (auto basisMultiIndex = BasisIndex::begin(); basisMultiIndex != BasisIndex::end(); basisMultiIndex++) {
					int const basisIndex = basisMultiIndex;
					stateDerivative[fieldIndex][basisIndex][cellFlatIndex] *= Type(0.5) * Type(2 * basisMultiIndex[dimension] + 1);
				}
			}
		}
	}
};

