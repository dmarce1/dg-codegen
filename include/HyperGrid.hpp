#pragma once
#include "Basis.hpp"
#include "ContainerArithmetic.hpp"
#include "Hdf5.hpp"
#include "Matrix.hpp"
#include "MultiIndex.hpp"
#include "Quadrature.hpp"
#include "dgTransforms.hpp"
#include "ValArray.hpp"

#include <functional>
#include <numeric>

template<typename Type>
inline constexpr Type minmod(Type const &a, Type const &b) {
	using namespace Math;
	Type const sgn = copysign(Type(0.5), a) + copysign(Type(0.5), b);
	Type const mag = min(abs(a), abs(b));
	return sgn * mag;
}

template<typename Type, int modeCount, TransformDirection transformDirection>
constexpr auto fourierLegendreTransform() {
	if constexpr (transformDirection == TransformDirection::forward) {
		return matrixInverse(fourierLegendreTransform<Type, modeCount, TransformDirection::backward>());
	} else {
		using namespace Math;
		SquareMatrix<Type, modeCount> transform;
		for (int basisIndex = 0; basisIndex < modeCount; basisIndex++) {
			for (int quadratureIndex = 0; quadratureIndex < modeCount; quadratureIndex++) {
				auto const [quadraturePosition, quadratureWeight] = gaussLegendreQuadraturePoint<Type, modeCount>(quadratureIndex);
				transform(quadratureIndex, basisIndex) = legendrePolynomial(basisIndex, quadraturePosition);
			}
		}
		return transform;
	}
}

template<typename State, int cellsAcrossInterior, int modeCount, typename RungeKutta>
class HyperGrid {
	static constexpr int dimensionCount = State::dimCount();
	static constexpr int ghostWidth = 2;
	static constexpr int fieldCount = State::fieldCount();
	static constexpr int rungeKuttaStageCount = RungeKutta::stageCount();
	static constexpr int modeVolume = BasisIndexType<modeCount, dimensionCount>::count();
	static constexpr int nodeVolume = ipow(modeCount, dimensionCount);
	static constexpr int cellsAcrossExterior = cellsAcrossInterior + 2 * ghostWidth;
	static constexpr int exteriorVolume = ipow(cellsAcrossExterior, dimensionCount);
	static constexpr RungeKutta butcherTable { };
	static constexpr Range<int, dimensionCount> exteriorBox { repeat<dimensionCount>(-ghostWidth), repeat<dimensionCount>(cellsAcrossInterior + ghostWidth) };
	static constexpr Range<int, dimensionCount> interiorBox { repeat<dimensionCount>(0), repeat<dimensionCount>(cellsAcrossInterior) };
	static constexpr auto inverseTransformMatrix = fourierLegendreTransform<typename State::value_type, modeCount, TransformDirection::backward>();
	static constexpr auto transformMatrix = fourierLegendreTransform<typename State::value_type, modeCount, TransformDirection::forward>();

	using Type = State::value_type;
	using BasisIndex = BasisIndexType<modeCount, dimensionCount>;
	using InteriorIndex = MultiIndex<exteriorBox, interiorBox>;

	std::array<std::array<ValArray<Type>, modeVolume>, rungeKuttaStageCount> stageDerivatives_;
	std::array<ValArray<Type>, modeVolume> currentState;
	std::array<ValArray<Type>, modeVolume> previousState;
	std::array<ValArray<Type>, dimensionCount> position;
	Type cellWidth;
	ValArray<int> gridSizes;
	ValArray<int> gridStrides;

	static constexpr auto interiorIndexMap(int i) {
		static constexpr auto map = createMultiIndexMap<exteriorBox, interiorBox>();
		return map[i];
	}
	static constexpr int stride(int d) {
		return ipow(cellsAcrossExterior, dimensionCount - 1 - d);
	}

public:
	HyperGrid(Type const &xNint = Type(1)) {
		cellWidth = xNint / Type(cellsAcrossInterior);
		for (int mode = 0; mode < modeVolume; mode++) {
			currentState[mode].resize(fieldCount * exteriorVolume);
		}
		gridSizes.resize(dimensionCount, cellsAcrossExterior);
		gridStrides.resize(dimensionCount);
		gridStrides[dimensionCount - 1] = 1;
		for (int dimension = dimensionCount - 1; dimension > 0; dimension--) {
			gridStrides[dimension - 1] = gridStrides[dimension] * gridSizes[dimension];
		}
		for (int dimension = 0; dimension < dimensionCount; dimension++) {
			position[dimension].resize(exteriorVolume);
			for (int i = 0; i < cellsAcrossExterior; i++) {
				auto const start = i * gridStrides[dimension];
				auto sizes = gridSizes;
				sizes[dimension] = 1;
				Type const x = Type(2 * (i - ghostWidth) + 1) * Type(0.5) * cellWidth;
				position[dimension][GSlice { start, sizes, gridStrides }] = x;
			}
		}
	}
	void initialize(std::function<State(std::array<Type, dimensionCount> const&)> const &initialState) {
		Type const halfCellWidth = Type(0.5) * cellWidth;
		std::array<ValArray<Type>, nodeVolume> nodalValues;
		nodalValues.fill(ValArray<Type> { exteriorVolume * fieldCount });
		for (int nodeIndex = 0; nodeIndex < nodeVolume; nodeIndex++) {
			auto quadraturePoint = getQuadraturePoint<Type, dimensionCount, modeCount>(nodeIndex);
			auto thisPosition = position;
			for (int dimension = 0; dimension < dimensionCount; dimension++) {
				thisPosition[dimension] += quadraturePoint[dimension] * halfCellWidth;
			}
			for (int i = 0; i < exteriorVolume; i++) {
				std::array<Type, dimensionCount> x;
				for (int dimension = 0; dimension < dimensionCount; dimension++) {
					x[dimension] = thisPosition[dimension][i];
				}
				auto const state = initialState(x);
				std::copy_n(state.begin(), fieldCount, nodalValues[nodeIndex].begin() + i * fieldCount);
			}
		}
		currentState = dgMassInverse<ValArray<Type>, dimensionCount, modeCount>(dgAnalyze<ValArray<Type>, dimensionCount, modeCount>(nodalValues));
//		GSlice const toSlice { 0, { fieldCount, exteriorVolume }, { exteriorVolume, 1 } };
//		GSlice const fromSlice { 0, { fieldCount, exteriorVolume }, { 1, fieldCount } };
	}
	void output(const char *filenameBase, int timeStepNumber, Type const &time) {
//		std::string filename = std::string(filenameBase) + "." + std::to_string(timeStepNumber) + ".h5";
//		writeHdf5<Type, dimensionCount, cellsAcrossInterior, modeCount, ghostWidth>(filename, cellWidth, currentState, State::getFieldNames());
//		writeList("X.visit", "!NBLOCKS 1\n", filename + ".xmf");
	}
	void enforceBoundaryConditions() {
//		using Index = MultiIndex<exteriorBox>;
//		for (auto ghostZoneIndex = Index::begin(); ghostZoneIndex != Index::end(); ghostZoneIndex++) {
//			Index interiorIndex;
//			bool isGhostZone = false;
//			for (int dimensionIndex = 0; dimensionIndex < dimensionCount; dimensionIndex++) {
//				if (ghostZoneIndex[dimensionIndex] < 0) {
//					isGhostZone = true;
//					interiorIndex[dimensionIndex] = ghostZoneIndex[dimensionIndex] + cellsAcrossInterior;
//				} else if (ghostZoneIndex[dimensionIndex] >= cellsAcrossInterior) {
//					isGhostZone = true;
//					interiorIndex[dimensionIndex] = ghostZoneIndex[dimensionIndex] - cellsAcrossInterior;
//				} else
//					interiorIndex[dimensionIndex] = ghostZoneIndex[dimensionIndex];
//			}
//			if (!isGhostZone) {
//				continue;
//			}
//			for (int fieldIndex = 0; fieldIndex < fieldCount; fieldIndex++) {
//				for (int basisIndex = 0; basisIndex < modeVolume; basisIndex++) {
//					currentState[fieldIndex][basisIndex][ghostZoneIndex] = currentState[fieldIndex][basisIndex][interiorIndex];
//				}
//			}
//		}
	}
	Type beginStep() {
//		applyLimiter();
//		using namespace Math;
//		previousState = currentState;
//		stageDerivatives_.fill(rungeKuttaStageCount, ValArray<Type> { fieldCount * modeVolume * exteriorVolume });
//		std::array<Type, dimensionCount> maximumEigenvalue;
//		maximumEigenvalue.fill(Type(0));
//		for (auto cellMultiIndex = InteriorIndex::begin(); cellMultiIndex != InteriorIndex::end(); cellMultiIndex++) {
//			int const cellFlatIndex = cellMultiIndex;
//			std::array<Type, modeVolume> modalState;
//			std::array<std::array<Type, nodeVolume>, fieldCount> nodalState;
//			for (int fieldIndex = 0; fieldIndex < fieldCount; fieldIndex++) {
//				for (int modeIndex = 0; modeIndex < modeVolume; modeIndex++) {
//					modalState[modeIndex] = currentState[fieldIndex][modeIndex][cellFlatIndex];
//				}
//				nodalState[fieldIndex] = dgSynthesize<Type, dimensionCount, modeCount>(modalState);
//			}
//			for (int nodeIndex = 0; nodeIndex < nodeVolume; nodeIndex++) {
//				State thisState;
//				for (int fieldIndex = 0; fieldIndex < fieldCount; fieldIndex++) {
//					thisState[fieldIndex] = nodalState[fieldIndex][nodeIndex];
//				}
//				for (int dimension = 0; dimension < dimensionCount; dimension++) {
//					auto const eigenvalues = thisState.eigenvalues(dimension);
//					for (auto thisEigenvalue : eigenvalues) {
//						maximumEigenvalue[dimension] = max(maximumEigenvalue[dimension], abs(thisEigenvalue));
//					}
//				}
//			}
//		}
//		Type maximumEigenvalueSum = Type(0);
//		for (int dimension = 0; dimension < dimensionCount; dimension++) {
//			maximumEigenvalueSum += maximumEigenvalue[dimension];
//		}
//		Type const timeStepSize = (cellWidth * butcherTable.cfl()) / (Type(2 * modeCount - 1) * maximumEigenvalueSum);
		return Type(0);
		//		return timeStepSize;
	}
	void subStep(Type const &timeStepSize, int stageIndex) {
//		currentState = previousState;
//		for (int fieldIndex = 0; fieldIndex < fieldCount; fieldIndex++) {
//			for (int modeIndex = 0; modeIndex < modeVolume; modeIndex++) {
//				for (int thisStage = 0; thisStage < stageIndex; thisStage++) {
//					for (auto cellMultiIndex = InteriorIndex::begin(); cellMultiIndex != InteriorIndex::end(); cellMultiIndex++) {
//						int const cellFlatIndex = cellMultiIndex;
//						currentState[fieldIndex][modeIndex][cellFlatIndex] +=
//								butcherTable.a(stageIndex, thisStage) * stageDerivatives_[thisStage][fieldIndex][modeIndex][cellFlatIndex];
//					}
//				}
//				stageDerivatives_[stageIndex][fieldIndex][modeIndex].fill(Type(0));
//			}
//		}
//		enforceBoundaryConditions();
//		applyLimiter();
//		computeDudt(timeStepSize, stageDerivatives_[stageIndex], std::make_integer_sequence<int, dimensionCount> { });
	}
	void endStep() {
//		for (int fieldIndex = 0; fieldIndex < fieldCount; fieldIndex++) {
//			for (auto modeMultiIndex = BasisIndex::begin(); modeMultiIndex != BasisIndex::end(); modeMultiIndex++) {
//				int const modeFlatIndex = modeMultiIndex;
//				for (auto cellMultiIndex = InteriorIndex::begin(); cellMultiIndex != InteriorIndex::end(); cellMultiIndex++) {
//					int const cellFlatIndex = cellMultiIndex;
//					currentState[fieldIndex][modeFlatIndex][cellFlatIndex] = previousState[fieldIndex][modeFlatIndex][cellFlatIndex];
//					for (int stageIndex = 0; stageIndex < rungeKuttaStageCount; stageIndex++) {
//						currentState[fieldIndex][modeFlatIndex][cellFlatIndex] +=
//								butcherTable.b(stageIndex) * stageDerivatives_[stageIndex][fieldIndex][modeFlatIndex][cellFlatIndex];
//					}
//				}
//			}
//		}
//		stageDerivatives_ = { };
//		previousState = { };
	}
	void applyLimiter() {
		constexpr int lowModeVolume = binco(dimensionCount + modeCount - 2, dimensionCount);
		std::array<std::array<ValArray<Type>, lowModeVolume>, dimensionCount> differenceState;
		ValArray<bool> mask { true, fieldCount * exteriorVolume };
		for (int modeIndex = 0; modeIndex < lowModeVolume; modeIndex++) {
			for (int dimension = 0; dimension < dimensionCount; dimension++) {
				differenceState[dimension][modeIndex].resize(fieldCount * exteriorVolume);
				auto const statePlus = currentState[modeIndex].shift(+fieldCount * gridStrides[dimension]);
				auto const &stateCentral = currentState[modeIndex];
				auto const stateMinus = currentState[modeIndex].shift(-fieldCount * gridStrides[dimension]);
				differenceState[dimension][modeIndex] = minmod(statePlus - stateCentral, stateCentral - stateMinus);
			}
		}
		for (int polynomialDegree = modeCount - 1; polynomialDegree > 0; polynomialDegree--) {
			int const begin = binco(polynomialDegree + dimensionCount - 1, dimensionCount);
			int const end = binco(polynomialDegree + dimensionCount, dimensionCount);
			auto const nextMask = mask;
			for (int modeIndex = begin; modeIndex < end; modeIndex++) {
				auto hiModeIndex = flatToTriangular<dimensionCount, modeCount>(modeIndex);
				for (int dimension = 0; dimension < dimensionCount; dimension++) {
					auto loModeIndex = hiModeIndex;
					loModeIndex[dimension]--;
					auto const cLimit = Type(1) / Type(2 * loModeIndex[dimension] + 1);
					auto const oldState = currentState[hiModeIndex];
					auto const nextState = minmod(oldState[mask], cLimit * differenceState[dimension][loModeIndex][mask]);
					currentState[hiModeIndex][mask] = nextState;
					nextMask = nextMask || (currentState[hiModeIndex] != oldState);
				}
			}
			mask = std::move(nextMask);
			if (mask.sum() == 0) {
				break;
			}
		}
	}

	template<int ... dimension>
	void computeDudt(Type timeStepSize, ValArray<ValArray<ValArray<Type>>> &stateDerivative, std::integer_sequence<int, dimension...>) {
		(computeDudtByDim<dimension>(timeStepSize, stateDerivative), ...);
	}
	template<int dimension>
	void computeDudtByDim(Type timeStepSize, ValArray<ValArray<ValArray<Type>>> &stateDerivative) {
		auto const nodalLeft = dgSynthesize<Type, dimensionCount - 1, modeCount>(dgTrace<Type, dimensionCount, modeCount>(2 * dimension + 1, currentState));
		auto const nodalRight = dgSynthesize<Type, dimensionCount - 1, modeCount>(dgTrace<Type, dimensionCount, modeCount>(2 * dimension - 1, currentState));
//		constexpr Range<int, dimensionCount> interiorBoxPlusOne { repeat<dimensionCount>(0), repeat<dimensionCount>(cellsAcrossInterior) + unit<dimensionCount>(
//				dimension) };
//		using FluxIndexType = MultiIndex<exteriorBox, interiorBoxPlusOne>;
//		Type const lambda = Type(2) * timeStepSize / cellWidth;
//		for (auto cellMultIndex = FluxIndexType::begin(); cellMultIndex != FluxIndexType::end(); cellMultIndex++) {
//			int const rightInterfaceIndex = cellMultIndex;
//			int const leftInterfaceIndex = rightInterfaceIndex - stride(dimension);
//
//			std::array<std::array<Type, triangleSize<dimensionCount, modeCount>>, fieldCount> volumeModesLeft, volumeModesRight;
//			std::array<std::array<Type, squareSize<dimensionCount - 1, modeCount>>, fieldCount> surfaceNodesLeft, surfaceNodesRight, surfaceNodesFlux;
//			for (int fieldIndex = 0; fieldIndex < fieldCount; fieldIndex++) {
//				for (int modeIndex = 0; modeIndex < modeVolume; modeIndex++) {
//					volumeModesLeft[fieldIndex][modeIndex] = currentState[fieldIndex][modeIndex][leftInterfaceIndex];
//					volumeModesRight[fieldIndex][modeIndex] = currentState[fieldIndex][modeIndex][rightInterfaceIndex];
//				}
//				surfaceNodesLeft[fieldIndex] = dgSynthesize<Type, dimensionCount - 1, modeCount>(
//						dgTrace<Type, dimensionCount, modeCount>(2 * dimension + 1, volumeModesLeft[fieldIndex]));
//				surfaceNodesRight[fieldIndex] = dgSynthesize<Type, dimensionCount - 1, modeCount>(
//						dgTrace<Type, dimensionCount, modeCount>(2 * dimension + 0, volumeModesRight[fieldIndex]));
//			}
//
//			for (int nodeIndex = 0; nodeIndex<squareSize < dimensionCount - 1, modeCount> ; nodeIndex++) {
//				State stateRight, stateLeft;
//				for (int fieldIndex = 0; fieldIndex < fieldCount; fieldIndex++) {
//					stateLeft[fieldIndex] = surfaceNodesLeft[fieldIndex][nodeIndex];
//					stateRight[fieldIndex] = surfaceNodesRight[fieldIndex][nodeIndex];
//				}
//				auto const riemannFlux = solveRiemannProblem(stateLeft, stateRight, dimension);
//				for (int fieldIndex = 0; fieldIndex < fieldCount; fieldIndex++) {
//					surfaceNodesFlux[fieldIndex][nodeIndex] = riemannFlux[fieldIndex];
//				}
//			}
//			for (int fieldIndex = 0; fieldIndex < fieldCount; fieldIndex++) {
//				auto const modes = dgAnalyze<Type, dimensionCount - 1, modeCount>(surfaceNodesFlux[fieldIndex]);
//				volumeModesLeft[fieldIndex] = dgTraceInverse<Type, dimensionCount, modeCount>(2 * dimension + 1, modes);
//				volumeModesLeft[fieldIndex] = dgMassInverse<Type, dimensionCount, modeCount>(volumeModesLeft[fieldIndex]);
//				volumeModesRight[fieldIndex] = dgTraceInverse<Type, dimensionCount, modeCount>(2 * dimension + 0, modes);
//				volumeModesRight[fieldIndex] = dgMassInverse<Type, dimensionCount, modeCount>(volumeModesRight[fieldIndex]);
//			}
//
//			if (cellMultIndex[dimension] > 0) {
//				for (int fieldIndex = 0; fieldIndex < fieldCount; fieldIndex++) {
//					for (int modeIndex = 0; modeIndex < modeVolume; modeIndex++) {
//						stateDerivative[fieldIndex][modeIndex][leftInterfaceIndex] -= lambda * volumeModesLeft[fieldIndex][modeIndex];
//					}
//				}
//			}
//			if (cellMultIndex[dimension] < cellsAcrossInterior) {
//				for (int fieldIndex = 0; fieldIndex < fieldCount; fieldIndex++) {
//					for (int modeIndex = 0; modeIndex < modeVolume; modeIndex++) {
//						stateDerivative[fieldIndex][modeIndex][rightInterfaceIndex] += lambda * volumeModesRight[fieldIndex][modeIndex];
//					}
//				}
//			}
//		}
//		for (auto cellMultiIndex = InteriorIndex::begin(); cellMultiIndex != InteriorIndex::end(); cellMultiIndex++) {
//			int const cellFlatIndex = cellMultiIndex;
//			std::array<State, nodeVolume> nodeStates;
//			for (int fieldIndex = 0; fieldIndex < fieldCount; fieldIndex++) {
//				std::array<Type, modeVolume> modalState;
//				for (int modeIndex = 0; modeIndex < modeVolume; modeIndex++) {
//					modalState[modeIndex] = currentState[fieldIndex][modeIndex][cellFlatIndex];
//				}
//				auto thisValue = dgSynthesize<Type, dimensionCount, modeCount>(modalState);
//				for (int nodeIndex = 0; nodeIndex < nodeVolume; nodeIndex++) {
//					nodeStates[nodeIndex][fieldIndex] = thisValue[nodeIndex];
//				}
//			}
//			for (int nodeIndex = 0; nodeIndex < nodeVolume; nodeIndex++) {
//				nodeStates[nodeIndex] = nodeStates[nodeIndex].flux(dimension);
//			}
//			for (int fieldIndex = 0; fieldIndex < fieldCount; fieldIndex++) {
//				std::array<Type, nodeVolume> nodes;
//				for (int nodeIndex = 0; nodeIndex < nodeVolume; nodeIndex++) {
//					nodes[nodeIndex] = nodeStates[nodeIndex][fieldIndex];
//				}
//				auto flux = dgAnalyze<Type, dimensionCount, modeCount>(nodes);
//				flux = dgMassInverse<Type, dimensionCount, modeCount>(flux);
//				flux = dgStiffness<Type, dimensionCount, modeCount>(dimension, flux);
//				flux = dgMassInverse<Type, dimensionCount, modeCount>(flux);
//				for (int modeIndex = 0; modeIndex < modeVolume; modeIndex++) {
//					stateDerivative[fieldIndex][modeIndex][cellFlatIndex] += lambda * flux[modeIndex];
//				}
//			}
//		}
	}
};

