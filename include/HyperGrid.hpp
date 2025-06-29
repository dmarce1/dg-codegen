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
#include <span>

template<typename Type>
inline constexpr Type minmod(Type const &a, Type const &b) {
	Type const sgn = copysign(Type(0.5), a) + copysign(Type(0.5), b);
	Type const mag = min(abs(a), abs(b));
	return sgn * mag;
}

template<typename Type, int dimensionCount, int cellsAcrossInterior, int modeCount, typename RungeKutta, template<typename, int> typename State>
class HyperGrid {
	static constexpr int ghostWidth = 2;
	static constexpr int fieldCount = State<Type, dimensionCount>::fieldCount();
	static constexpr int rungeKuttaStageCount = RungeKutta::stageCount();
	static constexpr int modeVolume = BasisIndexType<modeCount, dimensionCount>::count();
	static constexpr int nodeVolume = ipow(modeCount, dimensionCount);
	static constexpr int modeSurface = BasisIndexType<modeCount, dimensionCount - 1>::count();
	static constexpr int nodeSurface = ipow(modeCount, dimensionCount - 1);
	static constexpr int cellsAcrossExterior = cellsAcrossInterior + 2 * ghostWidth;
	static constexpr int exteriorVolume = ipow(cellsAcrossExterior, dimensionCount);
	static constexpr RungeKutta butcherTable { };
	static constexpr Range<int, dimensionCount> exteriorBox { repeat<dimensionCount>(-ghostWidth), repeat<dimensionCount>(cellsAcrossInterior + ghostWidth) };
	static constexpr Range<int, dimensionCount> interiorBox { repeat<dimensionCount>(0), repeat<dimensionCount>(cellsAcrossInterior) };
	template<typename V>
	static constexpr auto vAnalyze = dgAnalyze<V, dimensionCount, modeCount>;
	template<typename V>
	static constexpr auto vMassInverse = dgMassInverse<V, dimensionCount, modeCount>;
	template<typename V>
	static constexpr auto vSynthesize = dgSynthesize<V, dimensionCount, modeCount>;
	template<typename V>
	static constexpr auto vStiffness = dgStiffness<V, dimensionCount, modeCount>;
	template<typename V>
	static constexpr auto sAnalyze = dgAnalyze<V, dimensionCount - 1, modeCount>;
	template<typename V>
	static constexpr auto sSynthesize = dgSynthesize<V, dimensionCount - 1, modeCount>;
	template<typename V>
	static constexpr auto sTrace = dgTrace<V, dimensionCount, modeCount>;
	template<typename V>
	static constexpr auto sTraceInverse = dgTraceInverse<ValArray<Type>, dimensionCount, modeCount>;
	std::array<State<std::array<ValArray<Type>, modeVolume>, dimensionCount>, rungeKuttaStageCount> stageDerivatives_;
	State<std::array<ValArray<Type>, modeVolume>, dimensionCount> currentState;
	State<std::array<ValArray<Type>, modeVolume>, dimensionCount> previousState;
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
		for (int fieldIndex = 0; fieldIndex < fieldCount; fieldIndex++) {
			for (int modeIndex = 0; modeIndex < modeVolume; modeIndex++) {
				currentState[fieldIndex][modeIndex].resize(exteriorVolume);
			}
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
	void initialize(std::function<State<Type, dimensionCount>(std::array<Type, dimensionCount> const&)> const &initialState) {
		Type const halfCellWidth = Type(0.5) * cellWidth;
		State<std::array<ValArray<Type>, nodeVolume>, dimensionCount> nodalValues;
		for (int fieldIndex = 0; fieldIndex < fieldCount; fieldIndex++) {
			nodalValues[fieldIndex].fill(ValArray<Type> { exteriorVolume });
		}
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
				for (int fieldIndex = 0; fieldIndex < fieldCount; fieldIndex++) {
					nodalValues[fieldIndex][nodeIndex][i] = state[fieldIndex];
				}
			}
		}
		for (int fieldIndex = 0; fieldIndex < fieldCount; fieldIndex++) {
			currentState[fieldIndex] = vMassInverse<ValArray<Type>>(vAnalyze<ValArray<Type>>(nodalValues[fieldIndex]));
		}
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
		currentState = previousState;
		for (int fieldIndex = 0; fieldIndex < fieldCount; fieldIndex++) {
			for (int modeIndex = 0; modeIndex < modeVolume; modeIndex++) {
				for (int thisStage = 0; thisStage < stageIndex; thisStage++) {
					currentState[fieldIndex][modeIndex] += butcherTable.a(stageIndex, thisStage) * stageDerivatives_[thisStage][fieldIndex][modeIndex];
				}
				stageDerivatives_[stageIndex][fieldIndex][modeIndex].fill(Type(0));
			}
		}
		enforceBoundaryConditions();
		applyLimiter();
		computeTimeDerivative(timeStepSize, stageDerivatives_[stageIndex]);
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
		constexpr int hiModeVolume = modeVolume - 1;
		constexpr int loModeVolume = binco(dimensionCount + modeCount - 2, dimensionCount);
		State<ValArray<Type>, dimensionCount> meanState;
		for (int fieldIndex = 0; fieldIndex < fieldCount; fieldIndex++) {
			for (int i = 0; i < exteriorVolume; i++) {
				meanState[fieldIndex][i] = currentState[fieldIndex][0][i];
			}
		}
		State<ValArray<Type>, dimensionCount> differenceState;
		for (int fieldIndex = 0; fieldIndex < fieldCount; fieldIndex++) {
			differenceState[fieldIndex].resize(exteriorVolume);
		}
		for (int polynomialDegree = modeCount - 1; polynomialDegree > 0; polynomialDegree--) {
			int const begin = binco(polynomialDegree + dimensionCount - 1, dimensionCount);
			int const end = binco(polynomialDegree + dimensionCount, dimensionCount);
			for (int hiModeIndex = begin; hiModeIndex < end; hiModeIndex++) {
				std::array<State<ValArray<Type>, dimensionCount>, dimensionCount> transposedState;
				for (int dimension = 0; dimension < dimensionCount; dimension++) {
					auto hiModeIndices = flatToTriangular<dimensionCount, modeCount>(hiModeIndex);
					if (hiModeIndices[dimension] == 0) {
						continue;
					}
					auto loModeIndices = hiModeIndices;
					loModeIndices[dimension]--;
					auto const loModeIndex = triangularToFlat<dimensionCount, modeCount>(loModeIndices);
					for (int fieldIndex = 0; fieldIndex < fieldCount; fieldIndex++) {
						transposedState[dimension][fieldIndex] = currentState[fieldIndex][hiModeIndex];
					}
					for (int fieldIndex = 0; fieldIndex < fieldCount; fieldIndex++) {
						auto const statePlus = currentState[fieldIndex][loModeIndex].shift(+gridStrides[dimension]);
						auto const &stateCentral = currentState[fieldIndex][loModeIndex];
						auto const stateMinus = currentState[fieldIndex][loModeIndex].shift(-gridStrides[dimension]);
						differenceState[fieldIndex] = minmod(ValArray<Type>(statePlus - stateCentral), ValArray<Type>(stateCentral - stateMinus));
					}
					auto const rightEigenvectors = meanState.eigenSystem(dimension).second;
					auto const leftEigenvectors = matrixInverse(rightEigenvectors);
					differenceState = leftEigenvectors * differenceState;
					transposedState[dimension] = leftEigenvectors * transposedState[dimension];
					for (int fieldIndex = 0; fieldIndex < fieldCount; fieldIndex++) {
						auto const cLimit = Type(1) / Type(2 * loModeIndices[dimension] + 1);
						differenceState[fieldIndex] *= cLimit;
						transposedState[dimension][fieldIndex] = minmod(transposedState[dimension][fieldIndex], differenceState[fieldIndex]);
					}
				}
				for (int dimension = 0; dimension < dimensionCount; dimension++) {
					auto hiModeIndices = flatToTriangular<dimensionCount, modeCount>(hiModeIndex);
					if (hiModeIndices[dimension] == 0) {
						continue;
					}
					for (int fieldIndex = 0; fieldIndex < fieldCount; fieldIndex++) {
						auto& current = currentState[fieldIndex][hiModeIndex];
						current = minmod(current, transposedState[dimension][fieldIndex]);
					}
				}
			}
		}
		ValArray<Type> theta { Type(1), exteriorVolume };
		std::array<State<ValArray<Type>, dimensionCount>, nodeSurface> surfaceState;
		std::array<State<ValArray<Type>, dimensionCount>, nodeVolume> volumeState;
		for (int face = 0; face < 2 * dimensionCount; face++) {
			for (int fieldIndex = 0; fieldIndex < fieldCount; fieldIndex++) {
				auto const surfaceNodes = sSynthesize<ValArray<Type>>(sTrace<ValArray<Type>>(face, currentState[fieldIndex]));
				for (int nodeIndex = 0; nodeIndex < nodeSurface; nodeIndex++) {
					surfaceState[nodeIndex][fieldIndex] = surfaceNodes[nodeIndex];
				}
			}
			for (int nodeIndex = 0; nodeIndex < nodeSurface; nodeIndex++) {
				theta = min(Type(1), findPositivityPreservingTheta(meanState, surfaceState[nodeIndex]));
			}
		}
		for (int fieldIndex = 0; fieldIndex < fieldCount; fieldIndex++) {
			auto const volumeNodes = vSynthesize<ValArray<Type>>(currentState[fieldIndex]);
			for (int nodeIndex = 0; nodeIndex < nodeVolume; nodeIndex++) {
				volumeState[nodeIndex][fieldIndex] = volumeNodes[nodeIndex];
			}
		}
		for (int nodeIndex = 0; nodeIndex < nodeVolume; nodeIndex++) {
			theta = min(Type(1), findPositivityPreservingTheta(meanState, volumeState[nodeIndex]));
		}
		for (int fieldIndex = 0; fieldIndex < fieldCount; fieldIndex++) {
			for (int modeIndex = 1; modeIndex < modeVolume; modeIndex++) {
				currentState[fieldIndex][modeIndex] *= theta;
			}
		}
	}
	void computeTimeDerivative(Type timeStepSize, State<std::array<ValArray<Type>, modeVolume>, dimensionCount> &stateDerivative) {
		Type const lambda = Type(2) * timeStepSize / cellWidth;
		for (int dimension = 0; dimension < dimensionCount; dimension++) {
			auto const shift = gridStrides[dimension];
			std::array<State<ValArray<Type>, dimensionCount>, nodeSurface> leftState;
			std::array<State<ValArray<Type>, dimensionCount>, nodeSurface> rightState;
			State<std::array<ValArray<Type>, nodeSurface>, dimensionCount> riemannFlux;
			for (int fieldIndex = 0; fieldIndex < fieldCount; fieldIndex++) {
				auto leftModes = currentState[fieldIndex];
				for (int modeIndex = 0; modeIndex < modeCount; modeIndex++) {
					leftModes[modeIndex] = currentState[fieldIndex][modeIndex].cshift(-shift);
				}
				auto const left = sSynthesize<ValArray<Type>>(sTrace<ValArray<Type>>(2 * dimension + 1, leftModes));
				auto const right = sSynthesize<ValArray<Type>>(sTrace<ValArray<Type>>(2 * dimension + 0, currentState[fieldIndex]));
				for (int nodeIndex = 0; nodeIndex < nodeSurface; nodeIndex++) {
					leftState[nodeIndex][fieldIndex] = left[nodeIndex];
					rightState[nodeIndex][fieldIndex] = right[nodeIndex];
				}
			}
			for (int nodeIndex = 0; nodeIndex < nodeSurface; nodeIndex++) {
				auto const tmp = solveRiemannProblem(leftState[nodeIndex], rightState[nodeIndex], dimension);
				for (int fieldIndex = 0; fieldIndex < fieldCount; fieldIndex++) {
					riemannFlux[fieldIndex][nodeIndex] = tmp[fieldIndex];
				}
			}
			for (int fieldIndex = 0; fieldIndex < fieldCount; fieldIndex++) {
				auto const thisFlux = sAnalyze<ValArray<Type>>(riemannFlux[fieldIndex]);
				auto const leftFlux = vMassInverse<ValArray<Type>>(sTraceInverse<ValArray<Type>>(2 * dimension + 1, thisFlux));
				auto const rightFlux = vMassInverse<ValArray<Type>>(sTraceInverse<ValArray<Type>>(2 * dimension + 0, thisFlux));
				for (int modeIndex = 0; modeIndex < modeVolume; modeIndex++) {
					stateDerivative[fieldIndex][modeIndex] -= lambda * leftFlux[modeIndex].cshift(+shift);
					stateDerivative[fieldIndex][modeIndex] += lambda * rightFlux[modeIndex];
				}
			}
		}
		for (int dimension = 0; dimension < dimensionCount; dimension++) {
			std::array<State<ValArray<Type>, dimensionCount>, nodeVolume> volumeState;
			State<std::array<ValArray<Type>, nodeVolume>, dimensionCount> volumeFlux;
			for (int fieldIndex = 0; fieldIndex < fieldCount; fieldIndex++) {
				auto const nodes = vSynthesize<ValArray<Type>>(currentState[fieldIndex]);
				for (int nodeIndex = 0; nodeIndex < nodeVolume; nodeIndex++) {
					volumeState[nodeIndex][fieldIndex] = nodes[nodeIndex];
				}
			}
			for (int nodeIndex = 0; nodeIndex < nodeVolume; nodeIndex++) {
				auto tmp = volumeState[nodeIndex].flux(dimension);
				for (int fieldIndex = 0; fieldIndex < fieldCount; fieldIndex++) {
					volumeFlux[fieldIndex][nodeIndex] = tmp[fieldIndex];
				}
			}
			for (int fieldIndex = 0; fieldIndex < fieldCount; fieldIndex++) {
				auto const source = vMassInverse<ValArray<Type>>(
						vStiffness<ValArray<Type>>(dimension, vMassInverse<ValArray<Type>>(vAnalyze<ValArray<Type>>(volumeFlux[fieldIndex]))));
				for (int modeIndex = 0; modeIndex < modeVolume; modeIndex++) {
					stateDerivative[fieldIndex][modeIndex] += lambda * source[modeIndex];
				}
			}
		}
	}
};

