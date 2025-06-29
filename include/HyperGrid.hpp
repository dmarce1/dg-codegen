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
	static constexpr int gridVolume = ipow(cellsAcrossExterior, dimensionCount);
	static constexpr RungeKutta butcherTable { };
	static constexpr Range<int, dimensionCount> exteriorBox { repeat<dimensionCount>(-ghostWidth), repeat<dimensionCount>(cellsAcrossInterior + ghostWidth) };
	static constexpr Range<int, dimensionCount> interiorBox { repeat<dimensionCount>(0), repeat<dimensionCount>(cellsAcrossInterior) };
	static constexpr auto vAnalyze = dgAnalyze<ValArray<Type, gridVolume>, dimensionCount, modeCount>;
	static constexpr auto vMassInverse = dgMassInverse<ValArray<Type, gridVolume>, dimensionCount, modeCount>;
	static constexpr auto vSynthesize = dgSynthesize<ValArray<Type, gridVolume>, dimensionCount, modeCount>;
	static constexpr auto vStiffness = dgStiffness<ValArray<Type, gridVolume>, dimensionCount, modeCount>;
	static constexpr auto sAnalyze = dgAnalyze<ValArray<Type, gridVolume>, dimensionCount - 1, modeCount>;
	static constexpr auto sSynthesize = dgSynthesize<ValArray<Type, gridVolume>, dimensionCount - 1, modeCount>;
	static constexpr auto sTrace = dgTrace<ValArray<Type, gridVolume>, dimensionCount, modeCount>;
	static constexpr auto sTraceInverse = dgTraceInverse<ValArray<Type, gridVolume>, dimensionCount, modeCount>;
	using GridType = State<std::array<ValArray<Type, gridVolume>, modeVolume>, dimensionCount>;
	Type cellWidth;
	GridType currentState;
	GridType previousState;
	ValArray<int64_t, dimensionCount> gridSizes;
	ValArray<int64_t, dimensionCount> gridStrides;
	std::array<GridType, rungeKuttaStageCount> stageDerivatives_;
	std::array<ValArray<Type, gridVolume>, dimensionCount> position;

	static constexpr auto interiorIndexMap(int i) {
		static constexpr auto map = createMultiIndexMap<exteriorBox, interiorBox>();
		return map[i];
	}
	static constexpr int stride(int d) {
		return ipow(cellsAcrossExterior, dimensionCount - 1 - d);
	}

public:
	HyperGrid(Type const &xNint = Type(1)) :
			cellWidth { xNint / Type(cellsAcrossInterior) }, currentState { }, previousState { }, gridSizes { {
					cellsAcrossExterior,
					cellsAcrossExterior,
					cellsAcrossExterior } }, gridStrides { }, stageDerivatives_ { }, position { } {
		gridStrides[dimensionCount - 1] = 1;
		printf( "%i\n", dimensionCount);
		for (int dimension = dimensionCount - 1; dimension > 0; dimension--) {
			gridStrides[dimension - 1] = gridStrides[dimension] * gridSizes[dimension];
		}
		for (int dimension = 0; dimension < dimensionCount; dimension++) {
			for (int i = 0; i < cellsAcrossExterior; i++) {
				auto const start = i * gridStrides[dimension];
				auto sizes = gridSizes;
				sizes[dimension] = 1;
				Type const x = Type(2 * (i - ghostWidth) + 1) * Type(0.5) * cellWidth;
				position[dimension][GSlice<dimensionCount> { start, sizes, gridStrides }] = x;
			}
		}
	}
	void initialize(std::function<State<Type, dimensionCount>(std::array<Type, dimensionCount> const&)> const &initialState) {
		Type const halfCellWidth = Type(0.5) * cellWidth;
		State<std::array<ValArray<Type, gridVolume>, nodeVolume>, dimensionCount> nodalValues;
		for (int nodeIndex = 0; nodeIndex < nodeVolume; nodeIndex++) {
			auto quadraturePoint = getQuadraturePoint<Type, dimensionCount, modeCount>(nodeIndex);
			auto thisPosition = position;
			for (int dimension = 0; dimension < dimensionCount; dimension++) {
				thisPosition[dimension] += quadraturePoint[dimension] * halfCellWidth;
			}
			for (int i = 0; i < gridVolume; i++) {
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
			auto const tmp = vAnalyze(nodalValues[fieldIndex]);
			currentState[fieldIndex] = vMassInverse(tmp);
		}
	}
	void output(const char *filenameBase, int timeStepNumber, Type const &time) {
		std::string filename = std::string(filenameBase) + "." + std::to_string(timeStepNumber) + ".h5";
		writeHdf5<Type, dimensionCount, cellsAcrossInterior, modeCount, ghostWidth>(filename, cellWidth, currentState,
				State<Type, dimensionCount>::getFieldNames());
		writeList("X.visit", "!NBLOCKS 1\n", filename + ".xmf");
	}
	void enforceBoundaryConditions() {
		for (int fieldIndex = 0; fieldIndex < fieldCount; fieldIndex++) {
			for (int modeIndex = 0; modeIndex < modeVolume; modeIndex++) {
				for (int dimension = 0; dimension < dimensionCount; dimension++) {
					GSlice<dimensionCount> fromSlice, toSlice;
					int fromStart, toStart;
					auto sizes = gridSizes;
					sizes[dimension] = ghostWidth;
					fromStart = 0;
					for (int d = 0; d <= dimension; d++) {
						fromStart += ghostWidth * gridStrides[d];
					}
					toStart = fromStart + gridStrides[dimension] * cellsAcrossInterior;
					fromSlice = GSlice<dimensionCount> { fromStart, sizes, gridStrides };
					toSlice = GSlice<dimensionCount> { toStart, sizes, gridStrides };
					currentState[fieldIndex][modeIndex][toSlice] = currentState[fieldIndex][modeIndex][fromSlice];
					toStart = 0;
					for (int d = 0; d < dimension; d++) {
						toStart += ghostWidth * gridStrides[d];
					}
					fromStart = toStart - gridStrides[dimension] * cellsAcrossInterior;
					fromSlice = GSlice<dimensionCount> { fromStart, sizes, gridStrides };
					toSlice = GSlice<dimensionCount> { toStart, sizes, gridStrides };
					currentState[fieldIndex][modeIndex][toSlice] = currentState[fieldIndex][modeIndex][fromSlice];
				}
			}
		}
	}
	Type beginStep() {
		applyLimiter();
		using namespace Math;
		previousState = currentState;
		for (int stageIndex = 0; stageIndex < rungeKuttaStageCount; stageIndex++) {
			for (int modeIndex = 0; modeIndex < modeVolume; modeIndex++) {
				for (int fieldIndex = 0; fieldIndex < modeVolume; fieldIndex++) {
					stageDerivatives_[stageIndex][fieldIndex][modeIndex].fill(Type(0));
				}
			}
		}
		std::array<Type, dimensionCount> maximumEigenvalue;
		maximumEigenvalue.fill(Type(0));
		std::array<State<ValArray<Type, gridVolume>, dimensionCount>, nodeVolume> nodalState;
		for (int fieldIndex = 0; fieldIndex < fieldCount; fieldIndex++) {
			auto const tmp = vSynthesize(currentState[fieldIndex]);
			for (int nodeIndex = 0; nodeIndex < nodeVolume; nodeIndex++) {
				nodalState[nodeIndex][fieldIndex] = tmp[nodeIndex];
			}
		}
		for (int nodeIndex = 0; nodeIndex < nodeVolume; nodeIndex++) {
			for (int dimension = 0; dimension < dimensionCount; dimension++) {
				auto const eigenvalues = nodalState[nodeIndex].eigenvalues(dimension);
				for (int eigenIndex = 0; eigenIndex < int(eigenvalues.size()); eigenIndex++) {
					ValArray<Type, gridVolume> const lambda = abs(eigenvalues[eigenIndex]);
					maximumEigenvalue[dimension] = max(maximumEigenvalue[dimension], lambda.max());
				}
			}
		}
		Type maximumEigenvalueSum = Type(0);
		for (int dimension = 0; dimension < dimensionCount; dimension++) {
			maximumEigenvalueSum += maximumEigenvalue[dimension];
		}
		Type const timeStepSize = (cellWidth * butcherTable.cfl()) / (Type(2 * modeCount - 1) * maximumEigenvalueSum);
		return timeStepSize;
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
		for (int fieldIndex = 0; fieldIndex < fieldCount; fieldIndex++) {
			for (int modeIndex = 0; modeIndex < modeVolume; modeIndex++) {
				currentState[fieldIndex][modeIndex] = previousState[fieldIndex][modeIndex];
				for (int stageIndex = 0; stageIndex < rungeKuttaStageCount; stageIndex++) {
					currentState[fieldIndex][modeIndex] += butcherTable.b(stageIndex) * stageDerivatives_[stageIndex][fieldIndex][modeIndex];
				}
			}
		}
		stageDerivatives_ = { };
		previousState = { };
	}
	void applyLimiter() {
		State<ValArray<Type, gridVolume>, dimensionCount> meanState;
		for (int fieldIndex = 0; fieldIndex < fieldCount; fieldIndex++) {
			for (int i = 0; i < gridVolume; i++) {
				meanState[fieldIndex][i] = currentState[fieldIndex][0][i];
			}
		}
		State<ValArray<Type, gridVolume>, dimensionCount> differenceState;
		for (int polynomialDegree = modeCount - 1; polynomialDegree > 0; polynomialDegree--) {
			int const begin = binco(polynomialDegree + dimensionCount - 1, dimensionCount);
			int const end = binco(polynomialDegree + dimensionCount, dimensionCount);
			for (int hiModeIndex = begin; hiModeIndex < end; hiModeIndex++) {
				std::array<State<ValArray<Type, gridVolume>, dimensionCount>, dimensionCount> transposedState;
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
						differenceState[fieldIndex] = minmod(ValArray<Type, gridVolume>(statePlus - stateCentral),
								ValArray<Type, gridVolume>(stateCentral - stateMinus));
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
						auto &current = currentState[fieldIndex][hiModeIndex];
						current = minmod(current, transposedState[dimension][fieldIndex]);
					}
				}
			}
		}
		ValArray<Type, gridVolume> theta { Type(1), gridVolume };
		std::array<State<ValArray<Type, gridVolume>, dimensionCount>, nodeSurface> surfaceState;
		std::array<State<ValArray<Type, gridVolume>, dimensionCount>, nodeVolume> volumeState;
		for (int face = 0; face < 2 * dimensionCount; face++) {
			for (int fieldIndex = 0; fieldIndex < fieldCount; fieldIndex++) {
				auto const surfaceNodes = sSynthesize(sTrace(face, currentState[fieldIndex]));
				for (int nodeIndex = 0; nodeIndex < nodeSurface; nodeIndex++) {
					surfaceState[nodeIndex][fieldIndex] = surfaceNodes[nodeIndex];
				}
			}
			for (int nodeIndex = 0; nodeIndex < nodeSurface; nodeIndex++) {
				theta = min(Type(1), findPositivityPreservingTheta(meanState, surfaceState[nodeIndex]));
			}
		}
		for (int fieldIndex = 0; fieldIndex < fieldCount; fieldIndex++) {
			auto const volumeNodes = vSynthesize(currentState[fieldIndex]);
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
	void computeTimeDerivative(Type timeStepSize, State<std::array<ValArray<Type, gridVolume>, modeVolume>, dimensionCount> &stateDerivative) {
		Type const lambda = Type(2) * timeStepSize / cellWidth;
		for (int dimension = 0; dimension < dimensionCount; dimension++) {
			auto const shift = gridStrides[dimension];
			std::array<State<ValArray<Type, gridVolume>, dimensionCount>, nodeSurface> leftState;
			std::array<State<ValArray<Type, gridVolume>, dimensionCount>, nodeSurface> rightState;
			State<std::array<ValArray<Type, gridVolume>, nodeSurface>, dimensionCount> riemannFlux;
			for (int fieldIndex = 0; fieldIndex < fieldCount; fieldIndex++) {
				auto leftModes = currentState[fieldIndex];
				for (int modeIndex = 0; modeIndex < modeCount; modeIndex++) {
					leftModes[modeIndex] = currentState[fieldIndex][modeIndex].cshift(-shift);
				}
				auto const left = sSynthesize(sTrace(2 * dimension + 1, leftModes));
				auto const right = sSynthesize(sTrace(2 * dimension + 0, currentState[fieldIndex]));
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
				auto const thisFlux = sAnalyze(riemannFlux[fieldIndex]);
				auto const leftFlux = vMassInverse(sTraceInverse(2 * dimension + 1, thisFlux));
				auto const rightFlux = vMassInverse(sTraceInverse(2 * dimension + 0, thisFlux));
				for (int modeIndex = 0; modeIndex < modeVolume; modeIndex++) {
					stateDerivative[fieldIndex][modeIndex] -= lambda * leftFlux[modeIndex].cshift(+shift);
					stateDerivative[fieldIndex][modeIndex] += lambda * rightFlux[modeIndex];
				}
			}
		}
		for (int dimension = 0; dimension < dimensionCount; dimension++) {
			std::array<State<ValArray<Type, gridVolume>, dimensionCount>, nodeVolume> volumeState;
			State<std::array<ValArray<Type, gridVolume>, nodeVolume>, dimensionCount> volumeFlux;
			for (int fieldIndex = 0; fieldIndex < fieldCount; fieldIndex++) {
				auto const nodes = vSynthesize(currentState[fieldIndex]);
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
				auto const source = vMassInverse(vStiffness(dimension, vMassInverse(vAnalyze(volumeFlux[fieldIndex]))));
				for (int modeIndex = 0; modeIndex < modeVolume; modeIndex++) {
					stateDerivative[fieldIndex][modeIndex] += lambda * source[modeIndex];
				}
			}
		}
	}
}
;

