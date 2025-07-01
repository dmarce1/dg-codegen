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
	using EleType = typename ElementType<Type>::type;
	constexpr EleType half = EleType(1) / EleType(2);
	Type const sgn = copysign(half, a) + copysign(half, b);
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
	static constexpr int interiorVolume = ipow(cellsAcrossInterior, dimensionCount);
	static constexpr RungeKutta butcherTable { };
	static constexpr Range<int, dimensionCount> exteriorBox { repeat<dimensionCount>(-ghostWidth), repeat<dimensionCount>(cellsAcrossInterior + ghostWidth) };
	static constexpr Range<int, dimensionCount> interiorBox { repeat<dimensionCount>(0), repeat<dimensionCount>(cellsAcrossInterior) };
	template<int volume = exteriorVolume>
	static constexpr auto vAnalyze = dgAnalyze<ValArray<Type, volume>, dimensionCount, modeCount>;
	template<int volume = exteriorVolume>
	static constexpr auto vMassInverse = dgMassInverse<ValArray<Type, volume>, dimensionCount, modeCount>;
	template<int volume = exteriorVolume>
	static constexpr auto vSynthesize = dgSynthesize<ValArray<Type, volume>, dimensionCount, modeCount>;
	template<int volume = exteriorVolume>
	static constexpr auto vStiffness = dgStiffness<ValArray<Type, volume>, dimensionCount, modeCount>;
	template<int volume = exteriorVolume>
	static constexpr auto sAnalyze = dgAnalyze<ValArray<Type, volume>, dimensionCount - 1, modeCount>;
	template<int volume = exteriorVolume>
	static constexpr auto sSynthesize = dgSynthesize<ValArray<Type, volume>, dimensionCount - 1, modeCount>;
	template<int volume = exteriorVolume>
	static constexpr auto sTrace = dgTrace<ValArray<Type, volume>, dimensionCount, modeCount>;
	template<int volume = exteriorVolume>
	static constexpr auto sTraceInverse = dgTraceInverse<ValArray<Type, volume>, dimensionCount, modeCount>;
	using GridType = State<std::array<ValArray<Type, exteriorVolume>, modeVolume>, dimensionCount>;
	Type cellWidth;
	GridType currentState;
	GridType previousState;
	int64_t interiorStart;
	GSlice<dimensionCount, interiorVolume> interiorSlice;
	ValArray<int64_t, dimensionCount> interiorSizes;
	ValArray<int64_t, dimensionCount> exteriorSizes;
	ValArray<int64_t, dimensionCount> gridStrides;
	std::array<GridType, rungeKuttaStageCount> stageDerivatives_;
	std::array<ValArray<Type, exteriorVolume>, dimensionCount> position;

	static constexpr auto interiorIndexMap(int i) {
		static constexpr auto map = createMultiIndexMap<exteriorBox, interiorBox>();
		return map[i];
	}
	static constexpr int stride(int d) {
		return ipow(cellsAcrossExterior, dimensionCount - 1 - d);
	}

public:
	HyperGrid(Type const &xNint = Type(1)) :
			cellWidth { xNint / Type(cellsAcrossInterior) }, currentState { }, previousState { }, exteriorSizes { }, gridStrides { }, stageDerivatives_ { }, position { } {
		exteriorSizes.fill(cellsAcrossExterior);
		interiorSizes.fill(cellsAcrossInterior);
		gridStrides[dimensionCount - 1] = 1;
		for (int dimension = dimensionCount - 1; dimension > 0; dimension--) {
			gridStrides[dimension - 1] = gridStrides[dimension] * exteriorSizes[dimension];
		}
		interiorStart = 2 * std::accumulate(gridStrides.begin(), gridStrides.end(), int64_t(0));
		interiorSlice = GSlice<dimensionCount, interiorVolume>(interiorStart, interiorSizes, gridStrides);
		constexpr int64_t sliceSize = exteriorVolume / cellsAcrossExterior;
		for (int dimension = 0; dimension < dimensionCount; dimension++) {
			for (int i = 0; i < cellsAcrossExterior; i++) {
				auto const start = i * gridStrides[dimension];
				auto sizes = exteriorSizes;
				sizes[dimension] = 1;
				Type const x = Type(2 * (i - ghostWidth) + 1) * Type(0.5) * cellWidth;
				position[dimension][GSlice<dimensionCount, sliceSize> { start, sizes, gridStrides }] = x;
			}
		}
	}
	void initialize(std::function<State<Type, dimensionCount>(std::array<Type, dimensionCount> const&)> const &initialState) {
		Type const halfCellWidth = Type(0.5) * cellWidth;
		State<std::array<ValArray<Type, exteriorVolume>, nodeVolume>, dimensionCount> nodalValues;
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
			auto const tmp = vAnalyze<>(nodalValues[fieldIndex]);
			currentState[fieldIndex] = vMassInverse<>(tmp);
		}
	}
	void output(const char *filenameBase, int timeStepNumber, Type const &time) {
		std::string filename = std::string(filenameBase) + "." + std::to_string(timeStepNumber) + ".h5";
		writeHdf5<Type, dimensionCount, cellsAcrossInterior, modeCount, ghostWidth, State>(filename, cellWidth, currentState,
				State<Type, dimensionCount>::getFieldNames());
		writeList("X.visit", "!NBLOCKS 1\n", filename + ".xmf");
	}
	void enforceBoundaryConditions() {
		constexpr int boundaryVolume = exteriorVolume * ghostWidth / cellsAcrossExterior;
		using SliceType = GSlice<dimensionCount, boundaryVolume>;
		for (int fieldIndex = 0; fieldIndex < fieldCount; fieldIndex++) {
			for (int modeIndex = 0; modeIndex < modeVolume; modeIndex++) {
				for (int dimension = 0; dimension < dimensionCount; dimension++) {
					SliceType fromSlice, toSlice;
					int fromStart, toStart;
					auto sizes = exteriorSizes;
					sizes[dimension] = ghostWidth;
					fromStart = ghostWidth * gridStrides[dimension];
					toStart = fromStart + gridStrides[dimension] * cellsAcrossInterior;
					fromSlice = SliceType(fromStart, sizes, gridStrides);
					toSlice = SliceType(toStart, sizes, gridStrides);
					ValArray<Type, boundaryVolume> tmp = currentState[fieldIndex][modeIndex][fromSlice];
					currentState[fieldIndex][modeIndex][toSlice] = tmp;
					fromStart = cellsAcrossInterior * gridStrides[dimension];
					toStart = 0;
					fromSlice = SliceType(fromStart, sizes, gridStrides);
					toSlice = SliceType(toStart, sizes, gridStrides);
					tmp = currentState[fieldIndex][modeIndex][fromSlice];
					currentState[fieldIndex][modeIndex][toSlice] = tmp;
				}
			}
		}
	}
	Type beginStep() {
		applyLimiter();
		using namespace Math;
		previousState = currentState;
		std::array<Type, dimensionCount> maximumEigenvalue;
		maximumEigenvalue.fill(Type(0));
		std::array<State<ValArray<Type, exteriorVolume>, dimensionCount>, nodeVolume> nodalState;
		for (int fieldIndex = 0; fieldIndex < fieldCount; fieldIndex++) {
			auto const tmp = vSynthesize<>(currentState[fieldIndex]);
			for (int nodeIndex = 0; nodeIndex < nodeVolume; nodeIndex++) {
				nodalState[nodeIndex][fieldIndex] = tmp[nodeIndex];
			}
		}
		for (int nodeIndex = 0; nodeIndex < nodeVolume; nodeIndex++) {
			for (int dimension = 0; dimension < dimensionCount; dimension++) {
				auto const eigenvalues = nodalState[nodeIndex].eigenvalues(dimension);
				for (int eigenIndex = 0; eigenIndex < int(eigenvalues.size()); eigenIndex++) {
					ValArray<Type, exteriorVolume> const lambda = abs(eigenvalues[eigenIndex]);
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
	}
	void applyLimiter() {
		constexpr int limiterVolume = ipow(cellsAcrossInterior + 2, dimensionCount);
		using LimiterArray = ValArray<Type, limiterVolume>;
		ValArray<int64_t, dimensionCount> limiterSizes;
		limiterSizes.fill(cellsAcrossInterior + 2);
		int64_t limiterStart = std::accumulate(gridStrides.begin(), gridStrides.end(), int64_t(0));
		GSlice<dimensionCount, limiterVolume> const limiterSlice(limiterStart, limiterSizes, gridStrides);
		State<LimiterArray, dimensionCount> meanState;
		for (int fieldIndex = 0; fieldIndex < fieldCount; fieldIndex++) {
			meanState[fieldIndex] = currentState[fieldIndex][0][limiterSlice];
		}
		std::array<SquareMatrix<LimiterArray, fieldCount>, dimensionCount> leftEigenvectors;
		std::array<SquareMatrix<LimiterArray, fieldCount>, dimensionCount> rightEigenvectors;
		for (int dimension = 0; dimension < dimensionCount; dimension++) {
			rightEigenvectors[dimension] = meanState.eigenSystem(dimension).second;
			leftEigenvectors[dimension] = matrixInverse(rightEigenvectors[dimension]);
		}
		for (int polynomialDegree = modeCount - 1; polynomialDegree > 0; polynomialDegree--) {
			int const begin = binco(polynomialDegree + dimensionCount - 1, dimensionCount);
			int const end = binco(polynomialDegree + dimensionCount, dimensionCount);
			for (int hiModeIndex = begin; hiModeIndex < end; hiModeIndex++) {
				State<LimiterArray, dimensionCount> differenceState;
				State<LimiterArray, dimensionCount> transposedState;
				auto hiModeIndices = flatToTriangular<dimensionCount, modeCount>(hiModeIndex);
				for (int fieldIndex = 0; fieldIndex < fieldCount; fieldIndex++) {
					transposedState[fieldIndex] = currentState[fieldIndex][hiModeIndex][limiterSlice];
				}
				for (int dimension = 0; dimension < dimensionCount; dimension++) {
					GSlice<dimensionCount, limiterVolume> const limiterPlusSlice(limiterStart + gridStrides[dimension], limiterSizes, gridStrides);
					GSlice<dimensionCount, limiterVolume> const limiterMinusSlice(limiterStart - gridStrides[dimension], limiterSizes, gridStrides);
					if (hiModeIndices[dimension] == 0) {
						continue;
					}
					auto loModeIndices = hiModeIndices;
					loModeIndices[dimension]--;
					auto const loModeIndex = triangularToFlat<dimensionCount, modeCount>(loModeIndices);
					for (int fieldIndex = 0; fieldIndex < fieldCount; fieldIndex++) {
						LimiterArray const stateCentral = currentState[fieldIndex][loModeIndex][limiterSlice];
						LimiterArray const statePlus = currentState[fieldIndex][loModeIndex][limiterPlusSlice];
						LimiterArray const stateMinus = currentState[fieldIndex][loModeIndex][limiterMinusSlice];
						differenceState[fieldIndex] = minmod(LimiterArray(statePlus - stateCentral), LimiterArray(stateCentral - stateMinus));
					}
					//differenceState = leftEigenvectors[dimension] * differenceState;
					//transposedState = leftEigenvectors[dimension] * transposedState;
					auto const cLimit = Type(1) / Type(2 * loModeIndices[dimension] + 1);
					for (int fieldIndex = 0; fieldIndex < fieldCount; fieldIndex++) {
						differenceState[fieldIndex] *= cLimit;
						transposedState[fieldIndex] = minmod(transposedState[fieldIndex], differenceState[fieldIndex]);
					}
				//	transposedState = rightEigenvectors[dimension] * transposedState;
				}
				for (int fieldIndex = 0; fieldIndex < fieldCount; fieldIndex++) {
					currentState[fieldIndex][hiModeIndex][limiterSlice] = transposedState[fieldIndex];
				}
			}
		}

//		ValArray<Type, exteriorVolume> theta(Type(1));
//		std::array<State<ValArray<Type, exteriorVolume>, dimensionCount>, nodeSurface> surfaceState;
//		std::array<State<ValArray<Type, exteriorVolume>, dimensionCount>, nodeVolume> volumeState;
//		for (int face = 0; face < 2 * dimensionCount; face++) {
//			for (int fieldIndex = 0; fieldIndex < fieldCount; fieldIndex++) {
//				auto const surfaceNodes = sSynthesize(sTrace(face, currentState[fieldIndex]));
//				for (int nodeIndex = 0; nodeIndex < nodeSurface; nodeIndex++) {
//					surfaceState[nodeIndex][fieldIndex] = surfaceNodes[nodeIndex];
//				}
//			}
//			for (int nodeIndex = 0; nodeIndex < nodeSurface; nodeIndex++) {
//				theta = min(Type(1), findPositivityPreservingTheta(meanState, surfaceState[nodeIndex]));
//			}
//		}
//		for (int fieldIndex = 0; fieldIndex < fieldCount; fieldIndex++) {
//			auto const volumeNodes = vSynthesize(currentState[fieldIndex]);
//			for (int nodeIndex = 0; nodeIndex < nodeVolume; nodeIndex++) {
//				volumeState[nodeIndex][fieldIndex] = volumeNodes[nodeIndex];
//			}
//		}
//		for (int nodeIndex = 0; nodeIndex < nodeVolume; nodeIndex++) {
//			theta = min(Type(1), findPositivityPreservingTheta(meanState, volumeState[nodeIndex]));
//		}
//		for (int fieldIndex = 0; fieldIndex < fieldCount; fieldIndex++) {
//			for (int modeIndex = 1; modeIndex < modeVolume; modeIndex++) {
//				currentState[fieldIndex][modeIndex] *= theta;
//			}
//		}
	}
	void computeTimeDerivative(Type timeStepSize, State<std::array<ValArray<Type, exteriorVolume>, modeVolume>, dimensionCount> &stateDerivative) {
		constexpr int fluxVolume = interiorVolume + interiorVolume / cellsAcrossInterior;
		Type const lambda = Type(2) * timeStepSize / cellWidth;
		for (int dimension = 0; dimension < dimensionCount; dimension++) {
			auto const thisGridStride = gridStrides[dimension];
			auto fluxSizes = interiorSizes;
			fluxSizes[dimension]++;
			ValArray<int64_t, dimensionCount> fluxStrides;
			int stride = 1;
			for (int thisDimension = dimensionCount - 1; thisDimension >= 0; thisDimension--) {
				fluxStrides[thisDimension] = stride;
				stride *= fluxSizes[thisDimension];
			}
			auto const thisFluxStride = fluxStrides[dimension];
			GSlice<dimensionCount, fluxVolume> const fluxLeftSlice(interiorStart - thisGridStride, fluxSizes, gridStrides);
			GSlice<dimensionCount, fluxVolume> const fluxRightSlice(interiorStart, fluxSizes, gridStrides);
			std::array<State<ValArray<Type, fluxVolume>, dimensionCount>, nodeSurface> leftState;
			std::array<State<ValArray<Type, fluxVolume>, dimensionCount>, nodeSurface> rightState;
			for (int fieldIndex = 0; fieldIndex < fieldCount; fieldIndex++) {
				std::array<ValArray<Type, fluxVolume>, modeVolume> leftModes;
				std::array<ValArray<Type, fluxVolume>, modeVolume> rightModes;
				for (int modeIndex = 0; modeIndex < modeVolume; modeIndex++) {
					leftModes[modeIndex] = currentState[fieldIndex][modeIndex][fluxLeftSlice];
					rightModes[modeIndex] = currentState[fieldIndex][modeIndex][fluxRightSlice];
				}
				auto const left = sSynthesize<fluxVolume>(sTrace<fluxVolume>(2 * dimension + 1, leftModes));
				auto const right = sSynthesize<fluxVolume>(sTrace<fluxVolume>(2 * dimension + 0, rightModes));
				for (int nodeIndex = 0; nodeIndex < nodeSurface; nodeIndex++) {
					leftState[nodeIndex][fieldIndex] = left[nodeIndex];
					rightState[nodeIndex][fieldIndex] = right[nodeIndex];
				}
			}
			State<std::array<ValArray<Type, fluxVolume>, nodeSurface>, dimensionCount> nodalFlux;
			State<std::array<ValArray<Type, interiorVolume>, nodeSurface>, dimensionCount> leftFlux;
			State<std::array<ValArray<Type, interiorVolume>, nodeSurface>, dimensionCount> rightFlux;
			for (int nodeIndex = 0; nodeIndex < nodeSurface; nodeIndex++) {
				auto const tmp = solveRiemannProblem(leftState[nodeIndex], rightState[nodeIndex], dimension);
				for (int fieldIndex = 0; fieldIndex < fieldCount; fieldIndex++) {
					nodalFlux[fieldIndex][nodeIndex] = tmp[fieldIndex];
				}
			}
			GSlice<dimensionCount, interiorVolume> const fluxPlusSlice(thisFluxStride, interiorSizes, fluxStrides);
			GSlice<dimensionCount, interiorVolume> const fluxMinusSlice(0, interiorSizes, fluxStrides);
			for (int fieldIndex = 0; fieldIndex < fieldCount; fieldIndex++) {
				auto const modalFlux = sAnalyze<fluxVolume>(nodalFlux[fieldIndex]);
				std::array<ValArray<Type, interiorVolume>, modeSurface> modalSurfaceFluxPlus;
				std::array<ValArray<Type, interiorVolume>, modeSurface> modalSurfaceFluxMinus;
				std::array<ValArray<Type, interiorVolume>, modeVolume> modalVolumeFluxPlus;
				std::array<ValArray<Type, interiorVolume>, modeVolume> modalVolumeFluxMinus;
				for (int modeIndex = 0; modeIndex < modeSurface; modeIndex++) {
					modalSurfaceFluxPlus[modeIndex] = modalFlux[modeIndex][fluxPlusSlice];
					modalSurfaceFluxMinus[modeIndex] = modalFlux[modeIndex][fluxMinusSlice];
				}
				modalVolumeFluxPlus = sTraceInverse<interiorVolume>(2 * dimension + 1, modalSurfaceFluxPlus);
				modalVolumeFluxMinus = sTraceInverse<interiorVolume>(2 * dimension, modalSurfaceFluxMinus);
				modalVolumeFluxPlus = vMassInverse<interiorVolume>(modalVolumeFluxPlus);
				modalVolumeFluxMinus = vMassInverse<interiorVolume>(modalVolumeFluxMinus);
				for (int modeIndex = 0; modeIndex < modeVolume; modeIndex++) {
					stateDerivative[fieldIndex][modeIndex][interiorSlice] -= lambda * (modalVolumeFluxPlus[modeIndex] - modalVolumeFluxMinus[modeIndex]);
				}
			}
		}
		for (int dimension = 0; dimension < dimensionCount; dimension++) {
			std::array<State<ValArray<Type, interiorVolume>, dimensionCount>, nodeVolume> volumeState;
			State<std::array<ValArray<Type, interiorVolume>, nodeVolume>, dimensionCount> volumeFlux;
			for (int fieldIndex = 0; fieldIndex < fieldCount; fieldIndex++) {
				std::array<ValArray<Type, interiorVolume>, modeVolume> volumeModes;
				for (int modeIndex = 0; modeIndex < modeVolume; modeIndex++) {
					volumeModes[modeIndex] = currentState[fieldIndex][modeIndex][interiorSlice];
				}
				auto const nodes = vSynthesize<interiorVolume>(volumeModes);
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
				auto const tmp = vMassInverse<interiorVolume>(vAnalyze<interiorVolume>(volumeFlux[fieldIndex]));
				auto const source = vMassInverse<interiorVolume>(vStiffness<interiorVolume>(dimension, tmp));
				for (int modeIndex = 0; modeIndex < modeVolume; modeIndex++) {
					stateDerivative[fieldIndex][modeIndex][interiorSlice] += lambda * source[modeIndex];
				}
			}
		}
	}
}
;

