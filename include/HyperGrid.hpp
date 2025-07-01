#pragma once
#include "Basis.hpp"
#include "ContainerArithmetic.hpp"
#include "Hdf5.hpp"
#include "Matrix.hpp"
#include "MultiIndex.hpp"
#include "Quadrature.hpp"
#include "dgTransforms.hpp"
#include "ValArray.hpp"

#include <bitset>
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
	using EleType = typename ElementType<Type>::type;
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
		for (int node = 0; node < nodeVolume; node++) {
			auto quadraturePoint = getQuadraturePoint<Type, dimensionCount, modeCount>(node);
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
				for (int field = 0; field < fieldCount; field++) {
					nodalValues[field][node][i] = state[field];
				}
			}
		}
		for (int field = 0; field < fieldCount; field++) {
			auto const tmp = vAnalyze<>(nodalValues[field]);
			currentState[field] = vMassInverse<>(tmp);
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
		for (int field = 0; field < fieldCount; field++) {
			for (int mode = 0; mode < modeVolume; mode++) {
				for (int dimension = 0; dimension < dimensionCount; dimension++) {
					SliceType fromSlice, toSlice;
					int fromStart, toStart;
					auto sizes = exteriorSizes;
					sizes[dimension] = ghostWidth;
					fromStart = ghostWidth * gridStrides[dimension];
					toStart = fromStart + gridStrides[dimension] * cellsAcrossInterior;
					fromSlice = SliceType(fromStart, sizes, gridStrides);
					toSlice = SliceType(toStart, sizes, gridStrides);
					ValArray<Type, boundaryVolume> tmp = currentState[field][mode][fromSlice];
					currentState[field][mode][toSlice] = tmp;
					fromStart = cellsAcrossInterior * gridStrides[dimension];
					toStart = 0;
					fromSlice = SliceType(fromStart, sizes, gridStrides);
					toSlice = SliceType(toStart, sizes, gridStrides);
					tmp = currentState[field][mode][fromSlice];
					currentState[field][mode][toSlice] = tmp;
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
		for (int field = 0; field < fieldCount; field++) {
			auto const tmp = vSynthesize<>(currentState[field]);
			for (int node = 0; node < nodeVolume; node++) {
				nodalState[node][field] = tmp[node];
			}
		}
		for (int node = 0; node < nodeVolume; node++) {
			for (int dimension = 0; dimension < dimensionCount; dimension++) {
				auto const eigenvalues = nodalState[node].eigenvalues(dimension);
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
		for (int field = 0; field < fieldCount; field++) {
			for (int mode = 0; mode < modeVolume; mode++) {
				for (int thisStage = 0; thisStage < stageIndex; thisStage++) {
					currentState[field][mode] += butcherTable.a(stageIndex, thisStage) * stageDerivatives_[thisStage][field][mode];
				}
				stageDerivatives_[stageIndex][field][mode].fill(Type(0));
			}
		}
		enforceBoundaryConditions();
		applyLimiter();
		computeTimeDerivative(timeStepSize, stageDerivatives_[stageIndex]);
	}
	void endStep() {
		for (int field = 0; field < fieldCount; field++) {
			for (int mode = 0; mode < modeVolume; mode++) {
				currentState[field][mode] = previousState[field][mode];
				for (int stageIndex = 0; stageIndex < rungeKuttaStageCount; stageIndex++) {
					currentState[field][mode] += butcherTable.b(stageIndex) * stageDerivatives_[stageIndex][field][mode];
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
		std::array<std::array<State<LimiterArray, dimensionCount>, modeVolume>, dimensionCount> alpha;
		State<LimiterArray, dimensionCount> meanState;
		for (int field = 0; field < fieldCount; field++) {
			meanState[field] = currentState[field][0][limiterSlice];
		}
		std::array<SquareMatrix<LimiterArray, fieldCount>, dimensionCount> leftEigenvectors;
		std::array<SquareMatrix<LimiterArray, fieldCount>, dimensionCount> rightEigenvectors;
		std::array<std::bitset<modeVolume>, dimensionCount> wasLimited;
		for (int dimension = 0; dimension < dimensionCount; dimension++) {
			rightEigenvectors[dimension] = meanState.eigenSystem(dimension).second;
			leftEigenvectors[dimension] = matrixInverse(rightEigenvectors[dimension]);
		}
		for (int dimension = 0; dimension < dimensionCount; dimension++) {
			for (int mode = 0; mode < modeVolume; mode++) {
				for (int field = 0; field < fieldCount; field++) {
					alpha[dimension][mode][field] = Type(0.0);
				}
			}
		}
		for (int polynomialDegree = modeCount - 1; polynomialDegree > 0; polynomialDegree--) {
			int const begin = binco(polynomialDegree + dimensionCount - 1, dimensionCount);
			int const end = binco(polynomialDegree + dimensionCount, dimensionCount);
			for (int hiMode = begin; hiMode < end; hiMode++) {
				State<LimiterArray, dimensionCount> differenceState;
				State<LimiterArray, dimensionCount> transposedState;
				auto hiModeIndices = flatToTriangular<dimensionCount, modeCount>(hiMode);
				for (int dimension = 0; dimension < dimensionCount; dimension++) {
					if (hiModeIndices[dimension] == 0) {
						continue;
					}
					for (int field = 0; field < fieldCount; field++) {
						transposedState[field] = currentState[field][hiMode][limiterSlice];
					}
					GSlice<dimensionCount, limiterVolume> const limiterPlusSlice(limiterStart + gridStrides[dimension], limiterSizes, gridStrides);
					GSlice<dimensionCount, limiterVolume> const limiterMinusSlice(limiterStart - gridStrides[dimension], limiterSizes, gridStrides);
					auto loModeIndices = hiModeIndices;
					loModeIndices[dimension]--;
					auto const loMode = triangularToFlat<dimensionCount, modeCount>(loModeIndices);
					for (int field = 0; field < fieldCount; field++) {
						LimiterArray const stateCentral = currentState[field][loMode][limiterSlice];
						LimiterArray const statePlus = currentState[field][loMode][limiterPlusSlice];
						LimiterArray const stateMinus = currentState[field][loMode][limiterMinusSlice];
						differenceState[field] = minmod(LimiterArray(statePlus - stateCentral), LimiterArray(stateCentral - stateMinus));
					}
					constexpr auto tiny = EleType(std::numeric_limits<double>::min());
					State<LimiterArray, dimensionCount> initialStateInverse;
					for (int field = 0; field < fieldCount; field++) {
						initialStateInverse[field] = EleType(1) / (transposedState[field] + tiny);
					}
					differenceState = leftEigenvectors[dimension] * differenceState;
					transposedState = leftEigenvectors[dimension] * transposedState;
					auto const cLimit = Type(1) / Type(2 * loModeIndices[dimension] + 1);
					for (int field = 0; field < fieldCount; field++) {
						differenceState[field] *= cLimit;
						transposedState[field] = minmod(transposedState[field], differenceState[field]);
					}
					transposedState = rightEigenvectors[dimension] * transposedState;
					for (int field = 0; field < fieldCount; field++) {
						alpha[dimension][hiMode][field] =
						alpha[dimension][loMode][field] =
								max(max(EleType(0.0), transposedState[field] * initialStateInverse[field]), alpha[dimension][hiMode][field]);;
					}
					wasLimited[dimension][hiMode] = true;
				}
			}
		}
		for (int mode = 1; mode < modeVolume; mode++) {
			for (int field = 0; field < fieldCount; field++) {
				LimiterArray thisAlpha = EleType(0);
				for (int dimension = 0; dimension < dimensionCount; dimension++) {
					if(wasLimited[dimension][mode]) {
						thisAlpha = min(thisAlpha, alpha[dimension][mode][field]);
					}
				}
				currentState[field][mode][limiterSlice] *= thisAlpha;
			}
		}
//		ValArray<Type, exteriorVolume> theta(Type(1));
//		std::array<State<ValArray<Type, exteriorVolume>, dimensionCount>, nodeSurface> surfaceState;
//		std::array<State<ValArray<Type, exteriorVolume>, dimensionCount>, nodeVolume> volumeState;
//		for (int face = 0; face < 2 * dimensionCount; face++) {
//			for (int field = 0; field < fieldCount; field++) {
//				auto const surfaceNodes = sSynthesize(sTrace(face, currentState[field]));
//				for (int node = 0; node < nodeSurface; node++) {
//					surfaceState[node][field] = surfaceNodes[node];
//				}
//			}
//			for (int node = 0; node < nodeSurface; node++) {
//				theta = min(Type(1), findPositivityPreservingTheta(meanState, surfaceState[node]));
//			}
//		}
//		for (int field = 0; field < fieldCount; field++) {
//			auto const volumeNodes = vSynthesize(currentState[field]);
//			for (int node = 0; node < nodeVolume; node++) {
//				volumeState[node][field] = volumeNodes[node];
//			}
//		}
//		for (int node = 0; node < nodeVolume; node++) {
//			theta = min(Type(1), findPositivityPreservingTheta(meanState, volumeState[node]));
//		}
//		for (int field = 0; field < fieldCount; field++) {
//			for (int mode = 1; mode < modeVolume; mode++) {
//				currentState[field][mode] *= theta;
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
			for (int field = 0; field < fieldCount; field++) {
				std::array<ValArray<Type, fluxVolume>, modeVolume> leftModes;
				std::array<ValArray<Type, fluxVolume>, modeVolume> rightModes;
				for (int mode = 0; mode < modeVolume; mode++) {
					leftModes[mode] = currentState[field][mode][fluxLeftSlice];
					rightModes[mode] = currentState[field][mode][fluxRightSlice];
				}
				auto const left = sSynthesize<fluxVolume>(sTrace<fluxVolume>(2 * dimension + 1, leftModes));
				auto const right = sSynthesize<fluxVolume>(sTrace<fluxVolume>(2 * dimension + 0, rightModes));
				for (int node = 0; node < nodeSurface; node++) {
					leftState[node][field] = left[node];
					rightState[node][field] = right[node];
				}
			}
			State<std::array<ValArray<Type, fluxVolume>, nodeSurface>, dimensionCount> nodalFlux;
			State<std::array<ValArray<Type, interiorVolume>, nodeSurface>, dimensionCount> leftFlux;
			State<std::array<ValArray<Type, interiorVolume>, nodeSurface>, dimensionCount> rightFlux;
			for (int node = 0; node < nodeSurface; node++) {
				auto const tmp = solveRiemannProblem(leftState[node], rightState[node], dimension);
				for (int field = 0; field < fieldCount; field++) {
					nodalFlux[field][node] = tmp[field];
				}
			}
			GSlice<dimensionCount, interiorVolume> const fluxPlusSlice(thisFluxStride, interiorSizes, fluxStrides);
			GSlice<dimensionCount, interiorVolume> const fluxMinusSlice(0, interiorSizes, fluxStrides);
			for (int field = 0; field < fieldCount; field++) {
				auto const modalFlux = sAnalyze<fluxVolume>(nodalFlux[field]);
				std::array<ValArray<Type, interiorVolume>, modeSurface> modalSurfaceFluxPlus;
				std::array<ValArray<Type, interiorVolume>, modeSurface> modalSurfaceFluxMinus;
				std::array<ValArray<Type, interiorVolume>, modeVolume> modalVolumeFluxPlus;
				std::array<ValArray<Type, interiorVolume>, modeVolume> modalVolumeFluxMinus;
				for (int mode = 0; mode < modeSurface; mode++) {
					modalSurfaceFluxPlus[mode] = modalFlux[mode][fluxPlusSlice];
					modalSurfaceFluxMinus[mode] = modalFlux[mode][fluxMinusSlice];
				}
				modalVolumeFluxPlus = sTraceInverse<interiorVolume>(2 * dimension + 1, modalSurfaceFluxPlus);
				modalVolumeFluxMinus = sTraceInverse<interiorVolume>(2 * dimension, modalSurfaceFluxMinus);
				modalVolumeFluxPlus = vMassInverse<interiorVolume>(modalVolumeFluxPlus);
				modalVolumeFluxMinus = vMassInverse<interiorVolume>(modalVolumeFluxMinus);
				for (int mode = 0; mode < modeVolume; mode++) {
					stateDerivative[field][mode][interiorSlice] -= lambda * (modalVolumeFluxPlus[mode] - modalVolumeFluxMinus[mode]);
				}
			}
		}
		for (int dimension = 0; dimension < dimensionCount; dimension++) {
			std::array<State<ValArray<Type, interiorVolume>, dimensionCount>, nodeVolume> volumeState;
			State<std::array<ValArray<Type, interiorVolume>, nodeVolume>, dimensionCount> volumeFlux;
			for (int field = 0; field < fieldCount; field++) {
				std::array<ValArray<Type, interiorVolume>, modeVolume> volumeModes;
				for (int mode = 0; mode < modeVolume; mode++) {
					volumeModes[mode] = currentState[field][mode][interiorSlice];
				}
				auto const nodes = vSynthesize<interiorVolume>(volumeModes);
				for (int node = 0; node < nodeVolume; node++) {
					volumeState[node][field] = nodes[node];
				}
			}
			for (int node = 0; node < nodeVolume; node++) {
				auto tmp = volumeState[node].flux(dimension);
				for (int field = 0; field < fieldCount; field++) {
					volumeFlux[field][node] = tmp[field];
				}
			}
			for (int field = 0; field < fieldCount; field++) {
				auto const tmp = vMassInverse<interiorVolume>(vAnalyze<interiorVolume>(volumeFlux[field]));
				auto const source = vMassInverse<interiorVolume>(vStiffness<interiorVolume>(dimension, tmp));
				for (int mode = 0; mode < modeVolume; mode++) {
					stateDerivative[field][mode][interiorSlice] += lambda * source[mode];
				}
			}
		}
	}
}
;

