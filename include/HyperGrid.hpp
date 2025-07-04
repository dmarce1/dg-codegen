#pragma once
#include "ContainerArithmetic.hpp"
#include "Hdf5.hpp"
#include "EulerState.hpp"
#include "Matrix.hpp"
#include "MultiIndex.hpp"
#include "Quadrature.hpp"
#include "dgTransforms.hpp"
#include "ValArray.hpp"

#include <hpx/future.hpp>

#include <bitset>
#include <functional>
#include <memory>
#include <numeric>

template<typename Type>
inline constexpr Type minmod(Type const &a, Type const &b) {
	using EleType = typename ElementType<Type>::type;
	constexpr EleType half = EleType(1) / EleType(2);
	Type const sgn = copysign(half, a) + copysign(half, b);
	Type const mag = min(abs(a), abs(b));
	return sgn * mag;
}

using Face = int;

constexpr Face makeFace(int dim, int dir) {
	return Face(2 * dim + ((dir + 1) >> 1));
}

constexpr Face begFace() {
	return Face(0);
}

constexpr int faceCount(int dimCount) {
	return 2 * dimCount;
}

constexpr Face flipFace(Face face) {
	face ^= 1;
	return face;
}

constexpr int getFaceDim(Face face) {
	return face >> 1;
}

constexpr int getFaceDir(Face face) {
	return ((face & 1) << 1) - 1;
}

template<typename Type, int dimCount, int intWidth, int modeCount, typename RungeKutta, template<typename, int> typename State>
class HyperGrid {
	static constexpr int extWidth = intWidth + 2;
	static constexpr int fieldCount = State<Type, dimCount>::fieldCount();
	static constexpr int rkStageCount = RungeKutta::stageCount();
	static constexpr int modeVolume = binco(modeCount + dimCount - 1, dimCount);
	static constexpr int modeSurface = binco(modeCount + dimCount - 2, dimCount - 1);
	static constexpr int nodeVolume = ipow(modeCount, dimCount);
	static constexpr int nodeSurface = ipow(modeCount, dimCount - 1);
	static constexpr int intVolume = ipow(intWidth, dimCount);
	static constexpr int bndVolume = ipow(intWidth, dimCount - 1);
	static constexpr int pad1Volume = intVolume + bndVolume;
	static constexpr int pad2Volume = intVolume + 2 * bndVolume;
	static constexpr RungeKutta butcherTable { };
	template<int volume = intVolume>
	static constexpr auto& vAnalyze = dgAnalyze<ValArray<Type, volume>, dimCount, modeCount>;
	template<int volume = intVolume>
	static constexpr auto& vMassInverse = dgMassInverse<ValArray<Type, volume>, dimCount, modeCount>;
	template<int volume = intVolume>
	static constexpr auto& vSynthesize = dgSynthesize<ValArray<Type, volume>, dimCount, modeCount>;
	template<int volume = intVolume>
	static constexpr auto& vStiffness = dgStiffness<ValArray<Type, volume>, dimCount, modeCount>;
	template<int volume = intVolume>
	static constexpr auto& sAnalyze = dgAnalyze<ValArray<Type, volume>, dimCount - 1, modeCount>;
	template<int volume = intVolume>
	static constexpr auto& sSynthesize = dgSynthesize<ValArray<Type, volume>, dimCount - 1, modeCount>;
	template<int volume = intVolume>
	static constexpr auto& sTrace = dgTrace<ValArray<Type, volume>, dimCount, modeCount>;
	template<int volume = intVolume>
	static constexpr auto& sTraceInverse = dgTraceInverse<ValArray<Type, volume>, dimCount, modeCount>;

	using IntArray = ValArray<Type, intVolume>;
	using Pad1Array = ValArray<Type, pad1Volume>;
	using Pad2Array = ValArray<Type, pad2Volume>;
	using BndArray = ValArray<Type, bndVolume>;
	using IntGrid = State<std::array<IntArray, modeVolume>, dimCount>;
	using BndGrid = State<std::array<BndArray, modeVolume>, dimCount>;
	using Pad1Grid = State<std::array<Pad1Array, modeVolume>, dimCount>;
	using Pad2Grid = State<std::array<Pad2Array, modeVolume>, dimCount>;
	using BndHandle = std::function<hpx::future<BndGrid>()>;
	using BaseType = typename ElementType<Type>::type;

	Type cellWidth;
	IntGrid np1State;
	IntGrid nState;
	std::array<IntGrid, rkStageCount> jState;
	ValArray<int64_t, dimCount> intSizes;
	ValArray<int64_t, dimCount> intStrides;
	std::array<ValArray<Type, intVolume>, dimCount> position;
	std::array<BndHandle, faceCount(dimCount)> bndHandles;
	std::array<ValArray<int64_t, dimCount>, 2 * dimCount> pad1Strides;
	std::array<ValArray<int64_t, dimCount>, 2 * dimCount> pad1Sizes;
	std::array<GSlice<dimCount, intVolume>, 2 * dimCount> pad1Slices;
	std::array<ValArray<int64_t, dimCount>, dimCount> pad2Strides;
	std::array<ValArray<int64_t, dimCount>, dimCount> pad2Sizes;
	std::array<GSlice<dimCount, intVolume>, dimCount> pad2Slices;

	static constexpr BaseType zero = BaseType(0);
	static constexpr BaseType one = BaseType(1);
	static constexpr BaseType two = BaseType(2);
	static constexpr BaseType half = one / two;
public:
	HyperGrid(Type const &xNint = Type(1)) :
			cellWidth { xNint / Type(intWidth) } {
		intSizes.fill(intWidth);
		intStrides[dimCount - 1] = 1;
		for (int n = dimCount - 1; n > 0; n--) {
			intStrides[n - 1] = intStrides[n] * intSizes[n];
		}
		for (int dim = 0; dim < dimCount; dim++) {
			for (int i = 0; i < intWidth; i++) {
				auto const start = i * intStrides[dim];
				auto sizes = intSizes;
				sizes[dim] = 1;
				BaseType const x = (BaseType(i) + half) * half * cellWidth;
				GSlice<dimCount, bndVolume> slice(start, sizes, intStrides);
				position[dim][slice] = x;
			}
		}
		for (Face face = 0; face < faceCount(dimCount); face++) {
			bndHandles[face] = this->getBoundaryHandle(flipFace(face));
		}
		for (int dim = 0; dim < dimCount; dim++) {
			pad2Sizes[dim] = intSizes;
			pad2Sizes[dim][dim] += 2;
			pad2Strides[dim][dimCount - 1] = 1;
			for (int n = dimCount - 1; n > 0; n--) {
				pad2Strides[dim][n - 1] = pad2Strides[dim][n] * pad2Sizes[dim][n];
			}
			pad2Slices[dim] = GSlice<dimCount, intVolume>(pad2Strides[dim][dim], intSizes, pad2Strides[dim]);
		}
		for (Face face = 0; face < faceCount(dimCount); face++) {
			pad1Sizes[face] = intSizes;
			pad1Sizes[face][getFaceDim(face)]++;
			pad1Strides[face][dimCount - 1] = 1;
			for (int n = dimCount - 1; n > 0; n--) {
				pad1Strides[face][n - 1] = pad1Strides[face][n] * pad1Sizes[face][n];
			}
			pad1Slices[face] = GSlice<dimCount, intVolume>(pad1Strides[face][getFaceDim(face)], intSizes, pad1Strides[face]);
		}
	}
	void initialize(std::function<State<Type, dimCount>(std::array<Type, dimCount> const&)> const &initialState) {
		Type const halfCellWidth = Type(0.5) * cellWidth;
		State<std::array<ValArray<Type, intVolume>, nodeVolume>, dimCount> nodalValues;
		for (int node = 0; node < nodeVolume; node++) {
			auto quadraturePoint = getQuadraturePoint<Type, dimCount, modeCount>(node);
			auto thisPosition = position;
			for (int dim = 0; dim < dimCount; dim++) {
				thisPosition[dim] += quadraturePoint[dim] * halfCellWidth;
			}
			for (int i = 0; i < intVolume; i++) {
				std::array<Type, dimCount> x;
				for (int dim = 0; dim < dimCount; dim++) {
					x[dim] = thisPosition[dim][i];
				}
				auto const thisState = initialState(x);
				for (int field = 0; field < fieldCount; field++) {
					nodalValues[field][node][i] = thisState[field];
				}
			}
		}
		for (int field = 0; field < fieldCount; field++) {
			np1State[field] = vMassInverse<>(vAnalyze<>(nodalValues[field]));
		}
	}
	void output(const char *filenameBase, int timeStepNumber, Type const &time) {
		std::string filename = std::string(filenameBase) + "." + std::to_string(timeStepNumber) + ".h5";
		writeHdf5<Type, dimCount, intWidth, modeCount, State>(filename, cellWidth, np1State, State<Type, dimCount>::getFieldNames());
		writeList("X.visit", "!NBLOCKS 1\n", filename + ".xmf");
	}
	BndHandle getBoundaryHandle(Face face) const {
		return BndHandle([this, face]() {
			return hpx::async([this, face]() {
				auto const dim = getFaceDim(face);
				auto sizes = intSizes;
				sizes[dim] = 1;
				int64_t const start = (getFaceDir(face) < 0) ? int64_t(0) : ((intSizes[dim] - 1) * intStrides[dim]);
				GSlice<dimCount, bndVolume> slice(start, sizes, intStrides);
				BndGrid bndState;
				for (int field = 0; field < fieldCount; field++) {
					for (int mode = 0; mode < modeVolume; mode++) {
						bndState[field][mode] = np1State[field][mode][slice];
					}
				}
				return bndState;
			});
		});
	}
	Pad1Grid createPaddedGrid1(Face face) {
		Pad1Grid gState;
		hpx::future<BndGrid> bndFuture = bndHandles[face]();
		int64_t const dim = getFaceDim(face);
		int64_t const bndStart = (getFaceDir(face) > 0) ? (intSizes[dim] * pad1Strides[face][dim]) : int64_t(0);
		int64_t const intStart = (getFaceDir(face) < 0) ? pad1Strides[face][dim] : int64_t(0);
		auto bndSizes = intSizes;
		bndSizes[dim] = 1;
		auto const intSlice = GSlice<dimCount, intVolume>(intStart, intSizes, pad1Strides[face]);
		auto const bndSlice = GSlice<dimCount, bndVolume>(bndStart, bndSizes, pad1Strides[face]);
		for (int field = 0; field < fieldCount; field++) {
			for (int mode = 0; mode < modeVolume; mode++) {
				gState[field][mode][intSlice] = np1State[field][mode];
			}
		}
		auto const bndState = bndFuture.get();
		for (int field = 0; field < fieldCount; field++) {
			for (int mode = 0; mode < modeVolume; mode++) {
				gState[field][mode][bndSlice] = bndState[field][mode];
			}
		}
		return gState;
	}
	Pad2Grid createPaddedGrid2(int dim) {
		Pad2Grid gState;
		hpx::future<BndGrid> bndFuture1 = bndHandles[makeFace(dim, -1)]();
		hpx::future<BndGrid> bndFuture2 = bndHandles[makeFace(dim, +1)]();
		auto bndSizes = intSizes;
		bndSizes[dim] = 1;
		auto const intSlice = GSlice<dimCount, intVolume>(pad2Strides[dim][dim], intSizes, pad2Strides[dim]);
		auto const bndSlice1 = GSlice<dimCount, bndVolume>(0, bndSizes, pad2Strides[dim]);
		auto const bndSlice2 = GSlice<dimCount, bndVolume>((intSizes[dim] + 1) * pad2Strides[dim][dim], bndSizes, pad2Strides[dim]);
		for (int field = 0; field < fieldCount; field++) {
			for (int mode = 0; mode < modeVolume; mode++) {
				gState[field][mode][intSlice] = np1State[field][mode];
			}
		}
		auto const bndState1 = bndFuture1.get();
		auto const bndState2 = bndFuture2.get();
		for (int field = 0; field < fieldCount; field++) {
			for (int mode = 0; mode < modeVolume; mode++) {
				gState[field][mode][bndSlice1] = bndState1[field][mode];
				gState[field][mode][bndSlice2] = bndState2[field][mode];
			}
		}
		return gState;
	}
	ValArray<int64_t, dimCount> flxSizes(int dim) {
		ValArray<int64_t, dimCount> sizes;
		for (int n = 0; n < dimCount; n++) {
			sizes[n] = (intWidth + int64_t(n == dim));
		}
		return sizes;
	}
	ValArray<int64_t, dimCount> flxStrides(int dim) {
		ValArray<int64_t, dimCount> strides;
		strides[dimCount - 1] = 1;
		for (int n = dimCount - 1; n > 0; n--) {
			strides[n - 1] = strides[n] * (intWidth + int64_t(n == dim));
		}
		return strides;
	}
	Type beginStep() {
		using std::max;
		std::array<Type, dimCount> maximumEigenvalue;
		std::array<State<ValArray<Type, intVolume>, dimCount>, nodeVolume> stateAtNodes;
		nState = np1State;
		maximumEigenvalue.fill(Type(0));
		for (int field = 0; field < fieldCount; field++) {
			auto const tmp = vSynthesize<>(np1State[field]);
			for (int node = 0; node < nodeVolume; node++) {
				stateAtNodes[node][field] = tmp[node];
			}
		}
		for (int node = 0; node < nodeVolume; node++) {
			for (int dim = 0; dim < dimCount; dim++) {
				auto const eigenvalues = stateAtNodes[node].eigenvalues(dim);
				for (int i = 0; i < fieldCount; i++) {
					ValArray<Type, intVolume> const lambda = abs(eigenvalues[i]);
					maximumEigenvalue[dim] = max(maximumEigenvalue[dim], lambda.max());
				}
			}
		}
		Type maximumEigenvalueSum = Type(0);
		for (int dim = 0; dim < dimCount; dim++) {
			maximumEigenvalueSum += maximumEigenvalue[dim];
		}
		Type const timeStepSize = (cellWidth * butcherTable.cfl()) / (Type(2 * modeCount - 1) * maximumEigenvalueSum);
		return timeStepSize;
	}
	void subStep(Type const &timeStepSize, int stageIndex) {
		for (int field = 0; field < fieldCount; field++) {
			for (int mode = 0; mode < modeVolume; mode++) {
				np1State[field][mode] = nState[field][mode];
				for (int thisStage = 0; thisStage < stageIndex; thisStage++) {
					np1State[field][mode] += butcherTable.a(stageIndex, thisStage) * jState[thisStage][field][mode];
				}
				jState[stageIndex][field][mode].fill(Type(0));
			}
		}
		applyLimiter();
    	computeTimeDerivative(timeStepSize, jState[stageIndex]);
	}
	void endStep() {
		for (int field = 0; field < fieldCount; field++) {
			for (int mode = 0; mode < modeVolume; mode++) {
				np1State[field][mode] = nState[field][mode];
				for (int j = 0; j < rkStageCount; j++) {
					np1State[field][mode] += butcherTable.b(j) * jState[j][field][mode];
				}
			}
		}
		applyLimiter();
	}
	void applyLimiter() {
		{
			std::array<std::array<State<IntArray, dimCount>, modeVolume>, dimCount> alpha;
			std::array<SquareMatrix<IntArray, fieldCount>, dimCount> leftEigenvectors;
			std::array<SquareMatrix<IntArray, fieldCount>, dimCount> rightEigenvectors;
			State<IntArray, dimCount> meanState;
			std::array<std::bitset<modeCount>, dimCount> wasLimited;
			for (int field = 0; field < fieldCount; field++) {
				meanState[field] = np1State[field][0];
			}
			for (int dim = 0; dim < dimCount; dim++) {
				rightEigenvectors[dim] = meanState.eigenSystem(dim).second;
				leftEigenvectors[dim] = matrixInverse(rightEigenvectors[dim]);
			}
			for (int dim = 0; dim < dimCount; dim++) {
				for (int mode = 0; mode < modeVolume; mode++) {
					for (int field = 0; field < fieldCount; field++) {
						alpha[dim][mode][field] = Type(mode == 0);
					}
				}
			}
			State<Pad2Array, dimCount> loState;
			State<IntArray, dimCount> hiState;
			State<IntArray, dimCount> slopeState;
			for (int polynomialDegree = modeCount - 1; polynomialDegree > 0; polynomialDegree--) {
				constexpr auto tiny = BaseType(std::numeric_limits<double>::min());
				int const begin = binco(polynomialDegree + dimCount - 1, dimCount);
				int const end = binco(polynomialDegree + dimCount, dimCount);
				for (int hiMode = begin; hiMode < end; hiMode++) {
					auto hiModeIndices = flatToTriangular<dimCount, modeCount>(hiMode);
					for (int dim = 0; dim < dimCount; dim++) {
						if (hiModeIndices[dim] == 0) {
							continue;
						}
						auto loModeIndices = hiModeIndices;
						loModeIndices[dim]--;
						auto const loMode = triangularToFlat<dimCount, modeCount>(loModeIndices);
						auto gState = createPaddedGrid2(dim);
						int const stride = pad2Strides[dim][dim];
						GSlice<dimCount, intVolume> const intSlice(stride, intSizes, pad2Strides[dim]);
						GSlice<dimCount, intVolume> const intPlusSlice(2 * stride, intSizes, pad2Strides[dim]);
						GSlice<dimCount, intVolume> const intMinusSlice(0, intSizes, pad2Strides[dim]);
						for (int field = 0; field < fieldCount; field++) {
							hiState[field] = gState[field][hiMode][intSlice];
							loState[field] = gState[field][loMode];
						}
						for (int field = 0; field < fieldCount; field++) {
							IntArray const stateCentral = loState[field][intSlice];
							IntArray const statePlus = loState[field][intPlusSlice];
							IntArray const stateMinus = loState[field][intMinusSlice];
							slopeState[field] = minmod(IntArray(statePlus - stateCentral), IntArray(stateCentral - stateMinus));
						}
						State<IntArray, dimCount> initialStateInverse;
						for (int field = 0; field < fieldCount; field++) {
							initialStateInverse[field] = BaseType(1) / (hiState[field] + tiny);
						}
						slopeState = leftEigenvectors[dim] * slopeState;
						hiState = leftEigenvectors[dim] * hiState;
						auto const cLimit = Type(1) / Type(2 * loModeIndices[dim] + 1);
						for (int field = 0; field < fieldCount; field++) {
							slopeState[field] *= cLimit;
							hiState[field] = minmod(hiState[field], slopeState[field]);
						}
						hiState = rightEigenvectors[dim] * hiState;
						for (int field = 0; field < fieldCount; field++) {
							IntArray alphaHi, alphaLo;
							alphaHi = alpha[dim][hiMode][field];
							alphaLo = alpha[dim][loMode][field];
							alphaHi = max(alphaHi, max(BaseType(0.0), hiState[field] * initialStateInverse[field]));
							alpha[dim][hiMode][field] = alphaHi;
							alpha[dim][loMode][field] = alphaLo;
						}
						wasLimited[dim][hiMode] = true;
					}
				}
			}
			for (int mode = 1; mode < modeVolume; mode++) {
				for (int field = 0; field < fieldCount; field++) {
					IntArray thisAlpha = BaseType(0);
					int count = 0;
					for (int dim = 0; dim < dimCount; dim++) {
						if (wasLimited[dim][mode]) {
							thisAlpha += alpha[dim][mode][field];
							count++;
						}
					}
					np1State[field][mode] *= (thisAlpha / BaseType(count));
				}
			}
		}
		{
			ValArray<Type, intVolume> theta(Type(1));
			std::array<State<IntArray, dimCount>, nodeSurface> surfaceState;
			std::array<State<IntArray, dimCount>, nodeVolume> volumeState;
			State<IntArray, dimCount> meanState;
			std::array<IntArray, modeVolume> modes;
			for (int field = 0; field < fieldCount; field++) {
				meanState[field] = np1State[field][0];
			}
			for (Face face = 0; face < faceCount(dimCount); face++) {
				for (int field = 0; field < fieldCount; field++) {
					for (int mode = 0; mode < modeVolume; mode++) {
						modes[mode] = np1State[field][mode];
					}
					auto const surfaceNodes = sSynthesize<>(sTrace<>(face, modes));
					for (int node = 0; node < nodeSurface; node++) {
						surfaceState[node][field] = surfaceNodes[node];
					}
				}
				for (int node = 0; node < nodeSurface; node++) {
					theta = min(theta, findPositivityPreservingTheta(meanState, surfaceState[node]));
				}
			}
			for (int field = 0; field < fieldCount; field++) {
				for (int mode = 0; mode < modeVolume; mode++) {
					modes[mode] = np1State[field][mode];
				}
				auto const volumeNodes = vSynthesize<intVolume>(modes);
				for (int node = 0; node < nodeVolume; node++) {
					volumeState[node][field] = volumeNodes[node];
				}
			}
			for (int node = 0; node < nodeVolume; node++) {
				theta = min(theta, findPositivityPreservingTheta(meanState, volumeState[node]));
			}
			for (int field = 0; field < fieldCount; field++) {
				for (int mode = 1; mode < modeVolume; mode++) {
					np1State[field][mode] *= theta;
				}
			}
		}
	}
	IntGrid computeTimeDerivative(Type timeStepSize, IntGrid& kState) {
		Type const lambda = Type(2) * timeStepSize / cellWidth;
		{
			for (int dim = 0; dim < dimCount; dim++) {
				std::array<ValArray<Type, intVolume>, modeSurface> pModalFaceFlux;
				std::array<ValArray<Type, intVolume>, modeSurface> mModalFaceFlux;
				std::array<ValArray<Type, intVolume>, modeVolume> pModalVolumeFlux;
				std::array<ValArray<Type, intVolume>, modeVolume> mModalVolumeFlux;
				State<std::array<ValArray<Type, pad1Volume>, nodeSurface>, dimCount> nodalFlux;
				State<std::array<ValArray<Type, intVolume>, nodeSurface>, dimCount> lNodalFlux;
				State<std::array<ValArray<Type, intVolume>, nodeSurface>, dimCount> rNodalFlux;
				std::array<State<ValArray<Type, pad1Volume>, dimCount>, nodeSurface> lState;
				std::array<State<ValArray<Type, pad1Volume>, dimCount>, nodeSurface> rState;
				auto lModes = createPaddedGrid1(makeFace(dim, -1));
				auto rModes = createPaddedGrid1(makeFace(dim, +1));
				for (int field = 0; field < fieldCount; field++) {
					auto const lNodes = sTrace<pad1Volume>(makeFace(dim, +1), lModes[field]);
					auto const rNodes = sTrace<pad1Volume>(makeFace(dim, -1), rModes[field]);
					auto const left = sSynthesize<pad1Volume>(lNodes);
					auto const right = sSynthesize<pad1Volume>(rNodes);
					for (int node = 0; node < nodeSurface; node++) {
						lState[node][field] = left[node];
						rState[node][field] = right[node];
					}
				}
				for (int node = 0; node < nodeSurface; node++) {
					auto const tmp = solveRiemannProblem(lState[node], rState[node], dim);
					for (int field = 0; field < fieldCount; field++) {
						nodalFlux[field][node] = tmp[field];
					}
				}
				GSlice<dimCount, intVolume> const fluxPlusSlice(pad1Strides[dim][dim], intSizes, pad1Strides[2 * dim + 1]);
				GSlice<dimCount, intVolume> const fluxMinusSlice(0, intSizes, pad1Strides[2 * dim + 0]);
				for (int field = 0; field < fieldCount; field++) {
					auto const modalFlux = sAnalyze<pad1Volume>(nodalFlux[field]);
					for (int mode = 0; mode < modeSurface; mode++) {
						pModalFaceFlux[mode] = modalFlux[mode][fluxPlusSlice];
						mModalFaceFlux[mode] = modalFlux[mode][fluxMinusSlice];
					}
					pModalVolumeFlux = sTraceInverse<intVolume>(2 * dim + 1, pModalFaceFlux);
					mModalVolumeFlux = sTraceInverse<intVolume>(2 * dim + 0, mModalFaceFlux);
					pModalVolumeFlux = vMassInverse<intVolume>(pModalVolumeFlux);
					mModalVolumeFlux = vMassInverse<intVolume>(mModalVolumeFlux);
					for (int mode = 0; mode < modeVolume; mode++) {
						kState[field][mode] -= lambda * (pModalVolumeFlux[mode] - mModalVolumeFlux[mode]);
					}
				}
			}
		}
		{
			for (int dim = 0; dim < dimCount; dim++) {
				std::array<State<ValArray<Type, intVolume>, dimCount>, nodeVolume> volumeState;
				State<std::array<ValArray<Type, intVolume>, nodeVolume>, dimCount> volumeFlux;
				for (int field = 0; field < fieldCount; field++) {
					std::array<ValArray<Type, intVolume>, modeVolume> volumeModes;
					for (int mode = 0; mode < modeVolume; mode++) {
						volumeModes[mode] = np1State[field][mode];
					}
					auto const nodes = vSynthesize<>(volumeModes);
					for (int node = 0; node < nodeVolume; node++) {
						volumeState[node][field] = nodes[node];
					}
				}
				for (int node = 0; node < nodeVolume; node++) {
					auto tmp = volumeState[node].flux(dim);
					for (int field = 0; field < fieldCount; field++) {
						volumeFlux[field][node] = tmp[field];
					}
				}
				for (int field = 0; field < fieldCount; field++) {
					auto const tmp = vMassInverse<>(vAnalyze<>(volumeFlux[field]));
					auto const source = vMassInverse<>(vStiffness<>(dim, tmp));
					for (int mode = 0; mode < modeVolume; mode++) {
						kState[field][mode] += lambda * source[mode];
					}
				}
			}
		}
		return kState;
	}
}
;

