#pragma once
#include "ContainerArithmetic.hpp"
#include "Hdf5.hpp"
#include "EulerState.hpp"
#include "Matrix.hpp"
#include "MultiIndex.hpp"
#include "Quadrature.hpp"
#include "dgTransforms.hpp"
#include "Valarray.hpp"

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
	static constexpr auto vAnalyze = dgAnalyze<Valarray<Type>, dimCount, modeCount>;
	static constexpr auto vMassInverse = dgMassInverse<Valarray<Type>, dimCount, modeCount>;
	static constexpr auto vSynthesize = dgSynthesize<Valarray<Type>, dimCount, modeCount>;
	static constexpr auto vStiffness = dgStiffness<Valarray<Type>, dimCount, modeCount>;
	static constexpr auto sAnalyze = dgAnalyze<Valarray<Type>, dimCount - 1, modeCount>;
	static constexpr auto sSynthesize = dgSynthesize<Valarray<Type>, dimCount - 1, modeCount>;
	static constexpr auto sTrace = dgTrace<Valarray<Type>, dimCount, modeCount>;
	static constexpr auto sTraceInverse = dgTraceInverse<Valarray<Type>, dimCount, modeCount>;

//	using IntArray = Valarray<Type, intVolume>;
//	using Pad1Array = Valarray<Type, pad1Volume>;
//	using Pad2Array = Valarray<Type, pad2Volume>;
//	using BndArray = Valarray<Type, bndVolume>;
//	using IntGrid = State<std::array<IntArray, modeVolume>, dimCount>;
//	using BndGrid = State<std::array<BndArray, modeVolume>, dimCount>;
//	using Pad1Grid = State<std::array<Pad1Array, modeVolume>, dimCount>;
//	using Pad2Grid = State<std::array<Pad2Array, modeVolume>, dimCount>;
	using Grid = State<std::array<Valarray<Type>, modeVolume>, dimCount>;
	using BndHandle = std::function<hpx::future<Grid>()>;
	using BaseType = typename ElementType<Type>::type;

	Grid createGrid(size_t count) const {
		Grid grid;
		for (int field = 0; field < fieldCount; field++) {
			for (int mode = 0; mode < modeVolume; mode++) {
				grid[field][mode] = Valarray<Type>(count);
			}
		}
		return grid;
	}

	Type cellWidth;
	Grid np1State;
	Grid nState;
	std::array<Grid, rkStageCount> jState;
	Valarray<size_t> intSizes;
	Valarray<size_t> intStrides;
	std::array<Valarray<Type>, dimCount> position;
	std::array<BndHandle, faceCount(dimCount)> bndHandles;
	std::array<Valarray<size_t>, 2 * dimCount> pad1Strides;
	std::array<Valarray<size_t>, 2 * dimCount> pad1Sizes;
	std::array<GSlice, 2 * dimCount> pad1Slices;
	std::array<Valarray<size_t>, dimCount> pad2Strides;
	std::array<Valarray<size_t>, dimCount> pad2Sizes;
	std::array<GSlice, dimCount> pad2Slices;

	static constexpr BaseType zero = BaseType(0);
	static constexpr BaseType one = BaseType(1);
	static constexpr BaseType two = BaseType(2);
	static constexpr BaseType half = one / two;

public:
	HyperGrid(Type const &xNint = Type(1)) {
		cellWidth = (xNint / Type(intWidth));
		np1State = createGrid(intVolume);
		nState = createGrid(intVolume);
		jState = makeFilledArray<Grid, rkStageCount>(createGrid(intVolume));
		intSizes = Valarray<size_t>(dimCount);
		intStrides = Valarray<size_t>(dimCount);
		position = makeFilledArray<Valarray<Type>, dimCount>(Valarray<Type>(intVolume));
		pad1Strides = makeFilledArray<Valarray<size_t>, 2 * dimCount>(Valarray<size_t>(dimCount));
		pad1Sizes = makeFilledArray<Valarray<size_t>, 2 * dimCount>(Valarray<size_t>(dimCount));
		pad2Strides = makeFilledArray<Valarray<size_t>, dimCount>(Valarray<size_t>(dimCount));
		pad2Sizes = makeFilledArray<Valarray<size_t>, dimCount>(Valarray<size_t>(dimCount));
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
				GSlice slice(start, sizes, intStrides);
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
			pad2Slices[dim] = GSlice(pad2Strides[dim][dim], intSizes, pad2Strides[dim]);
		}
		for (Face face = 0; face < faceCount(dimCount); face++) {
			pad1Sizes[face] = intSizes;
			pad1Sizes[face][getFaceDim(face)]++;
			pad1Strides[face][dimCount - 1] = 1;
			for (int n = dimCount - 1; n > 0; n--) {
				pad1Strides[face][n - 1] = pad1Strides[face][n] * pad1Sizes[face][n];
			}
			size_t const thisStride = (getFaceDir(face) < 0) ? pad1Strides[face][getFaceDim(face)] : size_t(0);
			pad1Slices[face] = GSlice(thisStride, intSizes, pad1Strides[face]);
		}
	}
	void initialize(std::function<State<Type, dimCount>(std::array<Type, dimCount> const&)> const &initialState) {
		Type const halfCellWidth = Type(0.5) * cellWidth;
		State<std::array<Valarray<Type>, nodeVolume>, dimCount> nodalValues = makeFilledArray<std::array<Valarray<Type>, nodeVolume>, fieldCount>(
				makeFilledArray<Valarray<Type>, nodeVolume>(Valarray<Type>(intVolume)));
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
			np1State[field] = vMassInverse(vAnalyze(nodalValues[field]));
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
				size_t const start = (getFaceDir(face) < 0) ? size_t(0) : ((intSizes[dim] - 1) * intStrides[dim]);
				GSlice slice(start, sizes, intStrides);
				Grid bndState = createGrid(bndVolume);
				for (int field = 0; field < fieldCount; field++) {
					for (int mode = 0; mode < modeVolume; mode++) {
						bndState[field][mode] = np1State[field][mode][slice];
					}
				}
				return bndState;
			});
		});
	}
	Grid createPaddedGrid1(Face face) {
		Grid gState = createGrid(pad1Volume);
		hpx::future<Grid> bndFuture = bndHandles[face]();
		size_t const dim = getFaceDim(face);
		size_t const bndStart = (getFaceDir(face) > 0) ? (intSizes[dim] * pad1Strides[face][dim]) : size_t(0);
		size_t const intStart = (getFaceDir(face) < 0) ? pad1Strides[face][dim] : size_t(0);
		auto bndSizes = intSizes;
		bndSizes[dim] = 1;
		auto const intSlice = GSlice(intStart, intSizes, pad1Strides[face]);
		auto const bndSlice = GSlice(bndStart, bndSizes, pad1Strides[face]);
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
	Grid createPaddedGrid2(int dim) {
		Grid gState = createGrid(pad2Volume);
		;
		hpx::future<Grid> bndFuture1 = bndHandles[makeFace(dim, -1)]();
		hpx::future<Grid> bndFuture2 = bndHandles[makeFace(dim, +1)]();
		auto bndSizes = intSizes;
		bndSizes[dim] = 1;
		auto const intSlice = GSlice(pad2Strides[dim][dim], intSizes, pad2Strides[dim]);
		auto const bndSlice1 = GSlice(0, bndSizes, pad2Strides[dim]);
		auto const bndSlice2 = GSlice((intSizes[dim] + 1) * pad2Strides[dim][dim], bndSizes, pad2Strides[dim]);
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
	Type beginStep() {
		using std::max;
		std::array<Type, dimCount> maximumEigenvalue;
		std::array<State<Valarray<Type>, dimCount>, nodeVolume> stateAtNodes;
		nState = np1State;
		maximumEigenvalue.fill(Type(0));
		for (int field = 0; field < fieldCount; field++) {
			auto const tmp = vSynthesize(np1State[field]);
			for (int node = 0; node < nodeVolume; node++) {
				stateAtNodes[node][field] = tmp[node];
			}
		}
		for (int node = 0; node < nodeVolume; node++) {
			for (int dim = 0; dim < dimCount; dim++) {
				auto const eigenvalues = stateAtNodes[node].eigenvalues(dim);
				for (int i = 0; i < fieldCount; i++) {
					Valarray<Type> const lambda = abs(eigenvalues[i]);
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
			if (modeVolume == 1) {
				return;
			}
			State<std::array<std::array<Valarray<Type>, dimCount>, modeVolume>, dimCount> alpha = makeFilledArray<
					std::array<std::array<Valarray<Type>, dimCount>, modeVolume>, fieldCount>(
					makeFilledArray<std::array<Valarray<Type>, dimCount>, modeVolume>(makeFilledArray<Valarray<Type>, dimCount>(Valarray<Type>(intVolume))));
			std::array<SquareMatrix<Valarray<Type>, fieldCount>, dimCount> leftEigenvectors =
					makeFilledArray<SquareMatrix<Valarray<Type>, fieldCount>, dimCount>(Valarray<Type>(intVolume));
			std::array<SquareMatrix<Valarray<Type>, fieldCount>, dimCount> rightEigenvectors = makeFilledArray<SquareMatrix<Valarray<Type>, fieldCount>,
					dimCount>(Valarray<Type>(intVolume));
			State<Valarray<Type>, dimCount> meanState = makeFilledArray<Valarray<Type>, fieldCount>(Valarray<Type>(intVolume));
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
			State<Valarray<Type>, dimCount> loState = makeFilledArray<Valarray<Type>, fieldCount>(Valarray<Type>(pad2Volume));
			State<Valarray<Type>, dimCount> hiState = makeFilledArray<Valarray<Type>, fieldCount>(Valarray<Type>(intVolume));
			State<Valarray<Type>, dimCount> slopeState = makeFilledArray<Valarray<Type>, fieldCount>(Valarray<Type>(intVolume));
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
						GSlice const intSlice(stride, intSizes, pad2Strides[dim]);
						GSlice const intPlusSlice(2 * stride, intSizes, pad2Strides[dim]);
						GSlice const intMinusSlice(0, intSizes, pad2Strides[dim]);
						for (int field = 0; field < fieldCount; field++) {
							hiState[field] = gState[field][hiMode][intSlice];
							loState[field] = gState[field][loMode];
						}
						for (int field = 0; field < fieldCount; field++) {
							Valarray<Type> const stateCentral = loState[field][intSlice];
							Valarray<Type> const statePlus = loState[field][intPlusSlice];
							Valarray<Type> const stateMinus = loState[field][intMinusSlice];
							slopeState[field] = minmod(Valarray<Type>(statePlus - stateCentral), Valarray<Type>(stateCentral - stateMinus));
						}
						State<Valarray<Type>, dimCount> initialStateInverse = makeFilledArray<Valarray<Type>, fieldCount>(Valarray<Type>(intVolume));
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
							Valarray<Type> alphaHi(intVolume), alphaLo(intVolume);
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
					Valarray<Type> thisAlpha(BaseType(0), intVolume);
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
			Valarray<Type> theta(Type(1), intVolume);
			std::array<State<Valarray<Type>, dimCount>, nodeSurface> surfaceState = makeFilledArray<State<Valarray<Type>, dimCount>, nodeSurface>(
					makeFilledArray<Valarray<Type>, fieldCount>(Valarray<Type>(intVolume)));
			std::array<State<Valarray<Type>, dimCount>, nodeVolume> volumeState = makeFilledArray<State<Valarray<Type>, dimCount>, nodeVolume>(
					makeFilledArray<Valarray<Type>, fieldCount>(Valarray<Type>(intVolume)));
			State<Valarray<Type>, dimCount> meanState(makeFilledArray<Valarray<Type>, fieldCount>(Valarray<Type>(intVolume)));
			std::array<Valarray<Type>, modeVolume> modes;
			for (int field = 0; field < fieldCount; field++) {
				meanState[field] = np1State[field][0];
			}
			for (Face face = 0; face < faceCount(dimCount); face++) {
				for (int field = 0; field < fieldCount; field++) {
					for (int mode = 0; mode < modeVolume; mode++) {
						modes[mode] = np1State[field][mode];
					}
					auto const surfaceNodes = sSynthesize(sTrace(face, modes));
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
				auto const volumeNodes = vSynthesize(modes);
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
	Grid computeTimeDerivative(Type timeStepSize, Grid &kState) {
		Type const lambda = Type(2) * timeStepSize / cellWidth;
		{
			for (int dim = 0; dim < dimCount; dim++) {
				std::array<Valarray<Type>, modeSurface> pModalFaceFlux;
				std::array<Valarray<Type>, modeSurface> mModalFaceFlux;
				std::array<Valarray<Type>, modeVolume> pModalVolumeFlux;
				std::array<Valarray<Type>, modeVolume> mModalVolumeFlux;
				State<std::array<Valarray<Type>, nodeSurface>, dimCount> nodalFlux = makeFilledArray<std::array<Valarray<Type>, nodeSurface>, fieldCount>(
						makeFilledArray<Valarray<Type>, nodeSurface>(Valarray<Type>(pad1Volume)));
				std::array<State<Valarray<Type>, dimCount>, nodeSurface> lState = makeFilledArray<State<Valarray<Type>, dimCount>, nodeSurface>(
						makeFilledArray<Valarray<Type>, fieldCount>(Valarray<Type>(pad1Volume)));
				std::array<State<Valarray<Type>, dimCount>, nodeSurface> rState = makeFilledArray<State<Valarray<Type>, dimCount>, nodeSurface>(
						makeFilledArray<Valarray<Type>, fieldCount>(Valarray<Type>(pad1Volume)));
				auto lModes = createPaddedGrid1(makeFace(dim, -1));
				auto rModes = createPaddedGrid1(makeFace(dim, +1));
				for (int field = 0; field < fieldCount; field++) {
					auto const lNodes = sTrace(makeFace(dim, +1), lModes[field]);
					auto const rNodes = sTrace(makeFace(dim, -1), rModes[field]);
					auto const left = sSynthesize(lNodes);
					auto const right = sSynthesize(rNodes);
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
				GSlice const fluxPlusSlice(pad1Strides[2 * dim + 1][dim], intSizes, pad1Strides[2 * dim + 1]);
				GSlice const fluxMinusSlice(0, intSizes, pad1Strides[2 * dim + 0]);
				for (int field = 0; field < fieldCount; field++) {
					auto const modalFlux = sAnalyze(nodalFlux[field]);
					for (int mode = 0; mode < modeSurface; mode++) {
						pModalFaceFlux[mode] = modalFlux[mode][fluxPlusSlice];
						mModalFaceFlux[mode] = modalFlux[mode][fluxMinusSlice];
					}
					pModalVolumeFlux = sTraceInverse(2 * dim + 1, pModalFaceFlux);
					mModalVolumeFlux = sTraceInverse(2 * dim + 0, mModalFaceFlux);
					pModalVolumeFlux = vMassInverse(pModalVolumeFlux);
					mModalVolumeFlux = vMassInverse(mModalVolumeFlux);
					for (int mode = 0; mode < modeVolume; mode++) {
						kState[field][mode] -= lambda * (pModalVolumeFlux[mode] - mModalVolumeFlux[mode]);
					}
//					if(field==2){
//					for(int i = 0; i <= intWidth; i++) {
//						printf( "%i %e %e\n", i, mModalVolumeFlux[0][64*i]/lambda, pModalVolumeFlux[0][64*i]/lambda);
//					}}
				}
			}
		}
		{
			for (int dim = 0; dim < dimCount; dim++) {
				std::array<State<Valarray<Type>, dimCount>, nodeVolume> volumeState;
				State<std::array<Valarray<Type>, nodeVolume>, dimCount> volumeFlux;
				for (int field = 0; field < fieldCount; field++) {
					std::array<Valarray<Type>, modeVolume> volumeModes;
					for (int mode = 0; mode < modeVolume; mode++) {
						volumeModes[mode] = np1State[field][mode];
					}
					auto const nodes = vSynthesize(volumeModes);
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
					auto const tmp = vMassInverse(vAnalyze(volumeFlux[field]));
					auto const source = vMassInverse(vStiffness(dim, tmp));
					for (int mode = 0; mode < modeVolume; mode++) {
					//	kState[field][mode] += lambda * source[mode];
					}
				}
			}
		}
		return kState;
	}
}
;

