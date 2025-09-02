#pragma once
#include "ContainerArithmetic.hpp"
#include "Hdf5.hpp"
#include "EulerState.hpp"
#include "Face.hpp"
#include "Matrix.hpp"
#include "MultiIndex.hpp"
#include "dgTransforms.hpp"

#include <hpx/future.hpp>

#include <bitset>
#include <functional>
#include <memory>
#include <mutex>
#include <numeric>
#include <stack>
#include <valarray>

template<typename Type>
inline constexpr Type minmod(Type const &a, Type const &b) {
	using namespace Math;
	using EleType = typename ElementType<Type>::type;
	constexpr EleType half = EleType(1) / EleType(2);
	Type const sgn = copysign(half, a) + copysign(half, b);
	Type const aA = abs(a);
	Type const bB = abs(b);
	Type const mag = min(aA, bB);
	return Type(sgn * mag);
}

template<typename Type, int dimensionCount, int interiorWidth, int modeCount, template<typename > typename State>
struct LimiterWorkspace {
	static constexpr int modeVolume = binco(modeCount + dimensionCount - 1, dimensionCount);
	static constexpr int modeSurface = binco(modeCount + dimensionCount - 2, dimensionCount - 1);
	static constexpr int nodeVolume = ipow(modeCount, dimensionCount);
	static constexpr int nodeSurface = ipow(modeCount, dimensionCount - 1);
	static constexpr int interiorVolume = ipow(interiorWidth, dimensionCount);
	static constexpr int fieldCount = State<Type>::fieldCount();
	static hpx::mutex mutex;
	static std::stack<std::unique_ptr<LimiterWorkspace>> workspaces;
	static std::unique_ptr<LimiterWorkspace> getWorkspace() {
		std::unique_ptr<LimiterWorkspace> pointer;
		std::lock_guard<hpx::mutex> lock(mutex);
		if (workspaces.empty()) {
			pointer = std::make_unique<LimiterWorkspace>();
			workspaces.push(std::move(pointer));
		}
		pointer = std::move(workspaces.top());
		workspaces.pop();
		return pointer;

	}
	static void recycleWorkspace(std::unique_ptr<LimiterWorkspace> &&pointer) {
		std::lock_guard<hpx::mutex> lock(mutex);
		workspaces.push(std::move(pointer));
	}
	LimiterWorkspace() :
			thisAlpha(theta) {
		theta = std::valarray<Type>(interiorVolume);
		for (int field = 0; field < fieldCount; field++) {
			hiState[field] = std::valarray<Type>(interiorVolume);
		}
		for (int dim = 0; dim < dimensionCount; dim++) {
			for (int mode = 0; mode < modeVolume; mode++) {
				for (int field = 0; field < fieldCount; field++) {
					alpha[field][mode][dim] = std::valarray<Type>(interiorVolume);
				}
			}
		}
	}
	LimiterWorkspace(LimiterWorkspace const&) = delete;
	LimiterWorkspace(LimiterWorkspace&&) = default;
	LimiterWorkspace& operator=(LimiterWorkspace const&) = delete;
	LimiterWorkspace& operator=(LimiterWorkspace&&) = default;
	State<std::array<std::array<std::valarray<Type>, dimensionCount>, modeVolume>> alpha;
	std::array<SquareMatrix<std::valarray<Type>, fieldCount>, dimensionCount> leftEigenvectors;
	std::array<SquareMatrix<std::valarray<Type>, fieldCount>, dimensionCount> rightEigenvectors;
	State<std::valarray<Type>> loState;
	State<std::valarray<Type>> hiState;
	State<std::valarray<Type>> meanState;
	State<std::valarray<Type>> slopeState;
	std::array<std::bitset<modeCount>, dimensionCount> wasLimited;
	std::valarray<Type> theta;
	std::valarray<Type> &thisAlpha;
	std::array<std::valarray<Type>, nodeSurface> surfaceNodes;
	std::array<std::valarray<Type>, nodeVolume> volumeNodes;
	std::array<State<std::valarray<Type>>, nodeSurface> surfaceState;
	std::array<State<std::valarray<Type>>, nodeVolume> volumeState;
	std::array<std::valarray<Type>, modeVolume> modes;
	State<std::valarray<Type>> invState;
};

template<typename Type, int dimensionCount, int interiorWidth, int modeCount, template<typename > typename State>
hpx::mutex LimiterWorkspace<Type, dimensionCount, interiorWidth, modeCount, State>::mutex { };

template<typename Type, int dimensionCount, int interiorWidth, int modeCount, template<typename > typename State>
std::stack<std::unique_ptr<LimiterWorkspace<Type, dimensionCount, interiorWidth, modeCount, State>>> LimiterWorkspace<Type, dimensionCount, interiorWidth,
		modeCount, State>::workspaces { };

template<typename Type_, int dimensionCount_, int interiorWidth_, int modeCount_, typename RungeKutta_, template<typename, int> typename State_>
struct HyperSubgrid {
	using Type = Type_;
	using RungeKutta = RungeKutta_;
	template<typename T>
	using State = State_<T, dimensionCount_>;
	static constexpr int dimensionCount = dimensionCount_;
	static constexpr int interiorWidth = interiorWidth_;
	static constexpr int modeCount = modeCount_;
private:
	static constexpr Type zero = Type(0);
	static constexpr Type one = Type(1);
	static constexpr Type two = Type(2);
	static constexpr Type half = one / two;
	static constexpr int faceCount = Face<dimensionCount>::count();
	static constexpr int exteriorWidth = interiorWidth + 2;
	static constexpr int fieldCount = State<Type>::fieldCount();
	static constexpr int rkStageCount = RungeKutta::stageCount();
	static constexpr int modeVolume = binco(modeCount + dimensionCount - 1, dimensionCount);
	static constexpr int modeSurface = binco(modeCount + dimensionCount - 2, dimensionCount - 1);
	static constexpr int nodeVolume = ipow(modeCount, dimensionCount);
	static constexpr int nodeSurface = ipow(modeCount, dimensionCount - 1);
	static constexpr int interiorVolume = ipow(interiorWidth, dimensionCount);
	static constexpr int boundaryVolume = ipow(interiorWidth, dimensionCount - 1);
	static constexpr int singlePaddedVolume = interiorVolume + boundaryVolume;
	static constexpr int doublePaddedVolume = interiorVolume + 2 * boundaryVolume;
	static constexpr RungeKutta butcherTable { };
	static constexpr auto volumeAnalyze = dgAnalyze<std::valarray<Type>, dimensionCount, modeCount>;
	static constexpr auto volumeSynthesize = dgSynthesize<std::valarray<Type>, dimensionCount, modeCount>;
	static constexpr auto volumeAnalyzeDerivative = dgAnalyzeDerivative<std::valarray<Type>, dimensionCount, modeCount>;
	static constexpr auto surfaceAnalyze = dgAnalyze<std::valarray<Type>, dimensionCount - 1, modeCount>;
	static constexpr auto surfaceSynthesize = dgSynthesize<std::valarray<Type>, dimensionCount - 1, modeCount>;
	static constexpr auto surfaceTrace = dgTrace<std::valarray<Type>, dimensionCount, modeCount>;
	static constexpr auto surfaceTraceInverse = dgTraceInverse<std::valarray<Type>, dimensionCount, modeCount>;

	using Subgrid = State<std::array<std::valarray<Type>, modeVolume>>;
	using BndHandle = std::function<hpx::future<Subgrid>()>;
	using LimiterWorkspaceType = LimiterWorkspace<Type, dimensionCount, interiorWidth, modeCount, State>;
	using FaceType = Face<dimensionCount>;

	static std::valarray<size_t> const& interiorSizes() {
		static std::valarray<size_t> const a(interiorWidth, dimensionCount);
		return a;
	}
	static std::valarray<size_t> const& interiorStrides() {
		static auto const a = []() {
			std::valarray<size_t> v(dimensionCount);
			v[dimensionCount - 1] = 1;
			for (int n = dimensionCount - 1; n > 0; n--) {
				v[n - 1] = v[n] * interiorSizes()[n];
			}
			return v;
		}();
		return a;
	}
	static std::valarray<size_t> const& singlePaddedSizes(Face<dimensionCount> face) {
		static auto const a = []() {
			std::array<std::valarray<size_t>, faceCount> v;
			for (FaceType face = FaceType::begin(); face != FaceType::end(); face++) {
				v[face] = interiorSizes();
				v[face][face.dimension()]++;
			}
			return v;
		}();
		return a[face];
	}
	static std::valarray<size_t> const& doublePaddedSizes(int dim) {
		static auto const a = []() {
			std::array<std::valarray<size_t>, dimensionCount> v;
			for (int dim = 0; dim < dimensionCount; dim++) {
				v[dim] = interiorSizes();
				v[dim][dim] += 2;
			}
			return v;
		}();
		return a[dim];
	}
	static std::valarray<size_t> const& singlePaddedStrides(FaceType face) {
		static auto const a = []() {
			std::array<std::valarray<size_t>, faceCount> v;
			for (FaceType face = FaceType::begin(); face != FaceType::end(); face++) {
				v[face] = std::valarray<size_t>(dimensionCount);
				v[face][dimensionCount - 1] = 1;
				for (int dim = dimensionCount - 1; dim > 0; dim--) {
					v[face][dim - 1] = v[face][dim] * singlePaddedSizes(face)[dim];
				}
			}
			return v;
		}();
		return a[face];
	}
	;
	static std::valarray<size_t> const& doublePaddedStrides(int dim) {
		static auto const a = []() {
			std::array<std::valarray<size_t>, dimensionCount> v;
			for (int dim1 = 0; dim1 < dimensionCount; dim1++) {
				v[dim1] = std::valarray<size_t>(dimensionCount);
				v[dim1][dimensionCount - 1] = 1;
				for (int dim2 = dimensionCount - 1; dim2 > 0; dim2--) {
					v[dim1][dim2 - 1] = v[dim1][dim2] * doublePaddedSizes(dim1)[dim2];
				}
			}
			return v;
		}();
		return a[dim];
	}
	static std::gslice const& fluxPlusSlice(int dim) {
		static auto const a = []() {
			std::array<std::gslice, dimensionCount> v;
			for (int dim = 0; dim < dimensionCount; dim++) {
				v[dim] = std::gslice(singlePaddedStrides(FaceType(dim, +1))[dim], interiorSizes(), singlePaddedStrides(FaceType(dim, +1)));
			}
			return v;
		}();
		return a[dim];
	}
	static std::gslice const& fluxMinusSlice(int dim) {
		static auto const a = []() {
			std::array<std::gslice, dimensionCount> v;
			for (int dim = 0; dim < dimensionCount; dim++) {
				v[dim] = std::gslice(0, interiorSizes(), singlePaddedStrides(FaceType(dim, -1)));
			}
			return v;
		}();
		return a[dim];
	}
	static Subgrid createSubgrid(size_t count) {
		Subgrid grid;
		for (int field = 0; field < fieldCount; field++) {
			for (int mode = 0; mode < modeVolume; mode++) {
				grid[field][mode] = std::valarray<Type>(count);
			}
		}
		return grid;
	}

	Type cellWidth_;
	Subgrid np1State_;
	Subgrid nState_;
	std::array<Subgrid, rkStageCount> stageStates_;
	std::array<std::valarray<Type>, dimensionCount> position_;
	std::array<BndHandle, faceCount> boundaryHandles;
public:
	auto localSolve() const {
	}
	HyperSubgrid(Type const &xNint = Type(1)) {
		cellWidth_ = (xNint / Type(interiorWidth));
		np1State_ = createSubgrid(interiorVolume);
		nState_ = createSubgrid(interiorVolume);
		stageStates_ = makeFilledArray<Subgrid, rkStageCount>(createSubgrid(interiorVolume));
		position_ = makeFilledArray<std::valarray<Type>, dimensionCount>(std::valarray<Type>(interiorVolume));
		for (int dim = 0; dim < dimensionCount; dim++) {
			for (int i = 0; i < interiorWidth; i++) {
				auto const start = i * interiorStrides()[dim];
				auto sizes = interiorSizes();
				sizes[dim] = 1;
				Type const x = (Type(i) + half) * half * cellWidth_;
				std::gslice slice(start, sizes, interiorStrides());
				position_[dim][slice] = x;
			}
		}
		for (FaceType face = FaceType::begin(); face != FaceType::end(); face++) {
			boundaryHandles[face] = this->getBoundaryHandle(face.flip());
		}
	}
	Subgrid createSinglePaddedSubgrid(FaceType face) {
		Subgrid gState = createSubgrid(singlePaddedVolume);
		hpx::future<Subgrid> bndFuture = boundaryHandles[face]();
		size_t const dim = face.dimension();
		size_t const bndStart = (face.direction() > 0) ? (interiorSizes()[dim] * singlePaddedStrides(face)[dim]) : size_t(0);
		size_t const intStart = (face.direction() < 0) ? singlePaddedStrides(face)[dim] : size_t(0);
		auto bndSizes = interiorSizes();
		bndSizes[dim] = 1;
		auto const intSlice = std::gslice(intStart, interiorSizes(), singlePaddedStrides(face));
		auto const bndSlice = std::gslice(bndStart, bndSizes, singlePaddedStrides(face));
		for (int field = 0; field < fieldCount; field++) {
			for (int mode = 0; mode < modeVolume; mode++) {
				gState[field][mode][intSlice] = np1State_[field][mode];
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
	Subgrid createDoublePaddedSubgrid(int dim) {
		Subgrid gState = createSubgrid(doublePaddedVolume);
		hpx::future<Subgrid> bndFuture1 = boundaryHandles[FaceType(dim, -1)]();
		hpx::future<Subgrid> bndFuture2 = boundaryHandles[FaceType(dim, +1)]();
		auto bndSizes = interiorSizes();
		bndSizes[dim] = 1;
		auto const intSlice = std::gslice(doublePaddedStrides(dim)[dim], interiorSizes(), doublePaddedStrides(dim));
		auto const bndSlice1 = std::gslice(0, bndSizes, doublePaddedStrides(dim));
		auto const bndSlice2 = std::gslice((interiorSizes()[dim] + 1) * doublePaddedStrides(dim)[dim], bndSizes, doublePaddedStrides(dim));
		for (int field = 0; field < fieldCount; field++) {
			for (int mode = 0; mode < modeVolume; mode++) {
				gState[field][mode][intSlice] = np1State_[field][mode];
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
	void initialize(std::function<State<Type>(std::array<Type, dimensionCount> const&)> const &initialState) {
		Type const halfCellWidth = Type(0.5) * cellWidth_;
		State<std::array<std::valarray<Type>, nodeVolume>> nodalValues = makeFilledArray<std::array<std::valarray<Type>, nodeVolume>, fieldCount>(
				makeFilledArray<std::valarray<Type>, nodeVolume>(std::valarray<Type>(interiorVolume)));
		for (int node = 0; node < nodeVolume; node++) {
			auto quadraturePoint = getQuadraturePoint<Type, dimensionCount, modeCount>(node);
			auto thisPosition = position_;
			for (int dim = 0; dim < dimensionCount; dim++) {
				thisPosition[dim] += quadraturePoint[dim] * halfCellWidth;
			}
			for (int i = 0; i < interiorVolume; i++) {
				std::array<Type, dimensionCount> x;
				for (int dim = 0; dim < dimensionCount; dim++) {
					x[dim] = thisPosition[dim][i];
				}
				auto const thisState = initialState(x);
				for (int field = 0; field < fieldCount; field++) {
					nodalValues[field][node][i] = thisState[field];
				}
			}
		}
		for (int field = 0; field < fieldCount; field++) {
			np1State_[field] = volumeAnalyze(nodalValues[field]);
		}
	}
	void output(const char *filenameBase, int timeStepNumber, Type const &time) {
		std::string filename = std::string(filenameBase) + "." + std::to_string(timeStepNumber) + ".h5";
		writeHdf5<Type, dimensionCount, interiorWidth, modeCount, State>(filename, cellWidth_, np1State_, State<Type>::getFieldNames());
		writeList("X.visit", "!NBLOCKS 1\n", filename + ".xmf");
	}
	BndHandle getBoundaryHandle(FaceType face) const {
		return BndHandle([this, face]() {
			return hpx::async([this, face]() {
				auto const dim = face.dimension();
				auto sizes = interiorSizes();
				sizes[dim] = 1;
				size_t const start = (face.direction() < 0) ? size_t(0) : ((interiorSizes()[dim] - 1) * interiorStrides()[dim]);
				std::gslice slice(start, sizes, interiorStrides());
				Subgrid bndState = createSubgrid(boundaryVolume);
				for (int field = 0; field < fieldCount; field++) {
					for (int mode = 0; mode < modeVolume; mode++) {
						bndState[field][mode] = np1State_[field][mode][slice];
					}
				}
				return bndState;
			});
		});
	}
	Type beginStep() {
		using std::max;
		std::array<Type, dimensionCount> maximumEigenvalue;
		std::array<State<std::valarray<Type>>, nodeVolume> stateAtNodes;
		nState_ = np1State_;
		maximumEigenvalue.fill(Type(0));
		for (int field = 0; field < fieldCount; field++) {
			auto const tmp = volumeSynthesize(np1State_[field]);
			for (int node = 0; node < nodeVolume; node++) {
				stateAtNodes[node][field] = tmp[node];
			}
		}
		for (int node = 0; node < nodeVolume; node++) {
			for (int dim = 0; dim < dimensionCount; dim++) {
				auto const eigenvalues = stateAtNodes[node].eigenvalues(dim);
				for (int i = 0; i < fieldCount; i++) {
					std::valarray<Type> const lambda = abs(eigenvalues[i]);
					maximumEigenvalue[dim] = max(maximumEigenvalue[dim], lambda.max());
				}
			}
		}
		Type maximumEigenvalueSum = Type(0);
		for (int dim = 0; dim < dimensionCount; dim++) {
			maximumEigenvalueSum += maximumEigenvalue[dim];
		}
		Type const timeStepSize = (cellWidth_ * butcherTable.cfl()) / (Type(2 * modeCount - 1) * maximumEigenvalueSum);
		return timeStepSize;
	}
	void subStep(Type const &timeStepSize, int stageIndex) {
		for (int field = 0; field < fieldCount; field++) {
			for (int mode = 0; mode < modeVolume; mode++) {
				np1State_[field][mode] = nState_[field][mode];
				for (int thisStage = 0; thisStage < stageIndex; thisStage++) {
					np1State_[field][mode] += butcherTable.a(stageIndex, thisStage) * stageStates_[thisStage][field][mode];
				}
				stageStates_[stageIndex][field][mode] = std::valarray<Type>(Type(0), interiorVolume);
			}
		}
		applyLimiter();
		computeTimeDerivative(timeStepSize, stageStates_[stageIndex]);
	}
	void endStep() {
		for (int field = 0; field < fieldCount; field++) {
			for (int mode = 0; mode < modeVolume; mode++) {
				np1State_[field][mode] = nState_[field][mode];
				for (int j = 0; j < rkStageCount; j++) {
					np1State_[field][mode] += butcherTable.b(j) * stageStates_[j][field][mode];
				}
			}
		}
		applyLimiter();
	}
	void applyLimiter() {
		if constexpr (modeVolume > 1) {
			auto workspace = LimiterWorkspaceType::getWorkspace();
			auto &alpha = workspace->alpha;
			auto &leftEigenvectors = workspace->leftEigenvectors;
			auto &rightEigenvectors = workspace->rightEigenvectors;
			auto &meanState = workspace->meanState;
			auto &wasLimited = workspace->wasLimited;
			auto &loState = workspace->loState;
			auto &hiState = workspace->hiState;
			auto &slopeState = workspace->slopeState;
			auto &theta = workspace->theta;
			auto &surfaceNodes = workspace->surfaceNodes;
			auto &volumeNodes = workspace->volumeNodes;
			auto &surfaceState = workspace->surfaceState;
			auto &volumeState = workspace->volumeState;
			auto &modes = workspace->modes;
			auto &thisAlpha = workspace->thisAlpha;
			auto &invState = workspace->invState;
			for (int field = 0; field < fieldCount; field++) {
				meanState[field] = np1State_[field][0];
			}
			for (int dim = 0; dim < dimensionCount; dim++) {
				rightEigenvectors[dim] = meanState.eigenSystem(dim).second;
				leftEigenvectors[dim] = matrixInverse(rightEigenvectors[dim]);
			}
			for (int dim = 0; dim < dimensionCount; dim++) {
				for (int mode = 0; mode < modeVolume; mode++) {
					for (int field = 0; field < fieldCount; field++) {
						alpha[field][mode][dim] = std::valarray<Type>(Type(0), interiorVolume);
					}
				}
			}
			for (int polynomialDegree = modeCount - 1; polynomialDegree > 0; polynomialDegree--) {
				constexpr auto tiny = Type(std::numeric_limits<double>::min());
				int const begin = binco(polynomialDegree + dimensionCount - 1, dimensionCount);
				int const end = binco(polynomialDegree + dimensionCount, dimensionCount);
				for (int hiMode = begin; hiMode < end; hiMode++) {
					auto hiModeIndices = flatToTriangular<dimensionCount, modeCount>(hiMode);
					for (int dim = 0; dim < dimensionCount; dim++) {
						if (hiModeIndices[dim] == 0) {
							continue;
						}
						auto loModeIndices = hiModeIndices;
						loModeIndices[dim]--;
						auto const loMode = triangularToFlat<dimensionCount, modeCount>(loModeIndices);
						auto gState = createDoublePaddedSubgrid(dim);
						int const stride = doublePaddedStrides(dim)[dim];
						std::gslice const intSlice(stride, interiorSizes(), doublePaddedStrides(dim));
						std::gslice const intPlusSlice(2 * stride, interiorSizes(), doublePaddedStrides(dim));
						std::gslice const intMinusSlice(0, interiorSizes(), doublePaddedStrides(dim));
						for (int field = 0; field < fieldCount; field++) {
							hiState[field] = gState[field][hiMode][intSlice];
							loState[field] = gState[field][loMode];
						}
						for (int field = 0; field < fieldCount; field++) {
							std::valarray<Type> const stateCentral = loState[field][intSlice];
							std::valarray<Type> const statePlus = loState[field][intPlusSlice];
							std::valarray<Type> const stateMinus = loState[field][intMinusSlice];
							slopeState[field] = minmod(std::valarray<Type>(statePlus - stateCentral), std::valarray<Type>(stateCentral - stateMinus));
						}
						for (int field = 0; field < fieldCount; field++) {
							invState[field] = Type(1) / (hiState[field] + tiny);
						}
						slopeState = leftEigenvectors[dim] * slopeState;
						hiState = leftEigenvectors[dim] * hiState;
						auto const cLimit = Type(1) / Type(2 * loModeIndices[dim] + 1);
						for (int field = 0; field < fieldCount; field++) {
							hiState[field] = minmod(hiState[field], std::valarray<Type>(slopeState[field] * cLimit));
						}
						hiState = rightEigenvectors[dim] * hiState;
						for (int field = 0; field < fieldCount; field++) {
							auto &ref = alpha[field][hiMode][dim];
							ref = max(ref,
									std::valarray<Type>(max(std::valarray<Type>(zero, interiorVolume), std::valarray<Type>(hiState[field] * invState[field]))));
						}
						wasLimited[dim][hiMode] = true;
					}
				}
			}
			for (int mode = 1; mode < modeVolume; mode++) {
				for (int field = 0; field < fieldCount; field++) {
					thisAlpha = Type(0);
					int count = 0;
					for (int dim = 0; dim < dimensionCount; dim++) {
						if (wasLimited[dim][mode]) {
							thisAlpha += alpha[field][mode][dim];
							count++;
						}
					}
					thisAlpha /= Type(count);
					np1State_[field][mode] *= thisAlpha;
				}
			}
			theta = Type(1);
			for (int field = 0; field < fieldCount; field++) {
				meanState[field] = np1State_[field][0];
			}
			for (FaceType face = FaceType::begin(); face != FaceType::end(); face++) {
				for (int field = 0; field < fieldCount; field++) {
					for (int mode = 0; mode < modeVolume; mode++) {
						modes[mode] = np1State_[field][mode];
					}
					surfaceNodes = surfaceSynthesize(surfaceTrace(face, modes));
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
					modes[mode] = np1State_[field][mode];
				}
				volumeNodes = volumeSynthesize(modes);
				for (int node = 0; node < nodeVolume; node++) {
					volumeState[node][field] = volumeNodes[node];
				}
			}
			for (int node = 0; node < nodeVolume; node++) {
				theta = min(theta, findPositivityPreservingTheta(meanState, volumeState[node]));
			}
			for (int field = 0; field < fieldCount; field++) {
				for (int mode = 1; mode < modeVolume; mode++) {
					np1State_[field][mode] *= theta;
				}
			}
			LimiterWorkspaceType::recycleWorkspace(std::move(workspace));
		}
	}
	struct FluxWorkspace {
		static hpx::mutex mutex;
		static std::stack<std::unique_ptr<FluxWorkspace>> workspaces;
		static std::unique_ptr<FluxWorkspace> getWorkspace() {
			std::unique_ptr<FluxWorkspace> pointer;
			std::lock_guard<hpx::mutex> lock(mutex);
			if (workspaces.empty()) {
				pointer = std::make_unique<FluxWorkspace>();
				workspaces.push(std::move(pointer));
			}
			pointer = std::move(workspaces.top());
			workspaces.pop();
			return pointer;

		}
		static void recycleWorkspace(std::unique_ptr<FluxWorkspace> &&pointer) {
			std::lock_guard<hpx::mutex> lock(mutex);
			workspaces.push(std::move(pointer));
		}
		State<std::valarray<Type>> tmp;
		State<std::array<std::valarray<Type>, nodeSurface>> nodalFlux;
		State<std::array<std::valarray<Type>, nodeVolume>> volumeFlux;
		std::array<std::valarray<Type>, modeSurface> pModalFaceFlux;
		std::array<std::valarray<Type>, modeSurface> mModalFaceFlux;
		std::array<std::valarray<Type>, modeSurface> modalFlux;
		std::array<std::valarray<Type>, nodeSurface> left;
		std::array<std::valarray<Type>, nodeSurface> right;
		std::array<std::valarray<Type>, modeVolume> pModalVolumeFlux;
		std::array<std::valarray<Type>, modeVolume> mModalVolumeFlux;
		std::array<std::valarray<Type>, modeVolume> volumeModes;
		std::array<std::valarray<Type>, nodeVolume> volumeNodes;
		std::array<std::valarray<Type>, binco(modeCount + dimensionCount - 2, dimensionCount)> source;
		std::array<State<std::valarray<Type>>, nodeSurface> lState;
		std::array<State<std::valarray<Type>>, nodeSurface> rState;
		std::array<State<std::valarray<Type>>, nodeVolume> volumeState;
		Subgrid lModes;
		Subgrid rModes;
	};
	Subgrid computeTimeDerivative(Type timeStepSize, Subgrid &stageState) {
		auto workspace = FluxWorkspace::getWorkspace();
		auto &pModalFaceFlux = workspace->pModalFaceFlux;
		auto &mModalFaceFlux = workspace->mModalFaceFlux;
		auto &pModalVolumeFlux = workspace->pModalVolumeFlux;
		auto &mModalVolumeFlux = workspace->mModalVolumeFlux;
		auto &nodalFlux = workspace->nodalFlux;
		auto &modalFlux = workspace->modalFlux;
		auto &lState = workspace->lState;
		auto &rState = workspace->rState;
		auto &left = workspace->left;
		auto &right = workspace->right;
		auto &tmpState = workspace->tmp;
		auto &volumeNodes = workspace->volumeNodes;
		auto &volumeModes = workspace->volumeModes;
		auto &volumeState = workspace->volumeState;
		auto &volumeFlux = workspace->volumeFlux;
		auto &lModes = workspace->lModes;
		auto &rModes = workspace->rModes;
		auto &source = workspace->source;
		Type const lambda = Type(2) * timeStepSize / cellWidth_;
		{
			for (int dim = 0; dim < dimensionCount; dim++) {
				lModes = createSinglePaddedSubgrid(FaceType(dim, -1));
				rModes = createSinglePaddedSubgrid(FaceType(dim, +1));
				for (int field = 0; field < fieldCount; field++) {
					left = surfaceSynthesize(surfaceTrace(FaceType(dim, +1), lModes[field]));
					right = surfaceSynthesize(surfaceTrace(FaceType(dim, -1), rModes[field]));
					for (int node = 0; node < nodeSurface; node++) {
						lState[node][field] = left[node];
						rState[node][field] = right[node];
					}
				}
				for (int node = 0; node < nodeSurface; node++) {
					tmpState = solveRiemannProblem(lState[node], rState[node], dim);
					for (int field = 0; field < fieldCount; field++) {
						nodalFlux[field][node] = tmpState[field];
					}
				}
				for (int field = 0; field < fieldCount; field++) {
					modalFlux = surfaceAnalyze(nodalFlux[field]);
					for (int mode = 0; mode < modeSurface; mode++) {
						pModalFaceFlux[mode] = std::valarray<Type>(modalFlux[mode][fluxPlusSlice(dim)]);
						mModalFaceFlux[mode] = std::valarray<Type>(modalFlux[mode][fluxMinusSlice(dim)]);
					}
					pModalVolumeFlux = surfaceTraceInverse(FaceType(dim, +1), pModalFaceFlux);
					mModalVolumeFlux = surfaceTraceInverse(FaceType(dim, -1), mModalFaceFlux);
					for (int mode = 0; mode < modeVolume; mode++) {
						stageState[field][mode] -= lambda * (pModalVolumeFlux[mode] - mModalVolumeFlux[mode]);
					}
				}
			}
		}
		{
			for (int dim = 0; dim < dimensionCount; dim++) {
				for (int field = 0; field < fieldCount; field++) {
					for (int mode = 0; mode < modeVolume; mode++) {
						volumeModes[mode] = np1State_[field][mode];
					}
					volumeNodes = volumeSynthesize(volumeModes);
					for (int node = 0; node < nodeVolume; node++) {
						volumeState[node][field] = std::move(volumeNodes[node]);
					}
				}
				for (int node = 0; node < nodeVolume; node++) {
					tmpState = volumeState[node].flux(dim);
					for (int field = 0; field < fieldCount; field++) {
						volumeFlux[field][node] = std::move(tmpState[field]);
					}
				}
				for (int field = 0; field < fieldCount; field++) {
					int const loModeVolume = binco(modeCount + dimensionCount - 2, dimensionCount);
					source = volumeAnalyzeDerivative(dim, volumeFlux[field]);
					for (int loMode = 0; loMode < loModeVolume; loMode++) {
						auto modeIndices = flatToTriangular<dimensionCount, modeCount - 1>(loMode);
						modeIndices[dim]++;
						auto const hiMode = triangularToFlat<dimensionCount, modeCount>(modeIndices);
						stageState[field][hiMode] += lambda * source[loMode];
					}
				}
			}
		}
		FluxWorkspace::recycleWorkspace(std::move(workspace));
		return stageState;
	}
};

template<typename T, int D, int W, int M, typename RK, template<typename, int> typename S>
hpx::mutex HyperSubgrid<T, D, W, M, RK, S>::FluxWorkspace::mutex;

template<typename T, int D, int W, int M, typename RK, template<typename, int> typename S>
std::stack<std::unique_ptr<typename HyperSubgrid<T, D, W, M, RK, S>::FluxWorkspace>> HyperSubgrid<T, D, W, M, RK, S>::FluxWorkspace::workspaces;
