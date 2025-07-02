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
	static constexpr int extVolume = ipow(intWidth + 1, dimCount);
	static constexpr int bndVolume = ipow(intWidth, dimCount - 1);
	static constexpr int flxVolume = ipow(intWidth, dimCount - 1) * (intWidth + 1);
	static constexpr int intSurface = ipow(intWidth, dimCount - 1);
	static constexpr int extSurface = ipow(intWidth + 1, dimCount - 1);
	static constexpr RungeKutta butcherTable { };
	template<int volume = intVolume>
	static constexpr auto vAnalyze = dgAnalyze<ValArray<Type, volume>, dimCount, modeCount>;
	template<int volume = intVolume>
	static constexpr auto vMassInverse = dgMassInverse<ValArray<Type, volume>, dimCount, modeCount>;
	template<int volume = intVolume>
	static constexpr auto vSynthesize = dgSynthesize<ValArray<Type, volume>, dimCount, modeCount>;
	template<int volume = intVolume>
	static constexpr auto vStiffness = dgStiffness<ValArray<Type, volume>, dimCount, modeCount>;
	template<int volume = intVolume>
	static constexpr auto sAnalyze = dgAnalyze<ValArray<Type, volume>, dimCount - 1, modeCount>;
	template<int volume = intVolume>
	static constexpr auto sSynthesize = dgSynthesize<ValArray<Type, volume>, dimCount - 1, modeCount>;
	template<int volume = intVolume>
	static constexpr auto sTrace = dgTrace<ValArray<Type, volume>, dimCount, modeCount>;
	template<int volume = intVolume>
	static constexpr auto sTraceInverse = dgTraceInverse<ValArray<Type, volume>, dimCount, modeCount>;

	using IntArray = ValArray<Type, intVolume>;
	using ExtArray = ValArray<Type, extVolume>;
	using BndArray = ValArray<Type, bndVolume>;
	using FlxArray = ValArray<Type, flxVolume>;
	using IntGridType = State<std::array<IntArray, modeVolume>, dimCount>;
	using ExtGridType = State<std::array<ExtArray, modeVolume>, dimCount>;
	using BndGridType = State<std::array<ValArray<Type, intSurface>, modeVolume>, dimCount>;
	using FlxGridType = State<std::array<ValArray<Type, flxVolume>, modeVolume>, dimCount>;
	using BaseType = typename ElementType<Type>::type;
	using BoundaryGrid = State<std::array<ValArray<Type, intSurface>, modeVolume>, dimCount>;
	using BoundaryHandle = std::function<hpx::future<BoundaryGrid>(Face)>;

	Type cellWidth;
	IntGridType np1State;
	IntGridType nState;
	std::array<IntGridType, rkStageCount> jState;
	int64_t corner;
	ValArray<int64_t, dimCount> intSizes;
	ValArray<int64_t, dimCount> intStrides;
	ValArray<int64_t, dimCount> extSizes;
	ValArray<int64_t, dimCount> extStrides;
	std::array<ValArray<Type, intVolume>, dimCount> position;
	std::array<BoundaryHandle, faceCount(dimCount)> boundaryHandles;
	GSlice<dimCount, intVolume> intSlice;

	static constexpr BaseType zero = BaseType(0);
	static constexpr BaseType one = BaseType(1);
	static constexpr BaseType two = BaseType(2);
	static constexpr BaseType half = one / two;
public:
	HyperGrid(Type const &xNint = Type(1)) :
			cellWidth { xNint / Type(intWidth) } {
		extSizes.fill(intWidth + 2);
		intSizes.fill(intWidth);
		extStrides[dimCount - 1] = 1;
		intStrides[dimCount - 1] = 1;
		for (int n = dimCount - 1; n > 0; n--) {
			extStrides[n - 1] = extStrides[n] * extSizes[n];
			intStrides[n - 1] = intStrides[n] * intSizes[n];
		}
		corner = std::accumulate(extSizes.begin(), extSizes.end(), int64_t(1), std::multiplies<int64_t>());
		intSlice = GSlice<dimCount, intVolume>(corner, intSizes, extStrides);
		for (int dim = 0; dim < dimCount; dim++) {
			for (int i = 0; i < intWidth; i++) {
				auto const start = i * intStrides[dim];
				auto sizes = intSizes;
				sizes[dim] = 1;
				BaseType const x = (BaseType(i) + half) * half * cellWidth;
				GSlice<dimCount, intSurface> slice(start, sizes, intStrides);
				position[dim][slice] = x;
			}
		}
		for (Face face = 0; face < faceCount(dimCount); face++) {
			boundaryHandles[face] = this->getBoundaryHandle(flipFace(face));
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
			auto const tmp = vAnalyze<>(nodalValues[field]);
			np1State[field] = vMassInverse<>(tmp);
		}
	}
	void output(const char *filenameBase, int timeStepNumber, Type const &time) {
		std::string filename = std::string(filenameBase) + "." + std::to_string(timeStepNumber) + ".h5";
		writeHdf5<Type, dimCount, intWidth, modeCount, State>(filename, cellWidth, np1State, State<Type, dimCount>::getFieldNames());
		writeList("X.visit", "!NBLOCKS 1\n", filename + ".xmf");
	}
	BoundaryHandle getBoundaryHandle(Face face) const {
		return BoundaryHandle([this](Face face) {
			return hpx::async([this, face]() {
				auto const dim = getFaceDim(face);
				auto sizes = intSizes;
				sizes[dim] = 1;
				GSlice<dimCount, intSurface> slice(corner, sizes, extStrides);
				BndGridType bndState;
				for (int field = 0; field < fieldCount; field++) {
					for (int mode = 0; mode < modeVolume; mode++) {
						bndState[field][mode] = np1State[field][mode][slice];
					}
				}
				return bndState;
			});
		});
	}
	ExtGridType createGhostState() {
		constexpr int nFaces = faceCount(dimCount);
		std::array<GSlice<dimCount, intSurface>, nFaces> bndSlices;
		ExtGridType gState;
		std::vector<hpx::future<BndGridType>> bndFutures(nFaces);
		for (Face face = 0; face < nFaces; face++) {
			bndFutures[face] = boundaryHandles[face](nFaces);
		}
		for (Face face = 0; face < faceCount(dimCount); face++) {
			int64_t const dim = getFaceDim(face);
			int64_t const start = (getFaceDir(face) < 0) ? int64_t(0) : (extSizes[dim] * extStrides[dim]);
			auto sizes = intSizes;
			sizes[dim] = 1;
			bndSlices[face] = GSlice<dimCount, intSurface>(start, sizes, extStrides);
		}
		for (int field = 0; field < fieldCount; field++) {
			for (int mode = 0; mode < modeVolume; mode++) {
				gState[field][mode][intSlice] = np1State[field][mode];
			}
		}
		for (Face face = 0; face < nFaces; face++) {
			auto const bndState = bndFutures[face].get();
			for (int field = 0; field < fieldCount; field++) {
				for (int mode = 0; mode < modeVolume; mode++) {
					gState[field][mode][bndSlices[face]] = bndState[field][mode];
				}
			}
		}
		return gState;
	}
	void retireGhostState(ExtGridType &&gState) {
		for (int field = 0; field < fieldCount; field++) {
			for (int mode = 0; mode < modeVolume; mode++) {
				np1State[field][mode] = gState[field][mode][intSlice];
			}
		}
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
		GSlice<dimCount, intVolume> const intSlice(corner, intSizes, extStrides);
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
		jState[stageIndex] = computeTimeDerivative(timeStepSize);
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
	}
	void applyLimiter() {
		{
			auto gState = createGhostState();
			std::array<std::array<State<IntArray, dimCount>, modeVolume>, dimCount> alpha;
			std::array<SquareMatrix<IntArray, fieldCount>, dimCount> leftEigenvectors;
			std::array<SquareMatrix<IntArray, fieldCount>, dimCount> rightEigenvectors;
			State<IntArray, dimCount> meanState;
			State<IntArray, dimCount> transposedState;
			State<IntArray, dimCount> differenceState;
			std::array<std::bitset<modeCount>, dimCount> wasLimited;
			for (int field = 0; field < fieldCount; field++) {
				meanState[field] = gState[field][0][intSlice];
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
						GSlice<dimCount, intVolume> const intPlusSlice(corner + extStrides[dim], intSizes, extStrides);
						GSlice<dimCount, intVolume> const intMinusSlice(corner - extStrides[dim], intSizes, extStrides);
						for (int field = 0; field < fieldCount; field++) {
							transposedState[field] = gState[field][hiMode][intSlice];
						}
						auto loModeIndices = hiModeIndices;
						loModeIndices[dim]--;
						auto const loMode = triangularToFlat<dimCount, modeCount>(loModeIndices);
						for (int field = 0; field < fieldCount; field++) {
							IntArray const stateCentral = gState[field][loMode][intSlice];
							IntArray const statePlus = gState[field][loMode][intPlusSlice];
							IntArray const stateMinus = gState[field][loMode][intMinusSlice];
							differenceState[field] = minmod(IntArray(statePlus - stateCentral), IntArray(stateCentral - stateMinus));
						}
						State<IntArray, dimCount> initialStateInverse;
						for (int field = 0; field < fieldCount; field++) {
							initialStateInverse[field] = BaseType(1) / (transposedState[field] + tiny);
						}
						differenceState = leftEigenvectors[dim] * differenceState;
						transposedState = leftEigenvectors[dim] * transposedState;
						auto const cLimit = Type(1) / Type(2 * loModeIndices[dim] + 1);
						for (int field = 0; field < fieldCount; field++) {
							differenceState[field] *= cLimit;
							transposedState[field] = minmod(transposedState[field], differenceState[field]);
						}
						transposedState = rightEigenvectors[dim] * transposedState;
						for (int field = 0; field < fieldCount; field++) {
							IntArray alphaHi, alphaLo;
							alphaHi = alpha[dim][hiMode][field];
							alphaLo = alpha[dim][loMode][field];
							alphaHi = max(alphaHi, max(BaseType(0.0), transposedState[field] * initialStateInverse[field]));
							alphaLo = max(alphaLo, alphaHi);
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
					gState[field][mode][intSlice] *= (thisAlpha / BaseType(count));
				}
			}
			retireGhostState(std::move(gState));
		}
		{
			ValArray<Type, intVolume> theta(Type(1));
			std::array<State<ValArray<Type, intVolume>, dimCount>, nodeSurface> surfaceState;
			std::array<State<ValArray<Type, intVolume>, dimCount>, nodeVolume> volumeState;
			State<ValArray<Type, intVolume>, dimCount> meanState;
			std::array<IntArray, modeVolume> modes;
			for (int field = 0; field < fieldCount; field++) {
				meanState[field] = np1State[field][0];
			}
			for (Face face = 0; face < faceCount(dimCount); face++) {
				for (int field = 0; field < fieldCount; field++) {
					for (int mode = 0; mode < modeVolume; mode++) {
						modes[mode] = np1State[field][mode][intSlice];
					}
					auto const surfaceNodes = sSynthesize<>(sTrace<>(face, modes));
					for (int node = 0; node < nodeSurface; node++) {
						surfaceState[node][field] = surfaceNodes[node];
					}
				}
				for (int node = 0; node < nodeSurface; node++) {
					theta = min(BaseType(1), findPositivityPreservingTheta(meanState, surfaceState[node]));
				}
			}
			for (int field = 0; field < fieldCount; field++) {
				for (int mode = 0; mode < modeVolume; mode++) {
					modes[mode] = np1State[field][mode][intSlice];
				}
				auto const volumeNodes = vSynthesize<intVolume>(modes);
				for (int node = 0; node < nodeVolume; node++) {
					volumeState[node][field] = volumeNodes[node];
				}
			}
			for (int node = 0; node < nodeVolume; node++) {
				theta = min(BaseType(1), findPositivityPreservingTheta(meanState, volumeState[node]));
			}
			for (int field = 0; field < fieldCount; field++) {
				for (int mode = 1; mode < modeVolume; mode++) {
					np1State[field][mode][intSlice] *= theta;
				}
			}
		}
	}
	IntGridType computeTimeDerivative(Type timeStepSize) {
		IntGridType kState;
		Type const lambda = Type(2) * timeStepSize / cellWidth;
		{
			auto gState = createGhostState();
			for (int dim = 0; dim < dimCount; dim++) {
				GSlice<dimCount, flxVolume> const fluxLeftSlice(corner - extStrides[dim], flxSizes(dim), extStrides[dim]);
				GSlice<dimCount, flxVolume> const fluxRightSlice(corner, flxSizes(dim), extStrides[dim]);
				std::array<ValArray<Type, intVolume>, modeSurface> pModalFaceFlux;
				std::array<ValArray<Type, intVolume>, modeSurface> mModalFaceFlux;
				std::array<ValArray<Type, flxVolume>, modeVolume> lModes;
				std::array<ValArray<Type, flxVolume>, modeVolume> rModes;
				std::array<ValArray<Type, intVolume>, modeVolume> pModalVolumeFlux;
				std::array<ValArray<Type, intVolume>, modeVolume> mModalVolumeFlux;
				State<std::array<ValArray<Type, flxVolume>, nodeSurface>, dimCount> nodalFlux;
				State<std::array<ValArray<Type, intVolume>, nodeSurface>, dimCount> lNodalFlux;
				State<std::array<ValArray<Type, intVolume>, nodeSurface>, dimCount> rNodalFlux;
				std::array<State<ValArray<Type, flxVolume>, dimCount>, nodeSurface> lState;
				std::array<State<ValArray<Type, flxVolume>, dimCount>, nodeSurface> rState;
				for (int field = 0; field < fieldCount; field++) {
					for (int mode = 0; mode < modeVolume; mode++) {
						lModes[mode] = gState[field][mode][fluxLeftSlice];
						rModes[mode] = gState[field][mode][fluxRightSlice];
					}
					auto const left = sSynthesize<flxVolume>(sTrace<flxVolume>(2 * dim + 1, lModes));
					auto const right = sSynthesize<flxVolume>(sTrace<flxVolume>(2 * dim + 0, rModes));
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
				GSlice<dimCount, intVolume> const fluxPlusSlice(flxStrides(dim)[dim], intSizes, flxStrides(dim));
				GSlice<dimCount, intVolume> const fluxMinusSlice(0, intSizes, flxStrides(dim));
				for (int field = 0; field < fieldCount; field++) {
					auto const modalFlux = sAnalyze<flxVolume>(nodalFlux[field]);
					for (int mode = 0; mode < modeSurface; mode++) {
						pModalFaceFlux[mode] = modalFlux[mode][fluxPlusSlice];
						mModalFaceFlux[mode] = modalFlux[mode][fluxMinusSlice];
					}
					pModalVolumeFlux = sTraceInverse<intVolume>(2 * dim + 1, pModalFaceFlux);
					mModalVolumeFlux = sTraceInverse<intVolume>(2 * dim, mModalFaceFlux);
					pModalVolumeFlux = vMassInverse<intVolume>(pModalVolumeFlux);
					mModalVolumeFlux = vMassInverse<intVolume>(mModalVolumeFlux);
					for (int mode = 0; mode < modeVolume; mode++) {
						kState[field][mode] -= lambda * (pModalVolumeFlux[mode] - mModalVolumeFlux[mode]);
					}
				}
			}
			retireGhostState(std::move(gState));
		}
		{
			for (int dim = 0; dim < dimCount; dim++) {
				std::array<State<ValArray<Type, intVolume>, dimCount>, nodeVolume> volumeState;
				State<std::array<ValArray<Type, intVolume>, nodeVolume>, dimCount> volumeFlux;
				for (int field = 0; field < fieldCount; field++) {
					std::array<ValArray<Type, intVolume>, modeVolume> volumeModes;
					for (int mode = 0; mode < modeVolume; mode++) {
						volumeModes[mode] = np1State[field][mode][intSlice];
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
						kState[field][mode][intSlice] += lambda * source[mode];
					}
				}
			}
		}
		return kState;
	}
}
;

