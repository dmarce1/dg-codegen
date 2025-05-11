#include "IndexTuple.hpp"

#include <numbers>
#include <type_traits>
#include <utility>

template<int M>
struct Octant {
	constexpr Octant(int j) :
			i(j) {
	}
	constexpr Octant() {
	}
	inline constexpr Octant& operator++() {
		i++;
		return *this;
	}
	inline constexpr Octant& operator--() {
		i--;
		return *this;
	}
	inline constexpr Octant operator++(int) {
		auto const rc = *this;
		operator++();
		return rc;
	}
	inline constexpr Octant operator--(int) {
		auto const rc = *this;
		operator--();
		return rc;
	}
	inline constexpr bool operator[](int j) const {
		return (1 & (i >> j));
	}
	inline static constexpr Octant begin() {
		Octant o;
		o.i = 0;
		return o;
	}
	inline static constexpr Octant end() {
		Octant o;
		o.i = 1 << M;
		return o;
	}
	inline constexpr bool operator==(Octant const &other) const {
		return i == other.i;
	}
	inline constexpr bool operator!=(Octant const &other) const {
		return i != other.i;
	}
	inline constexpr bool operator<(Octant const &other) const {
		return i < other.i;
	}
	inline constexpr bool operator>(Octant const &other) const {
		return i > other.i;
	}
	inline constexpr bool operator<=(Octant const &other) const {
		return i <= other.i;
	}
	inline constexpr bool operator>=(Octant const &other) const {
		return i >= other.i;
	}
	inline constexpr static int size() {
		return 1 << M;
	}
	inline constexpr operator int() const {
		return i;
	}
private:
	int i;
};

template<int D, int N>
TriangularIndex<N, D> legendreToTriangular(std::array<int, D> const &orders) {
	TriangularIndex<N, D> tri { };
	for (int i = 0; i < D; i++) {
		for (int j = i; j < D; j++) {
			tri[i] += orders[j];
		}
	}
	return tri;
}

template<int D, int N>
inline constexpr std::array<int, D> triangularToLegendre(TriangularIndex<N, D> p) {
	for (int d = D - 1; d > 0; d--) {
		p[d - 1] -= p[d];
	}
	return p;
}

template<typename T, int D, int N>
struct Basis {
	static constexpr int Nh = Math::integerPower(N, D);
	static constexpr int Np = Math::binCo<int>(D + N - 1, D);
	static constexpr int Nf = Math::integerPower(N, D - 1);
	static constexpr int Nc = 1 << D;
	using hindex_type = IndexTuple<N, D>;
	using pindex_type = TriangularIndex<N, D>;
	using findex_type = IndexTuple<N, D - 1>;
	using basis_type = std::array<T, Np>;
	using face_type = std::array<T, Nf>;
	using volume_type = std::array<T, Nh>;
	using refined_type = std::array<basis_type, Nc>;
	template<int M>
	using vector_type = std::array<T, M>;
	static T constexpr zero = T(0);
	static T constexpr half = T(0.5);
	static T constexpr one = T(1);
	static T constexpr two = T(2);
	static T constexpr Pi = std::numbers::pi;
	static inline constexpr T legendreP(int n, T x) {
		if (n == 0) {
			return one;
		} else if (n == 1) {
			return x;
		} else {
			T Pnm1, Pn, Pnp1;
			Pnm1 = one;
			Pn = x;
			for (int k = 1; k < n; k++) {
				Pnp1 = (T(2 * k + 1) * Pn * x - T(k) * Pnm1) / T(k + 1);
				Pnm1 = Pn;
				Pn = Pnp1;
			}
			return Pn;
		}
	}
	static inline constexpr T dLegendrePdX(int n, T x) {
		if (n == 0) {
			return zero;
		} else if (n == 1) {
			return one;
		} else {
			T Pnm1, Pn, Pnp1, dPndX, dPnp1dX;
			Pnm1 = dPndX = one;
			Pn = x;
			for (int k = 1; k < n - 1; k++) {
				Pnp1 = (T(2 * k + 1) * Pn * x - T(k) * Pnm1) / T(k + 1);
				dPnp1dX = T(k + 1) * Pn + x * dPndX;
				Pnm1 = Pn;
				Pn = Pnp1;
				dPndX = dPnp1dX;
			}
			return T(n) * Pn + x * dPndX;
		}
	}
	template<int D1>
	static inline constexpr T legendreP(TriangularIndex<N, D1> p, Math::Vector<T, D1> x) {
		auto const n = triangularToLegendre(p);
		T value = one;
		for (int d = 0; d < D1; d++) {
			value *= legendreP(n[d], x[d]);
		}
		return value;
	}
	template<typename U, typename F>
	static constexpr std::array<U, Nh> toQuadPoints(F const &f) {
		std::array<U, Nh> hValues;
		for (hindex_type hIdx = hindex_type::begin(); hIdx != hindex_type::end(); hIdx++) {
			int const hi = hIdx.flatIndex();
			hValues[hi] = f(volumePoint[hi]);
		}
		return hValues;
	}
	template<typename U, typename F>
	static constexpr std::array<U, Np> volumeTransform(F const &f) {
		std::array<U, Np> pCo;
		auto const hValues = toQuadPoints<U>(f);
		for (pindex_type pIdx = pindex_type::begin(); pIdx != pindex_type::end(); pIdx++) {
			int const pi = pIdx.flatIndex();
			pCo[pi] = U(0);
			assert(pi < (int ) pCo.size());
			for (hindex_type hIdx = hindex_type::begin(); hIdx != hindex_type::end(); hIdx++) {
				int const hi = hIdx.flatIndex();
				assert(hi < (int ) hValues.size());
				assert(hi < (int ) phi[pi].size());
				assert(hi < (int ) volumeWeight.size());
				assert(pi < (int ) norm.size());
				pCo[pi] += norm[pi] * hValues[hi] * phi[pi][hi] * volumeWeight[hi];
			}
		}
		return pCo;
	}
	template<typename F>
	static constexpr basis_type volumeTransform(F const &f) {
		return volumeTransform<T>(f);
	}
	template<typename F>
	static constexpr volume_type inverseVolumeTransform(basis_type const &pC) {
		volume_type hValues;
		for (hindex_type hIdx = hindex_type::begin(); hIdx != hindex_type::end(); hIdx++) {
			int const hi = hIdx.flatIndex();
			hValues[hi] = zero;
			for (pindex_type pIdx = pindex_type::begin(); pIdx != pindex_type::end(); pIdx++) {
				int const pi = pIdx.flatIndex();
				hValues[hi] += pC[pi] * phi[pi][hi];
			}
		}
		return hValues;
	}
	static constexpr face_type inverseFaceTransform(basis_type const &pC, int face) {
		face_type fValues;
		for (findex_type fIdx = findex_type::begin(); fIdx != findex_type::end(); fIdx++) {
			int const fi = fIdx.flatIndex();
			fValues[fi] = zero;
			for (pindex_type pIdx = pindex_type::begin(); pIdx != pindex_type::end(); pIdx++) {
				int const pi = pIdx.flatIndex();
				fValues[fi] += pC[pi] * xi[face][pi][fi];
			}
		}
		return fValues;
	}
	static constexpr auto& getQuadraturePoints() {
		return volumeRules.second;
	}
	static constexpr T faceIntegrate(face_type const &phi) {
		T value = zero;
		for (int i = 0; i < Nf; i++) {
			value += half * faceWeight[i] * phi[i];
		}
		return value;
	}
	static constexpr basis_type weakDivergence(std::array<volume_type, D> const &f) {
		basis_type coeffs;
		static constexpr T C = one / T(1 << D);
		for (int pi = 0; pi < Np; pi++) {
			coeffs[pi] = zero;
			for (int hi = 0; hi < Nh; hi++) {
				for (int d = 0; d < D; d++) {
					coeffs[pi] += C * volumeWeight[hi] * f[d][hi] * dphi_dx[d][pi][hi];
				}
			}
		}
		return coeffs;
	}
	static constexpr basis_type restrict_(refined_type const &pFine) {
		static constexpr auto chi = []() {
			std::array<std::array<std::array<T, Np>, Np>, Octant<D>::size()> chi;
			for (pindex_type pIdx = pindex_type::begin(); pIdx != pindex_type::end(); pIdx++) {
				int const pi = pIdx.flatIndex();
				for (pindex_type qIdx = pindex_type::begin(); qIdx != pindex_type::end(); qIdx++) {
					int const qi = qIdx.flatIndex();
					for (Octant<D> o = Octant<D>::begin(); o != Octant<D>::end(); o++) {
						auto &ref = chi[o][pi][qi];
						ref = zero;
						for (hindex_type hIdx = hindex_type::begin(); hIdx != hindex_type::end(); hIdx++) {
							int const hi = hIdx.flatIndex();
							auto x = volumePoint[hi];
							T w = one;
							for (int d = 0; d < D; d++) {
								int const i = (o[d] ? +1 : -1) * ((x[d] > zero) ? +1 : ((x[d] < zero) ? -1 : 0));
								x[d] = two * x[d] + T(1 - 2 * int(o[d]));
								if (i == 0) {
									w *= half;
								} else if (i < 0) {
									w = zero;
									break;
								}
							}
							ref += w * norm[pi] * legendreP<D>(qIdx, x) * volumeWeight[hi] * phi[pi][hi];
						}
					}
				}
			}
			return chi;
		}();
		basis_type pCoarse;
		for (pindex_type qIdx = pindex_type::begin(); qIdx != pindex_type::end(); qIdx++) {
			int const qi = qIdx.flatIndex();
			pCoarse[qi] = zero;
			for (Octant<D> oct = Octant<D>::begin(); oct != Octant<D>::end(); oct++) {
				for (pindex_type pIdx = pindex_type::begin(); pIdx != pindex_type::end(); pIdx++) {
					int const pi = pIdx.flatIndex();
					pCoarse[qi] += chi[oct][qi][pi] * pFine[oct][pi];
				}
			}
		}
		return pCoarse;
	}
	static constexpr refined_type prolong(basis_type const &pCoarse) {
		static constexpr auto psi = []() {
			std::array<std::array<std::array<T, Np>, Np>, Octant<D>::size()> psi;
			for (Octant<D> o = Octant<D>::begin(); o != Octant<D>::end(); o++) {
				for (pindex_type pIdx = pindex_type::begin(); pIdx != pindex_type::end(); pIdx++) {
					int const pi = pIdx.flatIndex();
					for (pindex_type qIdx = pindex_type::begin(); qIdx != pindex_type::end(); qIdx++) {
						int const qi = qIdx.flatIndex();
						auto &ref = psi[o][pi][qi];
						ref = zero;
						for (hindex_type hIdx = hindex_type::begin(); hIdx != hindex_type::end(); hIdx++) {
							int const hi = hIdx.flatIndex();
							auto x = volumePoint[hi];
							for (int d = 0; d < D; d++) {
								x[d] = half * (x[d] + T(2 * int(o[d]) - 1));
							}
							ref += norm[pi] * legendreP<D>(qIdx, x) * volumeWeight[hi] * phi[pi][hi];
						}
					}
				}
			}
			return psi;
		}();
		refined_type pFine;
		for (Octant<D> oct = Octant<D>::begin(); oct != Octant<D>::end(); oct++) {
			for (pindex_type pIdx = pindex_type::begin(); pIdx != pindex_type::end(); pIdx++) {
				int const pi = pIdx.flatIndex();
				pFine[oct][pi] = zero;
				for (pindex_type qIdx = pindex_type::begin(); qIdx != pindex_type::end(); qIdx++) {
					int const qi = qIdx.flatIndex();
					pFine[oct][pi] += psi[oct][pi][qi] * pCoarse[qi];
				}
			}
		}
		return pFine;
	}
private:
	template<int D1>
	static constexpr auto computeRules() {
		static T constexpr toler = T(1e-12);
		static T constexpr dx = one / T(N);
		static constexpr int Nh = Math::integerPower(N, D1);
		using hindex_type = IndexTuple<N, D1>;
		std::pair<volume_type, std::array<Math::Vector<T, D1>, Nh>> rules;
		auto &wt = rules.first;
		auto &pt = rules.second;
		std::array<T, N> weight;
		std::array<T, N> pos;
		for (int l = 0; l < N; l++) {
			T dTheta;
			T theta = Pi * (one - (T(l) + half) * dx);
			do {
				T const cosTheta = cos(theta);
				T const sinTheta = sin(theta);
				T const Pn = legendreP(N, cosTheta);
				T const dPnDx = dLegendrePdX(N, cosTheta);
				T const iDen = one / (sinTheta * dPnDx);
				dTheta = Pn * iDen;
				theta += dTheta;
			} while (abs(dTheta) > toler);
			T const x = cos(theta);
			T const q = dLegendrePdX(N, x);
			pos[l] = x;
			weight[l] = two / (q * q * (one - x * x));
		}
		for (hindex_type hIdx = hindex_type::begin(); hIdx != hindex_type::end(); hIdx++) {
			int const hi = hIdx.flatIndex();
			wt[hi] = one;
			for (int d = 0; d < D1; d++) {
				int const dhi = hIdx[d];
				wt[hi] *= weight[dhi];
				pt[hi][d] = pos[dhi];
			}
		}
		return rules;
	}
	static std::pair<volume_type, std::array<Math::Vector<T, D>, Nh>> const volumeRules;
	static std::pair<volume_type, std::array<Math::Vector<T, D - 1>, Nh>> const faceRules;
	static volume_type const volumeWeight;
	static std::array<Math::Vector<T, D>, Nh> const volumePoint;
	static volume_type const faceWeight;
	static std::array<Math::Vector<T, D - 1>, Nh> const facePoint;
	static constexpr auto norm = []() {
		std::array<T, Np> norm;
		for (pindex_type pIdx = pindex_type::begin(); pIdx != pindex_type::end(); pIdx++) {
			int const pi = pIdx.flatIndex();
			auto const qIdx = triangularToLegendre<D>(pIdx);
			T n = one;
			for (int d = 0; d < D; d++) {
				n *= T(half + qIdx[d]);
			}
			norm[pi] = n;
		}
		return norm;
	}();
	static std::array<volume_type, Np> const phi;
	static std::array<std::array<volume_type, Np>, D> const dphi_dx;
	static std::array<std::array<face_type, Np>, 2 * D> const xi;
}
;

template<typename T, int D, int N>
std::pair<typename Basis<T, D, N>::volume_type, std::array<Math::Vector<T, D>, Basis<T, D, N>::Nh>> const Basis<T, D, N>::volumeRules =
		Basis<T, D, N>::computeRules<D>();

template<typename T, int D, int N>
std::pair<typename Basis<T, D, N>::volume_type, std::array<Math::Vector<T, D - 1>, Basis<T, D, N>::Nh>> const Basis<T, D, N>::faceRules =
		Basis<T, D, N>::computeRules<D - 1>();

template<typename T, int D, int N>
Basis<T, D, N>::volume_type const Basis<T, D, N>::volumeWeight = volumeRules.first;

template<typename T, int D, int N>
std::array<Math::Vector<T, D>, Basis<T, D, N>::Nh> const Basis<T, D, N>::volumePoint = volumeRules.second;

template<typename T, int D, int N>
Basis<T, D, N>::volume_type const Basis<T, D, N>::faceWeight = faceRules.first;

template<typename T, int D, int N>
std::array<Math::Vector<T, D - 1>, Basis<T, D, N>::Nh> const Basis<T, D, N>::facePoint = faceRules.second;

template<typename T, int D, int N>
std::array<typename Basis<T, D, N>::volume_type, Basis<T, D, N>::Np> const Basis<T, D, N>::phi = []() {
	std::array<volume_type, Np> phi;
	for (hindex_type hIdx = hindex_type::begin(); hIdx != hindex_type::end(); hIdx++) {
		int const hi = hIdx.flatIndex();
		for (pindex_type pIdx = pindex_type::begin(); pIdx != pindex_type::end(); pIdx++) {
			int const pi = pIdx.flatIndex();
			phi[pi][hi] = legendreP<D>(pIdx, volumePoint[hi]);
		}
	}
	return phi;
}();

template<typename T, int D, int N>
std::array<std::array<typename Basis<T, D, N>::face_type, Basis<T, D, N>::Np>, 2 * D> const Basis<T, D, N>::xi = []() {
	std::array<std::array<face_type, Np>, 2 * D> xi;
	for (int face = 0; face < 2 * D; face++) {
		for (findex_type fIdx = findex_type::begin(); fIdx != findex_type::end(); fIdx++) {
			int const fi = fIdx.flatIndex();
			int const fDim = face >> 1;
			vector_type<D> x;
			x[fDim] = T(2 * (face & 1) - 1);
			for (int d = 0; d < fDim; d++) {
				x[d] = facePoint[fi][d];
			}
			for (int d = fDim + 1; d < D; d++) {
				x[d] = facePoint[fi][d - 1];
			}
			for (pindex_type pIdx = pindex_type::begin(); pIdx != pindex_type::end(); pIdx++) {
				int const pi = pIdx.flatIndex();
				xi[face][pi][fi] = legendreP<D>(pIdx, x);
			}
		}
	}
	return xi;
}();

template<typename T, int D, int N>
std::array<std::array<typename Basis<T, D, N>::volume_type, Basis<T, D, N>::Np>, D> const Basis<T, D, N>::dphi_dx = []() {
	std::array<std::array<volume_type, Np>, D> dphi_dx;
	for (int dim = 0; dim < D; dim++) {
		for (hindex_type hIdx = hindex_type::begin(); hIdx != hindex_type::end(); hIdx++) {
			int const hi = hIdx.flatIndex();
			for (pindex_type pIdx = pindex_type::begin(); pIdx != pindex_type::end(); pIdx++) {
				int const pi = pIdx.flatIndex();
				auto const Q = triangularToLegendre(pIdx);
				T dphi = one;
				for (int d = 0; d < D; d++) {
					if (d == dim) {
						dphi *= dLegendrePdX(Q[d], volumePoint[hi][d]);
					} else {
						dphi *= legendreP(Q[d], volumePoint[hi][d]);
					}
				}
				dphi_dx[dim][pi][hi] = dphi;
			}
		}
	}
	return dphi_dx;
}();

