#pragma once
#include <tuple>
#include <cmath>
#include <algorithm>

struct SodSolver {
	SodSolver(double rhoL, double uL, double pL, double rhoR, double uR, double pR, double gamma = 1.4) :
			gamma_ { gamma }, rhoL_ { rhoL }, uL_ { uL }, pL_ { pL }, rhoR_ { rhoR }, uR_ { uR }, pR_ { pR } {
		computeStarRegion();
	}

	std::tuple<double, double, double> operator()(double x, double t) const {
		if (t <= 0.0) {
			return x < 0.0 ? std::tie(rhoL_, uL_, pL_) : std::tie(rhoR_, uR_, pR_);
		}
		double const xi = x / t;
		double rhoLStar { }, uLStar { }, pLStar { };
		if (isLeftRarefaction_) {
			auto const aL = soundSpeed(rhoL_, pL_);
			double const aStarL = aL * std::pow(pStar_ / pL_, (gamma_ - 1.0) / (2.0 * gamma_));
			double const sHead = uL_ - aL;
			double const sTail = uStar_ - aStarL;
			if (xi <= sHead) {
				return std::tie(rhoL_, uL_, pL_);
			}
			if (xi >= sTail) {
				rhoLStar = rhoStarL_;
				uLStar = uStar_;
				pLStar = pStar_;
			} else {
				double const u = 2.0 / (gamma_ + 1.0) * (aL + 0.5 * (gamma_ - 1.0) * uL_ + xi);
				double const a = u - xi;
				double const rho = rhoL_ * std::pow(a / aL, 2.0 / (gamma_ - 1.0));
				double const p = pL_ * std::pow(a / aL, 2.0 * gamma_ / (gamma_ - 1.0));
				return std::tie(rho, u, p);
			}
		} else {
			double const sL = shockSpeed(rhoL_, pL_);
			if (xi <= sL) {
				return std::tie(rhoL_, uL_, pL_);
			}
			rhoLStar = rhoStarL_;
			uLStar = uStar_;
			pLStar = pStar_;
		}
		double rhoRStar { }, uRStar { }, pRStar { };
		if (isRightRarefaction_) {
			auto const aR = soundSpeed(rhoR_, pR_);
			double const aStarR = aR * std::pow(pStar_ / pR_, (gamma_ - 1.0) / (2.0 * gamma_));
			double const sHead = uR_ + aR;
			double const sTail = uStar_ + aStarR;
			if (xi >= sHead) {
				return std::tie(rhoR_, uR_, pR_);
			}
			if (xi <= sTail) {
				rhoRStar = rhoStarR_;
				uRStar = uStar_;
				pRStar = pStar_;
			} else {
				double const u = 2.0 / (gamma_ + 1.0) * (-aR + 0.5 * (gamma_ - 1.0) * uR_ + xi);
				double const a = -u + xi;
				double const rho = rhoR_ * std::pow(a / aR, 2.0 / (gamma_ - 1.0));
				double const p = pR_ * std::pow(a / aR, 2.0 * gamma_ / (gamma_ - 1.0));
				return std::tie(rho, u, p);
			}
		} else {
			double const sR = shockSpeed(rhoR_, pR_);
			if (xi >= sR) {
				return std::tie(rhoR_, uR_, pR_);
			}
			rhoRStar = rhoStarR_;
			uRStar = uStar_;
			pRStar = pStar_;
		}
		return std::tie(rhoStarR_, uStar_, pStar_);
	}

private:
	double soundSpeed(double rho, double p) const {
		return std::sqrt(gamma_ * p / rho);
	}
	double shockSpeed(double rho, double p) const {
		double const factor = (pStar_ / p - 1.0);
		double const a = soundSpeed(rho, p);
		return (p == pL_ ? uL_ : uR_) - a * std::sqrt(0.5 * (gamma_ + 1.0) / gamma_ * factor + 1.0);
	}
	void computeStarRegion() {
		auto const aL = soundSpeed(rhoL_, pL_);
		auto const aR = soundSpeed(rhoR_, pR_);
		double pGuess = std::max(0.0, 0.5 * (pL_ + pR_) - 0.125 * (uR_ - uL_) * (rhoL_ + rhoR_) * (aL + aR));
		double p = pGuess;
		for (int i = 0; i < 20; ++i) {
			double fL { }, dfdpL { };
			fSide(p, rhoL_, pL_, aL, fL, dfdpL);
			double fR { }, dfdpR { };
			fSide(p, rhoR_, pR_, aR, fR, dfdpR);
			double f = fL + fR + uR_ - uL_;
			double dfdp = dfdpL + dfdpR;
			double pNew = p - f / dfdp;
			if (std::abs(pNew - p) / (pNew + 1.0e-14) < 1.0e-6) {
				p = pNew;
				break;
			}
			p = pNew;
		}
		pStar_ = p;
		uStar_ = 0.5 * (uL_ + uR_) + 0.5 * (fSideValue(p, rhoR_, pR_, aR) - fSideValue(p, rhoL_, pL_, aL));
		if (pStar_ > pL_) {                 // shock
			double const numerator = pStar_ / pL_ + (gamma_ - 1.0) / (gamma_ + 1.0);
			double const denominator = (gamma_ - 1.0) / (gamma_ + 1.0) * pStar_ / pL_ + 1.0;
			rhoStarL_ = rhoL_ * numerator / denominator;
			isLeftRarefaction_ = false;
		} else {                            // rarefaction
			rhoStarL_ = rhoL_ * std::pow(pStar_ / pL_, 1.0 / gamma_);
			isLeftRarefaction_ = true;
		}
		if (pStar_ > pR_) {
			double const numerator = pStar_ / pR_ + (gamma_ - 1.0) / (gamma_ + 1.0);
			double const denominator = (gamma_ - 1.0) / (gamma_ + 1.0) * pStar_ / pR_ + 1.0;
			rhoStarR_ = rhoR_ * numerator / denominator;
			isRightRarefaction_ = false;
		} else {
			rhoStarR_ = rhoR_ * std::pow(pStar_ / pR_, 1.0 / gamma_);
			isRightRarefaction_ = true;
		}
	}
	void fSide(double p, double rho, double p0, double a0, double &f, double &dfdp) const {
		if (p > p0) {
			double const A = 2.0 / ((gamma_ + 1.0) * rho);
			double const B = (gamma_ - 1.0) / (gamma_ + 1.0) * p0;
			double const sqrtTerm = std::sqrt(A / (p + B));
			f = (p - p0) * sqrtTerm;
			dfdp = sqrtTerm * (1.0 - 0.5 * (p - p0) / (p + B));
		} else {
			double const exponent = (gamma_ - 1.0) / (2.0 * gamma_);
			f = 2.0 * a0 / (gamma_ - 1.0) * (std::pow(p / p0, exponent) - 1.0);
			dfdp = (1.0 / (rho * a0)) * std::pow(p / p0, -exponent);
		}
	}
	double fSideValue(double p, double rho, double p0, double a0) const {
		double f { }, dummy;
		fSide(p, rho, p0, a0, f, dummy);
		return f;
	}
	double gamma_;
	double rhoL_, uL_, pL_;
	double rhoR_, uR_, pR_;
	double pStar_ { }, uStar_ { };
	double rhoStarL_ { }, rhoStarR_ { };
	bool isLeftRarefaction_ { };
	bool isRightRarefaction_ { };
};
