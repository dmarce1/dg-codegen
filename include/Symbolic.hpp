#pragma once

#include <symengine/basic.h>
#include <symengine/symbol.h>
#include <symengine/constants.h>
#include <symengine/functions.h>
#include <symengine/rational.h>
#include <symengine/add.h>
#include <symengine/mul.h>
#include <symengine/pow.h>
#include <symengine/real_double.h>
#include <symengine/visitor.h>
#include <symengine/derivative.h>
#include <symengine/printers/strprinter.h>
#include <symengine/symengine_rcp.h>
#include <iomanip>

using SymbolicExpression = SymEngine::RCP<const SymEngine::Basic>;
using SymbolicSymbol = SymEngine::RCP<const SymEngine::Symbol>;

inline auto operator+(SymbolicExpression const &a) {
	using namespace SymEngine;
	return a;
}

inline auto operator-(SymbolicExpression const &a) {
	using namespace SymEngine;
	return neg(a);
}

inline auto operator+(SymbolicExpression const &a, SymbolicExpression const &b) {
	using namespace SymEngine;
	return add(a, b);
}

inline auto operator-(SymbolicExpression const &a, SymbolicExpression const &b) {
	using namespace SymEngine;
	return sub(a, b);
}

inline auto operator*(SymbolicExpression const &a, SymbolicExpression const &b) {
	using namespace SymEngine;
	return mul(a, b);
}

inline auto operator/(SymbolicExpression const &a, SymbolicExpression const &b) {
	using namespace SymEngine;
	return div(a, b);
}

inline SymbolicExpression sqrt(SymbolicExpression &a) {
	using namespace SymEngine;
	return sqrt(a);
}

inline SymbolicExpression pow(SymbolicExpression &a, SymbolicExpression &b) {
	return SymEngine::pow(a, b);
}

inline SymbolicExpression& operator+=(SymbolicExpression &a, SymbolicExpression const &b) {
	a = a + b;
	return a;
}

inline SymbolicExpression& operator-=(SymbolicExpression &a, SymbolicExpression const &b) {
	a = a - b;
	return a;
}

inline SymbolicExpression& operator*=(SymbolicExpression &a, SymbolicExpression const &b) {
	a = a * b;
	return a;
}

inline SymbolicExpression& operator/=(SymbolicExpression &a, SymbolicExpression const &b) {
	a = a / b;
	return a;
}

inline auto derivative(SymbolicExpression &f, SymbolicSymbol const &x) {
	using namespace SymEngine;
	return diff(f, x);
}

static inline std::string toCxxCode(SymEngine::RCP<const SymEngine::Basic> expr) {
	using namespace SymEngine;
	if (is_a<Integer>(*expr)) {
		std::ostringstream oss;
		oss << down_cast<const Integer&>(*expr).as_integer_class();
		return oss.str();
	} else if (is_a<RealDouble>(*expr)) {
		const double eps = 1e-10;
		double val = down_cast<const RealDouble&>(*expr).as_double();
		if (std::abs(val) < eps) {
			val = 0.0;
		}
		std::ostringstream oss;
		oss << std::setprecision(17) << std::fixed;
		oss << "T(" << val << ")";
		return oss.str();
	} else if (is_a<Rational>(*expr)) {
		const auto &rat = down_cast<const Rational&>(*expr);
		std::ostringstream oss;
		oss << "T(";
		oss << toCxxCode(rat.get_num()) << ".";
		oss << "/";
		oss << toCxxCode(rat.get_den()) << ".";
		oss << ")";
		return oss.str();
	} else if (is_a<Symbol>(*expr)) {
		return down_cast<const Symbol&>(*expr).get_name();
	} else if (is_a<Add>(*expr)) {
		auto args = down_cast<const Add&>(*expr).get_args();
		std::string result;
		bool first = true;
		for (const auto &arg : args) {
			std::string term;
			bool isNegative = false;
			if (is_a<Mul>(*arg)) {
				auto mul_args = down_cast<const Mul&>(*arg).get_args();
				if (mul_args.size() == 2 && is_a<Integer>(*mul_args[0])) {
					const auto &val = down_cast<const Integer&>(*mul_args[0]).as_integer_class();
					if (val == -1) {
						term = toCxxCode(mul_args[1]);
						isNegative = true;
					}
				}
			}
			if (term.empty()) {
				term = toCxxCode(arg);
				if (term.size() >= 1 && term[0] == '-') {
					isNegative = true;
					term = term.substr(1);
				}
			}
			if (first) {
				result += (isNegative ? "-" : "") + term;
				first = false;
			} else {
				result += (isNegative ? " - " : " + ") + term;
			}
		}
		return result;
	} else if (is_a<Mul>(*expr)) {
		auto args = down_cast<const Mul&>(*expr).get_args();
		std::string result;
		bool negative = false;
		bool first = true;
		for (const auto &arg : args) {
			if (is_a<Integer>(*arg)) {
				const auto &val = down_cast<const Integer&>(*arg).as_integer_class();
				if (val == 1) {
					continue;
				}
				if (val == -1) {
					negative = true;
					continue;
				}
			}
			if (!first) {
				result += " * ";
			}
			result += toCxxCode(arg);
			first = false;
		}
		if (result.empty()) {
			return negative ? "-1" : "1";
		}
		if (negative) {
			return "-" + result;
		} else {
			return result;
		}
	} else if (is_a<Pow>(*expr)) {
		auto base = down_cast<const Pow&>(*expr).get_base();
		auto exp = down_cast<const Pow&>(*expr).get_exp();
		if (is_a<Integer>(*exp)) {
			const auto &ival = down_cast<const Integer&>(*exp).as_integer_class();
			if (ival == 2) {
				return "sqr(" + toCxxCode(base) + ")";
			} else if (ival == 3) {
				return "cube(" + toCxxCode(base) + ")";
			} else if (ival == -1) {
				return "(1. / (" + toCxxCode(base) + "))";
			} else if (ival == -2) {
				return "(1. / sqr(" + toCxxCode(base) + "))";
			} else if (ival == -3) {
				return "(1. / cube(" + toCxxCode(base) + "))";
			} else {
				return "ipow(" + toCxxCode(base) + ", " + std::to_string(ival.get_si()) + ")";
			}
		}
		if (is_a<Rational>(*exp)) {
			if (rational(+1, 2)->compare(*exp) == 0) {
				return "sqrt(" + toCxxCode(base) + ")";
			}
			if (rational(+1, 3)->compare(*exp) == 0) {
				return "cbrt(" + toCxxCode(base) + ")";
			}
			if (rational(-1, 2)->compare(*exp) == 0) {
				return "T(1) / sqrt(" + toCxxCode(base) + ")";
			}
			if (rational(-1, 3)->compare(*exp) == 0) {
				return "T(1) / cbrt(" + toCxxCode(base) + ")";
			}
		}
		return "pow(" + toCxxCode(base) + ", " + toCxxCode(exp) + ")";
	}
	return "/* unsupported */";
}

using namespace SymEngine;

inline bool isOneOverPow(RCP<const Basic> expr, RCP<const Basic> &base, int &exponent) {
	if (!is_a<Pow>(*expr))
		return false;

	const auto &p = down_cast<const Pow&>(*expr);
	if (!is_a<Integer>(*p.get_exp()))
		return false;

	int exp = down_cast<const Integer&>(*p.get_exp()).as_int();
	if (exp >= 0)
		return false;

	base = p.get_base();  // no need to cast
	exponent = -exp;
	return true;
}

static inline void factorPowersFromSubstitutions(vec_pair &subs) {
	// Map known: x0 = 1 / rho
	{
		std::map<RCP<const Basic>, RCP<const Basic>, RCPBasicKeyLess> invMap;
		std::map<RCP<const Basic>, RCP<const Symbol>, RCPBasicKeyLess> exprToVar;

		// First pass: find things like x0 = 1 / y
		for (auto const& [sym, expr] : subs) {
			if (is_a<Pow>(*expr)) {
				const auto &p = down_cast<const Pow&>(*expr);
				if (is_a<Integer>(*p.get_exp()) && down_cast<const Integer&>(*p.get_exp()).as_int() == -1) {
					invMap[p.get_base()] = sym;
				}
			}
		}

		// Second pass: rewrite things like x1 = 1 / y^n using x0 = 1 / y
		for (auto& [sym, expr] : subs) {
			RCP<const Basic> base;
			int exponent;
			if (isOneOverPow(expr, base, exponent)) {
				auto it = invMap.find(base);
				if (it != invMap.end()) {
					auto xj = it->second;
					RCP<const Basic> newExpr = xj;
					for (int i = 1; i < exponent; ++i) {
						newExpr = mul(newExpr, xj);
					}
					if (!has_symbol(*newExpr, *sym)) {
						expr = newExpr;
					}
				}
			}
		}
	}
	{
		vec_pair invFirst;
		vec_pair everythingElse;

		for (auto const& [lhs, rhs] : subs) {
			if (is_a<Pow>(*rhs)) {
				const auto &p = down_cast<const Pow&>(*rhs);
				if (is_a<Integer>(*p.get_exp()) && down_cast<const Integer&>(*p.get_exp()).as_int() == -1 && is_a<Symbol>(*p.get_base())) {
					invFirst.emplace_back(lhs, rhs);
					continue;
				}
			}
			everythingElse.emplace_back(lhs, rhs);
		}

		// Concatenate: inverses first
		subs.clear();
		subs.insert(subs.end(), invFirst.begin(), invFirst.end());
		subs.insert(subs.end(), everythingElse.begin(), everythingElse.end());
	}
}

static SymEngine::RCP<const SymEngine::Basic> zeroSmallConstants(SymEngine::RCP<const SymEngine::Basic> expr, double eps = 1e-10) {
	if (is_a<RealDouble>(*expr)) {
		double val = down_cast<const RealDouble&>(*expr).as_double();
		if (std::abs(val) < eps) {
			return zero;
		} else {
			return expr;
		}
	}

	vec_basic newArgs;
	bool changed = false;
	for (const auto &arg : expr->get_args()) {
		RCP<const Basic> newArg = zeroSmallConstants(arg, eps);
		changed |= (newArg != arg);
		newArgs.push_back(newArg);
	}

	if (!changed) {
		return expr;
	}

	if (is_a<Add>(*expr)) {
		RCP<const Basic> result = zero;
		for (const auto &term : newArgs) {
			result = add(result, term);
		}
		return result;
	} else if (is_a<Mul>(*expr)) {
		RCP<const Basic> result = one;
		for (const auto &factor : newArgs) {
			result = mul(result, factor);
		}
		return result;
	} else if (is_a<Pow>(*expr)) {
		return pow(newArgs[0], newArgs[1]);
	} else if (is_a<Sin>(*expr)) {
		return sin(newArgs[0]);
	} else if (is_a<Cos>(*expr)) {
		return cos(newArgs[0]);
	} else if (is_a<Tan>(*expr)) {
		return tan(newArgs[0]);
	} else if (is_a<Log>(*expr)) {
		return log(newArgs[0]);
	} else if (is_a<Abs>(*expr)) {
		return abs(newArgs[0]);
	}

	// Default: unhandled expression
	return expr;
}
