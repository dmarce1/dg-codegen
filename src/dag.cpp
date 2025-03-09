/*
 * dag.cpp
 *
 *  Created on: Mar 7, 2025
 *      Author: dmarce1
 */

#include <cctype>
#include <iostream>
#include <map>
#include <memory>
#include <numeric>
#include <regex>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#define TOKEN_PATTERN R"((\s+)|(_[A-Za-z()]+)|([A-Za-z0-9]+)|([()+\-*/=]))"

struct DagNode {
	struct term_t {
		double constant;
		std::map<DagNode, int> factors;
		friend bool operator==(term_t const &A, term_t const &B);
		friend bool operator<(term_t const &A, term_t const &B);
	};
	struct state_t {
		std::string varName;
		std::vector<term_t> terms;
		friend bool operator==(state_t const &A, state_t const &B);
		friend bool operator<(state_t const &A, state_t const &B);
	};
	friend bool operator==(DagNode const &A, DagNode const &B);
	friend bool operator<(DagNode const &A, DagNode const &B);
	void addTerm(std::vector<DagNode> const &factors, double constant = 1.0);
	DagNode();
	DagNode(DagNode const &other);
	DagNode(DagNode &&other);
	DagNode(std::string const&);
	DagNode& operator=(std::string const&);
	DagNode& operator=(DagNode const &other);
	DagNode& operator=(DagNode &&other);
	void setTag(std::string const &t);
	std::string getTag() const;
	std::string getVariableName() const;
	int getTermCount() const;
	friend std::ostream& operator<<(std::ostream &outStream, DagNode const&);
	static std::vector<DagNode> getInputNodes();
private:
	std::shared_ptr<state_t> state;
	std::string tag;
	static std::unordered_map<std::string, std::shared_ptr<state_t>> varNameDatabase;
	DagNode(std::shared_ptr<state_t> st);
};

std::unordered_map<std::string, std::shared_ptr<DagNode::state_t>> DagNode::varNameDatabase;

std::ostream& operator<<(std::ostream &outStream, DagNode const &dag) {
	std::string code;
	auto const tag = dag.getTag();
	if (tag.size()) {
		outStream << tag << " = ";
	}
	for (unsigned ti = 0; ti < dag.state->terms.size(); ti++) {
		auto const &term = dag.state->terms[ti];
		if (term.constant != 1.0) {
			outStream << " T(" << std::to_string(term.constant) << ")";
		}
		auto fit = term.factors.begin();
		for (int fi = 0; fi < term.factors.size(); fi++) {
			auto const &factor = *fit;
			for (int m = 0; m < factor.second; m++) {
				outStream << factor.first.getVariableName();
				outStream << ((m < factor.second - 1) ? " * " : "");
			}
			fit++;
			outStream << ((fi < term.factors.size() - 1) ? " * " : "");
		}
		outStream << ((ti == dag.state->terms.size() - 1) ? ";\n" : " + ");
	}
	return outStream;
}

int DagNode::getTermCount() const {
	return state->terms.size();
}

void DagNode::setTag(std::string const &t) {
	tag = t;
}

std::string DagNode::getTag() const {
	return tag;
}

std::string DagNode::getVariableName() const {
	return state->varName;
}

std::vector<DagNode> DagNode::getInputNodes() {
	std::vector<DagNode> nodes;
	for (auto it = varNameDatabase.begin(); it != varNameDatabase.end(); it++) {
		nodes.push_back(it->second);
	}
	return nodes;
}

DagNode::DagNode(std::shared_ptr<state_t> st) {
	state = st;
}

bool operator==(DagNode::term_t const &A, DagNode::term_t const &B) {
	if (A.constant == B.constant) {
		if (A.factors == B.factors) {
			return true;
		}
	}
	return false;
}

bool operator<(DagNode::term_t const &A, DagNode::term_t const &B) {
	if (A.constant < B.constant) {
		return true;
	} else if (A.constant > B.constant) {
		return false;
	} else {
		auto itA = A.factors.begin();
		auto itB = B.factors.begin();
		while (1) {
			if ((itA == A.factors.end()) && (itB == B.factors.end())) {
				return false;
			} else if (itA == A.factors.end()) {
				return true;
			} else if (itB == B.factors.end()) {
				return false;
			} else if (*itA < *itB) {
				return true;
			} else if (*itB < *itA) {
				return false;
			} else {
				itA++;
				itB++;
			}
		}
	}
}

bool operator==(DagNode::state_t const &A, DagNode::state_t const &B) {
	return (A.terms == B.terms) && (A.varName == B.varName);
}

bool operator<(DagNode::state_t const &A, DagNode::state_t const &B) {
	if (A.terms < B.terms) {
		return true;
	} else if (A.terms < B.terms) {
		return false;
	} else {
		return A.varName < B.varName;
	}
}

bool operator==(DagNode const &A, DagNode const &B) {
	return *A.state == *B.state;
}

bool operator<(DagNode const &A, DagNode const &B) {
	return *A.state < *B.state;
}

void DagNode::addTerm(std::vector<DagNode> const &factors, double constant) {
	term_t term;
	term.constant = constant;
	for (auto f : factors) {
		if (term.factors.find(f) == term.factors.end()) {
			term.factors[f] = 0;
		}
		term.factors[f]++;
	}
	auto del = std::move(state->varName);
	state->terms.push_back(std::move(term));
}

DagNode::DagNode() :
		state(std::make_shared<state_t>()) {
}

DagNode::DagNode(DagNode const &other) {
	state = other.state;
	tag = other.tag;
}

DagNode::DagNode(std::string const &varName) :
		state(std::make_shared<state_t>()) {
	*this = varName;
}

DagNode& DagNode::operator=(std::string const &varName) {
	if (varNameDatabase.find(varName) == varNameDatabase.end()) {
		state->varName = varName;
		auto del = std::move(state->terms);
		varNameDatabase[varName] = state;
	}
	state = varNameDatabase[varName];
	return *this;
}

DagNode::DagNode(DagNode &&other) {
	*this = std::move(other);
}

DagNode& DagNode::operator=(DagNode const &other) {
	state = other.state;
	tag = other.tag;
	return *this;
}

DagNode& DagNode::operator=(DagNode &&other) {
	state = other.state;
	tag = other.tag;
	other.tag = "";
	other.state = nullptr;
	return *this;
}

struct token_type {
	static constexpr int variable = 0;
	struct index {
		static constexpr int covariant = 2;
		static constexpr int contravariant = 3;
	};
	struct operation {
		static constexpr int add = 4;
		static constexpr int sub = 5;
		static constexpr int mul = 6;
		static constexpr int div = 7;
		static constexpr int equals = 8;
		static constexpr int open = 9;
		static constexpr int close = 10;
	};
	token_type() :
			value(0) {
	}
	token_type(int i) :
			value(i) {
	}
	operator int&() {
		return value;
	}
	operator int() const {
		return value;
	}
	bool isOperation() const {
		if (value >= operation::add) {
			if (value <= operation::close) {
				return true;
			}
		}
		return false;
	}
	bool isAdditive() const {
		if (value == operation::add) {
			return true;
		}
		if (value == operation::sub) {
			return true;
		}
		return false;
	}
	bool isIndex() const {
		if (value == index::covariant) {
			return true;
		}
		if (value == index::contravariant) {
			return true;
		}
		return false;
	}
private:
	int value;
};

struct Token {
	void setVariableName(std::string const &str) {
		name = str;
		dummyValue = value = 0;
		type = token_type::variable;
		symmetryGroup = 0;
	}
	void setCovariantIndex(char c) {
		value = c;
		type = token_type::index::covariant;
	}
	void setContravariantIndex(char c) {
		value = std::tolower(c);
		type = token_type::index::contravariant;
	}
	void setOperation(char c) {
		value = c;
		switch (c) {
		case '(':
			type = token_type::operation::open;
			break;
		case ')':
			type = token_type::operation::close;
			break;
		case '+':
			type = token_type::operation::add;
			break;
		case '-':
			type = token_type::operation::sub;
			break;
		case '*':
			type = token_type::operation::mul;
			break;
		case '/':
			type = token_type::operation::div;
			break;
		case '=':
			type = token_type::operation::equals;
			break;
		default:
			throw std::invalid_argument("Unrecognized character\n");
		}
	}
	friend std::string toString(Token const &token) {
		std::string str;
		auto const type = token.type;
		if (type == token_type::variable) {
			str += "Variable: " + token.name + "\n";
		} else if (token.isIndex()) {
			std::string indexType;
			if (token.isDummy()) {
				indexType = "dummy";
			} else {
				indexType = "free";
			}
			if (type == token_type::index::covariant) {
				str += "Covariant " + indexType + " index: ";
				str.push_back(token.value);
			} else if (type == token_type::index::contravariant) {
				str += "Contravariant " + indexType + " index: ";
				str.push_back(token.value);
			} else {
				throw std::invalid_argument("Unrecognized token\n");
			}
			if (token.isSymmetricIndex()) {
				str += " symmetry: group # "
						+ std::to_string(token.getSymmetryGroup());
			} else {
				str += " symmetry:  N/A "
						+ std::to_string(token.getSymmetryGroup());
			}
			str += "\n";
		} else if (type.isOperation()) {
			str += "Operation: ";
			str.push_back(token.value);
			str += "\n";
		} else {
			throw std::invalid_argument("Unrecognized token\n");
		}
		return str;
	}
	void setDummy(int i) {
		dummyValue = abs(i) + 1;
	}
	bool isAdditive() const {
		return type.isAdditive();
	}
	bool isEquals() const {
		return type == token_type::operation::equals;
	}
	bool isIndex() const {
		return type.isIndex();
	}
	bool isFree() const {
		return type.isIndex() && !isDummy();
	}
	bool isDummy() const {
		return dummyValue;
	}
	bool isVar() const {
		return type == token_type::variable;
	}
	char getChar() const {
		return value;
	}
	std::string getName() const {
		return name;
	}
	void setSymmetryGroup(int g) {
		symmetryGroup = g;
	}
	bool isSymmetricIndex() const {
		return symmetryGroup > 0;
	}
	int getSymmetryGroup() const {
		return symmetryGroup;
	}
private:
	token_type type;
	char value;
	int dummyValue;
	int symmetryGroup;
	std::string name;
};

auto processString(std::string const &input) {
	std::vector<Token> theseTokens;
	std::vector<Token> tokens;
	Token token;

	std::regex tokenRegex(TOKEN_PATTERN);

	auto const flush = [&tokens, &theseTokens]() {
		int nextDummy = 0;
		for (unsigned i = 0; i < theseTokens.size(); i++) {
			auto &token1 = theseTokens[i];
			for (unsigned m = i + 1; m < theseTokens.size(); m++) {
				auto &token2 = theseTokens[m];
				if (token1.isIndex() && token2.isIndex()) {
					if (token1.getChar() == token2.getChar()) {
						token1.setDummy(nextDummy);
						token2.setDummy(nextDummy);
						nextDummy++;
						break;
					}
				}
			}
		}
		tokens.insert(tokens.end(), theseTokens.begin(), theseTokens.end());
		theseTokens.clear();
	};
	for (auto it = std::sregex_iterator(input.begin(), input.end(), tokenRegex);
			it != std::sregex_iterator(); ++it) {
		std::smatch match = *it;
		if (match[1].matched) {
			continue;
		} else if (match[2].matched) {
			static int nextSymmetryGroup = 1;
			std::string indices = match[2].str();
			bool isSymmetryGroup = false;
			int groupNumber;
			for (size_t i = 1; i < indices.size(); ++i) {
				char letter = indices[i];
				Token token;
				if (letter == '(') {
					isSymmetryGroup = true;
					groupNumber = nextSymmetryGroup++;
				} else if (letter == ')') {
					isSymmetryGroup = false;
				} else {
					if (std::islower(letter)) {
						token.setCovariantIndex(letter);
					} else {
						token.setContravariantIndex(letter);
					}
					if (isSymmetryGroup) {
						token.setSymmetryGroup(groupNumber);
					}
					theseTokens.push_back(token);
				}
			}
		} else if (match[3].matched) {
			std::string varName = match[3].str();
			Token token;
			token.setVariableName(varName);
			theseTokens.push_back(token);
		} else if (match[4].matched) {
			char op = match[4].str()[0];
			Token token;
			token.setOperation(op);
			if (token.isAdditive() || token.isEquals()) {
				flush();
			}
			theseTokens.push_back(token);
		}
	}
	flush();
	return tokens;
}

template<>
struct std::hash<std::vector<int>> {
	size_t operator()(std::vector<int> I) const {
		size_t key = 1;
		for (unsigned i = 0; i < I.size(); i++) {
			key = (1664525 * key + I[i]) & 0xFFFFFFFF;
		}
		return key;
	}
};

std::vector<int> reduceIndices(std::vector<int> indices,
		std::vector<int> const &symmetries) {
	std::unordered_set<int> done;
	for (unsigned n = 0; n < symmetries.size(); n++) {
		if (symmetries[n] && (done.find(n) == done.end())) {
			std::vector<int> values;
			for (unsigned m = n; m < symmetries.size(); m++) {
				if (symmetries[m] == symmetries[n]) {
					values.push_back(indices[m]);
					done.insert(m);
				}
			}
			std::sort(values.begin(), values.end());
			for (unsigned m = n; m < symmetries.size(); m++) {
				if (symmetries[m] == symmetries[n]) {
					indices[m] = values.back();
					values.pop_back();
				}
			}
		}
	}
	return indices;
}

int computeIndex(std::vector<int> indices, std::vector<int> const &symmetries) {
	static std::unordered_map<std::vector<int>,
			std::unordered_map<std::vector<int>, int>> indexMaps;
	if (indices.size() == 0) {
		return 0;
	} else {
		if (indexMaps.find(symmetries) == indexMaps.end()) {
			auto &thisMap = indexMaps[symmetries];
			std::vector<int> I(indices.size(), 0);
			int nextIndex = 0;
			bool flag = false;
			while (!flag) {
				auto reducedI = reduceIndices(I, symmetries);
				if (thisMap.find(reducedI) == thisMap.end()) {
					thisMap[reducedI] = nextIndex++;
				}
				thisMap[I] = thisMap[reducedI];
				unsigned dim = 0;
				while (++I[dim] == NDIM) {
					I[dim++] = 0;
					if (dim == I.size()) {
						flag = true;
						break;
					}
				}
			}
		}
		return indexMaps[symmetries][indices];
	}
}

std::string genVarName(std::string const &base, int index) {
	std::string name = base;
	name += "[";
	name += std::to_string(index);
	name += "]";
	return name;
}

std::pair<std::vector<DagNode>, std::vector<DagNode>> processTokenGroup(
		std::vector<Token> const &tokens) {
	std::pair<std::vector<DagNode>, std::vector<DagNode>> rc;
	auto &inputVariables = rc.first;
	auto &outputVariables = rc.second;
	std::string lhsName;
	std::map<char, int> lhsGroups;
	std::vector<char> lhsIndexNames;
	if (!tokens[0].isVar()) {
		throw std::invalid_argument(
				"First token should be a variable name\currentTokenIndex");
	}
	lhsName = tokens[0].getName();
	unsigned tokenIndex = 1;
	while (tokens[tokenIndex].isIndex()) {
		char const c = tokens[tokenIndex].getChar();
		lhsIndexNames.push_back(c);
		lhsGroups[c] = tokens[tokenIndex].getSymmetryGroup();
		tokenIndex++;
	}
	if (!tokens[tokenIndex].isEquals()) {
		throw std::invalid_argument("Missing equals sign\currentTokenIndex");
	}
	int const Rank = lhsIndexNames.size();
	std::map<char, int> lhsIndices;
	for (unsigned i = 0; i < lhsIndexNames.size(); i++) {
		lhsIndices[lhsIndexNames[i]] = 0;
	}
	tokenIndex++;
	bool outerFlag = false;
	while (!outerFlag) {
		int currentTokenIndex = tokenIndex;
		std::vector<int> tmp1;
		std::vector<int> tmp2;
		for (auto i : lhsIndices) {
			tmp1.push_back(i.second);
			tmp2.push_back(lhsGroups[i.first]);
		}
		int const lhsIndex = computeIndex(std::move(tmp1), std::move(tmp2));
		auto const name = genVarName(lhsName.c_str(), lhsIndex);
		DagNode thisDagNode;
		thisDagNode.setTag(name);
		while (currentTokenIndex < tokens.size()) {
			std::vector<std::string> rhsNames;
			std::vector<std::map<char, int>> rhsGroups;
			std::vector<std::vector<char>> rhsIndexNames;
			std::vector<char> dummyIndexNames;
			while (currentTokenIndex < tokens.size()) {
				if (tokens[currentTokenIndex].isVar()) {
					rhsNames.push_back(tokens[currentTokenIndex].getName());
					currentTokenIndex++;
					std::vector<char> theseIndices;
					std::map<char, int> theseGroups;
					while ((currentTokenIndex < tokens.size())
							&& tokens[currentTokenIndex].isIndex()) {
						char const c = tokens[currentTokenIndex].getChar();
						theseIndices.push_back(c);
						theseGroups[c] =
								tokens[currentTokenIndex].getSymmetryGroup();
						if (tokens[currentTokenIndex].isDummy()) {
							dummyIndexNames.push_back(c);
						}
						currentTokenIndex++;
					}
					rhsIndexNames.push_back(std::move(theseIndices));
					rhsGroups.push_back(std::move(theseGroups));
				} else if (tokens[currentTokenIndex].isAdditive()) {
					currentTokenIndex++;
					break;
				} else {
					currentTokenIndex++;
				}
			}
			std::sort(dummyIndexNames.begin(), dummyIndexNames.end());
			for (unsigned i = 1; i < dummyIndexNames.size(); i += 2) {
				if (i < dummyIndexNames.size() - 1) {
					dummyIndexNames[i] = dummyIndexNames.back();
					dummyIndexNames.pop_back();
					dummyIndexNames.pop_back();
				} else {
					dummyIndexNames.pop_back();
				}
			}
			std::unordered_set<char> dummySet(dummyIndexNames.begin(),
					dummyIndexNames.end());

			std::map<char, int> dummyIndices;
			for (unsigned i = 0; i < dummyIndexNames.size(); i++) {
				dummyIndices[dummyIndexNames[i]] = 0;
			}
			bool innerFlag = false;
			while (!innerFlag) {
				std::vector<DagNode> factors;
				for (unsigned n = 0; n < rhsNames.size(); n++) {
					std::vector<int> theseIndices;
					std::vector<int> theseSymmetries;
					for (char index : rhsIndexNames[n]) {
						if (dummySet.find(index) == dummySet.end()) {
							theseIndices.push_back(lhsIndices[index]);
						} else {
							theseIndices.push_back(dummyIndices[index]);
						}
						theseSymmetries.push_back(rhsGroups[n][index]);
					}
					int const rhsIndex = computeIndex(theseIndices,
							theseSymmetries);
					auto const name = genVarName(rhsNames[n].c_str(), rhsIndex);
					factors.push_back(DagNode(name));
				}
				thisDagNode.addTerm(factors);
				if (dummyIndices.size() == 0) {
					innerFlag = true;
					break;
				} else {
					unsigned s = 0;
					while (++dummyIndices[dummyIndexNames[s]] == NDIM) {
						dummyIndices[dummyIndexNames[s++]] = 0;
						if (s == dummyIndices.size()) {
							innerFlag = true;
							break;
						}
					}
				}
			}
		}
		outputVariables.push_back(thisDagNode);
		if (lhsIndices.size() == 0) {
			outerFlag = true;
			break;
		} else {
			int r = 0;
			while (++lhsIndices[lhsIndexNames[r]] == NDIM) {
				lhsIndices[lhsIndexNames[r++]] = 0;
				if (r == Rank) {
					outerFlag = true;
					break;
				}
			}
		}
	}
	inputVariables = DagNode::getInputNodes();
	return rc;
}

void testStrings() {
	int N = 3;
	int index = 0;
	for (int i = 0; i < N; i++) {
		for (int k = 0; k < N; k++) {
			for (int j = 0; j <= i; j++) {
				int i2 = 0;
				i2 = j + (i + 1) * k + (i + 1) * N * i / 2;
				//(i + 1) * (N*i/2 + k) + j
				printf("%i %i | %i %i %i\n", index++, i2, i, k, j);
			}
			printf("\n");
		}
		printf("\n");
	}
	const char *A = "D1_ = g_(IJ) * A_(ij) + C_";
	auto tokens = processString(A);
	for (auto t : tokens) {
		std::cout << toString(t);
	}
	auto rc = processTokenGroup(tokens);
	auto outputs = std::move(rc.second);
	for (int i = 0; i < outputs.size(); i++) {
		std::cout << outputs[i];
	}

}
