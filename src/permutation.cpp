#include <iostream>
#include <vector>
#include <experimental/simd>
#include <numeric>
#include <cassert>
#include <algorithm>
#include <set>

using Real = long double;

using Matrix = std::vector<std::vector<Real>>;

static std::vector<int> buildPermutation(Matrix const& A, int simdWidth) {
    int ncols = A[0].size();
    std::set<int> used;
    std::vector<int> perm;

    while ((int)perm.size() < ncols) {
        for (int i = 0; i < ncols; ++i) {
            if (used.count(i)) continue;
            std::vector<int> group;
            for (int j = i; j < ncols && (int)group.size() < simdWidth; ++j) {
                if (!used.count(j)) group.push_back(j);
            }
            while ((int)group.size() < simdWidth)
                group.push_back(group.back()); // pad if needed
            for (int idx : group) used.insert(idx);
            perm.insert(perm.end(), group.begin(), group.end());
            break;
        }
    }
    return perm;
}

static Matrix buildTransformedMatrix(Matrix const& A, std::vector<int> const& perm, int simdWidth) {
    int rows = A.size();
    int cols = perm.size();
    Matrix B(rows, std::vector<Real>(cols, 0.0));
    for (int i = 0; i < rows; ++i) {
        for (int j = 0; j < (int)perm.size(); ++j) {
            B[i][j] = A[i][perm[j]];
        }
    }
    return B;
}


Matrix buildSimdPermutation(Matrix const& A, int simdWidth) {
	return buildTransformedMatrix(A, buildPermutation(A, simdWidth), simdWidth);
}


