/*
 * SparseMatrix.hpp
 *
 *  Created on: Mar 10, 2025
 *      Author: dmarce1
 */

#ifndef INCLUDE_SPARSEMATRIX_HPP_
#define INCLUDE_SPARSEMATRIX_HPP_

#include <map>
#include <memory>
#include <iostream>
#include <vector>
#include <hpx/hpx.hpp>

template<typename T>
struct SparseMatrix {
	SparseMatrix(int rows_, int cols_) :
			nRow(rows_), nCol(cols_), elements() {
	}
	T operator[](int r, int c) const {
		std::pair<int, int> key { r, c };
		auto it = std::find_first_of(elements.begin(), elements.end(), [r, c](std::tuple<int, int, T> tup) {
			return (std::get < 0 > (tup) == r) && (std::get < 1 > (tup) == c);
		});
		if (it == elements.end()) {
			return T(0);
		} else {
			return *it;
		}
	}
	void assign(int row, int col, T content) {
		elements.push_back(std::tuple<int, int, T>(row, col, content));
	}
	SparseMatrix subMatrix(int r, int c) const {
		SparseMatrix sub(nRow - 1, nCol - 1);
		for (auto it = elements.begin(); it != elements.end(); it++) {
			int const row = std::get < 0 > (*it);
			if (std::get < 0 > (*it) != r && std::get < 1 > (*it) != c) {
				int const thisRow = row - int(row > r);
				int const thisCol = std::get < 1 > (*it) - int(std::get < 1 > (*it) > c);
				sub.assign(thisRow, thisCol, std::get < 2 > (*it));
			}
		}
		return sub;
	}
	T determinant(int level = 0) const {
		if (nRow == 1) {
			if (elements.empty()) {
				return T(0);
			} else {
				return (std::get < 2 > (elements[0]));
			}
		}
		std::vector < hpx::future < T >> futs;
		T sum = T(0);
		for (auto it = elements.begin(); it != elements.end(); it++) {
			if (std::get < 0 > (*it) == 0) {
				if (level > 3) {
					T const sgn = T(1 - 2 * int(1 & std::get < 1 > (*it)));
					sum += sgn * std::get < 2 > (*it) * subMatrix(0, std::get < 1 > (*it)).determinant(level + 1);

				} else {
					T const sgn = T(1 - 2 * int(1 & std::get < 1 > (*it)));
					auto fut = hpx::async([sgn, it, this, level]() {
						return sgn * std::get < 2 > (*it) * subMatrix(0, std::get < 1 > (*it)).determinant(level + 1);
					});
					futs.push_back(std::move(fut));
				}
			}
		}
		for (auto &f : futs) {
			sum += f.get();
		}
		return sum;
	}
private:
	int nRow;
	int nCol;
	std::vector<std::tuple<int, int, T>> elements;
};

#endif /* INCLUDE_SPARSEMATRIX_HPP_ */
