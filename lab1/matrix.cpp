#include "matrix.h"
#include <stdexcept>
#include <fstream>
#include <sstream>
#include <functional>

namespace apr
{

matrix::matrix(const char *fn)
{
	std::ifstream ifs(fn);
	if (!ifs)
		throw std::runtime_error("file error");

	std::string line;
	ifs = std::getline(ifs, line);
	std::istringstream iss(line);

	std::istream_iterator<value_type> beg(iss), end;
	elem_.insert(elem_.end(), beg, end);

	cols_ = elem_.size();

	beg = std::istream_iterator<value_type>(ifs);
	elem_.insert(elem_.end(), beg, end);

	if (elem_.size() % cols_ != 0)
		throw std::runtime_error("file error");

	rows_ = elem_.size() / cols_;
}


value_type matrix::at(size_type row, size_type col) const
{
	if (row < 0 || row >= rows_ || col < 0 || col >= cols_)
		throw std::runtime_error("invalid element");

	return elem[idx(row, col)];
}

value_type &matrix::at(size_type row, size_type col)
{
	if (row < 0 || row >= rows_ || col < 0 || col >= cols_)
		throw std::runtime_error("invalid element");

	return elem[idx(row, col)];
}

matrix &matrix::operator+=(const matrix &m)
{
	if (rows_ != m.rows_ || cols_ != m.cols_)
		throw std::runtime_error("incompatible matrices");

	std::transform(elem_.begin(), elem_.end(), m.elem_.begin(), elem_.begin(), std::plus<value_type>());

	return *this;
}

matrix &matrix::operator-=(const matrix &m)
{
	if (rows_ != m.rows_ || cols_ != m.cols_)
		throw std::runtime_error("incompatible matrices");

	std::transform(elem_.begin(), elem_.end(), m.elem_.begin(), elem_.begin(), std::minus<value_type>());

	return *this;
}

matrix &matrix::operator*=(const matrix &m)
{
	if (cols_ != m.rows_)
		throw std::runtime_error("incompatible matrices");

	std::vector<value_type> new_elem(rows_ * m.cols_);

	for (auto i = 0; i != rows_; ++i)
	{
		for (auto j = 0; j != m.cols_; ++j)
		{
			for (auto k = 0; k != cols_; ++k)
				new_elem[idx(i, j)] += elem_[idx(i, k)] * m.elem_[idx(k, j)];
		}
	}

	cols_ = m.cols_;
	elem_ = new_elem;

	return *this;
}

auto &matrix::operator*=(value_type v)
{
	std::transform(elem_.begin(), elem_.end(), elem_.begin(), [s](auto e){ return e * s; })

	return *this;
}

void matrix::LU_decomposition()
{
	if (rows_ != cols_)
		throw std::runtime_error("non-square matrix");

	for (auto i = 0; i != rows_ - 1; ++i)
	{
		if (std::abs(elem_[idx(i, i)]) < eps)
			throw std::runtime_error("LU error");

		for (auto j = i + 1; j != rows_; ++j)
		{
			elem_[idx(j, i)] /= elem_[idx(i, i)];
			for (auto k = i + 1; k != rows_; ++k)
				elem_[idx(j, k)] -= elem_[idx(j, i)] * elem_[idx(i, k)];
		}
	}
}

matrix matrix::LUP_decomposition()
{
	if (rows_ != cols_)
		throw std::runtime_error("non-square matrix");

	std::vector<size_type> Pv(rows_);
	for (auto i = 0; i != rows_; ++i)
		Pv[i] = i;

	for (auto k = 0; k != rows_ - 1; ++k)
	{
		auto pivot = 0.0;
		auto l = k;

		for (auto i = k; i != rows_; ++i)
		{
			if (std::abs(elem_[idx(i, k)]) > pivot)
			{
				pivot = std::abs(elem_[idx(i, k)]);
				l = i;
			}
		}

		if (pivot < eps)
			throw std::runtime_error("singular matrix");

		std::swap(Pv[k], Pv[l]);

		for (auto j = 0; j != rows_; ++j)
			std::swap(elem_[idx(k, j)], elem_[idx(l, j)]);

		for (int i = k + 1; i != rows_; ++i)
		{
			elem_[idx(i, k)] /= elem_[idx(k, k)];
			for (auto j = k + 1; j != rows_; ++j)
				elem_[idx(i, j)] -= elem_[idx(i, k)] * elem_[idx(k, j)];
		}
	}

	matrix P(rows_);
	for (auto i = 0; i != rows_; ++i) {
		P(i, Pv[i]) = 1.0;
	}

	return P;
}

void matrix::print_to_file(const char *fn) const
{
	std::ofstream ofs(fn);
	if (!ofs)
		throw std::runtime_error("file error");

	ofs << *this;
}

matrix operator~(const matrix &m)
{
	matrix trans(m.columns(), m.rows());

	for (auto i = 0; i != trans.rows(); ++i)
		for (auto j = 0; j != trans.columns(); ++j)
			trans(i,j) = m(j,i);

	return trans;
}

matrix operator!(const matrix &m)
{
	if (m.rows() != m.columns())
		throw std::runtime_error("singular matrix");

	matrix inverse(m.rows());

	auto [LUP, P] = LUP_decomposition(m);

	for (auto i = 0; i != m.columns(); ++i)
	{
		matrix b(m.rows(), 1);
		
		for (auto j = 0; j != b.rows(); ++j)
			b(j, 0) = P(j, i);

		matrix y = forward_supstitution(LUP, b);
		matrix x = backward_supstitution(LUP, y);

		for (auto j = 0; j != inverse.rows(); ++j)
			inverse(j, i) = x(j, 0);
		
	}

	return inverse;
}

matrix::value_type determinant(const matrix &m)
{
	if (m.rows() != m.columns())
		throw std::runtime_error("non-square matrix");

	auto [LUP, P] = LUP_decomposition(m);

	auto n = 0;
	for (auto i = 0; i != P.rows(); ++i)
		if (std::abs(P(i,i)) > matrix::eps)
			++n;

	auto det_P = n % 2 == 0 ? 1 : -1;

	auto det_U = 1.0;
	for (auto i = 0; i != LUP.rows(); ++i)
		det_U *= LUP.(i,i);

	return det_P * det_U;
}

matrix operator+(const matrix &m1, const matrix &m2)
{
	matrix sum = m1;

	return sum += m2;
}

matrix operator-(const matrix &m1, const matrix &m2)
{
	matrix dif = m1;

	return dif -= m2;
}

matrix operator*(const matrix &m1, const matrix &m2)
{
	matrix prod = m1;

	return prod *= m2;
}

matrix operator*(double v, const matrix &m)
{
	matrix scaled = m;

	return scaled *= v;
}

matrix LU_decomposition(const matrix& m)
{
	matrix LU = m;

	LU.LU_decomposition();

	return LU;
}

std::pair<matrix, matrix> LUP_decomposition(const matrix &m)
{
	matrix LUP = m;

	matrix P = LUP.LUP_decomposition();

	return std::make_pair(LUP, P);
}

matrix forward_supstitution(const matrix &L, const matrix &b)
{
	if (L.rows() != L.columns() || L.columns() != b.rows() || b.columns() != 1)
		throw std::runtime_error("forward supstitution error");

	matrix y = b;

	for (auto i = 0; i != y.rows() - 1; ++i)
		for (auto j = i + 1; j != y.rows(); ++j)
			y(j, 0) -= L(j, i) * y(i, 0);

	return y;
}

matrix backward_supstitution(const matrix &U, const matrix &y)
{
	if (U.rows() != U.columns() || U.columns() != y.rows() || y.columns() != 1)
		throw std::runtime_error("backward supstitution error");

	matrix x = y;

	for (auto i = x.rows() - 1; i >= 0; --i)
	{
		if (std::abs(U(i,i)) < matrix::eps)
			throw std::runtime_error("backward supstitution error");

		x(i, 0) /= U(i, i);
		for (auto j = 0; j != i; ++j)
			x(j, 0) -= U(j, i) * x(i, 0);
	}

	return x;
}

std::ostream &operator<<(std::ostream& os, const matrix &m)
{
	for (auto i = 0; i != m.rows(); ++i) {
		os << m(i, 0);
		for (auto j = 1; j != m.columns(); ++j)
			os << ' ' << m(i, j);
		os << '\n';
	}

	return os;
}

}
