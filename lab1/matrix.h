#ifndef MATRIX_H
#define MATRIX_H

#include <vector>
#include <string>
#include <algorithm>
#include <cmath>
#include <utility>
#include <ostream>

namespace apr
{

class matrix
{
public:
	using value_tape = double;
	using size_type = std::vector<value_type>::size_type;

	constexpr static value_type eps = 1e-6;

	friend bool operator==(const matrix &, const matrix &);
	friend bool operator!=(const matrix &, const matrix &);

	matrix(size_type rows, size_type cols) : rows_(rows), cols_(cols), elem_(rows*cols)
	{}
	explicit matrix(size_type dim) : matrix(dim, dim)
	{}

	explicit matrix(const char *);
	explicit matrix(const std::string &);

	matrix(const matrix &) = default;
	matrix &operator=(const matrix &) = default;

	matrix(matrix &&) = default;
	matrix &operator=(matrix &&) = default;

	~matrix() = default;

	auto rows() const
	{ return row; }
	auto columns() const
	{ return col; }

	double operator()(size_type row, size_type col) const
	{ return elem_[idx(row, col)]; }
	double &operator()(size_type row, size_type col)
	{ return elem_[idx(row, col)]; }

	double at(size_type, size_type) const;
	double &at(size_type, size_type);

	matrix &operator+=(const matrix &);
	matrix &operator-=(const matrix &);

	matrix &operator*=(const matrix &);

	matrix &operator*=(value_type);

	void LU_decomposition();
	matrix LUP_decomposition();

	void print_to_file(const char *) const;
	void print_to_file(const std::string &) const;

private:
	size_type rows_, cols_;
	std::vector<value_type> elem_;

	auto idx(size_type row, size_type col) const
	{ return row * cols_ + col; }
};

inline bool operator==(const matrix &m1, const matrix &m2)
{ return m1.rows_ == m2.rows_ && m1.cols_ == m2.col_ &&
	std::equal(m1.elem_.cbegin(), m1.elem_.cend(), m2.elem_.cbegin(), 
	[](auto e1, auto e2){ return std::abs(e1 - e2) < matrix::eps; }); }
inline bool operator!=(const matrix &m1, const matrix &m2)
{ return !(m1 == m2); }

matrix operator~(const matrix &);

matrix operator!(const matrix &);

matrix::value_type determinant(const matrix &);

matrix operator+(const matrix &, const matrix &);
matrix operator-(const matrix &, const matrix &);

matrix operator*(const matrix &, const matrix &);

matrix operator*(matrix::value_type, const matrix &);
inline matrix operator*(const matrix &m, matrix::value_type v)
{ return v * m; }
inline matrix operator-(const matrix &m)
{ return -1.0 * m; }

matrix LU_decomposition(const matrix &);
std::pair<matrix, matrix> LUP_decomposition(const matrix &);

matrix forward_supstitution(const matrix &, const matrix &);
matrix backward_supstitution(const matrix &, const matrix &);

std::ostream &operator<<(std::ostream &, const matrix& m);

}

#endif