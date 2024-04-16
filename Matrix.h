#pragma once
#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <iomanip>
#include <initializer_list>
using namespace std;

class Matrix
{
	size_t rows;
	size_t columns;
	vector<vector<double> >matrix;
public:
	Matrix();
	Matrix(size_t m_rows, size_t m_columns, double fill = 0.0);
	Matrix(size_t m_rows, size_t m_columns, initializer_list<initializer_list<double>> element, double fill = 0.0);
	Matrix(const Matrix& m) = default;

	size_t Rows() const;
	size_t Columns() const;
	size_t size(size_t dim) const;
	Matrix T();
	double det();
	Matrix inv();
	Matrix Degree(size_t num);
	double Trace_of_matrix();
	vector<double> CharacteristicPolynomial();

	static Matrix IdentityMatrix(size_t dimension);
	static Matrix mult_by_elment(Matrix A, Matrix B);

	//повертає число за координатами
	double get(int row, int col);
	//заміна числа матриці за координатами на вказане число
	Matrix rep(int row, int col, double replacement);

	Matrix operator+(Matrix& other)const;
	Matrix operator+(const double& num);
	Matrix operator-(Matrix& other)const;
	Matrix operator-(const double& num);
	Matrix operator*(Matrix& other)const;
	Matrix operator*(const double& num);
	Matrix operator/(const double& scalar);

	friend istream& operator>>(std::istream& in, Matrix& m);
	friend ostream& operator<<(std::ostream& out, const Matrix& m);
	~Matrix();
};