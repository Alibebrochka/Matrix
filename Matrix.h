#pragma once
#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <iomanip>
#include <initializer_list>
#include <cmath>
#include <algorithm>

using namespace std;

const double PI = 3.141592653589793238;

enum EMethods_For_Resolving_Eigenvalues {
	EMFRE_Jacobi = 0,
	EMFRE_LU
};

	
class Matrix
{
	size_t rows;
	size_t columns;
	vector<vector<double> >matrix;

	vector<vector<double> > L;
	vector<vector<double> > U;
public:
	Matrix();
	Matrix(size_t m_rows, size_t m_columns, double fill = 0.0);
	Matrix(size_t m_rows, size_t m_columns, initializer_list<initializer_list<double>> element, double fill = 0.0);
	Matrix( vector<vector<double> > m);
	Matrix(const Matrix& m) = default;
	size_t Rows() const;
	size_t Columns() const;
	size_t size(size_t dim) const;
	Matrix T();
	double det() const;
	Matrix inv();
	Matrix Degree(size_t num);
	double Trace_of_matrix();
	// метод Леверье
	vector<double> CharacteristicPolynomial();

	void LU_Schedule(Matrix A);
	EMethods_For_Resolving_Eigenvalues egvalues{ EMFRE_Jacobi };
	Matrix Eigenvalues(EMethods_For_Resolving_Eigenvalues egvalues);
	const vector<vector<double>>& GetL() const;
	const vector<vector<double>>& GetU() const;

	static Matrix IdentityMatrix(size_t dimension);
	static Matrix mult_by_elment(Matrix A, Matrix B);

	//повертає число за координатами
	double get(size_t row, size_t col);
	//заміна числа матриці за координатами на вказане число
	Matrix rep(size_t row, size_t col, double replacement) const;

	Matrix operator+(const double& num);
	Matrix operator-(const double& num);
	Matrix operator*(const double& num);
	Matrix operator/(const double& scalar);
	friend Matrix operator+(const Matrix& m1, const Matrix& m2);
	friend Matrix operator-(const Matrix& m1, const Matrix& m2);
	friend Matrix operator*(const Matrix& m1, const Matrix& m2);
	friend bool operator==(const Matrix& m1, const Matrix& m2);

	friend istream& operator>>(std::istream& in, Matrix& m);
	friend ostream& operator<<(std::ostream& out, const Matrix& m);
	~Matrix();
};