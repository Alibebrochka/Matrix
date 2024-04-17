#include "Matrix.h"

Matrix::Matrix() :rows(1),  columns(1), matrix(2, vector<double>(2,0.0)){}

Matrix::Matrix(size_t m_rows, size_t m_columns, double fill)
{
	rows = m_rows;
	columns = m_columns;
	matrix = vector<vector<double> >(rows, vector<double>(columns, fill));
}

Matrix::Matrix(size_t m_rows, size_t m_columns, initializer_list<initializer_list<double>> element, double fill)
{
	rows = m_rows;
	columns = m_columns;
	matrix = vector<vector<double> >(rows, vector<double>(columns, fill));
	for (int i = 0; i < rows; i++) {
		int j = 0;
		for (double x : *(element.begin() + i)) {
			matrix[i][j] = x;
			j++; 
		}
	}
}

size_t Matrix::Rows() const
{
	return rows;
}

size_t Matrix::Columns() const
{
	return columns;
}

size_t Matrix::size(size_t dim) const
{
	try {
		if (dim == 0 || dim > 2)throw exception("Dimension must be 1 or 2!");
		return dim == 1 ? rows : columns;
	}
	catch (exception& e) {
		cout << e.what() << endl;
		exit(EXIT_FAILURE);
	}
}

Matrix Matrix::T()
{
	try {
		if (rows != columns)
			throw exception("The matrix must be square");
		Matrix transpose(matrix[0].size(), matrix.size());
		for (size_t i = 0; i < matrix.size(); ++i)
			for (size_t j = 0; j < matrix[0].size(); ++j)
				transpose.matrix[j][i] = matrix[i][j];
		return transpose;
	}
	catch (exception& e) {
		cout << e.what() << endl;
		exit(EXIT_FAILURE);
	}

}

double Matrix::det()
{
	try {
		if (rows != columns)
			throw exception("The matrix must be square");

		//копіювання матриці в triangl_view
		Matrix triangl_view(*this);

		//зведеня triangl_view до трикутного вигляду
		for (int i = 0; i < columns; ++i)
			for (int j = i + 1; j < rows; ++j) {
				double x = (triangl_view.matrix[j][i] / triangl_view.matrix[i][i]);
				for (int k = i; k < columns; ++k)
					triangl_view.matrix[j][k] = triangl_view.matrix[j][k] -(triangl_view.matrix[i][k] * x);
			}
		//перемноження елементів головної діагоналі
		double Det = 1;
		for (short i = 0; i < columns; ++i)
			Det *= triangl_view.matrix[i][i];
		return Det;
	}
	catch (exception& e) {
		cout << e.what() << endl;
		exit(EXIT_FAILURE);
	}
}

Matrix Matrix::inv()
{
	try {
		if (rows != columns)
			throw exception("The matrix must be square");
		double d = this->det();
		if (d == 0.0)
			throw exception("The matrix is singular");
		if (rows == 1)
			return *this;
		else {
			Matrix res(rows, columns);
			Matrix subm(rows - 1, columns - 1);
			for (size_t i = 0; i < rows; ++i) {
				for (size_t j = 0; j < columns; ++j) {
					size_t subi = 0;
					for (size_t k = 0; k < rows; ++k) {
						size_t subj = 0;
						if (k == i) continue;
						for (size_t n = 0; n < columns; ++n) {
							if (n == j)continue;
							subm.matrix[subi][subj] = matrix[k][n];
							++subj;
						}
						++subi;
					}
					res.matrix[j][i] = pow(-1.0, i + j) * subm.det() / this->det();
				}
			}
			return res;
		}
	}
	catch (exception& e) {
		cout << e.what() << endl;
		exit(EXIT_FAILURE);
	}
}

Matrix Matrix::Degree(size_t num)
{
	Matrix res(matrix[0].size(), matrix.size());

	Matrix copy(*this);

	res = copy;
	for (int i = 1; i < num; ++i) 
		res = res * copy;

	return res;
}

double Matrix::Trace_of_matrix()
{
	double trace_of_matrix{};
	for (int i = 0; i < matrix.size(); ++i)
		trace_of_matrix += matrix[i][i];

	return trace_of_matrix;
}

vector<double> Matrix::CharacteristicPolynomial()
{
	////цикл розмірності
	//for (int i = 0; i < matrix.size(); ++i) {
	//	double sum{};
	//	//цикл мінорів
	//	for (int ii = 0; ii <= i; ++ii) {
	//		Matrix minor(matrix.size() - i, matrix.size() - i);
	//
	//		//цикли копіювання
	//		for (int j = 0; j < minor.Columns(); ++j)
	//			for (int k = 0; k < minor.Rows(); ++k) 
	//				minor.matrix[j][k] = matrix[j + ii][k + ii];
	//
	//		sum += minor.det();
	//	}
	//	result.push_back(sum);
	//}

	Matrix copy(*this);

	vector<double> result{};
	result.push_back(0);
	for (size_t i = 1; i < matrix.size() + 1; ++i) {
		double differences = copy.Degree(i).Trace_of_matrix();
		for (size_t j = 1; j < i; ++j)
			differences -= (result[j] * copy.Degree(i - j).Trace_of_matrix());
		result.push_back((1.0 / i) * (differences));
	}

	return result;
}

Matrix Matrix::IdentityMatrix(size_t dimension)
{
	Matrix identity_matrix(dimension, dimension);
	for (size_t i = 0; i < dimension; ++i)
		identity_matrix.matrix[i][i] = 1;
	return identity_matrix;
}

Matrix Matrix::mult_by_elment(Matrix A, Matrix B)
{
	for (int i = 0; i < A.Rows(); ++i)
		for (int j = 0; j < A.Columns(); ++j)
			A = A.rep(i, j, A.get(i, j) * B.get(i, j));
	return A;
}

std::istream& operator>>(std::istream& in, Matrix& m)
{
	for ( auto& i : m.matrix) 
		for ( auto& j : i) 
			in >> j;
	return in;
}

std::ostream& operator<<(std::ostream& out, const Matrix& m)
{
	for (const auto& i : m.matrix) {
		for (const auto& j : i) 
			out << setw(13) << j;
		cout << endl;
	}
	return out;
}

double Matrix::get(int row, int col)
{
	try {
		if (row > rows || col > columns)
			throw exception("Going beyond the scope of the array");
		return matrix[row][col];
	}
	catch (exception& e) {
		cout << e.what() << endl;
		exit(EXIT_FAILURE);
	}
}

Matrix Matrix::rep(int row, int col, double replacement)
{
	try {
		if (row > rows || col > columns)
			throw exception("Going beyond the scope of the array");
		Matrix result(rows, columns);
		for (int i = 0; i < rows; ++i) {
			for (int j = 0; j < columns; ++j) {
				if (i == row && j == col)
					result.matrix[i][j] = replacement;
				else
					result.matrix[i][j] = matrix[i][j];
			}
		}
		return result;
	}
	catch (exception& e) {
		cout << e.what() << endl;
		exit(EXIT_FAILURE);
	}
}

Matrix Matrix::operator+(Matrix& other) const
{
	try {
		if (rows != other.rows || columns != other.columns)
			throw exception("Matrix dimensions do not match for addition");
		Matrix result(rows, columns);
		for (size_t i = 0; i < rows; ++i)
			for (size_t j = 0; j < columns; ++j)
				result.matrix[i][j] = matrix[i][j] + other.matrix[i][j];
		return result;
	}
	catch (exception& e) {
		cout << e.what() << endl;
		exit(EXIT_FAILURE);
	}
}

Matrix Matrix::operator+(const double& num)
{
	Matrix result(rows, columns);
	for (int i = 0; i < rows; ++i)
		for (int j = 0; j < columns; ++j)
			result.matrix[i][j] = matrix[i][j] + num;
	return result;
}

Matrix Matrix::operator-(Matrix& other) const
{
	try {
		if (rows != other.rows||columns !=other.columns)
			throw exception("Matrix dimensions do not match for subtraction");
		Matrix result(rows, columns);
		for (size_t i = 0; i < rows; ++i)
			for (size_t j = 0; j < columns; ++j)
				result.matrix[i][j] = matrix[i][j] - other.matrix[i][j];
		return result;
	}
	catch (exception& e) {
		cout << e.what() << endl;
		exit(EXIT_FAILURE);
	}
}

Matrix Matrix::operator-(const double& num)
{
	Matrix result(rows, columns);
	for (int i = 0; i < rows; ++i)
		for (int j = 0; j < columns; ++j)
			result.matrix[i][j] = matrix[i][j] - num;
	return result;
}

Matrix Matrix::operator*(Matrix& other) const
{
	try {
		if (columns != other.rows)
			throw exception("Dimension must be 1 or 2!");
		Matrix result(rows, other.columns);
		for (size_t i = 0; i < rows; ++i)
			for (size_t j = 0; j < other.columns; ++j)
				for (size_t k = 0; k < columns; ++k)
					result.matrix[i][j] += matrix[i][k] * other.matrix[k][j];
		return result;
	}
	catch (exception& e) {
		cout << e.what() << endl;
		exit(EXIT_FAILURE);
	}
}

Matrix Matrix::operator*(const double& num)
{
	Matrix result(rows, columns);
	for (int i = 0; i < rows; ++i)
		for (int j = 0; j < columns; ++j)
			result.matrix[i][j] = matrix[i][j] * num;
	return result;
}

Matrix Matrix::operator/(const double& scalar)
{
	try {
		if (scalar == 0.0)throw exception("Dimension by zero");
		Matrix result(rows, columns);
		for (size_t i = 0; i < rows; ++i)
			for (size_t j = 0; j < columns; ++j)
				result.matrix[i][j] = matrix[i][j] / scalar;
		return result;
	}
	catch (exception& e) {
		cout << e.what() << endl;
		exit(EXIT_FAILURE);
	}
}

Matrix::~Matrix(){}