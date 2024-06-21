#include "Matrix.h"

void Matrix::LU_Schedule(Matrix A)
{
	Matrix t_L(A.rows, A.columns);
	Matrix t_U = Matrix::IdentityMatrix(A.matrix.size());
	L = t_L.matrix;
	U = t_U.matrix;
	for (int i = 0; i < A.matrix.size(); ++i)
		for (int j = 0; j < A.matrix.size(); ++j) {
			if (i >= j) {
				for (int k = 0; k <= (j ? j - 1 : 0); ++k)
					L[i][j] += L[i][k] * U[k][j];
				L[i][j] = A.matrix[i][j] - L[i][j];
			}
			if (i < j) {
				for (int k = 0; k <= (i ? i - 1 : 0); ++k)
					U[i][j] += L[i][k] * U[k][j];
				U[i][j] = (A.matrix[i][j] - U[i][j]) * (1 / L[i][i]);
			}
		}
}

Matrix::Matrix() :rows(1),  columns(1), matrix(2, vector<double>(2,0.0)){}

Matrix::Matrix(size_t m_rows, size_t m_columns, double fill)
{
	rows = m_rows;
	columns = m_columns;
	matrix = vector<vector<double> >(rows, vector<double>(columns, fill));
	L = matrix;
	U = matrix;
}

Matrix::Matrix(size_t m_rows, size_t m_columns, initializer_list<initializer_list<double>> element, double fill)
{
	rows = m_rows;
	columns = m_columns;
	matrix = vector<vector<double> >(rows, vector<double>(columns, fill));
	for (size_t i = 0; i < rows; i++) {
		size_t j = 0;
		for (double x : *(element.begin() + i)) {
			matrix[i][j] = x;
			j++; 
		}
	}
}

Matrix::Matrix( vector<vector<double>> m)
{
	rows = m.size();
	columns = m.size();
	matrix = vector<vector<double> >(rows, vector<double>(columns, 0));
	for (int i = 0; i < m.size(); ++i) {
		for (int j = 0; j < m.size(); ++j)
			matrix[i][j] = m[i][j];
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

double Matrix::det() const
{
	try {
		if (rows != columns)
			throw exception("The matrix must be square");
		//копіювання матриці в triangl_view
		Matrix triangl_view(*this);
		//зведеня triangl_view до трикутного вигляду
		for (size_t  i = 0; i < columns; ++i)
			for (size_t  j = i + 1; j < rows; ++j) {
				double x = (triangl_view.matrix[j][i] / triangl_view.matrix[i][i]);
				for (size_t  k = i; k < columns; ++k)
					triangl_view.matrix[j][k] = triangl_view.matrix[j][k] -(triangl_view.matrix[i][k] * x);
			}
		//перемноження елементів головної діагоналі
		double Det = 1;
		for (size_t i = 0; i < columns; ++i)
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
	for (size_t i = 1; i < num; ++i)
		res = res * copy;

	return res;
}

double Matrix::Trace_of_matrix()
{
	double trace_of_matrix{};
	for (size_t i = 0; i < matrix.size(); ++i)
		trace_of_matrix += matrix[i][i];

	return trace_of_matrix;
}

vector<double> Matrix::CharacteristicPolynomial()
{ 
	Matrix copy(*this);
	vector<double> result{};
	result.push_back(0);
	for (size_t i = 1; i < matrix.size() + 1; ++i) {
		double differences = copy.Degree(i).Trace_of_matrix();
		for (size_t j = 1; j < i; ++j)
			differences -= (result[j] * copy.Degree(i - j).Trace_of_matrix());
		result.push_back((1.0 / i) * (differences));
	}
	result[0] = 1;
	return result;
}

Matrix Matrix::Eigenvalues(EMethods_For_Resolving_Eigenvalues egvalues)
{
	auto sgn = [](double a) { return a > 0 ? 1 : a < 0 ? -1 : 0; };
	switch (egvalues)
	{
	case EMFRE_Jacobi:
		try {
			Matrix A(*this), B = A;
			if (!(A.T() == A))throw exception("Matrix must be symmetric");

			double p{}, q{}, s{}, r{}, d{}, c{};

			for (int it = 0; it < 100; ++it) {
				size_t  i = 1, j = 0;
				double max_elem = A.get(1, 0);
				//Знаходження максимального внедіагонального елемента
				for (size_t ii = 0; ii < A.Rows(); ++ii)
					for (size_t jj = 0; jj < A.Columns(); ++jj)
						if (ii != jj)
							if (max_elem < A.get(ii, jj)) {
								i = ii;
								j = jj;
								max_elem = A.get(ii, jj);
							}

				p = 2 * A.get(i, j);
				q = A.get(i, i) - A.get(j, j);
				d = sqrt(pow(p, 2) + pow(q, 2));

				if (q != 0) {
					r = abs(q) / (2 * d);
					c = sqrt(0.5 + r);
					s = sqrt(0.5 - r) * sgn(p * q);
				}
				else
					c = s = (sqrt(2) / 2);

				//елементи нової матриці
				B.matrix[i][i] = (pow(c, 2) * A.get(i, i)) + (pow(s, 2) * A.get(j, j)) + (2 * c * s * A.get(i, j));
				B.matrix[j][j] = (pow(s, 2) * A.get(i, i)) + (pow(c, 2) * A.get(j, j)) - (2 * c * s * A.get(i, j));
				B.matrix[i][j] = B.matrix[j][i] = (pow(c, 2) - pow(s, 2)) * A.get(1, j) + c * s * (A.get(j, j) - A.get(i, i));
				// змінювані позадіагональні елементи
				for (size_t m = 0; m < A.matrix.size(); ++m)
					if (m != i && m != j) {
						B.matrix[i][m] = B.matrix[m][i] = (c * A.get(m, i)) + (s * A.get(m, j));
						B.matrix[j][m] = B.matrix[m][j] = (-s * A.get(m, i)) + (c * A.get(m, j));
					}

				double fin{};
				for (int ii = 0; ii < A.matrix.size(); ++ii)
					fin += pow(A.matrix[ii][ii] - B.matrix[ii][ii], 2);
				double f = sqrt(fin);

				//оновлення
				A = B;

				//cout << fixed << A << endl;
				if (f == 0)
					return A;
			}
		}
		catch (exception& e) {
			cout << e.what() << endl;
			exit(EXIT_FAILURE);
		}
		break;
	case EMFRE_LU:
		LU_Schedule(*this);
		Matrix AL(L), AU(U);
		for (int i = 1;; ++i) {
			LU_Schedule(AU*AL);
			Matrix mL(L), mU(U);
			if ((abs((mL * mU).matrix[1][1] - (AL * AU).matrix[1][1])) < 0.1) {
				cout << "iteration:\t" << i << endl;
				return  AU * AL;
			}
			AL = L;
			AU = U;
		}
	}
	return Matrix();
}

const vector<vector<double>>& Matrix::GetL() const
{
	return L;
}

const vector<vector<double>>& Matrix::GetU() const
{
	return U;
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
	for (size_t  i = 0; i < A.Rows(); ++i)
		for (size_t  j = 0; j < A.Columns(); ++j)
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

double Matrix::get(size_t  row, size_t  col)
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

Matrix Matrix::rep(size_t  row, size_t  col, double replacement) const
{
	try {
		if (row > rows || col > columns)
			throw exception("Going beyond the scope of the array");
		Matrix result(*this);
		result.matrix[row][col] = replacement;
		return result;
	}
	catch (exception& e) {
		cout << e.what() << endl;
		exit(EXIT_FAILURE);
	}
}

Matrix operator+(const Matrix& m1, const Matrix& m2) 
{
	try {
		if (m2.rows != m1.rows || m2.columns != m1.columns)
			throw exception("Matrix dimensions do not match for addition");
		Matrix result(m2.rows, m2.columns);
		for (size_t i = 0; i < m2.rows; ++i)
			for (size_t j = 0; j < m2.columns; ++j)
				result.matrix[i][j] = m2.matrix[i][j] + m1.matrix[i][j];
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
	for (size_t  i = 0; i < rows; ++i)
		for (size_t  j = 0; j < columns; ++j)
			result.matrix[i][j] = matrix[i][j] + num;
	return result;
}

Matrix operator-(const Matrix& m1, const Matrix& m2)
{
	try {
		if (m2.rows != m1.rows|| m2.columns !=m1.columns)
			throw exception("Matrix dimensions do not match for subtraction");
		Matrix result(m2.rows, m2.columns);
		for (size_t i = 0; i < m2.rows; ++i)
			for (size_t j = 0; j < m2.columns; ++j)
				result.matrix[i][j] = m2.matrix[i][j] - m1.matrix[i][j];
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
	for (size_t  i = 0; i < rows; ++i)
		for (size_t  j = 0; j < columns; ++j)
			result.matrix[i][j] = matrix[i][j] - num;
	return result;
}

Matrix operator*(const Matrix& m1,const Matrix& m2) 
{
	try {
		if (m1.columns != m2.rows) throw exception("Dimension must be 1 or 2!");
		Matrix result(m1.rows, m2.columns);
		for (size_t i = 0; i < m1.rows; ++i)
			for (size_t j = 0; j < m2.columns; ++j)
				for (size_t k = 0; k < m1.columns; ++k)
					result.matrix[i][j] += m1.matrix[i][k] * m2.matrix[k][j];
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
	for (size_t  i = 0; i < rows; ++i)
		for (size_t  j = 0; j < columns; ++j)
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

bool operator==(const Matrix& m1, const Matrix& m2) 
{
	try {
		if (m2.rows != m1.rows) throw exception("Size must be the same");
		if (m2.columns != m1.columns) throw exception("Size must be the same");
		for (size_t  i = 0; i < m2.matrix.size(); ++i)
			for (size_t  j = 0; j < m2.matrix.size(); ++j)
				if (m2.matrix[i][j] != m1.matrix[i][j])
					return false;
		return true;
	}
	catch (exception& e) {
		cout << e.what() << endl;
		exit(EXIT_FAILURE);
	}
}

Matrix::~Matrix(){}
