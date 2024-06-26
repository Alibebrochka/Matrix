﻿#include <cmath>
#include <Windows.h>
#include"Matrix.h"

void Neural_Network_Learning();
void Nonlinear_Equations();
void linear_Equations();
void Finding_eigenvalues();

int main() {
	cout << "\033[92m";
	SetConsoleCP(1251);
	SetConsoleOutputCP(1251);

	Matrix A(3, 3, {
		{5,1,2},
		{5,4,3},
		{2,1,3} }),
		B(4, 4, {
		{33,10,-10,55},
		{-21,5,88,-15},
		{13,31,-18,52},
		{13,17,51,22 } }),
		C(5, 5, {
			{1,1,1,2,3} ,
			{5,8,7,1,9} ,
			{2,6,7,9,5} ,
			{6,8,6,4,1} ,
			{6,8,5,2,9} });
	Matrix D(4, 4, {
		{33,10,-10,55},
		{-21,5,88,-15},
		{13,31,-18,52},
		{13,17,51,22 } }),
		E(4, 1, {
			{-12},
			{29 },
			{66 },
			{-19} });
	//cout << "fin:\n" << C.Eigenvalues(EMFRE_LU);
	//cout << D.SLE(E, EMFRS_Matrix) << endl;
	//cout << D.SLE(E, EMFRS_LU) << endl;
	cout << D.Characteristic_Polynomial() << endl;
	cout << D << endl;

	//for (int x = 0; x <= 2; ++x) 
		//cout << exp(x) - exp(-x) - 2 << endl;
	//const double a = 0.5, b = 1.5;
	//double x = a, Fa{}, Fb{};
	//
	//Fa = (pow(x, 3) - 2 * x) * (x + 1);
	//x = b;
	//Fb = (pow(x, 3) - 2 * x)* (x + 1);
	//
	//cout << "Fa = " << Fa << '\n' << "Fb = " << Fb << endl;
	//cout << 0.1e-2;
	//Matrix A(3, 3, { 
	//	{1,    0,       0      },
	//	{0,    3,       sqrt(5)},
	//	{0,    sqrt(5), 4      } }),
	//	B = Matrix::IdentityMatrix(3);
	//cout << A << endl;
	//cout << B << endl;
	//A = A - B;
	//cout << A << endl;
	//cout << A.det() << endl;
	return 0;
}

void Neural_Network_Learning()
{
	auto sigmoid = [](Matrix* A) {
		Matrix& temp = *A;
		for (size_t i = 0; i < A->Rows(); ++i)
			for (size_t j = 0; j < A->Columns(); ++j)
				temp = A->rep(i, j, 1 / (1 + exp((A->get(i, j)) * (-1))));
		return A;
		};

	Matrix W12(5, 4, {
		{ 0.3,0.8,0.7,0.2},
		{ 0.6,0.9,0.1,0.4},
		{ 0.4,0.5,0.6,0.9},
		{ 0.8,0.2,0.3,0.7},
		{ 0.5,0.7,0.8,0.1} }),
		W23(6, 5, {
			{ 0.7, 0.3, 0.5, 0.4, 0.9 },
			{ 0.5, 0.2, 0.8, 0.3, 0.1 },
			{ 0.8, 0.6, 0.3, 0.5, 0.2 },
			{ 0.3, 0.7, 0.4, 0.1, 0.9 },
			{ 0.4, 0.8, 0.6, 0.2, 0.7 },
			{ 0.9, 0.4, 0.2, 0.6, 0.1 } }),
			W34(3, 6, {
				{ 0.5, 0.2, 0.7, 0.1, 0.8, 0.4 },
				{ 0.7, 0.3, 0.8, 0.9, 0.1, 0.6 },
				{ 0.4, 0.5, 0.2, 0.7, 0.9, 0.3 } }),
				I(4, 1, {
					{ 0.6 },
					{ 0.2 },
					{ 0.5 },
					{ 0.3 } }),
					T(3, 1, {
						{ 0.1 },
						{ 0.9 },
						{ 0.5 } });
	
	Matrix O1 = I;	//вихідні сигнали із шару 1
	Matrix X2 = W12 * O1;	//комбіновані згладжені вхідні сигнали, які надходять на вузли шару 2 
	Matrix* O2 = sigmoid(&X2);	//обчислення сигналів які виходять із вузлів шару 2
	Matrix X3 = W23 * (*O2);	//комбіновані згладжені вхідні сигнали, які надходять на вузли шару 3
	Matrix* O3 = sigmoid(&X3);	//обчислення сигналів які виходять із вузлів шару 3
	Matrix X4 = W34 * (*O3);	//комбіновані згладжені вхідні сигнали, які надходять на вузли шару 4
	Matrix* O4 = sigmoid(&X4);	//обчислення сигналів які виходять із вузлів шару 4
	cout << "output" << '\n' << *O4 << endl;


	//розрахунок помилок які утворилися на кожному шарі
	Matrix E4 = T - *O4;		//шар 4
	Matrix E3 = W34.T() * E4;	//шар 3
	Matrix E2 = W23.T() * E3;	//шар 2
	const double lr = 0.5;	//коефіцієнт швидкості навчання


	//обчислення поправок до вагових коефіцієнтів шару 3 та 4
	Matrix OxT = (*O3).T();	//збереження транспонованої матриці значень вихідних сигналів
	Matrix Exlr = E4 * lr;	//збереження добутку помилки на коефіцієнт швидкості навчання
	Matrix Oxi = (*O4 - 1) * (-1);	//(1 - Ох)
	Matrix dWxx = Matrix::mult_by_elment(Matrix::mult_by_elment(Exlr, *O4), Oxi) * OxT;	//обчислення поправок 
	W34 = W34 + dWxx;	//оновлення значень вагових коефіціентів зв'язків між вузлами шару 3 і 4
	std::cout << W34 << endl;


	//обчислення поправок до вагових коефіцієнтів шару 2 та 3
	OxT = (*O2).T();	
	Exlr = E3 * lr;	
	Oxi = (*O3 - 1) * (-1);
	dWxx = Matrix::mult_by_elment(Matrix::mult_by_elment(Exlr, *O3), Oxi) * OxT;
	W23 = W23 + dWxx;	//оновлення значень вагових коефіціентів зв'язків між вузлами шару 2 і 3
	cout << W23 << endl;

	
	//обчислення поправок до вагових коефіцієнтів шару 1 та 2 
	// 
	OxT = O1.T();	
	Exlr = E2 * lr;	
	Oxi = (*O2 - 1) * (-1);
	dWxx = Matrix::mult_by_elment(Matrix::mult_by_elment(Exlr, *O2), Oxi) * OxT;
	W12 = W12 + dWxx;	//оновлення значень вагових коефіціентів зв'язків між вузлами шару 1 і 2
	cout << W12 << endl;


	O1 = I;	
	X2 = W12 * O1;	 
	O2 = sigmoid(&X2);	
	X3 = W23 * (*O2);	
	O3 = sigmoid(&X3);	
	X4 = W34 * (*O3);	
	O4 = sigmoid(&X4);	
	cout << "output" << '\n' << *O4 << endl;

}

void Nonlinear_Equations() {
	double x{}, y{}, x1{}, y1{}, er1{}, er2{};

	//метод простих ітерацій
	cout << "<<<<<<<<<<<<<<<метод простих ітерацій>>>>>>>>>>>>>>>" << endl;
	x = 1;
	y = 1;
	cout << "x = " << x << '\n' << "y = " << y << endl;

	for (int i = 1; i <= 7; ++i) {
		x1 = (sin(y) - 1.6) / 2;
		y1 = 0.8 - cos(x + 0.5);

		y = y1;
		x = x1;
		cout	<< "x1 = " << x1 << '\n' 
				<< "y1 = " << y1 << '\n' 
				<< "x = " << x << '\n' 
				<< "y = " << y << '\n' 
			<< '\t' << i  << endl;
	}
	er1 = cos(x + 0.5) + y - 0.8;
	er2 = sin(y) - 2 * x - 1.6;
	cout << "er1 = " << er1 << '\n' << "er2 = " << er2 << endl;
	cout << endl;

	//метод ітерації по Зейделю
	cout << "<<<<<<<<<<<<<<<метод ітерації по Зейделю>>>>>>>>>>>>>>>" << endl;
	x = 1;
	y = 1;
	cout << "x = " << x << '\n' << "y = " << y << endl;
	for (int i = 1; i <= 7; ++i) {
		x = (sin(y) - 1.6) / 2;
		y = 0.8 - cos(x + 0.5);
		cout << "x = " << x << '\n' << "y = " << y << endl;
	}
	er1 = cos(x + 0.5) + y - 0.8;
	er2 = sin(y) - 2 * x - 1.6;
	cout << "er1 = " << er1 << '\n' << "er2 = " << er2 << endl;
	cout << endl;

	//МЕТОД НЬЮТОНА
	cout << "<<<<<<<<<<<<<<<МЕТОД НЬЮТОНА>>>>>>>>>>>>>>>" << endl;
	double f{}, g{}, fx{}, fy{}, gx{}, gy{}, det{}, detx{}, dety{}, hx{}, hy{};
	x = 9;
	y = -9;
	cout << "x = " << x << '\n' << "y = " << y << endl;
	for (int i = 1; i <= 12; ++i) {
		f = cos(x + 0.5) + y - 0.8;
		g = sin(y) - 2 * x - 1.6;
		fx = -sin(x + 0.5);
		fy = 1.;
		gx = -2.;
		gy = cos(y);
		det = fx * gy - fy * gx;
		detx = -f * gy + fy * g;
		dety = f * gx - fx * g;
		hx = detx / det;
		hy = dety / det;
		x = x + hx;
		y = y + hy;
		cout << "f = " << f << '\n'
			<< "g = " << g << '\n'
			<< "fx = " << fx << '\n'
			<< "fy = " << fy << '\n'
			<< "gx = " << gx << '\n'
			<< "gy = " << gy << '\n'
			<< "det = " << det << '\n'
			<< "detx = " << detx << '\n'
			<< "dety = " << dety << '\n'
			<< "hx = " << hx << '\n'
			<< "hy = " << hy << '\n'
			<< "x = " << x << '\n'
			<< "y = " << y << '\n'
			<< '\t' << i << '\n' << endl;
	}
	er1 = cos(x + .5) + y - .8;
	er2 = sin(y) - 2 * x - 1.6;
	cout << "er1 = " << er1 << '\n' << "er2 = " << er2 << endl;
}
void linear_Equations()
{
	Matrix A(4, 4, {
		{33,10,-10,55},
		{-21,5,88,-15},
		{13,31,-18,52},
		{13,17,51,22 } }),
		B(4, 1, {
			{-12},
			{29},
			{66},
			{-19} });
	cout << "A:\n" << A;
	cout << "A.det(): " << A.det() << endl;
	cout << "B:\n" << B << endl;

	Matrix A1 = A.inv();
	cout << "A1:\n" << A1 << endl;

	Matrix X = A1 * B ;
	cout <<"X:\n" << X << endl;

	Matrix X1 = A * X;
	cout << "X1:\n" << X1 << endl;

	cout << "A * A1:\n" << fixed << A * A1 << endl;
}

void Finding_eigenvalues()
{
	Matrix A(3, 3, {
		{1, 1, 1 },
		{1, 3, 1 },
		{1, 1, 10} 
		}),
		B_pro(5, 5, {
			{1,1,1,2,3} ,
			{5,8,7,1,9} ,
			{2,6,7,9,5} ,
			{6,8,6,4,1} ,
			{6,8,5,2,9}
			}),
		C(3, 3, {
			{ 5, -2, 1 } ,
			{ 3, -2, 3 } ,
			{ 3, -6, 7 } 
			}),
			E(3, 3, {
				{ 1, 4, 3 } ,
				{ 4, 5, 4 } ,
				{ 3, 4, 1 } 
				}),
				M_LU0(3, 3, {
					{ 1, 2, 3 } ,
					{ 4, 2, 1 } ,
					{ 5, 2, 2 } 
					}),
					M_LU1(2, 2, {
						{ 1, 5 } ,
						{ 4, 2 } 
						}),
						D(6, 6, {
							{ 1, 1, 1, 1, 1, 5 } ,
							{ 1, 5, 1, 1, 1, 1 } ,
							{ 1, 1, 3, 1, 1, 1 } ,
							{ 1, 1, 1, 7, 1, 1 } , 
							{ 1, 1, 1, 1, 5, 1 } , 
							{ 5, 1, 1, 1, 1, 10} });
	cout << "CharacteristicPolynomial:\n";
		cout << A.Characteristic_Polynomial() << endl;
	cout << endl;

	cout << "Eigenvalues:\n";
	Matrix A_egnvls = A.Eigenvalues(EMFRE_Jacobi);
	for (int i = 0; i < A_egnvls.Columns(); ++i)
		cout << A_egnvls.get(i, i) << endl;
	cout << endl;

	cout<<"x9:\n";
	Matrix A9, A1 = A, x9;
	for (int i = 0; i < 9; ++i) {
		A9 = A * A1;
		A1 = A9;
	}
	Matrix X0(3, 1, {
		{0},
		{1},
		{0} });
	x9 = A9 * X0;
	cout << x9 << endl;
	cout << endl;

	cout << "x10:\n";
	Matrix A10, x10;
	A1 = A;
	for (int i = 0; i < 10; ++i) {
		A10 = A * A1;
		A1 = A10;
	}
	x10 = A10 * X0;
	cout << x10 << endl;
	cout << endl;

	cout << "vmax:\n";
	double vmax{};
	for (int i = 0; i < 3; ++i) {
		vmax = x10.get(i, 0) / x9.get(i, 0);
		cout << vmax << endl;
	}
	cout << endl;

	Matrix B = A.inv();
	cout << "B:\n" << B << endl;
	cout << endl;

	Matrix B1 = B, B9, B10;
	for (int i = 0; i < 9; ++i) {
		B9 = B1 * B;
		B1 = B9;
	}
	for (int i = 0; i < 10; ++i) {
		B10 = B1 * B;
		B1 = B10;
	}
	Matrix y9, y10;
	y9 = B9 * X0;
	cout << "y9:\n" << y9 << endl;
	y10 = B10 * X0;
	cout << "y10:\n" << y10 << endl;

	double vmin{};
	for (int i = 0; i < 3; ++i) {
		vmin = y9.get(i, 0) / y10.get(i, 0);
		cout << vmin << endl;
	}
	cout << endl;

	cout << "vmidl:\n" << A.det() / (vmin * vmax) << endl;
}
