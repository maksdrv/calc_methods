#define _USE_MATH_DEFINES
#include <cmath>
#include <iostream>
#include <typeinfo>
#include <vector>
#include <fstream>
#include <iomanip>
#define type double
using namespace std;


void set_boundaries(int num, type& a, type& b) {
	if (num == 3) {
		a = -3; b = 3;
	}
	else {
		a = -1; b = 1;
	}
}

type f(type x, int num = 1) { // num - номер теста
	switch (num) {
	case 1:
		return pow(x, 2);
		//return 1;
	case 2:
		return 1 / (1 + 25 * pow(x, 2));
	case 3:
		return 1 / (atan(1 + 10 * pow(x, 2)));
	case 4:
		return pow(pow(4 * x, 3) + pow(2 * x, 2) - 4 * x + 75, sqrt(2)) + asin(1 / (5 + x - pow(x, 2))) - 5;
	}
}

vector<type> regular_grid(type a, type b, int n) {
	vector<type> x(n + 1);
	for (int i = 0; i < n + 1; i++) {
		x[i] = a + (b - a) / n * i;
	}
	return x;
}

vector<type> Chebyshev_grid(type a, type b, int n) {
	vector<type> x(n + 1);
	for (int i = n; i >= 0; i--) {
		x[n - i] = (a + b) / 2 + (b - a) / 2 * cos((2 * i + 1) * M_PI / (2 * (n + 1)));
	}
	return x;
}

vector<type> func_grid(vector<type> x) {
	vector<type> y(x.size());
	for (int i = 0; i < x.size(); i++) {
		y[i] = f(x[i]);
	}
	return y;
}

void print(vector<type> b) {
	for (int i = 0; i < b.size(); i++) cout << setw(14) << b[i] << setw(14);
}

vector<type> output(vector<type> b, ofstream& fs) {
	for (int i = 0; i < b.size(); i++) {
		fs << b[i] << " ";
	}
	return b;
}

type c(type x, int k, vector<type> x_grid) {
	type result = 1;
	for (int j = 0; j < x_grid.size(); j++) {
		if (j != k) {
			result *= (x - x_grid[j]) / (x_grid[k] - x_grid[j]);
		}
	}
	return result;
}

type Lagrange(type x, vector<type> y_grid, vector<type> x_grid) {
	type result = 0;
	for (int k = 0; k < y_grid.size(); k++) {
		result += c(x, k, x_grid) * y_grid[k];
	}
	return result;
}

vector<type> Lagrange_interpolation(vector<type> y_grid, vector<type> x_grid, type a, type b) {
	type step = 0.01;
	vector<type> result;
	for (type x = a; x <= b; x += step) {
		result.push_back(Lagrange(x, y_grid, x_grid));
	}
	return result;
}

type Lagrange_eps(vector<type> lagrange, type a, type b) {
	vector<type> func;
	for (type x = a; x <= b; x += 0.01) {
		func.push_back(f(x));
	}
	type eps = 0;
	for (int i = 0; i < lagrange.size(); i++) {
		if (abs(lagrange[i] - func[i]) > eps) {
			eps = abs(lagrange[i] - func[i]);
		}
	}
	return eps;
}

vector<type>run_through_method(vector<type> a, vector<type> b, vector<type> c, vector<type> d, int n) { // Метод прогонки
	vector<type> result(n, 0);
	vector<type> alpha(n-1, 0), beta(n, 0);
	
	for (int i = 0; i < n; i++) { //прямой ход
		if (i == 0) {
			alpha[0] = -1.0 * c[0] / b[0];
			beta[0] = d[0] / b[0];
		}
		else if (i == n-1) {
			beta[i] = (d[i] - beta[i - 1] * a[i]) / (b[i] + alpha[i - 1] * a[i]);
		}
		else {
			alpha[i] = -1.0 * c[i] / (b[i] + alpha[i - 1] * a[i]);
			beta[i] = (d[i] - beta[i - 1] * a[i]) / (b[i] + alpha[i - 1] * a[i]);
		}
	}

	for (int i = n - 1; i >= 0; i--) { // обратный ход
		if (i == n - 1) {
			result[i] = beta[i];
		}
		else {
			result[i] = alpha[i] * result[i + 1] + beta[i];
		}
	}
	
	return result;
}

vector<vector<type>>Cubic_splain_coefficients(vector<type> y_grid, vector<type> x_grid)
{
	int n = x_grid.size() - 1;
	vector<type>a(n, 0), b(n, 0), c(n + 1, 0), d(n, 0), h(n, 0), g(n, 0);
	vector<type> A(n-1, 0), B(n-1, 0), C(n-1, 0), D(n-1, 0);
	vector<type> others_c(n - 1, 0);
	vector<vector<type>> coeff;

	for (int i = 0; i < n; i++) {
		h[i] = x_grid[i + 1] - x_grid[i];
		g[i] = (y_grid[i + 1] - y_grid[i]) / h[i];
	}

	if (n == 2) {
		c[n - 1] = (3 * (g[1] - g[0])) / (2 * (h[0] + h[1]));
	}
	else {
		for (int i = 0; i < n-1; i++) {
			if (i == 0) {
				B[i] = 2 * (h[i] + h[i + 1]);
				C[i] = h[i+1];
				D[i] = 3 * (g[i + 1] - g[i]);
			}
			else if (i == n - 2) {
				B[i]= 2 * (h[i] + h[i + 1]);
				A[i] = h[i];
				D[i] = 3 * (g[i + 1] - g[i]);
			}
			else {
				A[i] = h[i];
				B[i] = 2 * (h[i] + h[i + 1]);
				C[i] = h[i + 1];
				D[i] = 3 * (g[i + 1] - g[i]);
			}

		}

		others_c = run_through_method(A, B, C, D, n - 1);
		for (int i = 0; i < others_c.size(); i++) {
			c[i + 1] = others_c[i];
		}
	}

	for (int i = 0; i < n; i++) {
		a[i] = y_grid[i];
		b[i] = g[i] - (((c[i + 1] + 2 * c[i]) * h[i]) / 3);
		d[i] = (c[i + 1] - c[i]) / (3 * h[i]);
	}

	c.pop_back();

	coeff.push_back(a);
	coeff.push_back(b);
	coeff.push_back(c);
	coeff.push_back(d);
	coeff.push_back(h);
	coeff.push_back(g);

	return coeff;
}

type Cubic_splain(type x, vector<type> x_grid, vector<vector<type>> coeff, int i)
{
	int n = x_grid.size() - 1;
	vector<type>a(n, 0), b(n, 0), c(n, 0), d(n, 0), h(n, 0), g(n, 0);
	a = coeff[0];
	b = coeff[1];
	c = coeff[2];
	d = coeff[3];
	h = coeff[4];
	g = coeff[5];
	type s;
	s = a[i] + b[i] * (x - x_grid[i]) + c[i] * (pow(x - x_grid[i], 2)) + d[i] * (pow(x - x_grid[i], 3));

	return s;
}

vector<type> Cubic_splain_interpolation(vector<type> x_grid, vector<vector<type>> coeff)
{
	int s = 0;
	type step = 0.01;
	vector<type> result;
	int n = x_grid.size() - 1;
	
	for (type x = x_grid[0]; x <= x_grid[n]; x += step) {
		if (x > x_grid[s+1]) {
			s += 1;
		}
		result.push_back(Cubic_splain(x, x_grid, coeff, s));
	}

	return result;
}

type Cubic_splain_eps(vector<type> splain, vector<type> x_grid) {
	type step = 0.01;
	vector<type> func;
	int n = x_grid.size() - 1;

	for (int i = 0; i < n; i++) {
		if (i != n - 1) {
			for (type x = x_grid[i]; x < x_grid[i + 1]; x += step) {
				func.push_back(f(x));
			}
		}
		else {
			for (type x = x_grid[i]; x <= x_grid[i + 1]; x += step) {
				func.push_back(f(x));
			}
		}
	}

	type eps = 0;
	for (int i = 0; i < splain.size(); i++) {
		if (abs(splain[i] - func[i]) > eps) {
			eps = abs(splain[i] - func[i]);
		}
	}

	return eps;
}


int main() {
	setlocale(LC_ALL, "Russian");
	ofstream fs1_l("output1_l.txt");
	ofstream fs2_l("output2_l.txt");
	ofstream fs1_s("output1_s.txt");
	ofstream fs2_s("output2_s.txt");
	type a, b;
	vector<type> x_grid, y_grid, result_Lagrange, result_splain;
	vector<vector<type>> coeff;
	coeff.resize(6);
	set_boundaries(1, a, b); // первый параметр - номер теста, для которого задаются границы
	int n1 = 10; // узлов в равномерной сетке
	int	n2 = 10; // узлов в чебышевской сетке

	cout << " ======= Табличное задание функции f(x) (равномерная сетка) =======\n\n";
	x_grid = regular_grid(a, b, n1);
	y_grid = func_grid(x_grid);
	cout << "x: "; print(x_grid);
	cout << "\nf: "; print(y_grid);
	result_Lagrange = Lagrange_interpolation(y_grid, x_grid, a, b);
	cout << "\n\nПогрешность интерполяции полиномом лагранжа с помощью равномерной сетки = " << Lagrange_eps(result_Lagrange, a, b) << endl;
	output(result_Lagrange, fs1_l);

	coeff = Cubic_splain_coefficients(y_grid, x_grid);
	result_splain = Cubic_splain_interpolation(x_grid, coeff);
	cout << "\nПогрешность интерполяции кубическим сплайном с помощью равномерной сетки = " << Cubic_splain_eps(result_splain, x_grid) << endl;
	fs1_s << x_grid[0] << " " << x_grid[x_grid.size() - 1] << " ";
	output(result_splain, fs1_s);

	cout << "\n\n\n ======= Табличное задание функции f(x) (чебышевская сетка) =======\n\n";
	x_grid = Chebyshev_grid(a, b, n2);
	y_grid = func_grid(x_grid);
	cout << "x: "; print(x_grid);
	cout << "\nf: "; print(y_grid);
	result_Lagrange = Lagrange_interpolation(y_grid, x_grid, a, b);
	cout << "\n\nПогрешность интерполяции полиномом лагранжа с помощью чебs2ышевской сетки = " << Lagrange_eps(result_Lagrange, a, b) << endl;
	output(result_Lagrange, fs2_l);

	coeff = Cubic_splain_coefficients(y_grid, x_grid);
	result_splain = Cubic_splain_interpolation(x_grid, coeff);
	cout << "\nПогрешность интерполяции кубическим сплайном с помощью чебышевской сетки = " << Cubic_splain_eps(result_splain, x_grid) << endl;
	fs2_s << x_grid[0] << " " << x_grid[x_grid.size() - 1] << " ";
	output(result_splain, fs2_s);

}