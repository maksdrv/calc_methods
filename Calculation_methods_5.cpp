#include <iostream>
#include <cmath>
#include <typeinfo>
#include <vector>
#include <fstream>
#include <iomanip>
#include "Matrix.h"
#define type double
#define eps 1e-6
using namespace std;


type f(type x, int inp) {
    switch (inp) {
    case 1: return (x - 0.1) * (x - 0.22) * (x - 0.55) * (x - 0.7) * (x - 0.75);
    case 2: return sqrt(x + 1) - 1;
    case 3: return 35 * pow(x, 3) - 67 * pow(x, 2) - 3 * x + 3;
    }
}

type df(type x, int inp) {
    return (f(x + eps, inp) - f(x, inp)) / eps;
}

void setup(int inp, type& a, type& b, type& step) {
    switch (inp) {
    case 1: {
        a = 0;
        b = 1;
        step = 0.04;
        break;
    };
    case 2: {
        a = -1;
        b = 10;
        step = 0.3;
        break;
    };
    case 3: {
        a = 0;
        b = 1;
        step = 0.6;
        break;
    }
    }
}


void localization(int inp, vector<type>& locals) {
    type a, b, step;
    setup(inp, a, b, step);
    vector<type> x;
    vector<type> F;

    x.push_back(a);
    F.push_back(f(a, inp));

    // таблица значений
    while (true) {
        a += step;
        x.push_back(a);
        F.push_back(f(a, inp));
        if (a > b) break;
    }

    // локализация
    for (size_t i = 1; i < x.size(); i++) {
        if (F[i] * F[i - 1] < 0) {
            locals.push_back(x[i - 1]);
            locals.push_back(x[i]);
        }
    }
}
type bisection_method(int inp, type a, type b, int& iterations) {
    type x = (a + b) / 2;
    while (abs(a - b) >= 2 * eps) {
        iterations++;
        if (f(a, inp) * f(x, inp) < 0) b = x;
        else  a = x;
        x = (a + b) / 2.;
    }
    if (abs(x) < eps) x = 0;
    return x;
}
type newton_method(int inp, type a, type b, int& iterations) {
    type x0 = (f(a, inp) * b - f(b, inp) * a) / (f(a, inp) - f(b, inp));
    type x = x0;
    type xprev;
    do {
        iterations++;
        xprev = x;
        x = xprev - f(xprev, inp) / df(xprev, inp);
        if (x < a or x > b) {
            cout << endl << "Приближение вышло за границы отрезка. a = " << a << " b = " << b << endl;
            return 1e-10;
        }
    } while (abs(x - xprev) > eps);
    if (abs(x) < eps) x = 0;
    return x;
}
type newton_method_modificated(int inp, type a, type b, int& iterations) {
    type x0 = (f(a, inp) * b - f(b, inp) * a) / (f(a, inp) - f(b, inp));
    type x = x0;
    type xprev;
    do {
        iterations++;
        xprev = x;
        x = xprev - f(xprev, inp) / df(xprev, inp);
        if (x < a or x > b) {
            type alpha = 0.9;
            while (x < a or x > b) {
                x = xprev - alpha * f(xprev, inp) / df(xprev, inp);
                alpha += -0.01;
            }
        }
    } while (abs(x - xprev) > eps);
    if (abs(x) < eps) x = 0;
    return x;
}


void setup_system(int inp_system, vector<vector<type>>& X_0, type& L1, type& L2, type& h1, type& h2) {
    switch (inp_system)
    {
    case 1: {
        X_0 = { { 2, 1} ,{-2, -1} };
        L1 = 10;
        L2 = 10;
        h1 = 1;
        h2 = 1;
        cout << "Система имеет 2 решения";
        break;
    };
    case 2: {
        X_0 = { { -5, 2} ,{2, -5}, {1, 5}, {5, 1} };
        L1 = 10;
        L2 = 10;
        h1 = 1;
        h2 = 1;
        cout << "Система имеет 4 решения";
        break;
    }
    }
}

void convergence_area(int inp_system, type L1, type L2, type h1, type h2, vector<type>& w1, vector<type>& w2 ) {
    int N1 = (2 * L1) / h1;
    int N2 = (2 * L2) / h2;

    for (int i = 0; i <= N1; i++) {
        w1.push_back(-1.0 * L1 + i * h1);
    }
    for (int i = 0; i <= N2; i++) {
        w2.push_back(-1.0 * L2 + i * h2);
    }
}

vector<type> F(vector<type> X, int inp_system) {
    type eps_f = 1e-12;
    type x1, x2;
    x1 = X[0];
    x2 = X[1];
    vector<type> result(2);

    switch (inp_system)
    {
    case 1: {
        result[0] = pow(x1, 2) - pow(x2, 2) - 15;
        result[1] = x1 * x2 + 4;
        break;
    };
    case 2: {
        result[0] = pow(x1, 2) + pow(x2, 2) + x1 + x2 - 8;
        result[1] = pow(x1, 2) + pow(x2, 2) + x1 * x2 - 7;
        break;
    }
    }

    if (abs(result[0]) < eps_f) {
        result[0] = 0;
    }

    if (abs(result[1]) < eps_f) {
        result[1] = 0;
    }

    return result;
}


matrix Jacobi_analitical(vector<type> X, int inp_system) {
    type x1, x2;
    x1 = X[0];
    x2 = X[1];
    matrix J(2);

    switch (inp_system)
    {
    case 1: {
        J[0][0] = 2 * x1;
        J[0][1] = -2 * x2;
        J[1][0] = x2;
        J[1][1] = x1;
        break;
    };
    case 2: {
        J[0][0] = 2 * x1 + 1;
        J[0][1] = 2 * x2 + 1;
        J[1][0] = 2 * x1 + x2;
        J[1][1] = 2 * x2 + x1;
        break;
    }
    }
    return J;
}

matrix Jacobi_numerical(vector<type> X, int inp_system) {
    type x1, x2;
    x1 = X[0];
    x2 = X[1];
    matrix J(2);

    J[0][0] = (F({ x1 + eps, x2 }, inp_system)[0] - F(X, inp_system)[0]) / eps;
    J[0][1] = (F({ x1, x2 + eps }, inp_system)[0] - F(X, inp_system)[0]) / eps;
    J[1][0] = (F({ x1 + eps, x2 }, inp_system)[1] - F(X, inp_system)[1]) / eps;
    J[1][1] = (F({ x1, x2 + eps }, inp_system)[1] - F(X, inp_system)[1]) / eps;

    return J;
}

vector<type> Newton_method_system(int inp_system, vector<type> X_0, int L1, int L2, int& iterations, bool analitical) {
    vector<type> X_prev(2), X(2), zero_v(2);
    matrix J(2), J_inv(2), zero_m(2);
    X = X_0;
    do {
        iterations++;
        X_prev = X;

        if (analitical) {
            J = Jacobi_analitical(X_prev, inp_system);
        }
        else {
            J = Jacobi_numerical(X_prev, inp_system);
        }
    
        J_inv = inverse(J);
        if (J_inv[0] == zero_m[0] and J_inv[1] == zero_m[1]) {
            return zero_v;
            break;
        }
        
        X = X_prev - (J_inv * F(X_prev, inp_system));
        if (abs(X[0]) > L1 or abs(X[1]) > L2) {
            type alpha = 0.9;
            while (abs(X[0]) > L1 or abs(X[1]) > L2) {
                X = X_prev - alpha * (J_inv * F(X_prev, inp_system));
                alpha += -0.01;
            }
        }
    } while (abs(vnorm2((X - X_prev))) > eps);

    if (abs(X[0]) < eps) {
        X[0] = 0;
    }
    if (abs(X[1]) < eps) {
        X[1] = 0;
    }

    return X;
}


int main()
{
    setlocale(LC_ALL, "Russian");

    cout << "\n=========== Решение одного уравнения ===========\n";
    int inp = 1;
    cout << "\nВведён тестовый пример " << inp << endl;

    vector<type> solutions, locals;
    localization(inp, locals);
    cout << "\nКорней уравнения: " << locals.size() / 2;
    cout << "\n\n\n=== Метод бисекции ===\n";
    int iterations = 0;
    for (size_t i = 1; i <= locals.size(); i += 2) {
        solutions.push_back(bisection_method(inp, locals[i - 1], locals[i], iterations));
    }

    cout << "\nКоличество итераций = " << iterations;
    cout << "\nКорни уравнения: ";
    print(solutions);


    solutions.clear();

    cout << "\n\n\n=== Метод Ньютона ===\n";
    iterations = 0;
    for (size_t i = 1; i <= locals.size(); i += 2) {
        type x = newton_method(inp, locals[i - 1], locals[i], iterations);
        if (x != 1e-10) solutions.push_back(x);
    }

    cout << "\nКоличество итераций = " << iterations;
    cout << "\nКорни уравнения: ";
    print(solutions);

    solutions.clear();

    cout << "\n\n\n=== Метод Ньютона с модификацией ===\n";
    iterations = 0;
    for (size_t i = 1; i <= locals.size(); i += 2) {
        solutions.push_back(newton_method_modificated(inp, locals[i - 1], locals[i], iterations));
    }

    cout << "\nКоличество итераций = " << iterations;
    cout << "\nКорни уравнения: ";
    print(solutions);


    cout << "\n\n\n\n\n=========== Решение системы уравнений методом Ньютона ===========\n\n\n";

    int inp_system = 1;
    cout << "Введён тестовый пример " << inp + 3 << "\n\n";
    type L1, L2, h1, h2;
    vector<vector<type>> X_0;
    setup_system(inp_system, X_0, L1, L2, h1, h2);

    vector<type> solutions_system, zero_v(2);
    cout << "\n\n";
    cout << "\n\n===== Матрица Якоби вычисляется аналитически =====\n\n";

    for (int i = 0; i < X_0.size(); i++) {
        iterations = 0;
        solutions_system = Newton_method_system(inp_system, X_0[i], L1, L2, iterations, true);

        if (solutions_system == zero_v and (F(zero_v, inp_system)[0] != 0 or F(zero_v, inp_system)[1] != 0)) {
            cout << "При заданном начальном приближении метод не сошёлся\n\n";
        }

        else if (solutions_system == zero_v and F(zero_v, inp_system)[0] == 0 or F(zero_v, inp_system)[1] == 0){
            cout << "Решение " << i + 1 << ":   x1 = " << solutions_system[0] << ",   x2 = " << solutions_system[1] << "\n\n";
            }

        else {
            cout << "Решение " << i + 1 << ":   x1 = " << solutions_system[0] << ",   x2 = " << solutions_system[1] << "\t   " << "получено за " << iterations << " итераций" << "\n\n";
        }
    }

    cout << "\n\n==== Матрица Якоби вычисляется численно ====\n\n";

    for (int i = 0; i < X_0.size(); i++) {
        iterations = 0;
        solutions_system = Newton_method_system(inp_system, X_0[i], L1, L2, iterations, false);

        if (solutions_system == zero_v and (F(zero_v, inp_system)[0] != 0 or F(zero_v, inp_system)[1] != 0)) {
            cout << "При заданном начальном приближении метод не сошёлся\n\n";
        }

        else if (solutions_system == zero_v and F(zero_v, inp_system)[0] == 0 or F(zero_v, inp_system)[1] == 0) {
            cout << "Решение " << i + 1 << ":   x1 = " << solutions_system[0] << ",   x2 = " << solutions_system[1] << "\n\n";
        }

        else {
            cout << "Решение " << i + 1 << ":   x1 = " << solutions_system[0] << ",   x2 = " << solutions_system[1] << "\t   " << "получено за " << iterations << " итераций" << "\n\n";
        }
    }




    // Ддя построения диаграммы области сходимости

    /*vector<type> w1, w2;
    vector<type> x_coord, y_coord, iter;
    ofstream c_x("coordinats_x.txt");
    ofstream c_y("coordinats_y.txt");
    ofstream it("iterations.txt");


    convergence_area(inp_system, L1, L2, h1, h2, w1, w2);

    for (int i = 0; i < w1.size(); i++) {
        for (int j = 0; j < w2.size(); j++) {
            iterations = 0;
            solutions_system = Newton_method_system(inp_system, { w1[i], w2[j] }, L1, L2, iterations, true);
            if (F(solutions_system, inp_system)[0] != 0 or F(solutions_system, inp_system)[1] != 0) {
                //cout << "При " << w1[i] << " и " << w2[j] <<  "  метод не сошёлся" << endl;
                iterations = 31;
            }

            else if (solutions_system == zero_v and F(zero_v, inp_system)[0] == 0 and F(zero_v, inp_system)[1] == 0) {
                //cout << "Решение " << i + 1 << ":   x1 = " << solutions_system[0] << ",   x2 = " << solutions_system[1] << endl;
            }

            else {
                //cout << "При " << w1[i] << " и " << w2[j] << "  найдено решение " << solutions_system[0] << " " << solutions_system[1] << "  за " << iterations << " итераций" << endl;
            }

           
            x_coord.push_back(w1[i] - h1 / 2);
            x_coord.push_back(w1[i] + h1 / 2);
            x_coord.push_back(w1[i] + h1 / 2);
            x_coord.push_back(w1[i] - h1 / 2);
            y_coord.push_back(w2[j] - h2 / 2);
            y_coord.push_back(w2[j] - h2 / 2);
            y_coord.push_back(w2[j] + h2 / 2);
            y_coord.push_back(w2[j] + h2 / 2);
            iter.push_back(iterations);
        }
    }

    output(x_coord, c_x);
    output(y_coord, c_y);
    output(iter, it);*/

}
