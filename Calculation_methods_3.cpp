#include <iostream>
#include <cmath>
#include <typeinfo>
#include <vector>
#include <fstream>
#include <iomanip>

#define type double
#define EPS 1e-08
using namespace std;

class matrix {
private:
    int n;
    vector<vector<type>> data;
public:
    matrix() : n(0) {};
    matrix(int n) {
        this->n = n;
        data.resize(n, vector<type>(n));
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                data[i][j] = 0;
            }
        }
    }

    void E() {
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (i == j) data[i][j] = 1;
            }
        }
    }

    int get_size() {
        return n;
    }

    void input(ifstream& fs) { // метод ввода данных из файла
        fs >> n;
        data.resize(n, vector<type>(n));
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                fs >> data[i][j];
            }
        }
    }
    void print() { // метод вывода матрицы
        cout << endl;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                cout << setw(12) << data[i][j] << setw(12);
            }
            cout << endl;
        }
    }


    vector<type>& operator [] (int i) {
        return data[i];
    }
    matrix transpose() {
        matrix result(n);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                result.data[i][j] = data[j][i];
            }
        }
        return result;
    }
    matrix operator /(type h) {
        matrix result(n);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                result.data[i][j] = data[i][j] / h;
            }
        }
        return result;
    }
    void swap_rows(int r1, int r2) {
        data[r1].swap(data[r2]);
    }
    type norm1() {
        type max = 0;
        for (int j = 0; j < n; ++j) {
            type sum = 0;
            for (int i = 0; i < n; ++i)
                sum += abs(data[i][j]);
            if (sum > max)
                max = sum;
        }
        return max;
    }
    type norm_inf() {
        type max = 0;
        for (int i = 0; i < n; ++i) {
            type sum = 0;
            for (int j = 0; j < n; ++j)
                sum += abs(data[i][j]);
            if (sum > max)
                max = sum;
        }
        return max;
    }
    friend matrix operator *(matrix a, matrix b);
    friend matrix operator -(matrix a, matrix b);
    friend matrix operator *(matrix a, type h);
    matrix inverse();
};

matrix operator *(matrix a, matrix b) {
    if (a.n == b.n) {
        int n = a.get_size();
        matrix result(n);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                result[i][j] = 0;
                for (int t = 0; t < b.n; t++) {
                    result[i][j] += a[i][t] * b[t][j];
                }
            }
        }
        return result;
    }
}

matrix operator +(matrix a, matrix b) {
    {
        int n = a.get_size();
        matrix result(n);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                result[i][j] = a[i][j] + b[i][j];
            }
        }
        return result;
    }
}

matrix operator -(matrix a, matrix b) {
    {
        int n = a.get_size();
        matrix result(n);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                result[i][j] = a[i][j] - b[i][j];
            }
        }
        return result;
    }
}

vector<type> operator *(matrix a, vector<type> b) {
    int n = b.size();
    vector<type> result(n);
    for (int i = 0; i < b.size(); i++) {
        double c = 0;
        for (int j = 0; j < a.get_size(); j++) {
            c += a[i][j] * b[j];
        }
        result[i] = c;
    }
    return result;
}

matrix operator *(matrix a, type h) {
    int n = a.get_size();
    matrix result(n);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            result[i][j] = a[i][j] * h;
        }
    }
    return result;
}

vector<type> operator -(vector<type> a, vector<type> b) {
    vector<type> result(a.size(), 0);
    for (int i = 0; i < a.size(); i++) {
        result[i] = a[i] - b[i];
    }
    return result;
}

vector<type> operator +(vector<type> a, vector<type> b) {
    vector<type> result(a.size(), 0);
    for (int i = 0; i < a.size(); i++) {
        result[i] = a[i] + b[i];
    }
    return result;
}

vector<type> operator *(vector<type> a, type h) {
    vector<type> result(a.size(), 0);
    for (int i = 0; i < a.size(); i++) {
        result[i] = a[i] * h;
    }
    return result;
}

vector<type> operator /(vector<type> a, type h) {
    vector<type> result(a.size(), 0);
    for (int i = 0; i < a.size(); i++) {
        result[i] = a[i] / h;
    }
    return result;
}

vector<type> operator *(vector<type> a, vector<type> b) {
    vector<type> result(a.size(), 0);
    for (int i = 0; i < a.size(); i++) {
        result[i] = a[i] * b[i];
    }
    return result;
}

void print(vector<type> b) {
    cout << endl;
    for (int i = 0; i < b.size(); i++) cout << b[i] << " ";
    cout << endl;
}

vector<type> input(vector<type> b, ifstream& fs) {
    for (int i = 0; i < b.size(); i++) {
        fs >> b[i];
    }
    return b;
}
type vnorm1(vector<type> b) {
    type temp = 0;
    for (int i = 0; i < b.size(); i++) {
        temp += abs(b[i]);
    }
    return temp;
}
type vnorm_inf(vector<type> b) {
    type max = 0;
    for (int i = 0; i < b.size(); i++) {
        if (max < abs(b[i])) max = abs(b[i]);
    }
    return max;
}

type vnorm2(vector<type> b) {
    type temp = 0;
    for (int i = 0; i < b.size(); i++) {
        temp += pow(b[i], 2);
    }
    return sqrt(temp);
}

vector<type> gauss_algorythm(matrix A, vector<type> b, bool flag) {
    // cout << "Метод Гаусса (частичный выбор главного элемента):\n";
    type eps = 1e-14;
    {
        bool flag1 = false;
        int n = A.get_size();
        vector<type> solution(n, 0.);
        if (flag) {
            for (int j = 0; j < n; j++) {
                for (int i = 0; i < n; i++) {
                    if (A[i][j] != 0) {
                        flag1 = true;
                    }
                }
                if (!flag1) {
                    solution[j] = 1;
                }
                flag1 = false;
            }
        }
        for (int k = 0; k < n; k++) {

            // частичный выбор главного элемента
            int max_i = k;
            for (int i = k; i < n; i++) {
                if (abs(A[i][k]) > abs(A[max_i][k])) max_i = i;
            }
            A.swap_rows(k, max_i);
            type temp = b[k];
            b[k] = b[max_i];
            b[max_i] = temp;
            if (flag) {
                if (abs(A[k][k]) < eps) {
                    return solution;
                }
            }
            // прямой ход Гаусса
            for (int i = k + 1; i < n; i++) {
                type c = A[i][k] / A[k][k];
                for (int j = k; j < n; j++) {
                    A[i][j] = A[i][j] - c * A[k][j];
                }
                b[i] = b[i] - c * b[k];
            }
        }
        // обратный ход Гаусса
        for (int i = n - 1; i > -1; i--) {
            type temp = 0;
            for (int j = i; j < n; j++) {
                temp += solution[j] * A[i][j];
            }
            solution[i] = (b[i] - temp) / A[i][i];
        }


        for (int i = 0; i < n; i++) {
            if (solution[i] == -0) {
                solution[i] = 0;
            }
        }

        return solution;

    }
}

matrix matrix::inverse() {
    matrix result(n);
    matrix matrix_(*this);
    bool flag = false;
    vector<type> A(n), e(n, 0.);
    for (int j = 0; j < n; j++) {
        for (int i = 0; i < n; i++) {
            A[i] = data[i][j];
            e[i] = 0;
        }
        e[j] = 1.;
        A = gauss_algorythm(matrix_, e, flag);
        for (int i = 0; i < n; i++) {
            result[i][j] = A[i];
        }

    }
   
    return result;
}

matrix Hessenberg_form(matrix A) {
    type eps = 1e-14;
    int n = A.get_size();
    type a, b;
    for (int k = 1; k < n - 1; k++) {
        for (int l = k + 1; l < n; l++) {
            if (A[l][k - 1] == 0) {
                continue;
            }
            a = A[k][k - 1] / (sqrt(pow(A[k][k - 1], 2) + pow(A[l][k - 1], 2)));
            b = A[l][k - 1] / (sqrt(pow(A[k][k - 1], 2) + pow(A[l][k - 1], 2)));

            A[k][k - 1] = a * A[k][k - 1] + b * A[l][k - 1];
            A[l][k - 1] = 0;

            for (int i = k; i < n; i++) {
                type buf = A[k][i];
                A[k][i] = a * buf + b * A[l][i];
                A[l][i] = -b * buf + a * A[l][i];
            }

            for (int i = 0; i < n; i++) {
                type buf = A[i][k];
                A[i][k] = a * buf + b * A[i][l];
                A[i][l] = -b * buf + a * A[i][l];
            }

        }

    }
    return A;
}

vector<matrix> QR_decomposion(matrix A) {
    type eps = 1e-14;
    {
        int n = A.get_size();
        vector<matrix> Q_R(2);
        matrix R_inv(n);
        matrix A_temp(n);
        A_temp = A;
        type c, s;
        for (int i = 0; i < n - 1; i++) {
            for (int j = i + 1; j < n; j++) {
                if (A[j][i] == 0) {
                    continue;
                }
                c = A[i][i] / (sqrt(pow(A[i][i], 2) + pow(A[j][i], 2)));
                s = A[j][i] / (sqrt(pow(A[i][i], 2) + pow(A[j][i], 2)));

                for (int k = i; k < n; k++) {
                    type buf = A[i][k];
                    A[i][k] = c * buf + s * A[j][k];
                    A[j][k] = -s * buf + c * A[j][k];
                }

                if (abs(A[j][i]) < eps) {
                    A[j][i] = 0;
                }


            }
        }

        matrix R(n);
        R = A;

       for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (R[i][j] == -0) {
                    R[i][j] = 0;
                }
            }
        }

       R_inv = R.inverse();
       for (int i = 0; i < n; i++) {
           for (int j = 0; j < n; j++) {
               if (R_inv[i][j] == -0) {
                   R_inv[i][j] = 0;
                }
            }
        }
        matrix Q(n);
        Q = A_temp * R_inv;
        Q_R[0] = Q;
        Q_R[1] = R;

        return Q_R;
    }
}

matrix rank_reduction(matrix A) {
    int n = A.get_size();
    matrix result(n - 1);
    for (int i = 0; i < n - 1; i++) {
        for (int j = 0; j < n - 1; j++) {
            result[i][j] = A[i][j];
        }
    }
    return result;
}

matrix diagonal_shift(matrix A, type sigma, bool sign) {
    int n = A.get_size();
    matrix result(n);
    for (int i = 0; i < n; i++) {
        if (sign) {
            A[i][i] += sigma;
        }
        else {
            A[i][i] -= sigma;
        }
    }
    return A;
}
type scalar_multiply(vector<type> a, vector<type> b) {
    type res = 0;
    for (int i = 0; i < a.size(); i++) {
        res += a[i] * b[i];
    }
    return res;
}
vector<type> eigenvalues_searching_by_QR_decomposion(matrix A, type eps, bool shift, int &counter_iterations_shift, int &counter_iterations) {
    matrix Q(A.get_size());
    matrix R(A.get_size());
    matrix Z(A.get_size());
    int n = A.get_size();
    bool zero = true;
    vector<type> eigenvalues(n, 0);
    matrix A_temp(n);
    type sigma;
    type sum;
    for (int r = n - 1; r > 0; --r) {
        A = Hessenberg_form(A);
        if (shift) {
            do {
                sum = 0;
                zero = true;
                A_temp = A;
                sigma = A[r][r];
                A = diagonal_shift(A, sigma, false);

                for (int i = 0; i < A.get_size(); i++) {
                    if (A[r][i] != 0) {
                        zero = false;
                        break;

                    }

                }
                if (!zero) {
                    Q = QR_decomposion(A)[0];
                    R = QR_decomposion(A)[1];
                    A = R * Q;

                    for (int i = 0; i < A.get_size(); i++) {
                        for (int j = 0; j < A.get_size(); j++) {
                            if (abs(A[i][j]) < 1e-14) {
                                A[i][j] = 0;
                            }
                        }
                    }
                }
                else {

                    A = Z;
                }
                A = diagonal_shift(A, sigma, true);

                for (int i = 0; i < r; i++) {
                    sum += abs(A[r][i]);
                }

                counter_iterations_shift += 1;

            } while (sum > eps);
        }
        else {
            do {
                sum = 0;
                zero = true;
                A_temp = A;
                sigma = A[r][r];

                for (int i = 0; i < A.get_size(); i++) {
                    if (A[r][i] != 0) {
                        zero = false;
                        break;

                    }

                }
                if (!zero) {
                    Q = QR_decomposion(A)[0];
                    R = QR_decomposion(A)[1];
                    A = R * Q;

                    for (int i = 0; i < A.get_size(); i++) {
                        for (int j = 0; j < A.get_size(); j++) {
                            if (abs(A[i][j]) < 1e-14) {
                                A[i][j] = 0;
                            }
                        }
                    }
                }
                else {

                    A = Z;
                }

                for (int i = 0; i < r; i++) {
                    sum += abs(A[r][i]);
                }

                counter_iterations += 1;

            } while (sum > eps);
        }

        eigenvalues[r] = A[r][r];

        if (r == 1) {
            eigenvalues[0] = A[0][0];
            break;
        }
        A_temp = rank_reduction(A_temp);
        A = rank_reduction(A);
        Z = rank_reduction(Z);
        A = A_temp;
    }

    return eigenvalues;
}

vector<type> reverse_iterations(matrix A, type lambda, vector<type> x, bool flag) {
    int n = A.get_size();
    bool flag1 = false;
    matrix E(n);
    E.E();
    vector<type> xprev(n, 0);
    vector<type> y(n, 0);
    vector<type> b(n, 0);
    x = x / vnorm2(x);
    matrix M(n);
    M = A - (E * lambda);
    for (int i = 0; i < n; i++) {
        if (M[i][i] == 0) {
            flag1 = true;
        }
    }
    if (!flag) {
        do {
            if (flag1) {
                y = gauss_algorythm(A - (E * lambda), b, true);
                return y;
            }
            else {
                y = gauss_algorythm(A - (E * lambda), x, false);

                xprev = x;
                x = y / vnorm2(y);
            }
        } while ((vnorm2(x - xprev)) > EPS && (vnorm2(x + xprev)) > EPS);
    }
    else {
        cout << "\n\n\nНачальное приближение:\n";
        print(x);
        do {
            lambda = scalar_multiply(A * x, x);
            y = gauss_algorythm(A - (E * lambda), x, false);
            xprev = x;
            x = y / abs(vnorm2(y));
        } while ((vnorm2((A - E * lambda) * x) > EPS) && ((vnorm2(x - xprev)) > EPS && (vnorm2(x + xprev)) > EPS));
        cout << "\n\nНайдено собственное число: " << lambda << endl;
    }
    return x;
}


int main()
{
    setlocale(LC_ALL, "Russian");
    cout << "Тип данных: " << typeid(type).name() << endl;
    string file = "inp1.txt";
    ifstream fs(file);
    matrix A;

    vector<vector<type>> tests;
    tests = { {-0.4, 0.25, -0.3, -0.2}, { 0.4, 0.3, -0.2, 0 }, { -0.05, 0.35, 0.2, -0.1 }, { 1, 0.1, -0.2, -0.5 } };
    vector<type> eigenvector;
    type eps = 1e-8;
    bool flag = false;
    int counter_iterations_shift = 0;
    int counter_iterations = 0;
    // Ввод данных
    A.input(fs);
    int n = A.get_size();
    cout << "\n================ Исходные данные ================\n";
    cout << "\nМатрица линейного оператора A:\n";
    A.print();
    cout << endl << endl;

    cout << "\n\n ======= Матрица линейного оператора A, приведённая к форме Хессенберга: =======\n";
    matrix A_Hessenberg(n);
    A_Hessenberg = Hessenberg_form(A);
    A_Hessenberg.print();

    vector<type> eigenvalues(n, 0);
    vector<type> eigenvalues_shift(n, 0);

    cout << "\n\n\n======= Собственные значения, найденные методом QR-разложения со сдвигом: ======= \n";
    eigenvalues_shift = eigenvalues_searching_by_QR_decomposion(A, eps, true, counter_iterations_shift, counter_iterations);
    print(eigenvalues_shift);
    cout << "\nКоличество итераций:" << counter_iterations_shift;

    cout << "\n\n\n======= Собственные значения, найденные методом QR-разложения без сдвига: ======= \n";
    eigenvalues = eigenvalues_searching_by_QR_decomposion(A, eps, false, counter_iterations_shift, counter_iterations);
    print(eigenvalues);
    cout << "\nКоличество итераций:" << counter_iterations;

    cout << "\n\n\n\n======= Метод обратной итерации =======\n";
    for (int i = 0; i < A.get_size(); i++) {
        eigenvector = reverse_iterations(A, eigenvalues_shift[i], tests[0], flag);
        cout << "\nСобственный вектор, соответствующий собственному числу " << eigenvalues_shift[i];
        print(eigenvector);
    }
    flag = true;
    cout << "\n\n\n========= Модификация метода обратной итерации с использованием отношения Рэлея =========";
    for (vector<type> x : tests) {
        eigenvector = reverse_iterations(A, eigenvalues_shift[0], x, flag);
        cout << "\n\nНайден собственный вектор:\n";
        print(eigenvector);
    }
}