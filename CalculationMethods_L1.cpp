
#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>
#include <typeinfo>

#define type double
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
                cout << data[i][j] << "\t";
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
    
    type cond1() {
        return this->norm1() * this->inverse().norm1();
    }
    type cond_inf() {
        return this->norm_inf() * this->inverse().norm_inf();
    }
    friend matrix operator *(matrix a, matrix b);
    matrix inverse();
};

matrix operator *(matrix a, matrix b) {
    if (a.n == b.n) {
        int n = a.get_size();
        matrix result(n);
        for (int i = 0; i < a.n; i++) {
            for (int j = 0; j < b.n; j++) {
                result.data[i][j] = 0;
                for (int t = 0; t < b.n; t++) {
                    result.data[i][j] += a.data[i][t] * b.data[t][j];
                }
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

void print(vector<type> b) {
    cout << endl;
    for (int i = 0; i < b.size(); i++) cout << b[i] << " ";
    cout << endl;
}

vector<type> gauss_algorythm(matrix A, vector<type> b) {
    // cout << "Метод Гаусса (частичный выбор главного элемента):\n";
    type eps = 1e-14;
    {
        int n = A.get_size();
        vector<type> solution(n, 0.);
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
            if (abs(A[k][k]) < eps) {
                cout << "Матрица вырождена. Решений нет.";
                return solution;
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

        return solution;

    }
}

vector<type> QR_decomposion(matrix A, vector<type> b) {
    bool flag = false;
    int n = A.get_size();
    vector<type> QR_solution(n, 0);
    type c, s;
    vector<matrix> T_vec;
    for (int i = 0; i < n - 1; i++) {
        for (int j = i + 1; j < n; j++) {
            if (A[j][i] == 0) {
                continue;
            }
            c = A[i][i] / (sqrt(pow(A[i][i], 2) + pow(A[j][i], 2)));
            s = A[j][i] / (sqrt(pow(A[i][i], 2) + pow(A[j][i], 2)));
            matrix T(n);
            T.E();
            T[i][i] = c;
            T[i][j] = s;
            T[j][i] = -1.0 * s;
            T[j][j] = c;
            T_vec.push_back(T);
            A = T * A;
            if (A[j][i] < 1e-14) {
                A[j][i] = 0;
            }
        }
    }
    matrix R(n);
    R = A;
    for (int i = 0; i < n; i++) {
        if (R[i][i] == 0) {
            cout << "\nМатрица вырождена. Решений нет.\n";
            return QR_solution;
        }
    }

    matrix T_rotation(n);
    T_rotation.E();
    for (int k = 0; k < T_vec.size(); k++) {
        T_rotation = T_vec[k] * T_rotation;
    }
    matrix Q(n);
    Q = T_rotation.transpose();

    cout << "\n\nQ-matrix:" << endl;
    Q.print();
    cout << "\n\nR-matrix:" << endl;
    R.print();

    vector<type> bt(n, 0);
    bt = T_rotation * b;

    for (int i = n - 1; i > -1; i--) {
        type temp = 0;
        for (int j = i; j < n; j++) {
            temp += QR_solution[j] * R[i][j];
        }
        QR_solution[i] = (bt[i] - temp) / R[i][i];
    }
    return QR_solution;
}

matrix matrix::inverse() {
    matrix result(n);
    matrix matrix_(*this);
    vector<type> A(n), e(n, 0.);
    for (int j = 0; j < n; j++) {
        for (int i = 0; i < n; i++) {
            A[i] = data[i][j];
            e[i] = 0;
        }
        e[j] = 1.;
        A = gauss_algorythm(matrix_, e);
        for (int i = 0; i < n; i++) {
            result[i][j] = A[i];
        }
    }
    return result;
}
vector<type> input(vector<type> b, ifstream &fs) {
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
void residual_norm(matrix A, vector<type> b, vector<type> solution, type& residual_norm1, type& residual_norm2, type& residual_norm_inf) {
    vector<type> b1(b);
    int n = A.get_size();
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            b1[i] -= A[i][j] * solution[j];
        }
    }
    
    cout << "\nВектор невязки:";
    print(b1);
    residual_norm2 = vnorm2(b1);
    residual_norm1 = vnorm1(b1);
    residual_norm_inf = vnorm_inf(b1);
}

type cond1_estimation(matrix A, vector<type> delta, vector<type> b, vector<type> solution) {
    vector<type> delta_solution(A.inverse() * delta);
    type delta_relative_solution = vnorm1(delta_solution) / vnorm1(solution);
    type delta_relative_b = vnorm1(delta) / vnorm1(b);
    return delta_relative_solution / delta_relative_b;
}

type condinf_estimation(matrix A, vector<type> delta, vector<type> b, vector<type> solution) {
    vector<type> delta_solution(A.inverse() * delta);
    type delta_relative_solution = vnorm_inf(delta_solution) / vnorm_inf(solution);
    type delta_relative_b = vnorm_inf(delta) / vnorm_inf(b);
    return delta_relative_solution / delta_relative_b;
}

int main()
{
    setlocale(LC_ALL, "Russian");
    cout << "Тип данных: " << typeid(type).name() << endl;
    string file = "inp5.txt";
    ifstream fs(file);
    matrix A;

    type residual_norm1;
    type residual_norm2;
    type residual_norm_inf;

    // Ввод данных

    A.input(fs);
    int n = A.get_size();

    matrix E(n);
    E.E();

    vector<type> b(n), solution(n), solution1(n);
    b = input(b, fs);
    cout << "\n================ Исходные данные ================\n";
    cout << "\nМатрица системы A:";
    A.print();

    cout << "\nСтолбец b:";
    print(b);
    cout << endl;
    cout << "\n================ Решение системы ================\n";
    //Решаем алгоритмом Гаусса
    solution = gauss_algorythm(A, b);
    cout << "\nРешение, полученное методом Гаусса (частичный выбор):";
    print(solution);

    //Решаем алгоритмом QR-разложения
    solution = QR_decomposion(A, b);
    cout << "\nРешение, полученное методом QR-разложения:";
    print(solution);

    cout << "\n================ Вектор невязки ================\n";
    //Посчитать норму вектора невязки
    residual_norm(A, b, solution, residual_norm1, residual_norm2, residual_norm_inf);
    cout << "\nНорма вектора невязки: \n";
    cout << "Октаэдрическая норма = " << residual_norm1 << "\n";
    cout << "Евклидова норма = " << residual_norm2 << "\n";
    cout << "Кубическая норма = " << residual_norm_inf << "\n";
    cout << "\n================ Решение системы с возмущением ================\n";
    //Вносим возмущение в систему
    vector<type> b1(b);
    vector<type> delta(n, 0.01);

    for (int i = 0; i < n; i++) { delta[i] *= pow(-1, 3); }
    for (int i = 0; i < n; i++) { b1[i] += delta[i]; }
    print(b1);
    //Решаем систему с возмущением
    solution1 = gauss_algorythm(A, b1);
    cout << "\nРешение, полученное методом Гаусса (частичный выбор):";
    print(solution1);

    solution1 = QR_decomposion(A, b1);
    cout << "\nРешение, полученное методом QR-разложения:";
    print(solution1);
    cout << "\n================ Оценка снизу числа обусловленности ================\n";
    //Оценка снизу числа обусловленности
    type max1 = 0;
    type max_inf = 0;
    for (int i = 0; i < n; i++) {
        b1 = b;
        for (int i = 0; i < n; i++) {
            delta[i] = 0.01;
        }
        int random_number = 1 + rand() % 9;
        for (int i = 0; i < n; i++) delta[i] *= pow(-1., rand()) * random_number;
        for (int i = 0; i < n; i++) b1[i] += delta[i];
        solution = gauss_algorythm(A, b1);
        if (cond1_estimation(A, delta, b1, solution) > max1) max1 = cond1_estimation(A, delta, b1, solution);
        if (condinf_estimation(A, delta, b1, solution) > max1) max_inf = condinf_estimation(A, delta, b1, solution);
    }
    cout << "\ncond1(A) >= " << max1;
    cout << "\ncond_inf(A) >= " << max_inf;

    cout << "\n\n================ Подсчет числа обусловленности ================\n";
    //Подсчет числа обусловленности
    cout << "\nОктаэдрическая обусловленность = " << A.cond1() << "\nКубическая обусловленность = " << A.cond_inf() << endl;

    cout << "\n================ Произведение обратной матрицы и исходной ================\n";
    cout << "\nA-1 * A: ";
    (A.inverse() * A).print();

    /*
    // вывод
    string file_out = "result_";
    file_out += typeid(type).name();
    file_out += ".txt";
    ofstream os(file_out);
    if (os.is_open()) {
        os << "Тип данных: " << typeid(type).name() << endl;
        os << "Решение с исходными данными:\n";
        for (int i = 0; i < n; i++) {
            os << solution[i] << ' ';
        }
        os << "\nНорма вектора невязки = ";
        os << residual_norm(A, b, solution);
        os << "\nРешение с отклонением:\n ";
        for (int i = 0; i < n; i++) {
            os << solution1[i] << ' ';
        }
        os << endl;
        os << "A-1 * A\n";
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                os << (A.inverse() * A)[i][j] << "\t";
            }
            os << endl;
        }
        os << "\nОктаэдрическая обсуловленность = " << A.cond1() << "\nКубическая обусловленность = " << A.cond_inf() << endl;
        os.close();*/
    
    

}


