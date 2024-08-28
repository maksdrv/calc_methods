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

    // type cond1() {
    //     return this->norm1() * this->inverse().norm1();
   //  }
    // type cond_inf() {
    //     return this->norm_inf() * this->inverse().norm_inf();
   //  }
    friend matrix operator *(matrix a, matrix b);
    friend matrix operator -(matrix a, matrix b);
    friend matrix operator *(matrix a, type h);
    //matrix inverse();
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

int iteration_estimate(matrix C, vector<type> r, bool norm_type) {
    int n = 1;
    matrix temp = C;
    while (true) {
        if (norm_type) {
            if (temp.norm1() * vnorm1(r) <= 1e-7) {
                break;
            }
            else n++;
        }
        else {
            if (temp.norm_inf() * vnorm_inf(r) <= 1e-7) {
                break;
            }
            else n++;
        }
    }
    return n;
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

vector<type> simple_iteration_method(matrix A, vector<type> b, type eps, bool norm_type, bool logs) {

    for (int i = 0; i < A.get_size(); i++) {
        if (A[i][i] < 0) {
            for (int j = 0; j < A.get_size(); j++) {
                A[i][j] = -A[i][j];
            }
            b[i] = -b[i];
        }
    }

    type tau = 0.;
    if (tau == 0.) {    // поиск наилучшего tau
        type maxA = (abs(A[0][0]));
        for (int i = 1; i < A.get_size(); i++)
            if (maxA < (abs(A[i][i])))
                maxA = (abs(A[i][i]));
        tau = 1. / maxA;
    }

    matrix E(A.get_size());
    
    E.E();
    A = (A * -1) * tau + E;
    
    
    vector<type> y = b * tau;
    vector<type> Xprev(y), Xnext = Xprev;
    /*
    if (logs) {
        vector <type> r = A * Xprev - b;
        cout << "Теоретическое число итераций = " << iteration_estimate(C, r, true);
    }*/
    type norm, vnorm;
    if (norm_type) {
        norm = A.norm1();
    }
    else {
        norm = A.norm_inf();
    }
    if (logs) {
        cout << "\n\nМатрица C:";
        A.print();
        cout << "\nСтолбец y:";
        print(y);
        cout << "\nОктаэдрическая норма C = ";
        cout << A.norm1();
        cout << "\nКубическая норма C = ";
        cout << A.norm_inf() << endl;
    }
    int counter_i = 0; // счётчик итераций
    do {
        Xprev = Xnext;
        Xnext = A * Xprev + y;
        if (norm_type) {
            vnorm = vnorm1(Xnext - Xprev);
        }
        else {
            vnorm = vnorm_inf(Xnext - Xprev);
        }
        counter_i++;
    } while (vnorm > ((1 - norm) / norm * eps));
    if (logs) cout << "\nМетод сошелся за " << counter_i << " итераций\n";
    vector<type> solution = Xnext;
    return solution;
}

vector<type> jacobi_method(matrix A, vector<type> b, type eps, bool norm_type, bool logs) {
    matrix C(A.get_size());
    vector<type> y(A.get_size());
    for (int i = 0; i < A.get_size(); i++) {
        for (int j = 0; j < A.get_size(); j++) {
            if (i != j) {
                C[i][j] = -A[i][j] / A[i][i];
            }
            else C[i][j] = 0;
        }
        y[i] = b[i] / A[i][i];
    }
    type norm, vnorm;
    if (norm_type) {
        norm = C.norm1();
    }
    else {
        norm = C.norm_inf();
    }
    if (logs) {
        cout << "\n\nМатрица C:";
        C.print();
        cout << "\nСтолбец y:";
        print(y);
        cout << "\nОктаэдрическая норма C = ";
        cout << C.norm1();
        cout << "\nКубическая норма C = ";
        cout << C.norm_inf() << endl;
    }
    vector<type> Xprev(A.get_size(), 0), Xnext = Xprev;
    int counter_i = 0;
    do {
        Xprev = Xnext;
        for (int i = 0; i < A.get_size(); i++) {
            type rowsum = 0;
            for (int j = 0; j < A.get_size(); j++) {
                if (i != j) {
                    rowsum += A[i][j] * Xprev[j];
                }
            }
            rowsum -= b[i];
            rowsum *= -1;
            Xnext[i] = rowsum / A[i][i];
        }
        if (norm_type) {
            vnorm = vnorm1(Xnext - Xprev);
        }
        else {
            vnorm = vnorm_inf(Xnext - Xprev);
        }
        counter_i++;
    } while (vnorm > ((1 - norm) / norm * eps));
    vector<type> solution = Xnext;
    if (logs) cout << "\nМетод сошелся за " << counter_i << " итераций\n";
    return solution;
}

int main()
{
    setlocale(LC_ALL, "Russian");
    cout << "Тип данных: " << typeid(type).name() << endl;
    string file = "inp2.txt";
    ifstream fs(file);
    matrix A;
    type residual_norm1;
    type residual_norm2;
    type residual_norm_inf;
    bool norm_type = true; // выбор нормы для матрицы C, true - октаэдрическая, false - кубическая
    // Ввод данных

    A.input(fs);
    int n = A.get_size();
    type eps = 1e-7;
    vector<type> b(n), solution(n), solution1(n), solution2(n);
    b = input(b, fs);
    cout << "\n================ Исходные данные ================\n";
    cout << "\nМатрица системы A:";
    A.print();

    cout << "\nСтолбец b:";
    print(b);
    cout << endl;
    cout << "\n================ Решение системы ================\n";

    cout << "\n ==> Метод простой итерации <== ";
    solution1 = simple_iteration_method(A, b, eps, norm_type, true);
    cout << "\nПолученное решение:";
    print(solution1);

    cout << "\n ==> Метод Якоби <==";
    solution2 = jacobi_method(A, b, eps, norm_type, true);
    cout << "\nПолученное решение:";
    print(solution2);

    cout << "\n================ Решение системы с возмущением ================\n";
    //Вносим возмущение в систему
    for (type d : {1e-4, 1e-7}) {
        vector<type> b1(b);
        cout << "\nПогрешность ~ " << d << "\n";
        vector<type> delta(n, d);

        int random_number = 1 + rand() % 9;
        for (int i = 0; i < n; i++) delta[i] *= pow(-1., rand()) * random_number;
        for (int i = 0; i < n; i++) b1[i] += delta[i];
        cout << "\nСтолбец b с погрешностью:";
        print(b1);

        //Решаем систему с возмущением
        cout << "\n ==> Метод простой итерации <==\n";
        solution1 = simple_iteration_method(A, b1, eps, norm_type, false);
        cout << "\nПолученное решение:";
        print(solution1);


        residual_norm(A, b1, solution1, residual_norm1, residual_norm2, residual_norm_inf);
        cout << "\nНормы вектора невязки: \n";
        cout << "Октаэдрическая норма = " << residual_norm1 << "\n";
        cout << "Евклидова норма = " << residual_norm2 << "\n";
        cout << "Кубическая норма = " << residual_norm_inf << "\n";


        cout << "\n ==> Метод Якоби <==\n";
        solution2 = jacobi_method(A, b1, eps, norm_type, false);
        cout << "\nПолученное решение:";
        print(solution2);

        residual_norm(A, b1, solution1, residual_norm1, residual_norm2, residual_norm_inf);
        cout << "\nНормы вектора невязки: \n";
        cout << "Октаэдрическая норма = " << residual_norm1 << "\n";
        cout << "Евклидова норма = " << residual_norm2 << "\n";
        cout << "Кубическая норма = " << residual_norm_inf << "\n";

    }

}
