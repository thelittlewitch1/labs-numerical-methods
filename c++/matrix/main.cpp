#include <QCoreApplication>
#include <iostream>
#include <vector>

using namespace std;

/*Gauss method*/
vector <double> Gauss (vector <vector <double>> A, vector <double> b, const double eps)
{
    double max;
    vector <double> x; x.resize(b.size());
    unsigned long long k, index;
    k = 0;
    while (k < A.size())
    {
        // Поиск строки с максимальным a[i][k]
        max = abs(A.at(k).at(k));
        index = k;
        for (unsigned long long i = k + 1; i < A.at(k).size(); i++)
        {
            if (abs(A[i][k]) > max)
            {
                max = abs(A[i][k]);
                index = i;
            }
        }
        // Перестановка строк
        if (max < eps)
        {
            // нет ненулевых диагональных элементов
            cout << "Решение получить невозможно из-за нулевого столбца ";
            cout << index << " матрицы A" << endl;
            exit(0);
        }
        for (unsigned long long j = 0; j < A.at(k).size(); j++)
        {
            double temp = A[k][j];
            A[k][j] = A[index][j];
            A[index][j] = temp;
        }
        double temp = b[k];
        b[k] = b[index];
        b[index] = temp;
        // Нормализация уравнений
        for (unsigned long long i = k; i < A.size(); i++)
        {
            double temp = A[i][k];
            if (abs(temp) < eps) continue; // для нулевого коэффициента пропустить
            for (unsigned long long j = 0; j < A.size(); j++)
                A[i][j] = A[i][j] / temp;
            b[i] = b[i] / temp;
            if (i == k)  continue; // уравнение не вычитать само из себя
            for (int j = 0; j < A.size(); j++)
                A[i][j] = A[i][j] - A[k][j];
            b[i] = b[i] - b[k];
        }
        k++;
    }
    // обратная подстановка
    for (int j = b.size() - 1; j >= 0; j--)
    {
        x[j] = b[j];
        for (int i = 0; i < j; i++)
            b[i] = b[i] - A[i][j] * x[j];
    }
    return x;
}


/*Gauss-Seidel method*/
bool isConvergent (vector< vector<double>> A)
{
    bool res = true;
    for (unsigned long long i = 0; i < A.size(); i++)
    {
        double sum = 0;
        for (unsigned long long j = 0; j < A.at(i).size(); j++)
        {
            if (j != i)
            {
                sum = sum + abs(A[i][j]);
            }
        }
        if (abs(A[i][i]) < sum)
            res = false;
    }
    return res;
}

double curdif (vector<vector<double>>A, vector<double>x, vector<double> p, unsigned long long i)
{
    double res = 0;
    for (unsigned long long j = 0; j < i; ++j)
    {
        res = res - A[i][j]*x[j];
    }
    for (unsigned long long j = i+1; j < A.at(i).size(); ++j)
    {
        res = res - A[i][j]*p[j];
    }
    return res;
}

bool isDiagContZero (vector<vector<double>> A, double eps)
{
    bool flag = true;
    for (unsigned long long i = 0; i < A.size(); ++i)
    {
        if (abs(A[i][i] < eps))
            flag = false;
    }
    return flag;
}

double curEps (vector<double>x, vector<double>p)
{
    double m = abs (x[0] - p[0]);
    for (unsigned long long i=1; i<x.size(); i++)
    {
        m = max (m, abs(x[i] - p[i]));
    }
    return m;
}

vector<double> GaussSeidel (vector<vector<double>>A, vector<double>b, double eps)
{
    vector<double> x; x.resize(b.size());
    vector<double> p; p.resize(b.size());
    for (unsigned long long i = 0; i < x.size(); ++i)
    {
        x[i] = 1; p[i]=0;
    }
    int counter = 0;
    double epsPrev = 0;
    while (isDiagContZero(A, eps))
    {
        for (unsigned long long i=0; i<A.size(); i++)
        {
            if (A[i][i] < eps)
            {
                if (i>1)
                {
                    vector<double>buf; buf.resize(A.at(i).size());
                    for (unsigned long long j = 0; j < A.at(i).size(); ++j)
                    {
                        buf[j] = A [(i-1)][j];
                        A [(i-1)][j] = A [i][j];
                        A [i][j] = buf[j];
                    }
                }
                else
                {
                    vector<double>buf; buf.resize(A.at(i).size());
                    for (unsigned long long j = 0; j < A.at(i).size(); ++j)
                    {
                        buf[j] = A [i+1][j];
                        A [i+1][j] = A [i][j];
                        A [i][j] = buf[j];
                    }
                }
            }
        }
    }
    if (!isConvergent(A))
    {
        cout << "I cannot solve this matrix equation, because there is no predominance of diagonal elements. Try another matrix." << endl;
        exit (0);
    }
    else
    {
        while (curEps(x,p)>eps)
        {
            if (epsPrev == curEps(x,p))
            {
                cout << "I cannot solve this matrix equation." << endl;
                exit (0);
            }
            else
            {
                epsPrev = curEps(x,p);
                for (unsigned long long i = 0; i < p.size(); ++i)
                {
                    p[i] = x[i];
                    x[i] = ( b[i] + curdif(A,x,p,i) )/A[i][i];
                }
            }
            counter++;
        }
    }
    cout << "number of iterations required to solve the equation: " << counter << endl;
    return x;
}


/*Jacoby method*/
double norm (vector<vector<double>> B)
{
    double m=-10000;
    for (unsigned long long i = 0; i < B.size(); ++i)
    {
        double s = 0;
        for (unsigned long long j = 0; j < B.at(i).size(); ++j)
        {
            s = s + abs(B[i][j]);
        }
        m = max(m, s);
    }
    return m;
}
double norm (vector<double> B)
{
    double m=-10000;
    for (unsigned long long i = 0; i < B.size(); ++i)
    {
        m = max(m, abs(B[i]));
    }
    return m;
}

vector<double> dif (vector<double> L, vector<double> R)
{
    vector<double>res; res.resize(L.size());
    for (unsigned long long i = 0; i < L.size(); ++i)
        {
            res[i] = L[i] - R[i];
        }
    return res;
}

vector<double> Jacobi (vector<vector<double>>A, vector<double>b, double eps)
{
    int counter = 0;
    vector<double>x; x.resize(b.size());
    vector<double>p; p.resize(b.size());
    vector<vector<double>>B; B.resize(A.size());
    for (unsigned long long i = 0; i < A.size(); ++i)
    {
        B.at(i).resize(A.at(i).size());
    }
    for (unsigned long long i = 0; i < A.size(); ++i)
    {
        x[i] = 1;
        p[i] = 0;
        for (unsigned long long j = 0; j < A.at(i).size(); j++)
        {
            if (i == j)
                {B[i][j] = 0;}
            else
                {B[i][j] = A[i][j]/A[i][i];}
        }
    }

    if (norm(B)>1)
    {
        cout << "Cannot solve this equaletion with Jacobi method." << endl;
        exit(-1);
    }
    else {
        double q = eps*((1 - norm(B))/norm(B));
        while ( norm(dif(x, p)) >= q )
        {
            for (unsigned long long i = 0; i < x.size(); ++i)
            {
                p[i] = x[i];
                double buf = b[i];
                for (unsigned long long j = 0; j < A.size(); j++)
                {
                    if (i != j)
                            buf = buf - A[i][j]*p[j];
                }
                x[i] = buf/A[i][i];
            }
            counter++;
        }

        cout << "number of iterations required to solve the equation: " << counter << endl;
        return x;
    }
}

/*main function*/
int main()
{
    char method;
    cout << "Please, chose method. Enter:" << endl << "1 for Gauss,\t2 for Gauss-Saidel,\t3 for Jacobi." << endl << "> ";
    cin >> method;
    double eps = 0.001;
    vector <vector <double>> A = {{7.77, 0.27, - 0.29}, {1.15, - 6.22, 1.77}, {1.05, 4.52, 9.544}};
    vector <double> b = {1.450, 1.050, - 1.310};
    vector <double> x;
    x.resize(b.size());
    //cout << "Please, enter error of the solution:" << endl << "> ";
    //cin >> eps;

    switch (method)
    {
    case '1':
        cout << "You choose Gaussian elimination for Ax=b." << endl;
        x = Gauss (A, b, eps);
        break;
    case '2':
        cout << "You choose Gauss-Seidel method." << endl;
        x = GaussSeidel(A,b,eps);
        break;
    case '3':
        cout << "You choose Jacoby method." << endl;
        x = Jacobi(A,b,eps);
        break;
    default:
        cout << "You enter something wrong" << endl;
        exit (0);
    }


    cout << "Your matrixes: " << endl << "A:" << endl;
    for (unsigned long long i = 0; i < A.size(); ++i) {
        for (unsigned long long j = 0; j < A.at(i).size(); ++j) {
            cout << A.at(i).at(j) << "\t";
        }
        cout << endl;
    }
    cout << "b: " << endl;
    for (unsigned long long i = 0; i < A.size(); ++i) {
            cout << b.at(i) << endl;
    }
    cout << "x:"<< endl;
    for (unsigned long long i = 0; i < A.size(); ++i) {
            cout << x.at(i) << endl;
    }

    cout << endl << endl << "---------" << endl << "Result: " << endl;
    for (int i = 0; i < A.size(); ++i)
    {
        cout << "| ";
        for (int j = 0; j < A.at(i).size(); ++j)
        {
            cout << A[i][j] << " ";
        }
        cout << " |";

        if (i == 0) {cout << " * ";}
            else {cout << "   ";}

        cout << "| " << x[i] << " |";

        if (i == 0) {cout << " = ";}
            else {cout << "   ";}

        double c = 0;
        for (int k=0; k < A.at(i).size(); ++k)
        {
            c = c + A[i][k]*x[k];
        }
        cout << "| " << c << " |";

        if (i == 0) {cout << " ~= ";}
        else {cout << "    ";}
        cout << "| " << b[i] << " |" << endl;
    }
}
