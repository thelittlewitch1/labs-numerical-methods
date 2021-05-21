#include <QCoreApplication>
#include <vector>
#include <iostream>

using namespace std;

vector<vector<double>> mult (vector<vector<double>> A, vector<vector<double>> B)
{
    if (A.at(1).size() != B.size())
    {
        cout << "eror";
        exit (-1);
    }
    unsigned long long l = A.size(), m = A.at(1).size(), n = B.at(1).size();

    vector<vector<double>> C; C.resize(l);


    for (unsigned long long i = 0; i < l; ++i)
    {
        C.at(i).resize(n);
        for (unsigned long long j = 0; j < n; ++j)
        {
            for (unsigned long long r = 0; r < m; ++r)
            {
                C.at(i).at(j) = C.at(i).at(j) + A.at(i).at(r)*B.at(r).at(j);
            }
        }
    }
    return C;
}

double max (vector<vector<double>> M)
{
    double max = -1000;
    for (unsigned long long i = 0; i < M.size(); ++i)
    {
        for (unsigned long long j =0; j < M.at(i).size(); ++j)
        {
            if (M.at(i).at(j) > max)
                max = M.at(i).at(j);
        }
    }
    return max;
}

vector<vector<double>> pow (vector<vector<double>> basis, int exponent)
{
    vector<vector<double>> result = basis;
    for (int i = 0; i < exponent-1; ++i)
    {
        result = mult (result, basis);
    }
    return result;
}

bool isSA (vector<vector<double>> M, double eps)
{
    bool flag = true;
    if (M.size() == M.at(1).size())
    {
        unsigned long long n = M.size();

        for (unsigned long long i = 0; i < n; ++i)
        {
            for (unsigned long long j = 0; j < n; ++j)
            {
                if (M.at(i).at(j) - M.at(j).at(i) > eps)
                {
                    flag = false;
                }
            }
        }
    }
    else
    {
        flag = false;
    }
    return flag;
}

double method (vector<vector<double>> A, vector<vector<double>> *u, double eps)
{
    double lambd0 = 0;
    double lambd1 = max ( mult(A, *u) );

    int k = 0;
    while ( abs(lambd0 - lambd1) > eps )
    {
        lambd0 = lambd1;
        k = k+1;
        lambd1 = max ( mult(pow(A, k+1), *u) )/max ( mult(pow(A, k), *u) );
    }

    *u = mult(pow(A, k+1), *u);
    return lambd1;
}

void print (vector<vector<double>> M)
{
    for (int i = 0; i < M.size(); ++i)
    {
        for (int j = 0; j < M.at(i).size(); ++j)
        {
            cout << M.at(i).at(j) << ", ";
        }
        cout << endl;
    }
}

int main()
{
    vector<vector<double>> A = {{4,1,1},{1,2,1},{1,1,3}}, u = {{1},{0},{1}};
    double eps = 0.001;

    if ( (isSA(A, eps)) && (A.at(1).size() == u.size()) )
    {
        double lambd = method (A, &u, eps);
        cout << "lambd = " << lambd << endl;
        print (u);
    }
}
