#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

void print (int i, double x1, double x2, double f1, double f2)
{
    int n = 13;
    cout << "|\t" << i
         << "\t| " << setw(n) << x1 << " | " << setw(n) << x2
         << "\t| " << setw(n) << f1 << " | " << setw(n) << f2 << " |" << endl;
}

int main()
{
    double const eps = pow(10, -9);
    double const pi = 	3.141592653589;
    int NIT = 20;
    int i = 0;
    double x1 = 0, x2 = 0;
    double y1 = 1, y2 = 0;
    double f1 = y1*y1 - y2 + 1, f2 = y1 - cos(pi*y2/2);

    cout << "Solving a nonlinear system of equations by the Newton method." << endl
         << "f1(x1,x2) = x1^2 - x2 + 1 = 0" << endl
         << "f2(x1,x2) = x1 - cos(pi/2 x2)" << endl
         << "Accuracy of calculations: eps = " << eps << "." << endl
         << "limit number of iterations: NIT = " << NIT << "." << endl
         << "Initial value: x = (1, 0)." << endl << "--------" << endl << endl;

    cout << "|\ti\t|\tx1\t|\tx2\t|\tf1\t|\tf2\t|" << endl;
    print (i, y1, y2, f1, f2);

    while ( sqrt( (y1-x1)*(y1-x1) + (y2-x2)*(y2-x2) ) >= eps )
    {
        x1 = y1; x2 = y2;

        i = i + 1;

        double h = ( x2 - 1 + x1*x1 - 2*x1*cos(pi*x2/2) )/( -pi*x1*sin(pi*x2/2) - 1 );
        double g = - x1 + cos(pi*x2/2) - pi*sin(pi*x2/2)*h/2;

        y1 = x1 + g;
        y2 = x2 + h;

        f1 = y1*y1 - y2 + 1;
        f2 = y1 - cos(pi*y2/2);


        print (i, y1, y2, f1, f2);

        if (i >= NIT)
        {
            cout << "IER = 2" << endl;
            exit(-1);
        }
    }

    cout << endl << "--------" << endl << "Result: " << endl
         << "x1 = " << y1 << " x1 = " << y2 << endl
         << "f1 = " << f1 << " f2 = " << f2 << endl;

}
