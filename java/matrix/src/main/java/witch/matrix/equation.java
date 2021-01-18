package witch.matrix;


import java.util.ArrayList;

public abstract class equation
{
    ArrayList<ArrayList<Double>> A;
    ArrayList<Double> b;
    ArrayList<Double> x;
    double eps;
    int n;
    int counter;

    public equation () {}
    public equation (int _n, ArrayList<ArrayList<Double>> _A, ArrayList<Double> _b, double _eps)
    {
        n = _n;
        A = _A;
        b = _b;
        eps = _eps;
    }

    abstract void solve();

    public ArrayList<ArrayList<Double>> getA ()
    { return A; }
    public ArrayList<Double> getB ()
    { return b; }
    public ArrayList<Double> getX ()
    { return x; }
    public int getCounter ()
    {return counter;}
}
