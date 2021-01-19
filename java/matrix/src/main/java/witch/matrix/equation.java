package witch.matrix;


import java.util.ArrayList;
import java.util.Arrays;

abstract class equation
{
    ArrayList<ArrayList<Double>> A;
    ArrayList<Double> b;
    ArrayList<Double> x;
    double eps;
    int n;
    int counter;
    boolean isSolved;

    public equation ()
    {
        n = 3; eps = 0.001;
        ArrayList<Double> buf1 = new ArrayList<>(Arrays.asList(7.77, 0.27, - 0.29));
        ArrayList<Double> buf2 = new ArrayList<>(Arrays.asList(1.15, - 6.22, 1.77));
        ArrayList<Double> buf3 = new ArrayList<>(Arrays.asList(1.05, 4.52, 9.544));
        A = new ArrayList<>(Arrays.asList(buf1, buf2, buf3));

        b = new ArrayList<>(Arrays.asList(1.450, 1.050, - 1.310));
        x = new ArrayList<>(); x.ensureCapacity(n);
        isSolved = solve();
    }
    public equation (int _n, ArrayList<ArrayList<Double>> _A, ArrayList<Double> _b, double _eps)
    {
        n = _n;
        A = _A;
        b = _b;
        eps = _eps;
        isSolved = solve();
    }

    abstract boolean solve();

    public ArrayList<ArrayList<Double>> getA ()
    { return A; }
    public ArrayList<Double> getB ()
    { return b; }
    public ArrayList<Double> getX ()
    { return x; }
    public int getCounter ()
    {return counter;}
    public boolean getIsSolved ()
    {return isSolved;}

    public void show ()
    {
        System.out.printf("Result.\nNumber of iterations required to solve the equation: %d\nEquation:\n", counter);
        for (int i = 0; i < n; i++)
        {
            System.out.print("| ");
            for (int j = 0; j < n; j++)
                {System.out.printf("%4f ", A.get(i).get(j));}
            System.out.print(" |");

            if (i == 0 ) {System.out.print(" * ");}
                    else {System.out.print("   ");}

            System.out.printf ("| %4f |", x.get(i));

            if (i == 0 ) {System.out.print(" = ");}
                    else {System.out.print("   ");}

            double c = 0;
            for (int j = 0; j < n; j++)
            {
                c = c + A.get(i).get(j)*x.get(j);
            }
            System.out.printf ("| %4f |", c);

            if (i == 0 ) {System.out.print(" ~= ");}
            else {System.out.print("    ");}
            System.out.printf ("| %4f |\n", b.get(i));
        }
        System.out.print("\n\n");
    }
}
