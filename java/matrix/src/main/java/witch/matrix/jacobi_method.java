package witch.matrix;

import java.util.ArrayList;

class jacobi_method extends equation
{
    double norm2 (ArrayList<ArrayList<Double>> C)
    {
        double m=-10000;
        for (int i = 0; i < C.size(); i++)
        {
            double s = 0;
            for (int j = 0; j < C.get(i).size(); j++)
            {
                s = s + Math.abs(C.get(i).get(j));
            }
            m = Math.max(m,s);
        }
        return m;
    }
    double norm1 (ArrayList<Double> C)
    {
        double m = -10000;
        for (int i = 0; i < C.size(); i++)
        {
            m = Math.max(m, Math.abs(C.get(i)));
        }
        return m;
    }

    ArrayList<Double> dif (ArrayList<Double> L, ArrayList<Double> R)
    {
        ArrayList<Double> res = new ArrayList<>();
        res.ensureCapacity(L.size());
        for (int i = 0; i < L.size(); i++)
        {res.add(L.get(i)-R.get(i));}
        return res;
    }

    @Override
    boolean solve()
    {
        counter = 0;
        ArrayList<Double> curX = new ArrayList<>();
        ArrayList<Double> prevX = new ArrayList<>();
        ArrayList<ArrayList<Double>> B = new ArrayList<>();
        for (int i = 0; i < n; i++)
        {
            ArrayList<Double> buf = new ArrayList<>();
            for (int j = 0; j < n; j++)
            {
                buf.add(0.0);
            }
            B.add(buf);
        }

        for (int i = 0; i < n; i++)
        {
            curX.add(1.0);
            prevX.add(0.0);
            for (int j = 0; j < n; j++)
            {
                if (i == j)
                {
                    ArrayList<Double> buf = B.get(i);
                    buf.set(j, 0.0);
                    B.set(i, buf);
                }
                else
                {
                    ArrayList<Double> buf = B.get(i);
                    buf.set(j, A.get(i).get(j)/A.get(i).get(i));
                    B.set(i, buf);
                }
            }
        }

        if (norm2(B) > 1)
        {
            System.out.print("\nCannot solve this equation with Jacobi method.\n");
            return false;
        }
        double norma = norm2(B);
        double q = eps*((1-norm2(B))/(norm2(B)));
        while (norm1(dif(curX, prevX)) >= q)
        {
            for (int i = 0; i < n; i++)
            {
                prevX.set(i, curX.get(i));
                double buf = b.get(i);
                for (int j = 0; j < n; j++)
                {
                    if (i != j)
                    {
                        buf = buf - A.get(i).get(j)*prevX.get(j);
                    }
                }
                curX.set(i, buf/A.get(i).get(i));
            }
            counter++;
        }
        x = curX;
        return true;
    }
}
