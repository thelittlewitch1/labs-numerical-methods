package witch.matrix;

import java.util.ArrayList;

class gauss_seidel_method extends equation
{
    boolean isCovergent (ArrayList<ArrayList<Double>> C)
    {
        boolean res = true;
        for (int i = 0; i < C.size(); i++)
        {
            double sum = 0;
            for (int j = 0; j < C.get(i).size(); j++)
            {
                if (j != i)
                {sum = sum + Math.abs(C.get(i).get(j));}
            }
            if ( Math.abs(C.get(i).get(i)) < sum )
            {res = false;}
        }
        return res;
    }

    double curdif (ArrayList<ArrayList<Double>> C, ArrayList<Double> q, ArrayList<Double> p, int i)
    {
        double res = 0;
        for (int j = 0; j < i; ++j)
            {res = res - C.get(i).get(j)*q.get(j); }
        for (int j = i+1; j < C.get(i).size(); ++j)
            {res = res - C.get(i).get(j)*p.get(j); }
        return res;
    }

    boolean isDiagContZero (ArrayList<ArrayList<Double>> C, double e)
    {
        boolean flag = false;
        for (int i = 0; i < C.size(); ++i)
        {
            if (Math.abs(C.get(i).get(i)) < e)
            {
                flag = true;
                break;
            }
        }
        return flag;
    }

    double curEps (ArrayList<Double> q, ArrayList<Double>p)
    {
        double m = Math.abs (q.get(0) - p.get(0));
        for (int i=1; i<q.size(); i++)
            {m = Math.max (m, Math.abs(q.get(i) - p.get(i))); }
        return m;
    }

    @Override
    boolean solve()
    {
        counter = 0;
        ArrayList<Double> curX = new ArrayList<>();
        curX.ensureCapacity(n);
        ArrayList<Double> prevX = new ArrayList<>();
        prevX.ensureCapacity(n);
        for (int i = 0; i < n; i++)
        {
            curX.add(1.0);
            prevX.add(0.0);
        }
        double epsPrev = 0;
        while (isDiagContZero(A,eps))
        {
            for (int i = 0; i < n; i++)
            {
                if (A.get(i).get(i) < eps)
                {
                    if (i>1)
                    {
                        ArrayList<Double> buf = A.get(i-1);
                        A.set(i-1, A.get(i));
                        A.set(i, buf);
                    }
                    else
                    {
                        ArrayList<Double> buf = A.get(i+1);
                        A.set(i+1, A.get(i));
                        A.set(i, buf);
                    }
                }
            }
        }
        if (!isCovergent(A))
        {
            System.out.print("I cannot solve this matrix equation, because there is no predominance of diagonal elements. Try another matrix.");
            return false;
        }
        while (curEps(curX, prevX)>eps)
        {
            if (epsPrev == curEps(curX, prevX))
            {
                System.out.print("I cannot solve this matrix equation.");
                return false;
            }
            else
            {
                epsPrev = curEps(curX,prevX);
                for (int i = 0; i < n; ++i)
                {
                    prevX.set(i, curX.get(i));
                    curX.set(i, (b.get(i)+curdif(A, curX, prevX, i))/A.get(i).get(i));
                }
            }
            counter++;
        }
        x = curX;
        return true;
    }
}
