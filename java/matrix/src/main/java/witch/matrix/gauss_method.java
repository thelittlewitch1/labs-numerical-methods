package witch.matrix;

import java.util.ArrayList;

class gauss_method extends equation
{
    @Override
    boolean solve()
    {
        double max;
        ArrayList<Double> curX = new ArrayList<>(n);
        ArrayList<ArrayList<Double>> A1 = new ArrayList<>(n);
        ArrayList<Double> b1 = new ArrayList<>(n);
        for (int i = 0; i < n; i++)
        {
            curX.add(0.0);
            b1.add(b.get(i));
            ArrayList<Double> buf = new ArrayList<>(n);
            for (int j = 0; j < n; j++)
                {buf.add(A.get(i).get(j));}
            A1.add(buf);
        }
        int k = 0; int index;

        while (k < n)
        {
            max = Math.abs(A1.get(k).get(k));
            index = k;
            for (int i = k+1; i < n; i++)
            {
                if (Math.abs(A1.get(i).get(k)) > max)
                {
                    max = Math.abs(A1.get(i).get(k));
                    index = i;
                }
            }

            if (max < eps)
            {
                System.out.printf("Решение получить невозможно из-за нулевого %d-того столбца матрицы A.\n", index);
                return false;
            }
            ArrayList<Double> tempAr = A1.get(k);
            A1.set(k, A1.get(index));
            A1.set(index, tempAr);
            double tempDo = b1.get(k);
            b1.set(k, b1.get(index));
            b1.set(index, tempDo);

            for (int i = k; i<n; i++)
            {
                tempDo = A1.get(i).get(k);
                if (Math.abs(tempDo) < eps) continue;
                for (int j = 0; j < n; j++)
                {
                    tempAr = A1.get(i);
                    tempAr.set(j, tempAr.get(j)/tempDo);
                    A1.set(i,tempAr);
                }
                b1.set(i, b1.get(i)/tempDo);
                if (i == k) continue;
                for (int j=0; j<n; j++)
                {
                    tempAr = A1.get(i);
                    tempAr.set(j, A1.get(i).get(j)-A1.get(k).get(j));
                    A1.set(i,tempAr);
                }
                b1.set(i, b1.get(i)-b1.get(k));
            }
            k++;
        }

        for (int j = n-1; j>=0; j--)
        {
            curX.set(j, b1.get(j));
            for (int i = 0; i < j; i++)
            {
                b1.set(i, b1.get(i)-A1.get(i).get(j)*curX.get(j));
            }
        }

        counter = 1;
        x = curX;
        return true;
    }
}
