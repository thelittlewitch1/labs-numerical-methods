package witch.matrix;

import java.util.ArrayList;

public class gauss_method extends equation
{
    @Override
    void solve()
    {
        double max;
        ArrayList<Double> curX = new ArrayList<>(n);
        int k = 0; int index;

        while (k < n)
        {
            max = Math.abs(A.get(k).get(k));
            index = k;
            for (int i = k+1; i < n; i++)
            {
                if (Math.abs(A.get(i).get(k)) > max)
                {
                    max = Math.abs(A.get(i).get(k));
                    index = i;
                }
            }

            if (max < eps)
            {
                System.out.printf("Решение получить невозможно из-за нулевого %d-того столбца матрицы A.\n", index);
                return;
            }
            ArrayList<Double> tempAr = A.get(k);
            A.set(k, A.get(index));
            A.set(index, tempAr);
            double tempDo = b.get(k);
            b.set(k, b.get(index));
            b.set(index, tempDo);

            for (int i = k; i<n; i++)
            {
                tempDo = A.get(i).get(k);
                if (Math.abs(tempDo) < eps) continue;
                for (int j = 0; j < n; j++)
                {
                    tempAr = A.get(i);
                    tempAr.set(j, tempAr.get(j)/tempDo);
                    A.set(i,tempAr);
                }
                b.set(i, b.get(i)/tempDo);
                if (i == k) continue;
                for (int j=0; j<n; j++)
                {
                    tempAr = A.get(i);
                    tempAr.set(j, A.get(i).get(j)-A.get(k).get(j));
                    A.set(i,tempAr);
                }
                b.set(i, b.get(i)-b.get(k));
            }
            k++;
        }

        for (int j = n-1; j>=0; j--)
        {
            curX.set(j, b.get(j));
            for (int i = 0; i < j; i++)
            {
                b.set(i, b.get(i)-A.get(i).get(j)*curX.get(j));
            }
        }

        counter = 1;
        x = curX;
    }
}
