package witch.matrix;

public class Main
{
    public static void main(String[] args)
    {
        gauss_method G = new gauss_method();
        if (G.getIsSolved())
        {
            System.out.print("The equation is solved by the Gauss method:\n");
            G.show();
        }

        gauss_seidel_method GS = new gauss_seidel_method();
        if (GS.getIsSolved())
        {
            System.out.print("The equation is solved by the Gauss-Seidel method:\n");
            GS.show();
        }

        jacobi_method J = new jacobi_method();
        if (J.getIsSolved())
        {
            System.out.print("The equation is solved by the Jacobi method:\n");
            J.show();
        }
    }
}
