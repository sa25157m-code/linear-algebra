namespace LinearAlgebra
{
    public class LinearSystemParameters
    {
        public double Tolerance { get; set; } = 1e-8;
        public int MaxIterations { get; set; } = 100;
    }

    public class LinearSystemResults
    {
        public double[] Solution { get; set; } = new double[0];
    }

    public class LinearSystem
    {
        // ベクトルの内積を計算するヘルパーメソッド
        private static double InnerProduct(double[] u, double[] v)
        {
            double sum = 0.0;
            for (int i = 0; i < u.Length; i++)
            {
                sum += u[i] * v[i];
            }
            return sum;
        }

        public static LinearSystemResults Solve(double[,] A, double[] b, double[] x0, LinearSystemParameters param)
        {
            int n = b.Length;
            double[] x = (double[])x0.Clone();

            // r = Ax - b
            double[] Ax = Matrix.Multiply(A, x);
            double[] r = new double[n];
            for (int i = 0; i < n; i++)
            {
                r[i] = Ax[i] - b[i];
            }

            // p = -r
            double[] p = new double[n];
            for (int i = 0; i < n; i++)
            {
                p[i] = -r[i];
            }

            for (int k = 0; k < param.MaxIterations; k++)
            {
                double rNormSq = InnerProduct(r, r);
                if (Math.Sqrt(rNormSq) < param.Tolerance)
                {
                    break;
                }

                double[] Ap = Matrix.Multiply(A, p);
                double pAp = InnerProduct(p, Ap);
                double alpha = rNormSq / pAp;

                for (int i = 0; i < n; i++)
                {
                    x[i] += alpha * p[i];
                }

                double[] rNew = new double[n];
                for (int i = 0; i < n; i++)
                {
                    rNew[i] = r[i] + alpha * Ap[i];
                }

                double rNewNormSq = InnerProduct(rNew, rNew);
                double beta = rNewNormSq / rNormSq;

                for (int i = 0; i < n; i++)
                {
                    p[i] = -rNew[i] + beta * p[i];
                }

                r = rNew;
            }

            return new LinearSystemResults
            {
                Solution = x
            };
        }
    }
}