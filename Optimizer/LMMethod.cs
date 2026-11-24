namespace Optimizer;

using LinearAlgebra;
using static LinearAlgebra.Matrix;

///<summary>LM法のパラメータ</summary>
public class LMMethodParameters
{
    private LMMethodParameters() { }

    public static LMMethodParameters Instance { get; } = new();
    public double Tolerance { get; set; } = 1e-6;
    public int MaxIterations { get; set; } = 1000;
    public double IncreaseRatio { get; set; } = 2.0;
    public double DecreaseRatio { get; set; } = 0.8;
    public double TrustLimit { get; set; } = 1e+6;
    // あとで何か追加する
}

///<summary>LM法の結果</summary>
public class LMMethodResults(double[] solution)
{
    public double[] Solution { get; init; } = Copy(solution);
    // あとで何か追加する
    public double[]? Residuals { get; set; }
}

/// <summary>
/// LevenbergMarquardt法（修正ガウス・ニュートン法）
/// </summary>
public static class LMMethod
{
    /// <summary>
    /// LevenbergMarquardt法の実行
    /// </summary>
    /// <param name="func">残差関数: Funcデリゲート，1番目の引数は入力ベクトル，2番目の引数は残差ベクトル</param>
    /// <param name="deriv">残差関数のヤコビ行列: Funcデリゲート，1番目の引数は入力ベクトル，2番目の引数は行が残差ベクトルで列が入力ベクトルのヤコビ行列</param>
    /// <param name="x">解の初期値を入力する</param>
    /// <returns>LevenbergMarquardt法で得られた解</returns>
    public static double[] Solve(Func<double[], double[]> func, Func<double[], double[,]> deriv, double[] x)
    {
        LMMethodParameters param = LMMethodParameters.Instance;
        var res = Solve(func, deriv, x, param);
        return res.Solution;
    }

    /// <summary>
    /// LevenbergMarquardt法の実行
    /// </summary>
    /// <param name="func">残差関数: Funcデリゲート，1番目の引数は入力ベクトル，2番目の引数は残差ベクトル</param>
    /// <param name="deriv">残差関数のヤコビ行列: Funcデリゲート，1番目の引数は入力ベクトル，2番目の引数は行が残差ベクトルで列が入力ベクトルのヤコビ行列</param>
    /// <param name="x">解の初期値を入力する</param>
    /// <param name="param">LevenbergMarquardt法のパラメータ</param>
    /// <returns>LevenbergMarquardt法の結果</returns>
    public static LMMethodResults Solve(Func<double[], double[]> func, Func<double[], double[,]> deriv, double[] x, LMMethodParameters param)
    {
        int n = x.Length;
        double[] currentX = Copy(x);
        double lambda = 1.0;

        double[] f = func(currentX);
        double currentCost = NormSquared(f);

        for (int k = 0; k < param.MaxIterations; k++)
        {
            double[,] J = deriv(currentX);
            double[,] JT = Transpose(J);
            double[,] JTJ = Multiply(JT, J);
            double[] JTf = Multiply(JT, f);

            // 終了条件: || 2 * J^T * f || < Tolerance
            double[] grad = Scale(JTf, 2.0);
            if (Norm(grad) < param.Tolerance)
            {
                break;
            }

            // (J^T * J + lambda * I) * d = -J^T * f
            double[,] A = AddIdentityScaled(JTJ, n, lambda);
            double[] rhs = Scale(JTf, -1.0);
            double[] x0 = new double[n];
            var sysRes = LinearSystem.Solve(A, rhs, x0, new LinearSystemParameters());
            double[] d = sysRes.Solution;
            double[] candidateX = Add(currentX, d);
            double[] candidateF = func(candidateX);
            double candidateCost = NormSquared(candidateF);

            if (candidateCost < currentCost)
            {
                currentX = candidateX;
                f = candidateF;
                currentCost = candidateCost;
                lambda *= param.DecreaseRatio;
            }
            else
            {
                lambda *= param.IncreaseRatio;
            }

            if (lambda > param.TrustLimit)
            {
                break;
            }
        }

        return new LMMethodResults(currentX)
        {
            Residuals = func(currentX)
        };
    }

    private static double NormSquared(double[] v)
    {
        double sum = 0.0;
        for (int i = 0; i < v.Length; i++)
        {
            sum += v[i] * v[i];
        }
        return sum;
    }

    private static double Norm(double[] v) => Math.Sqrt(NormSquared(v));

    private static double[] Scale(double[] v, double s)
    {
        double[] res = new double[v.Length];
        for (int i = 0; i < v.Length; i++)
        {
            res[i] = v[i] * s;
        }
        return res;
    }

    private static double[] Add(double[] a, double[] b)
    {
        double[] res = new double[a.Length];
        for (int i = 0; i < a.Length; i++)
        {
            res[i] = a[i] + b[i];
        }
        return res;
    }

    private static double[,] Transpose(double[,] m)
    {
        int rows = m.GetLength(0);
        int cols = m.GetLength(1);
        double[,] res = new double[cols, rows];
        for (int i = 0; i < rows; i++)
        {
            for (int j = 0; j < cols; j++)
            {
                res[j, i] = m[i, j];
            }
        }
        return res;
    }

    private static double[,] Multiply(double[,] a, double[,] b)
    {
        int rows = a.GetLength(0);
        int inner = a.GetLength(1);
        int cols = b.GetLength(1);
        double[,] res = new double[rows, cols];
        for (int i = 0; i < rows; i++)
        {
            for (int k = 0; k < inner; k++)
            {
                double aik = a[i, k];
                for (int j = 0; j < cols; j++)
                {
                    res[i, j] += aik * b[k, j];
                }
            }
        }
        return res;
    }

    private static double[] Multiply(double[,] m, double[] v)
    {
        int rows = m.GetLength(0);
        int cols = m.GetLength(1);
        double[] res = new double[rows];
        for (int i = 0; i < rows; i++)
        {
            double sum = 0.0;
            for (int j = 0; j < cols; j++)
            {
                sum += m[i, j] * v[j];
            }
            res[i] = sum;
        }
        return res;
    }

    private static double[,] AddIdentityScaled(double[,] m, int n, double lambda)
    {
        double[,] res = Copy(m);
        for (int i = 0; i < n; i++)
        {
            res[i, i] += lambda;
        }
        return res;
    }
}