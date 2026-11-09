using Xunit;

namespace LinearAlgebra.Tests
{
    public class LinearSystemTest
    {
        [Fact]
        public void TestSolve()
        {
            double[,] A = {
                { 4.0, 1.0 },
                { 1.0, 3.0 }
            };

            double[] b = { 1.0, 2.0 };
            double[] x0 = { 0.0, 0.0 };

            var param = new LinearSystemParameters
            {
                Tolerance = 1e-8,
                MaxIterations = 100
            };

            var result = LinearSystem.Solve(A, b, x0, param);

            Assert.NotNull(result.Solution);
            Assert.Equal(2, result.Solution.Length);
        }
    }
}