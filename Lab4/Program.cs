using System;
using System.Collections.Generic;
using System.Linq;

namespace Lab4
{
    class Program
    {
        static double[][] MatrixInverse(double[][] matrix)
        {
            int n = matrix.Length;
            double[][] result = MatrixDuplicate(matrix);

            int[] perm;
            int toggle;
            double[][] lum = MatrixDecompose(matrix, out perm,
              out toggle);
            if (lum == null)
                throw new Exception("Unable to compute inverse");

            double[] b = new double[n];
            for (int i = 0; i < n; ++i)
            {
                for (int j = 0; j < n; ++j)
                {
                    if (i == perm[j])
                        b[j] = 1.0;
                    else
                        b[j] = 0.0;
                }

                double[] x = HelperSolve(lum, b);

                for (int j = 0; j < n; ++j)
                    result[j][i] = x[j];
            }
            return result;
        }
        static double[][] MatrixDuplicate(double[][] matrix)
        {
            // allocates/creates a duplicate of a matrix.
            double[][] result = MatrixCreate(matrix.Length, matrix[0].Length);
            for (int i = 0; i < matrix.Length; ++i) // copy the values
                for (int j = 0; j < matrix[i].Length; ++j)
                    result[i][j] = matrix[i][j];
            return result;
        }
        static double[] HelperSolve(double[][] luMatrix, double[] b)
        {
            // before calling this helper, permute b using the perm array
            // from MatrixDecompose that generated luMatrix
            int n = luMatrix.Length;
            double[] x = new double[n];
            b.CopyTo(x, 0);

            for (int i = 1; i < n; ++i)
            {
                double sum = x[i];
                for (int j = 0; j < i; ++j)
                    sum -= luMatrix[i][j] * x[j];
                x[i] = sum;
            }

            x[n - 1] /= luMatrix[n - 1][n - 1];
            for (int i = n - 2; i >= 0; --i)
            {
                double sum = x[i];
                for (int j = i + 1; j < n; ++j)
                    sum -= luMatrix[i][j] * x[j];
                x[i] = sum / luMatrix[i][i];
            }

            return x;
        }

        static double vectxvect(double[] a, double[] b)
        {
            int l = a.Length;
            double res = 0;
            for (int i = 0; i < l; i++)
            {
                res += a[i] * b[i];
            }
            return res;
        }


        static double[][] MatrixDecompose(double[][] matrix, out int[] perm, out int toggle)
        {
            // Doolittle LUP decomposition with partial pivoting.
            // rerturns: result is L (with 1s on diagonal) and U;
            // perm holds row permutations; toggle is +1 or -1 (even or odd)
            int rows = matrix.Length;
            int cols = matrix[0].Length; // assume square
            if (rows != cols)
                throw new Exception("Attempt to decompose a non-square m");

            int n = rows; // convenience

            double[][] result = MatrixDuplicate(matrix);

            perm = new int[n]; // set up row permutation result
            for (int i = 0; i < n; ++i) { perm[i] = i; }

            toggle = 1; // toggle tracks row swaps.
                        // +1 -greater-than even, -1 -greater-than odd. used by MatrixDeterminant

            for (int j = 0; j < n - 1; ++j) // each column
            {
                double colMax = Math.Abs(result[j][j]); // find largest val in col
                int pRow = j;

                // reader Matt V needed this:
                for (int i = j + 1; i < n; ++i)
                {
                    if (Math.Abs(result[i][j]) > colMax)
                    {
                        colMax = Math.Abs(result[i][j]);
                        pRow = i;
                    }
                }
                // Not sure if this approach is needed always, or not.

                if (pRow != j) // if largest value not on pivot, swap rows
                {
                    double[] rowPtr = result[pRow];
                    result[pRow] = result[j];
                    result[j] = rowPtr;

                    int tmp = perm[pRow]; // and swap perm info
                    perm[pRow] = perm[j];
                    perm[j] = tmp;

                    toggle = -toggle; // adjust the row-swap toggle
                }

                if (result[j][j] == 0.0)
                {
                    // find a good row to swap
                    int goodRow = -1;
                    for (int row = j + 1; row < n; ++row)
                    {
                        if (result[row][j] != 0.0)
                            goodRow = row;
                    }

                    if (goodRow == -1)
                        throw new Exception("Cannot use Doolittle's method");

                    // swap rows so 0.0 no longer on diagonal
                    double[] rowPtr = result[goodRow];
                    result[goodRow] = result[j];
                    result[j] = rowPtr;

                    int tmp = perm[goodRow]; // and swap perm info
                    perm[goodRow] = perm[j];
                    perm[j] = tmp;

                    toggle = -toggle; // adjust the row-swap toggle
                }

                for (int i = j + 1; i < n; ++i)
                {
                    result[i][j] /= result[j][j];
                    for (int k = j + 1; k < n; ++k)
                    {
                        result[i][k] -= result[i][j] * result[j][k];
                    }
                }


            } // main j column loop

            return result;
        }
        static double[][] MatrixCreate(int rows, int cols)
        {
            double[][] result = new double[rows][];
            for (int i = 0; i < rows; ++i)
                result[i] = new double[cols];
            return result;
        }

        static double[] ArrXVect(double[][] A, double[] v)//multiply matrix on vector
        {
            int n = v.Length;
            double[] res = new double[n];
            for (int i = 0; i < n; i++)
            {
                double tmp = 0;
                for (int j = 0; j < n; j++)
                {
                    tmp += A[i][j] * v[j];
                }
                res[i] = tmp;
            }
            return res;
        }

        static double[] MArrXVect(double[][] A, double[] v)//multiply matrix on vector
        {

            double[][] AT = new double[A[0].Length][];
            for (int i = 0; i < AT.Length; i++)
            {
                AT[i] = new double[A.Length];
            }

            for (int i = 0; i < A.Length; i++)
            {
                for (int j = 0; j < A[0].Length; j++)
                {
                    AT[j][i] = A[i][j];
                }
            }
            int columns = AT[0].Length;
            int rows = AT.Length;

            double[] res = new double[rows];
            for (int i = 0; i < rows; i++)
            {
                double tmp = 0;
                for (int j = 0; j < columns; j++)
                {
                    tmp += AT[i][j] * v[j];
                }
                res[i] = tmp;
            }
            return res;
        }

        static double[][] Transpose(double[][] a)
        {
            double[][] at = new double[a[0].Length][];
            for (int i = 0; i < at.Length; i++)
            {
                at[i] = new double[a.Length];
            }
            for (int i = 0; i < a.Length; i++)
            {
                for (int j = 0; j < a[0].Length; j++)
                {
                    at[j][i] = a[i][j];
                }
            }
            return at;
        }

        static Tuple<double[], double[]> Lab4(double[][] a, double[] b, int[] basis, double[] c)
        {
            
                //1
                var Abasis = new double[a.Length][];
                for (int i = 0; i < a.Length; i++)
                    Abasis[i] = new double[basis.Length];
                var AbInv = new double[a.Length][];
                for (int i = 0; i < a.Length; i++)
                    AbInv[i] = new double[basis.Length];
                for (int i = 0; i < Abasis.Length; i++)
                {
                    for (int j = 0; j < Abasis[0].Length; j++)
                    {
                        Abasis[i][j] = a[j][basis[i] - 1];
                    }
                }
                AbInv = MatrixInverse(Abasis);
                //2
                var Cbasis = new double[basis.Length];
                for (int i = 0; i < basis.Length; i++)
                    Cbasis[i] = c[basis[i] - 1];
                var y = new double[Cbasis.Length];
                y = ArrXVect(AbInv, Cbasis);
                //3
                var tmpH = new double[b.Length];
                tmpH = MArrXVect(AbInv, b);
                var Kappa = new double[a[0].Length];

                Console.WriteLine();
                for (int i = 0; i < Kappa.Length; i++)
                {
                    for (int j = 0; j < basis.Length; j++)
                    {
                        if (i == basis[j] - 1)
                        {
                            Kappa[i] = tmpH[j];
                            //H[i] = 0;
                        }
                    }
                }
                var f = false;
                for (int i = 0; i < Kappa.Length; i++)
                    if (Kappa[i] < 0)
                        f = true;
                //4
                if (f)
                {
                    var j1 = 0;
                    for (int i = 0; i < Kappa.Length; i++)
                    {
                        if (Kappa[i] < 0)
                        {
                            j1 = i;
                        }
                    }
                    var delta = AbInv[0];
                    //5
                    var muA = new double[a[0].Length - b.Length][];
                    for (int i = 0; i < muA.Length; i++)
                    {
                        muA[i] = new double[a.Length];
                    }
                    var atranspose = Transpose(a);
                    var no = new List<int>();
                    for (int i = 0; i < a[0].Length; i++)
                    {
                        no.Add(i + 1);
                    }
                    for (int i = 0; i < basis.Length; i++)
                    {
                        for (int j = 0; j < no.Count; j++)
                        {
                            if (basis[i] == no[j])
                            {
                                no.Remove(basis[i]);
                            }
                        }
                    }
                    for (int i = 0; i < no.Count; i++)
                    {
                        if (i == no[i] - 1)
                        {
                            muA[i] = atranspose[i];
                        }

                    }
                    var mu = new double[a[0].Length - b.Length];
                    for (int i = 0; i < mu.Length; i++)
                    {
                        mu[i] = vectxvect(delta, muA[i]);
                        if (mu[i] >= 0)
                        {
                            throw new Exception("taks isnt compatible");
                        }
                    }
                    //6
                    var sigma = new double[a[0].Length - basis.Length];
                    for (int i = 0; i < mu.Length; i++)
                    {
                        sigma[i] = (c[i] - vectxvect(atranspose[i], y)) / mu[i];
                    }
                    //7
                    var sigma0 = sigma.Min();
                    var tmpj = -1;
                    for (int i = 0; i < sigma.Length; i++)
                    {
                        if (sigma[i] == sigma0)
                        {
                            tmpj = i + 1;
                        }
                    }
                    for (int i = 0; i < basis.Length; i++)
                    {
                        if (basis[i] == j1)
                        {
                            basis[i] = tmpj;
                        }
                    }
                    return Lab4(a, b, basis, c);

                }
                else return Tuple.Create(Kappa, y);

            
            //в пункте 6 я скорее всего объебался с тем, как идет цикл

            //я зациклил через вайл тру, не уверен, что правильно

        }
        static void Main(string[] args)
        {
            double[] c = new double[] { -4, -3, -7, 0, 0 };
            double[][] a = new double[2][]
            {
                new double[]{-2, -1, -4, 1, 0},
                new double[]{-2, -2, -2, 0, 1}
            };
            int[] B = new int[] { 4, 5 };
            double[] b = new double[] { -1, -1.5 };
            try
            {

                var (resH, resY) = Lab4(a, b, B, c);
                for (int i = 0; i < resH.Length; i++)
                {
                    Console.WriteLine(resH[i]);
                }
                for (int i = 0; i < resY.Length; i++)
                {
                    // Console.WriteLine(resY[i]);
                }
            }
            catch (Exception e)
            {
                Console.WriteLine(e.Message);
            }
        }
    }
}
