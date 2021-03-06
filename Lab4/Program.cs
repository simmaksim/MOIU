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
<<<<<<< HEAD

            //1
            var nonBasis = new List<int>();
            for (int i = 0; i < a[0].Length; i++)
            {
               
                    if (!basis.Contains(i+1))
                    {
                        nonBasis.Add(i + 1);
                        
                    }
                    
                
               
            }
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
            var kappa = new double[a[0].Length];

            Console.WriteLine();
            for (int i = 0; i < kappa.Length; i++)
            {
                for (int j = 0; j < basis.Length; j++)
                {
                    if (i == basis[j] - 1)
                    {
                        kappa[i] = tmpH[j];
                        //H[i] = 0;
                    }
                }
            }
            var f = false;
            for (int i = 0; i < kappa.Length; i++)
                if (kappa[i] < 0)
                    f = true;
            //4
            if (f)//MISTAKE!!!
            {
                var ji = 0;
                var j1 = 0;
                List<double> minusKappa = new List<double>();
                foreach (var item in kappa)
                {
                    if (item < 0)
                    {
                        minusKappa.Add(item);
                    }
                }
                double mKappa = minusKappa.Min();
                for (int i = 0; i < kappa.Length; i++)
                {
                    if (kappa[i] == mKappa)
                    {
                        j1 = i + 1;
                    }
                }
                foreach (var item in minusKappa)
                {
                    if (item == mKappa)
                    {
                        ji = minusKappa.IndexOf(item);//index of otric kappa
                        //j1 = item;
                    }
                }
                //bool check = false;
                
                /*bool f2 = false;
                for (int i = 0; !f2; i++)
                {
                    if (kappa[i] < 0)
=======
            
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
                var kappa = new double[a[0].Length];

                Console.WriteLine();
                for (int i = 0; i < kappa.Length; i++)
                {
                    for (int j = 0; j < basis.Length; j++)
                    {
                        if (i == basis[j] - 1)
                        {
                            kappa[i] = tmpH[j];
                            //H[i] = 0;
                        }
                    }
                }
                var f = false;
                for (int i = 0; i < kappa.Length; i++)
                    if (kappa[i] < 0)
                        f = true;
                //4
                if (f)
                {
                    var j1 = 0;
                    bool f2 = false;
                    for (int i = 0; !f2; i++)
                    {
                        if (kappa[i] < 0)
                        {
                            f2 = true;
                            j1 = i+1;
                        }
                    }
                    var delta = AbInv[0];
                    for (int i = 0; i < AbInv[0].Length; i++)
                    {
                        delta[i] = AbInv[i][0];
                    }
                    //5
                    var muA =new List<double[]>();
                    //for (int i = 0; i < muA.Length; i++)
                    //{
                    //    muA[i] = new double[a.Length];
                    //}
                    var atranspose = Transpose(a);
                    var no = new List<int>();
                    for (int i = 0; i < a[0].Length; i++)
>>>>>>> 76bae06e93872da1482759d7a4835a7ed1bb9dc3
                    {
                        f2 = true;
                        j1 = i + 1;
                    }
<<<<<<< HEAD
                }*/
                var delta = AbInv[ji];
                for (int i = 0; i < AbInv[0].Length; i++)
                {
                    delta[i] = AbInv[i][ji];
                }
                //5
                var muA = new List<double[]>();
                //for (int i = 0; i < muA.Length; i++)
                //{
                //    muA[i] = new double[a.Length];
                //}
                var atranspose = Transpose(a);
                var no = new List<int>();
                for (int i = 0; i < a[0].Length; i++)
                {
                    no.Add(i + 1);
                }
                for (int i = 0; i < basis.Length; i++)
                {
                    for (int j = 0; j < no.Count; j++)
=======
                    for (int i = 0; i < basis.Length; i++)
>>>>>>> 76bae06e93872da1482759d7a4835a7ed1bb9dc3
                    {
                        if (basis[i] == no[j])
                        {
<<<<<<< HEAD
                            no.Remove(basis[i]);
                        }
                    }
                }
                for (int i = 0, j = 0; i < 5; i++)
                {
                    if (basis.Contains(i + 1))
                    {
                        continue;
                    }
                    if (i == no[j] - 1)
                    {
                        j++;
                        // muA[i] = add(atranspose[i]);
                        muA.Add(atranspose[i]);
                    }

                }
                var mu = new double[a[0].Length - b.Length];
                var f3 = false;
                for (int i = 0; i < mu.Length; i++)
                    mu[i] = vectxvect(delta, muA[i]);
                for (int i = 0; i < mu.Length; i++)
                {
                    if (mu[i] < 0)
                    {
                        f3 = true;
                    }
                }
                if (!f3)
                {
                    throw new Exception("task incompatible");
                }

                //6
                var sigma = new double[a[0].Length - basis.Length];
                for (int i = 0; i < mu.Length; i++)
                {

                    sigma[i] = (c[nonBasis[i]-1] - vectxvect(atranspose[nonBasis[i]-1], y)) / mu[i];
                }
                //7
                List<double> sigmaZero = new List<double>();
                for (int i = 0; i < sigma.Length; i++)
                {
                    if (mu[i] < 0)
                    {
                        sigmaZero.Add(sigma[i]);
                    }
                }
                var sigma0 = sigmaZero.Min();
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
            else return Tuple.Create(kappa, y);
=======
                            if (basis[i] == no[j])
                            {
                                no.Remove(basis[i]);
                            }
                        }
                    }
                    for (int i = 0, j = 0; i < 5; i++)
                    {
                        if (basis.Contains(i+1))
                        {
                            continue;
                        }
                        if (i == no[j] - 1)
                        {
                            j++;
                            // muA[i] = add(atranspose[i]);
                            muA.Add(atranspose[i]);
                        }

                    }
                    var mu = new double[a[0].Length - b.Length];
                    var f3 = false;
                    for (int i = 0; i < mu.Length; i++)
                        mu[i] = vectxvect(delta, muA[i]);
                    for (int i = 0; i < mu.Length; i++)
                    {
                        if (mu[i] < 0)
                        {
                            f3 = true;
                        }
                    }
                    if (!f3)
                    {
                        throw new Exception("task incompatible");
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
                else return Tuple.Create(kappa, y);

            
            //?? ???????????? 6 ?? ???????????? ?????????? ?????????????????? ?? ??????, ?????? ???????? ????????

            //?? ???????????????? ?????????? ???????? ??????, ???? ????????????, ?????? ??????????????????
>>>>>>> 76bae06e93872da1482759d7a4835a7ed1bb9dc3

        }
        static void Main(string[] args)
        {
            double[] c = new double[] { -3, -3, 2, 1 };
            double[][] a = new double[2][]
            {
                new double[]{-3, 1, 2, 0},
                new double[]{1, -2, 0, 1}
            };
            int[] B = new int[] { 3, 4 };
            double[] b = new double[] { -3, -4 };
            try
            {

                var (resH, resY) = Lab4(a, b, B, c);
                Console.WriteLine("H:");
                for (int i = 0; i < resH.Length; i++)
                {
<<<<<<< HEAD

=======
>>>>>>> 76bae06e93872da1482759d7a4835a7ed1bb9dc3
                    Console.WriteLine(resH[i]);
                }
                Console.WriteLine("Y:");
                for (int i = 0; i < resY.Length; i++)
                {
                    Console.WriteLine(resY[i]);
                }
            }
            catch (Exception e)
            {
                Console.WriteLine(e.Message);
            }
        }
    }
}