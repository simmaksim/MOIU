using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;

namespace ConsoleApp1
{
    class Program
    {
        static bool isFisrtIteration = true;
        static int j0 = 0;
        const string aPath = @"C:\Users\Dell\Desktop\MOIU\Lab2\bin\Debug\net5.0\a.txt";//way to file
        const string bPath = @"C:\Users\Dell\Desktop\MOIU\Lab2\bin\Debug\net5.0\b.txt";//way to file
        const string cPath = @"C:\Users\Dell\Desktop\MOIU\Lab2\bin\Debug\net5.0\c.txt";//way to file
        const string mPath = @"C:\Users\Dell\Desktop\MOIU\Lab2\bin\Debug\net5.0\m.txt";//way to file
        const string xPath = @"C:\Users\Dell\Desktop\MOIU\Lab2\bin\Debug\net5.0\x.txt";//way to file


        static double[][] MatrixCreate(int rows, int cols)
        {
            double[][] result = new double[rows][];
            for (int i = 0; i < rows; ++i)
                result[i] = new double[cols];
            return result;
        }

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
        static double[][] CreateEmptyArr(int n)//a function that creates an array that is filled with zeros
        {
            double[][] tmp = new double[n][];
            for (int i = 0; i < n; i++)
            {
                tmp[i] = new double[n];
            }
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    tmp[i][j] = 0;
                }
            }
            return tmp;
        }

        static double[][] ReplaceColumn(double[][] A, double[] vect, int pos)//function in which the column (position which is equal to pos)
                                                                             //of the array is replaced by a column vector 
        {
            int n = vect.Length;
            for (int i = 0; i < n; i++)
            {
                A[i][pos] = vect[i];
            }
            return A;
        }

        static double[][] CreateIdentityArr(int n)//a function that creates an identity matrix
        {
            double[][] tmp = new double[n][];
            for (int i = 0; i < n; i++)
            {
                tmp[i] = new double[n];
            }
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    if (i == j)
                    {
                        tmp[i][j] = 1;
                    }
                    else
                        tmp[i][j] = 0;
                }
            }
            return tmp;
        }
        static double[] MArrXVect(double[][] A, double[] v)//multiply matrix on vector
        {
            int columns = A[0].Length;
            int rows = A.Length;
            double[] res = new double[rows];
            for (int i = 0; i < rows; i++)
            {
                double tmp = 0;
                for (int j = 0; j < columns; j++)
                {
                    tmp += A[i][j] * v[j];
                }
                res[i] = tmp;
            }
            return res;
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
        static double[][] ArrxArr(double[][] Q, int pos, double[][] A)
        {
            int n = Q[0].Length;
            double[][] res = CreateEmptyArr(n);
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    res[i][j] = Q[i][pos] * A[pos][j];
                    if (i != pos)
                    {
                        res[i][j] += A[i][j];
                    }
                }
            }
            return res;
        }

        static void PrintArray(double[][] tmp)//funtion for printing arrays
        {
            int n = tmp[0].Length;
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    Console.Write(tmp[i][j] + " ");
                }
                Console.Write("\n");
            }
        }
        static void PrintVector(double[] tmp)//function for printing vectors
        {
            int n = tmp.Length;
            for (int i = 0; i < n; i++)
            {
                Console.Write(tmp[i] + " ");
            }
            Console.WriteLine();
        }
        static double[][] Get_Colswap(double[][] A, double[] x, int i)//main funcionality of fist lab
        {
            int n = A[0].Length;
            double[] l = ArrXVect(A, x);
            Console.WriteLine("Vector L");
            PrintVector(l);
            Console.WriteLine();
            if (l[i] == 0)
            {
                throw new Exception("Smth goes wrong");
            }
            Console.WriteLine("Alright");
            double Li = l[i];
            l[i] = -1;
            for (int j = 0; j < n; j++)
            {
                l[j] = -l[j] / Li;
            }
            Console.WriteLine("Vector l-crown:");
            PrintVector(l);
            double[][] t = CreateIdentityArr(n);
            t = ReplaceColumn(t, l, i);
            Console.WriteLine();
            Console.WriteLine("Matrix T:");
            PrintArray(t);
            return ArrxArr(t, i, A);
        }

        static Tuple<double[], int[]> Simplex(double[] c, double[] x, double[][] a, int[] b)
        {
            double[][] abInv = new double[a.Length][];
            double[][] abInvT = new double[a.Length][];
            double[][] aT = new double[a[0].Length][];
            for (int i = 0; i < a.Length; i++)
            {
                abInv[i] = new double[a.Length];
                abInvT[i] = new double[a.Length];

            }
            for (int i = 0; i < a[0].Length; i++)
            {
                aT[i] = new double[a.Length];
            }

            while (true)
            {
                double[][] ab = CreateEmptyArr(a.Length);
                List<double> cb = new List<double>();
                for (int i = 0; i < a.Length; i++)
                {
                    int colID = b[i] - 1;
                    for (int j = 0; j < a.Length; j++)
                    {
                        ab[i][j] = a[j][colID];
                    }

                }
                foreach (var id in b)
                    cb.Add(c[id - 1]);
                double[] Cb = new double[cb.Count];
                for (int i = 0; i < cb.Count; i++)
                {
                    Cb[i] = cb[i];
                }

                if (isFisrtIteration)
                {
                    abInv = MatrixInverse(ab);
                    isFisrtIteration = false;
                }
                else
                {
                    int pos = 0;
                    for (int i = 1; i < a.Length; i++)
                    {
                        if (b[i] == j0)
                        {
                            pos = i;
                            break;
                        }
                    }
                    double[] swapCol = new double[a.Length];
                    for (int i = 0; i < a.Length; i++)
                        swapCol[i] = a[i][b[pos] - 1];
                    abInv = Get_Colswap(abInv, swapCol, pos);

                }
                for (int i = 0; i < a.Length; i++)
                {
                    for (int j = 0; j < a.Length; j++)
                    {
                        abInvT[j][i] = abInv[i][j];
                    }
                }
                for (int i = 0; i < a[0].Length; i++)
                {
                    for (int j = 0; j < a.Length; j++)
                    {
                        aT[i][j] = a[j][i];
                    }
                }
                double[] U = MArrXVect(abInvT, Cb);


                double[] tmp = MArrXVect(aT, U);
                double[] delta = new double[a[0].Length];
                for (int i = 0; i < delta.Length; i++)
                    delta[i] = 0;
                for (int i = 0; i < tmp.Length; i++)
                {
                    delta[i] += tmp[i];
                }
                for (int i = 0; i < delta.Length; i++)
                {
                    delta[i] -= c[i];
                }
                bool stop = true;
                foreach (var item in delta)
                {
                    for (int i = 0; i < b.Length; i++)
                        if (item == b[i])
                            continue;
                    if (item < 0)
                    {
                        stop = false;
                        break;
                    }
                }
                if (stop)
                {
                    Console.Write("\nAnswer:\nX: ");
                    for (int i = 0; i < x.Length; i++)
                    {
                        Console.Write(x[i] + " ");
                    }
                    Console.Write("\nb: ");
                    for (int i = 0; i < b.Length; i++)
                    {
                        Console.Write(b[i] + " ");
                    }
                    Console.WriteLine();
                    return Tuple.Create(x, b);
                }
                j0 = 0;
                bool temp = false;
                for (int i = 0; i < delta.Length; i++)
                {
                    for (int j = 0; j < b.Length; j++)
                    {
                        if (i + 1 != b[j] && delta[i] < 0)
                        {
                            j0 = i + 1;
                            temp = true;
                        }
                    }
                    if (temp)
                        break;
                }
                double[] col = new double[a.Length];
                for (int i = 0; i < a.Length; i++)
                {
                    col[i] = a[i][j0 - 1];
                }
                double[] z = ArrXVect(abInv, col);
                double[] teta = new double[a.Length];
                for (int i = 0; i < a.Length; i++)
                {
                    if (z[i] <= 0)
                        teta[i] = double.PositiveInfinity;
                    else
                        teta[i] = x[b[i] - 1] / z[i];
                }

                int minIndex = 1;

                for (int i = 1; i < a.Length; i++)
                    if (teta[i] < teta[minIndex - 1])
                        minIndex = i + 1;

                if (teta[minIndex - 1] == double.PositiveInfinity)
                    throw new Exception("The objective function isn't bounded from above on the set of admissible plans");

                var jStar = b[minIndex - 1];
                b[minIndex - 1] = j0;
                for (int i = 0; i < x.Length; i++)
                {
                    
                    if (i + 1 == b[0] || i + 1 == b[1])
                    {
                        continue;
                    }
                    else
                        x[i] = 0;

                }
                x[j0 - 1] = teta[minIndex - 1];

                for (int i = 0; i < b.Length; i++)
                {
                    var value = b[i];
                    if (value == j0)
                        continue;
                    x[value - 1] -= teta[minIndex - 1] * z[i];
                }
            }
            
            return Tuple.Create(x, b);
        }
        static int Input(string path)
        {
            StreamReader sr = new StreamReader(path);
            string vect = sr.ReadLine();
            int res = Convert.ToInt32(vect);
            return res;
        }
        static double[] Input(string path, int n)
        {
            StreamReader sr = new StreamReader(path);
            string vect = sr.ReadLine();
            double[] res = new double[n];
            string[] vArr = vect.Split(' ');
            for (int i = 0; i < n; i++)
            {
                res[i] = Convert.ToDouble(vArr[i]);
            }
            return res;
        }
        static int[] IInput(string path, int n)
        {
            StreamReader sr = new StreamReader(path);
            string vect = sr.ReadLine();
            int[] res = new int[n];
            string[] vArr = vect.Split(' ');
            for (int i = 0; i < n; i++)
            {
                res[i] = Convert.ToInt32(vArr[i]);
            }
            return res;
        }
        static double[][] InputMatr(string path, int n)
        {
            int k = 0, j = 0;
            String inp = File.ReadAllText(path);
            double[][] res = new double[n][];
            for (int i = 0; i < n; i++)
            {
                res[i] = new double[5];
            }
            foreach (var row in inp.Split('\n'))
            {
                j = 0;
                foreach (var col in row.Trim().Split(' '))
                {
                    res[k][j] = double.Parse(col.Trim());
                    j++;
                }
                k++;
            }
            return res;
        }

        static Tuple<ArraySegment<double>, List<int>> StartSimplex(double[][] A, double[] v, int[] b)
        {
            int cols = A[0].Length;
            int rows = A.Length;
            double[] zeroVect = new double[] { 0, 0, 0 };
            double[][] identity = CreateIdentityArr(2);
            double[] ones = new double[] { 1, 1 };
            for (int i = 0; i < v.Length; i++)
            {
                if (v[i] < 0)
                {
                    v[i] *= -1;
                    for (int j = 0; j < cols; j++)
                    {
                        A[i][j] *= -1;
                    }
                }
            }
            double[][] newArr = new double[rows][];
            int tmp1 = 0, tmp2 = 0;
            for (int i = 0; i < rows; i++)
            {
                newArr[i] = new double[cols + identity[0].Length];
            }
            double[][] newArrT = new double[cols+identity[0].Length][];
            for (int i = 0; i < newArrT.Length; i++)
            {
                newArrT[i] = new double[rows];
            }

            for (int i = 0; i < newArr.Length; i++)
            {
                for (int j = 0; j < newArr[0].Length; j++)
                {
                    
                    if (j < cols)
                    {
                        newArr[i][j] = A[i][j];
                    }
                    if (j>=cols)
                    {
                        newArr[i][j] = identity[tmp1][tmp2];
                        tmp2++;
                    }
                }
                tmp1++;
                tmp2 = 0;
            }
            for (int i = 0; i < newArr.Length; i++)
            {
                for (int j = 0; j < newArr[0].Length; j++)
                {
                    newArrT[j][i] = newArr[i][j];
                }
            }
            int[] jb = new int[newArr[0].Length - cols];
            for (int i = 0, j = newArr.Length; i < jb.Length; i++, j--)
            {
                jb[i] = newArr[0].Length - j + 1;
            }
            double[] vectC = new double[newArr[0].Length];
            for (int i = 0; i < vectC.Length; i++)
            {
                for (int j = 0; j < jb.Length; j++)
                {
                    if (i+1 == jb[j])
                    {
                        vectC[i] = -1;
                    }
                }
            }
            double[] x = new double[zeroVect.Length + v.Length];
            for (int i = 0; i < x.Length; i++)
            {
                for (int j = 0; j < v.Length; j++)
                {
                    if (i > (zeroVect.Length - 1))
                    {
                        x[i] = v[j];
                    }
                }
            }
            var (resX, resB) = Simplex(vectC, x, newArr, jb);
            for (int i = resX.Length -1; i > A[0].Length; i--)
            {
                if (resX[i] != 0)
                {
                    throw new Exception("Task incompatible");
                }
            }
            int k = 0;
            for (int i = 0; i < resB.Length; i++)
            {
                if (resB[i] > 0)
                {
                    k = i;
                }
            }

            double[] array = new double[cols - b.Length];
            for (int i = 0; i < cols - 1; i++)
            {
                array[i] = i + 1;
            }
            int[] vct = new int[rows * b.Length];
            for (int i = 0; i < vct.Length; i++)
            {
                for (int j = 0; j < b.Length; j++)
                {
                    vct[i] = (int)newArr[i][j];//!!!!!!!!!!!!!!!!!!!!!!!!!
                } 
            }
            double[][] test = new double[rows][];
            for (int i = 0; i < rows; i++)
            {
                test[i] = new double[resB.Length];
            }
            for (int i = 0; i < test.Length; i++)
            {
                for (int j = 0; j < test.Length; j++)
                {
                    test[j][i] = newArr[i][resB[j] - 1];
                }
            }
            double[][] testInv = new double[rows][];
            for (int i = 0; i < rows; i++)
            {
                testInv[i] = new double[resB.Length];
            }
            testInv = MatrixInverse(test);
            double[][] testInvT = new double[resB.Length][];
            for (int i = 0; i < resB.Length; i++)
            {
                testInvT[i] = new double[rows];
            }
            for (int i = 0; i < rows; i++)
            {
                for (int j = 0; j < resB.Length; j++)
                {
                    testInvT[j][i] = testInv[i][j];
                }
            }
            double[][] l = new double[rows][];
            for (int i = 0; i < rows; i++)
            {
                l[i] = new double[vct.Length];
            }

            for (int i = 0; i < l.Length; i++)
            {
                l[i] = ArrXVect(testInvT, newArrT[vct[i]]);
            }
            var newB = new List<double>();
            var newC = new List<double>();
            var newL = new List<double>();
            var newJb = new List<int>();
            for (int i = 0; i < l.Length; i++)
            {
                for (int j = 0; j < l[0].Length; j++)
                {
                    if (l[i][k] == 0)
                    {
                        newJb = vct.ToList();
                        newJb.RemoveAt(k);
                        newC = vectC.ToList();
                        newC.RemoveAt(k);
                        newB = v.ToList();
                        newB.RemoveAt(k);
                        
                    }
                }
            }

            var res = new ArraySegment<double>(resX, 0, cols);
            return Tuple.Create(res, newJb);
        }
        static void Main(string[] args)
        {
            
            int[] b = new int[] { 0 };
            double[] vect = new double[] { -10, 0 };
            double[][] A =
            {
                new double[] {1, 1, 1},
                new double[] {2, 2, 2}
            };

            try
            {
                var (resX, resB) = StartSimplex(A, vect, b);
                Console.WriteLine();
                Console.WriteLine("X:");
                foreach (var item in resX)
                {
                    Console.WriteLine(item);
                }
                Console.WriteLine();
                Console.WriteLine("B:");
                foreach (var item in resB)
                {
                    Console.WriteLine("");
                    Console.WriteLine(item);
                }
            }
            catch (Exception e)
            {
                Console.WriteLine(e.Message);
            }
        }
    }
}
