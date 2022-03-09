using System;
using System.IO;
using System.Collections.Generic;

namespace Lab1
{
    class Program
    {
        const string matrixPath = @"C:\Users\Dell\Desktop\MOIU\Lab1\bin\Debug\net5.0\Matrix.txt";//way to file
        const string vectPath = @"C:\Users\Dell\Desktop\MOIU\Lab1\bin\Debug\net5.0\Vect.txt";//way to file
        const string posPath = @"C:\Users\Dell\Desktop\MOIU\Lab1\bin\Debug\net5.0\Pos.txt";//way to file
        const string invPath = @"C:\Users\Dell\Desktop\MOIU\Lab1\bin\Debug\net5.0\Inv.txt";//way to file
        const string nPath = @"C:\Users\Dell\Desktop\MOIU\Lab1\bin\Debug\net5.0\N.txt";//way to file

        static void PrintArray(double[][] tmp)//funtion for printing arrays
        {
            int n = tmp[0].Length;
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    Console.Write(tmp[i][j]+" ");
                }
                Console.Write("\n");
            }
        }
        static void PrintVector(double[] tmp)//function for printing vectors
        {
            int n = tmp.Length;
            for (int i = 0; i < n; i++)
            {
                Console.Write(tmp[i]+" ");
            }
            Console.WriteLine();
        }
        static double[][] ReplaceColumn(double[][] A, double[]vect, int pos)//function in which the column (position which is equal to pos)
                                                                      //of the array is replaced by a column vector 
        {
            int n = vect.Length;
            for (int i = 0; i < n; i++)
            {
                A[i][pos] = vect[i];
            }
            return A;
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
                    if (i==j)
                    {
                        tmp[i][j] = 1;
                    }
                    else
                        tmp[i][j] = 0;
                }
            }
            return tmp;
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
            double[][] res= CreateEmptyArr(n);
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    res[i][j] = Q[i][pos] * A[pos][j];
                    if (i!= pos)
                    {
                        res[i][j] += A[i][j];
                    }
                }
            }
            return res;
        }

        static bool check(double[][] A, double[][] B)//checks is result matrix is correct
        {
            int n = A[0].Length;
            double[][] res = CreateEmptyArr(n);
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    for (int k = 0; k < n; k++)
                    {
                        res[i][j] += A[i][k] * B[k][j];
                    }
                }
            }
            double[][] e = CreateIdentityArr(n);
            for (int i = 0; i < n; i++)
            {
                for (int j = 0; j < n; j++)
                {
                    if (res[i][j] != e[i][j])
                    {
                        return false;
                    }
                }
            }
            return true;
        }

        static double[][] Get_Colswap(double[][] A, double[] x, int i)//main funcionality
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
            for (int  i = 0;  i < n;  i++)
            {
                res[i] = Convert.ToDouble(vArr[i]);
            }
            return res;
        }
        static double[][] InputMatr(string path, int n)
        {
            string[] lines = File.ReadAllLines(path);
            string tmp = "";
            double[][] res = new double [n][];
            for (int i = 0; i < n; i++)
            {
                res[i] = new double[n];
                tmp += lines[i];
                tmp += " ";
            }
            List<string> parsedArr = new List<string>();
            string[] strTmp = tmp.Split(' ');
            for (int i = 0; i < strTmp.Length; i++)
            {
                if (strTmp[i]!="")
                {
                    parsedArr.Add(strTmp[i]);
                }
            }
            int k = 0;
            while (k < 9)
            {
                for (int i = 0; i < n; i++)
                {
                    for (int j = 0; j < n; j++)
                    {
                        res[i][j] = Convert.ToDouble(parsedArr[k]);
                        k++;
                    }
                }
            }

            return res;           
        }

        static void Main(string[] args)
        {
            int n = Input(nPath);             
            double[][] A = InputMatr(matrixPath, n);
            double[][] Inv = InputMatr(invPath, n);           
            double[] x = Input(vectPath, n); ;
            int i = Input(posPath);
            Console.WriteLine("Source matrix");
            PrintArray(A);
            Console.WriteLine();
            Console.WriteLine("Inversed matrix");
            PrintArray(Inv);
            Console.WriteLine();
            Console.WriteLine("Source vector");
            PrintVector(x);
            A = ReplaceColumn(A, x, i);
            Console.WriteLine();
            Console.WriteLine("Source matrix after replacing");
            PrintArray(A);
            Console.WriteLine("Position of replacing");
            Console.WriteLine(i+1);
            Console.WriteLine();
            try
            {
                double[][] Res = Get_Colswap(Inv, x, i);
                Console.WriteLine();
                Console.WriteLine("Result matrix:");
                PrintArray(Res);
                if (check(Res, A))
                {
                    Console.WriteLine();
                    Console.WriteLine("Result matrix is correct");
                }
            }
            catch(Exception e)
            {
                Console.WriteLine();
                Console.WriteLine($"Mistake: {e.Message}");
            }
        }
    }
}