using System;
using System.Linq;
using MathNet.Numerics.LinearAlgebra;

namespace Lr6
{
    class Program
    {
        static void Main(string[] args)
        {
            var c = Vector<double>.Build.DenseOfArray(new double[] { -8, -6, -4, -6 });
            var b = Vector<double>.Build.DenseOfArray(new double[] { 2, 3 });
            var x = Vector<double>.Build.DenseOfArray(new double[] { 2, 3, 0, 0 });
            var a = Matrix<double>.Build.DenseOfArray(new double[,]
            {
                {1, 0, 2, 1},
                {0, 1, -1, 2}
            });
            var d = Matrix<double>.Build.DenseOfArray(new double[,]
            {
                {2, 1, 1, 0},
                {1, 1, 0, 0},
                {1, 0, 1, 0},
                {0, 0, 0, 0}
            });
            var supportConstraints = Vector<double>.Build.DenseOfArray(new double[] { 0, 1 });
            var supportConstraintsExtended = Vector<double>.Build.DenseOfArray(new double[] { 0, 1 });
            var res = QuadraticProgramming(c, x, a, d, supportConstraints, supportConstraintsExtended);
            foreach(var item in res){
                Console.WriteLine(item);
            }
        }
        
        private static Vector<double> QuadraticProgramming(Vector<double> vectorC, Vector<double> vectorX, Matrix<double> matrixA,
                                         Matrix<double> matrixD, Vector<double> supportConstraints, Vector<double> supportConstraintsExtended)
        {
            var aBasis = Matrix<double>.Build.DenseOfColumns(supportConstraints.Select(j => matrixA.Column((int)j)));
            var invABasis = aBasis.Inverse();
            var cX = vectorC + (vectorX * matrixD);
            var cXBasis = Vector<double>.Build.DenseOfEnumerable(supportConstraints.Select(j => -cX[(int)j]));
            var uX = cXBasis * invABasis;
            var deltaX = (uX * matrixA) + cX;
            var j0 = deltaX.ToList().FindIndex(x => x < 0 && !supportConstraintsExtended.Contains(x));
            if (j0 == -1)
                return vectorX;
            

            var l = Vector<double>.Build.DenseOfArray(new double[vectorX.Count]);
            l[j0] = 1;

            var test = matrixD.ToColumnArrays();
            var dX = Matrix<double>.Build.DenseOfColumns(supportConstraintsExtended.Select(j => test[(int)j].Where((_, i) => supportConstraintsExtended.Contains(i))));
            var aX = Matrix<double>.Build.DenseOfColumns(supportConstraintsExtended.Select(j => matrixA.Column((int)j)));
            var aXTr = aX.Transpose();

            var topMatrix = dX.Append(aXTr);
            int dimension = dX.RowCount + aX.RowCount - aXTr.RowCount;
            var bottomMatrix = aX.Append(Matrix<double>.Build.Dense(dimension, dimension));
            Matrix<double>[,] testH =
            {
                {topMatrix},
                {bottomMatrix }
            };
            var h = Matrix<double>.Build.DenseOfMatrixArray(testH);

            var topPart = Vector<double>.Build.DenseOfEnumerable(supportConstraintsExtended.Select(j => matrixD[(int)j, j0]).ToArray());
            var bottomPart = matrixA.Column(j0);
            var bXTest = Vector<double>.Build.DenseOfEnumerable(topPart.Concat(bottomPart));
            var bX = bXTest * -1;

            var hInv = h.Inverse();

            var res = hInv * bX;
            for (int i = 0; i < supportConstraintsExtended.Count; i++)
            {
                var j = (int)supportConstraintsExtended[i];
                l[j] = res[j];
            }

            var delta = (l * matrixD) * l;
            var tettaj0 = delta == 0 ? Double.PositiveInfinity : Math.Abs(deltaX[j0]) / delta;

            double[] tetta = new double[vectorX.Count];
            for (int i = 0; i < tetta.Length; i++)
                tetta[i] = Double.PositiveInfinity;
            

            for (int i = 0; i < supportConstraintsExtended.Count; i++)
            {
                var j = (int)supportConstraintsExtended[i];
                if(l[j] < 0)
                    tetta[j] = -vectorX[j] / l[j];
                
                else
                    tetta[j] = Double.PositiveInfinity;
                
            }
            tetta[j0] = tettaj0;

            var q = tetta.Min();
            if(q == Double.PositiveInfinity)
                throw new Exception("Function isnt ranged!");
            

            var jq = Array.IndexOf(tetta, q);
            var newX = vectorX + (l * q);
            var diff = supportConstraintsExtended.Except(supportConstraints);

            if(j0 == jq)
                supportConstraintsExtended = Vector<double>.Build.DenseOfEnumerable(supportConstraintsExtended.Append(jq));
            
            else if(diff.Contains(jq))
                supportConstraintsExtended = Vector<double>.Build.DenseOfEnumerable(supportConstraintsExtended.Where(val => val != jq));
            
            else
            {
                var sup = Array.IndexOf(supportConstraints.ToArray(), jq);
                int jp = 0;
                for (int i = 0; i < diff.Count(); i++)
                {
                    var j = (int)diff.ElementAt(i);
                    var t = (invABasis * matrixA.Column(j))[sup];
                    jp = Array.IndexOf(diff.ToArray(), t);
                }
                if(jp != -1)
                    supportConstraintsExtended = Vector<double>.Build.DenseOfEnumerable(supportConstraintsExtended.Where(val => val != jq));              
                else
                    supportConstraintsExtended[sup] = j0;

                if(jp != -1)
                    supportConstraints[sup] = jp;
                
                else
                    supportConstraints[sup] = j0;
                
            }
            return QuadraticProgramming(vectorC, newX, matrixA, matrixD, supportConstraints, supportConstraintsExtended);
        }
    }
}

