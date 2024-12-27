using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Numerics;
using System.Reflection;
using System.Reflection.Metadata;
using System.Text;
using static Program;
using static System.Runtime.InteropServices.JavaScript.JSType;

public class Program
{
    public static void Main()
    {
        Program program = new Program();
    }
    #region Level 1
    public long Task_1_1(int n, int k)
    {
        long answer = 0;

        // code here
        // create and use Combinations(n, k);
        // create and use Factorial(n);
        if (n < k || k < 0 || n < 0) return 0;

        answer = Combinations(n, k);
        // end

        return answer;
    }

    public long Combinations(int n, int k)
    {
        return Factorial(n) / (Factorial(k) * Factorial(n - k));
    }
    public long Factorial(int n)
    {
        long result = 1;

        for (int i = 1; i <= n; i++)
        {
            result *= i;
        }

        return result;
    }
    public int Task_1_2(double[] first, double[] second)
    {
        int answer = 0;

        // code here
        // create and use GeronArea(a, b, c);
        double a1 = first[0], b1 = first[1], c1 = first[2];
        double a2 = second[0], b2 = second[1], c2 = second[2];

        if (!((a1 < b1 + c1) && (b1 < a1 + c1) && (c1 < a1 + b1)) || !((a2 < b2 + c2) && (b2 < a2 + c2) && (c2 < a2 + b2))) return -1;

        double sq1 = GeronArea(a1, b1, c1);
        double sq2 = GeronArea(a2, b2, c2);

        if (sq1 > sq2)
            answer = 1;
        else if (sq2 > sq1)
            answer = 2;
        // end

        // first = 1, second = 2, equal = 0, error = -1
        return answer;
    }

    public double GeronArea(double a, double b, double c)
    {
        double p = (a + b + c) / 2;
        double S = Math.Sqrt(p * (p - a) * (p - b) * (p - c));

        return S;
    }

    public int Task_1_3a(double v1, double a1, double v2, double a2, int time)
    {
        int answer = 0;

        // code here
        // create and use GetDistance(v, a, t); t - hours
        double destination1 = GetDistance(v1, a1, time);
        double destination2 = GetDistance(v2, a2, time);

        if (destination1 > destination2)
            answer = 1;
        else if (destination2 > destination1)
            answer = 2;
        // end

        // first = 1, second = 2, equal = 0
        return answer;
    }

    public double GetDistance(double v, double a, double t)
    {
        double S = v * t + a * t * t / 2;

        return S;
    }

    public int Task_1_3b(double v1, double a1, double v2, double a2)
    {
        int answer = 0;

        // code here
        // use GetDistance(v, a, t); t - hours
        int time = 1;

        while (GetDistance(v1, a1, time) > GetDistance(v2, a2, time))
            time++;

        answer = time;
        // end

        return answer;
    }
    #endregion

    #region Level 2
    public void Task_2_1(int[,] A, int[,] B)
    {
        // code here

        // create and use FindMaxIndex(matrix, out row, out column);

        // end
    }

    public void Task_2_2(double[] A, double[] B)
    {
        //code here
        // create and use FindMaxIndex(array);
        // only 1 array has to be changed!
        int maxIndexA = FindMaxIndex(A);
        int maxIndexB = FindMaxIndex(B);

        double[] targetArray = (A.Length - maxIndexA > B.Length - maxIndexB) ? A : B;
        int targetMaxIndex = (targetArray == A) ? maxIndexA : maxIndexB;

        UpdateWithAverage(targetArray, targetMaxIndex);
        //end
    }

    int FindMaxIndex(double[] array)
    {
        int maxIndex = 0;

        for (int i = 1; i < array.Length; i++)
        {
            if (array[i] > array[maxIndex])
            {
                maxIndex = i;
            }
        }
        return maxIndex;
    }

    public void UpdateWithAverage(double[] array, int indexMax)
    {
        double sum = 0;
        int count = 0;

        for (int i = indexMax + 1; i < array.Length; i++)
        {
            sum += array[i];
            count++;
        }

        double average = count > 0 ? sum / count : 0; 
        double max = array[indexMax];

        for (int i = 0; i < array.Length; i++)
        {
            if (array[i] == max)
            {
                array[i] = average;
            }
        }
    }

    public void Task_2_3(ref int[,] B, ref int[,] C)
    {
        // code here

        //  create and use method FindDiagonalMaxIndex(matrix);

        // end
    }

    public void Task_2_4(int[,] A, int[,] B)
    {
        // code here
        //  create and use method FindDiagonalMaxIndex(matrix); like in Task_2_3
        int maxIndexA = FindDiagonalMaxIndex(A);
        int maxIndexB = FindDiagonalMaxIndex(B);
        int size = A.GetLength(0);

        for (int i = 0; i < size; i++)
        {
            (A[maxIndexA, i], B[i, maxIndexB]) = (B[i, maxIndexB], A[maxIndexA, i]);
        }
        // end
    }

    public int FindDiagonalMaxIndex(int[,] matrix)
    {
        int maxIndex = 0;
        int size = matrix.GetLength(0);

        for (int i = 1; i < size; i++)
        {
            if (matrix[maxIndex, maxIndex] < matrix[i, i])
            {
                maxIndex = i;
            }
        }

        return maxIndex;
    }

    public void Task_2_5(int[,] A, int[,] B)
    {
        // code here

        // create and use FindMaxInColumn(matrix, columnIndex, out rowIndex);

        // end
    }

    public void Task_2_6(ref int[] A, int[] B)
    {
        // code here
        // create and use FindMax(matrix, out row, out column); like in Task_2_1
        // create and use DeleteElement(array, index);
        int maxIndexA = FindMax(A);
        int maxIndexB = FindMax(B);

        DeleteElement(ref A, maxIndexA);
        DeleteElement(ref B, maxIndexB);

        int[] C = new int[A.Length + B.Length];

        for (int i = 0; i < A.Length; i++)
        {
            C[i] = A[i];
        }

        for (int i = 0; i < B.Length; i++)
        {
            C[A.Length + i] = B[i];
        }

        A = C;
        // end
    }

    public int[] DeleteElement(ref int[] arr, int index)
    {
        int[] result = new int[arr.Length - 1];

        for (int i = 0, j = 0; i < arr.Length; i++)
        {
            if (i != index)
            {
                result[j] = arr[i];
                j++;
            }
        }

        arr = result;
        return arr;
    }

    public int FindMax(int[] arr)
    {
        int maxIndex = 0;

        for (int i = 1; i < arr.Length; i++)
        {
            if (arr[maxIndex] < arr[i])
            {
                maxIndex = i;
            }
        }

        return maxIndex;
    }

    public void Task_2_7(ref int[,] B, int[,] C)
    {
        // code here

        // create and use CountRowPositive(matrix, rowIndex);
        // create and use CountColumnPositive(matrix, colIndex);

        // end
    }

    public void Task_2_8(int[] A, int[] B)
    {
        // code here
        // create and use SortArrayPart(array, startIndex);
        int maxIndexA = FindMax(A);
        int maxIndexB = FindMax(B);

        SortArrayPart(ref A, maxIndexA);
        SortArrayPart(ref B, maxIndexB);
        // end
    }

    public int[] SortArrayPart(ref int[] arr, int startIndex)
    {
        for (int i = startIndex + 1; i < arr.Length;)
        {
            if (i == startIndex + 1 || arr[i] >= arr[i - 1])
            {
                i++;
            }
            else
            {
                int temp = arr[i];
                arr[i] = arr[i - 1];
                arr[i - 1] = temp;
                i--;
            }
        }

        return arr;
    }

    public int[] Task_2_9(int[,] A, int[,] C)
    {
        int[] answer = default(int[]);

        // code here

        // create and use SumPositiveElementsInColumns(matrix);

        // end

        return answer;
    }

    public void Task_2_10(ref int[,] matrix)
    {
        // code here
        // create and use RemoveColumn(matrix, columnIndex);
        int n = matrix.GetLength(0);
        int max = int.MinValue, min = int.MaxValue;
        int maxIndex = -1, minIndex = -1;

        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j <= i; j++)
            {
                if (matrix[i, j] > max)
                {
                    max = matrix[i, j];
                    maxIndex = j;
                }
            }
        }

        for (int i = 0; i < n; i++)
        {
            for (int j = i + 1; j < n; j++)
            {
                if (matrix[i, j] < min)
                {
                    min = matrix[i, j];
                    minIndex = j;
                }
            }
        }

        if (maxIndex != minIndex)
        {
            int firstIndexToRemove = Math.Max(maxIndex, minIndex);
            int secondIndexToRemove = Math.Min(maxIndex, minIndex);

            RemoveColumn(ref matrix, firstIndexToRemove);
            RemoveColumn(ref matrix, secondIndexToRemove);
        }
        else
        {
            RemoveColumn(ref matrix, minIndex);
        }
        // end
    }

    public int[,] RemoveColumn(ref int[,] matrix, int columnIndex)
    {
        int row = matrix.GetLength(0);
        int col = matrix.GetLength(1);
        int[,] result = new int[row, col - 1];

        for (int j = 0; j < col; j++)
        {
            for (int i = 0; i < row; i++)
            {
                if (j < columnIndex)
                {
                    result[i, j] = matrix[i, j];
                }
                else if (j > columnIndex)
                {
                    result[i, j - 1] = matrix[i, j];
                }
            }
        }

        matrix = result;
        return matrix;
    }
    public void Task_2_11(int[,] A, int[,] B)
    {
        // code here

        // use FindMaxIndex(matrix, out row, out column); from Task_2_1

        // end
    }
    public void Task_2_12(int[,] A, int[,] B)
    {
        // code here
        // create and use FindMaxColumnIndex(matrix);
        int maxIndexA = FindMaxColumnIndex(A);
        int maxIndexB = FindMaxColumnIndex(B);

        for (int i = 0; i < A.GetLength(0); i++)
        {
            int temp = A[i, maxIndexA];
            A[i, maxIndexA] = B[i, maxIndexB];
            B[i, maxIndexB] = temp;
        }
        // end
    }

    public int FindMaxColumnIndex(int[,] matrix)
    {
        int maxIndex = 0;
        int max = matrix[0, 0];

        for (int j = 0; j < matrix.GetLength(1); j++)
        {
            for (int i = 0; i < matrix.GetLength(0); i++)
            {
                if (matrix[i, j] > max)
                {
                    max = matrix[i, j];
                    maxIndex = j;
                }
            }
        }
        return maxIndex;
    }

    public void Task_2_13(ref int[,] matrix)
    {
        // code here

        // create and use RemoveRow(matrix, rowIndex);

        // end
    }

    public void Task_2_14(int[,] matrix)
    {
        // code here
        // create and use SortRow(matrix, rowIndex);
        for (int i = 0; i < matrix.GetLength(0); i++)
        {
            SortRow(ref matrix, i);
        }
        // end
    }

    public int[,] SortRow(ref int[,] matrix, int rowIndex)
    {
        int n = matrix.GetLength(1);
        int i = 0;

        while (i < n)
        {
            if (i == 0 || matrix[rowIndex, i] >= matrix[rowIndex, i - 1])
            {
                i++;
            }
            else
            {
                int temp = matrix[rowIndex, i];
                matrix[rowIndex, i] = matrix[rowIndex, i - 1];
                matrix[rowIndex, i - 1] = temp;
                i--;
            }
        }
        return matrix;
    }

    public int Task_2_15(int[,] A, int[,] B, int[,] C)
    {
        int answer = 0;

        // code here

        // create and use GetAverageWithoutMinMax(matrix);

        // end

        // 1 - increasing   0 - no sequence   -1 - decreasing
        return answer;
    }

    public void Task_2_16(int[] A, int[] B)
    {
        // code here
        // create and use SortNegative(array);
        SortNegative(ref A);
        SortNegative(ref B);
        // end
    }

    public int[] SortNegative(ref int[] array)
    {
        int negCount = 0;

        for (int i = 0; i < array.Length; i++)
        {
            if (array[i] < 0)
                negCount++;
        }

        if (negCount == 0) return array;

        int[] negative = new int[negCount];

        for (int i = 0, j = 0; i < array.Length; i++)
        {
            if (array[i] < 0)
            {
                negative[j] = array[i];
                j++;
            }
        }

        for (int i = 1; i < negative.Length; i++)
        {
            int key = negative[i];
            int j = i - 1;

            while (j >= 0 && negative[j] < key)
            {
                negative[j + 1] = negative[j];
                j--;
            }
            negative[j + 1] = key;
        }

        int negIndex = 0;
        for (int i = 0; i < array.Length; i++)
        {
            if (array[i] < 0)
            {
                array[i] = negative[negIndex];
                negIndex++;
            }
        }

        return array;
    }

    public void Task_2_17(int[,] A, int[,] B)
    {
        // code here

        // create and use SortRowsByMaxElement(matrix);

        // end
    }

    public void Task_2_18(int[,] A, int[,] B)
    {
        // code here
        // create and use SortDiagonal(matrix);
        SortDiagonal(ref A);
        SortDiagonal(ref B);
        // end
    }

    public int[,] SortDiagonal(ref int[,] matrix)
    {
        int n = matrix.GetLength(0);
        int[] diagonal = new int[n];

        for (int i = 0; i < n; i++)
        {
            diagonal[i] = matrix[i, i];
        }

        for (int i = 1; i < n; i++)
        {
            int key = diagonal[i];
            int j = i - 1;

            while (j >= 0 && diagonal[j] > key)
            {
                diagonal[j + 1] = diagonal[j];
                j--;
            }
            diagonal[j + 1] = key;
        }

        for (int i = 0; i < n; i++)
        {
            matrix[i, i] = diagonal[i];
        }

        return matrix;
    }

    public void Task_2_19(ref int[,] matrix)
    {
        // code here

        // use RemoveRow(matrix, rowIndex); from 2_13

        // end
    }
    public void Task_2_20(ref int[,] A, ref int[,] B)
    {
        // code here
        // use RemoveColumn(matrix, columnIndex); from 2_10
        RemoveAllColumns(ref A);
        RemoveAllColumns(ref B);
        // end
    }

    public int[,] RemoveAllColumns(ref int[,] matrix)
    {
        int col = matrix.GetLength(1);

        for (int j = col - 1; j >= 0; j--)
        {
            bool hasZero = false;

            for (int i = 0; i < matrix.GetLength(0); i++)
            {
                if (matrix[i, j] == 0)
                {
                    hasZero = true;
                    break;
                }
            }

            if (hasZero == false)
            {
                RemoveColumn(ref matrix, j);
            }
        }

        return matrix;
    }

    public void Task_2_21(int[,] A, int[,] B, out int[] answerA, out int[] answerB)
    {
        answerA = null;
        answerB = null;

        // code here

        // create and use CreateArrayFromMins(matrix);

        // end
    }

    public void Task_2_22(int[,] matrix, out int[] rows, out int[] cols)
    {
        rows = null;
        cols = null;

        // code here
        // create and use CountNegativeInRow(matrix, rowIndex);
        // create and use FindMaxNegativePerColumn(matrix);
        rows = CountNegativeInRows(matrix);
        cols = FindMaxNegativePerColumn(matrix);
        // end
    }

    public int[] CountNegativeInRows(int[,] matrix)
    {
        int rowCount = matrix.GetLength(0);
        int[] negativeCounts = new int[rowCount];

        for (int i = 0; i < rowCount; i++)
        {
            int count = 0;
            for (int j = 0; j < matrix.GetLength(1); j++)
            {
                if (matrix[i, j] < 0)
                {
                    count++;
                }
            }
            negativeCounts[i] = count;
        }
        return negativeCounts;
    }

    public int[] FindMaxNegativePerColumn(int[,] matrix)
    {
        int colCount = matrix.GetLength(1);
        int[] maxNegatives = new int[colCount];

        for (int j = 0; j < colCount; j++)
        {
            int maxNegative = int.MinValue;
            bool hasNegative = false;

            for (int i = 0; i < matrix.GetLength(0); i++)
            {
                if (matrix[i, j] < 0)
                {
                    hasNegative = true;

                    if (matrix[i, j] > maxNegative)
                    {
                        maxNegative = matrix[i, j];
                    }
                }
            }

            maxNegatives[j] = hasNegative ? maxNegative : 0;
        }

        return maxNegatives;
    }

    public void Task_2_23(double[,] A, double[,] B)
    {
        // code here

        // create and use MatrixValuesChange(matrix);

        // end
    }

    public void Task_2_24(int[,] A, int[,] B)
    {
        // code here
        // use FindMaxIndex(matrix, out row, out column); like in 2_1
        // create and use SwapColumnDiagonal(matrix, columnIndex);
        int rowA, columnA, rowB, columnB;

        FindMaxIndex(A, out rowA, out columnA);
        FindMaxIndex(B, out rowB, out columnB);

        SwapColumnDiagonal(ref A, columnA);
        SwapColumnDiagonal(ref B, columnB);
        // end
    }

    public void FindMaxIndex(int[,] matrix, out int row, out int column)
    {
        row = 0;
        column = 0;
        int maxValue = matrix[0, 0];

        for (int i = 0; i < matrix.GetLength(0); i++)
        {
            for (int j = 0; j < matrix.GetLength(1); j++)
            {
                if (matrix[i, j] > maxValue)
                {
                    maxValue = matrix[i, j];
                    row = i;
                    column = j;
                }
            }
        }
    }

    public int[,] SwapColumnDiagonal(ref int[,] matrix, int columnIndex)
    {
        int size = matrix.GetLength(0);

        for (int i = 0; i < size; i++)
        {
            int temp = matrix[i, i];
            matrix[i, i] = matrix[i, columnIndex];
            matrix[i, columnIndex] = temp;
        }

        return matrix;
    }

    public void Task_2_25(int[,] A, int[,] B, out int maxA, out int maxB)
    {
        maxA = 0;
        maxB = 0;

        // code here

        // create and use FindRowWithMaxNegativeCount(matrix);
        // in FindRowWithMaxNegativeCount create and use CountNegativeInRow(matrix, rowIndex); like in 2_22

        // end
    }

    public void Task_2_26(int[,] A, int[,] B)
    {
        // code here
        // create and use FindRowWithMaxNegativeCount(matrix); like in 2_25
        // in FindRowWithMaxNegativeCount use CountNegativeInRow(matrix, rowIndex); from 2_22
        int[] negCountA = CountNegativeInRows(A);
        int[] negCountB = CountNegativeInRows(B);

        int indexA = FindRowWithMaxNegativeCount(negCountA);
        int indexB = FindRowWithMaxNegativeCount(negCountB);

        for (int j = 0; j < A.GetLength(1); j++)
        {
            int temp = A[indexA, j];
            A[indexA, j] = B[indexB, j];
            B[indexB, j] = temp;
        }
        // end
    }
    public int FindRowWithMaxNegativeCount(int[] arr)
    {
        int maxIndex = 0;

        for (int i = 1; i < arr.Length; i++)
        {
            if (arr[i] > arr[maxIndex])
                maxIndex = i;
        }

        return maxIndex;
    }

    public void Task_2_27(int[,] A, int[,] B)
    {
        // code here

        // create and use FindRowMaxIndex(matrix, rowIndex, out columnIndex);
        // create and use ReplaceMaxElementOdd(matrix, row, column);
        // create and use ReplaceMaxElementEven(matrix, row, column);

        // end
    }

    public void Task_2_28a(int[] first, int[] second, ref int answerFirst, ref int answerSecond)
    {
        // code here
        // create and use FindSequence(array, A, B); // 1 - increasing, 0 - no sequence,  -1 - decreasing
        // A and B - start and end indexes of elements from array for search
        answerFirst = FindSequence(first, 0, first.Length - 1);
        answerSecond = FindSequence(second, 0, second.Length - 1);
        // end
    }

    public int FindSequence(int[] arr, int A, int B)
    {
        if (A >= B) return 0;

        bool isIncreasing = arr[A] < arr[A + 1];

        for (int i = A; i < B; i++)
        {
            if (isIncreasing && arr[i] > arr[i + 1]) return 0;
            if (!isIncreasing && arr[i] < arr[i + 1]) return 0;
        }

        return isIncreasing ? 1 : -1;
    }

    public void Task_2_28b(int[] first, int[] second, ref int[,] answerFirst, ref int[,] answerSecond)
    {
        // code here
        // use FindSequence(array, A, B); from Task_2_28a or entirely Task_2_28a
        // A and B - start and end indexes of elements from array for search
        answerFirst = FindSequences(first, out int nFirst);
        answerSecond = FindSequences(second, out int nSecond);
        // end
    }

    public int[,] FindSequences(int[] array, out int count)
    {
        int maxIntervals = array.Length * (array.Length - 1) / 2;
        int[,] sequences = new int[maxIntervals, 2];
        count = 0;

        for (int i = 0; i < array.Length; i++)
        {
            for (int j = i + 1; j < array.Length; j++)
            {
                if (FindSequence(array, i, j) != 0)
                {
                    sequences[count, 0] = i;
                    sequences[count, 1] = j;
                    count++;
                }
            }
        }

        return ResizeArray(sequences, count);
    }

    private int[,] ResizeArray(int[,] original, int newSize)
    {
        int[,] resized = new int[newSize, 2];

        for (int i = 0; i < newSize; i++)
        {
            resized[i, 0] = original[i, 0];
            resized[i, 1] = original[i, 1];
        }
        return resized;
    }

    public void Task_2_28c(int[] first, int[] second, ref int[] answerFirst, ref int[] answerSecond)
    {
        // code here
        // use FindSequence(array, A, B); from Task_2_28a or entirely Task_2_28a or Task_2_28b
        // A and B - start and end indexes of elements from array for search
        answerFirst = FindLongestSequence(first);
        answerSecond = FindLongestSequence(second);
        // end
    }
    public int[] FindLongestSequence(int[] array)
    {
        int maxLength = 0;
        int start = 0, end = 0;

        for (int i = 0; i < array.Length; i++)
        {
            for (int j = i + 1; j < array.Length; j++)
            {
                if (FindSequence(array, i, j) != 0)
                {
                    int length = j - i + 1;

                    if (length > maxLength)
                    {
                        maxLength = length;
                        start = i;
                        end = j;
                    }
                }
            }
        }

        return new int[] { start, end };
    }
    #endregion

    #region Level 3
    public void Task_3_1(ref double[,] firstSumAndY, ref double[,] secondSumAndY)
    {
        // code here

        // create and use public delegate SumFunction(x) and public delegate YFunction(x)
        // create and use method GetSumAndY(sFunction, yFunction, a, b, h);
        // create and use 2 methods for both functions calculating at specific x

        // end
    }

    public void Task_3_2(int[,] matrix)
    {
        // code here
        // create and use public delegate SortRowStyle(matrix, rowIndex);
        // create and use methods SortAscending(matrix, rowIndex) and SortDescending(matrix, rowIndex)
        // change method in variable sortStyle in the loop here and use it for row sorting
        SortRowStyle sortStyle = default(SortRowStyle);

        for (int i = 0; i < matrix.GetLength(0); i++)
        {
            sortStyle = (i % 2 == 0) ? SortAscending : SortDescending;
            sortStyle(ref matrix, i);
        }
        // end
    }

    public delegate void SortRowStyle(ref int[,] matrix, int rowIndex);

    public void SortAscending(ref int[,] matrix, int rowIndex)
    {
        for (int j = 1; j < matrix.GetLength(1); j++)
        {
            int key = matrix[rowIndex, j];
            int jj = j - 1;

            while (jj >= 0 && matrix[rowIndex, jj] > key)
            {
                matrix[rowIndex, jj + 1] = matrix[rowIndex, jj];
                jj--;
            }

            matrix[rowIndex, jj + 1] = key;
        }
    }

    public void SortDescending(ref int[,] matrix, int rowIndex)
    {
        for (int j = 1; j < matrix.GetLength(1); j++)
        {
            int key = matrix[rowIndex, j];
            int jj = j - 1;

            while (jj >= 0 && matrix[rowIndex, jj] < key)
            {
                matrix[rowIndex, jj + 1] = matrix[rowIndex, jj];
                jj--;
            }

            matrix[rowIndex, jj + 1] = key;
        }
    }

    public double Task_3_3(double[] array)
    {
        double answer = 0;
        // SwapDirection swapper = default(SwapDirection); - uncomment me

        // code here

        // create and use public delegate SwapDirection(array);
        // create and use methods SwapRight(array) and SwapLeft(array)
        // create and use method GetSum(array, start, h) that sum elements with a negative indexes
        // change method in variable swapper in the if/else and than use swapper(matrix)

        // end

        return answer;
    }

    public int Task_3_4(int[,] matrix, bool isUpperTriangle)
    {
        int answer = 0;

        // code here
        // create and use public delegate GetTriangle(matrix);
        // create and use methods GetUpperTriange(array) and GetLowerTriange(array)
        // create and use GetSum(GetTriangle, matrix)

        // end
        GetTriangle getTriangle = default(GetTriangle);

        getTriangle = (isUpperTriangle) ? GetUpperTriangle : GetLowerTriangle;
        answer = GetSum(getTriangle, matrix);

        return answer;
    }

    public delegate int[] GetTriangle(int[,] matrix);

    public int[] GetUpperTriangle(int[,] arr)
    {
        int size = arr.GetLength(0);
        int[] upperTriangle = new int[size * (size + 1) / 2];
        int index = 0;

        for (int i = 0; i < size; i++)
        {
            for (int j = i; j < size; j++)
            {
                upperTriangle[index++] = arr[i, j];
            }
        }

        return upperTriangle;
    }
    
    public int[] GetLowerTriangle(int[,] arr)
    {
        int size = arr.GetLength(0);
        int[] lowerTriangle = new int[size * (size + 1) / 2];
        int index = 0;

        for (int i = 0; i < size; i++)
        {
            for (int j = 0; j <= i; j++)
            {
                lowerTriangle[index++] = arr[i, j];
            }
        }

        return lowerTriangle;
    }

   
    public int GetSum(GetTriangle getTriangle, int[,] arr)
    {
        int sum = 0;
        int [] array = getTriangle(arr);

        foreach (int value in array)
        {
            sum += value * value; 
        }

        return sum;
    }

    public void Task_3_5(out int func1, out int func2)
    {
        func1 = 0;
        func2 = 0;

        // code here

        // use public delegate YFunction(x, a, b, h) from Task_3_1
        // create and use method CountSignFlips(YFunction, a, b, h);
        // create and use 2 methods for both functions

        // end
    }

    public void Task_3_6(int[,] matrix)
    {
        // code here
        // create and use public delegate FindElementDelegate(matrix);
        // use method FindDiagonalMaxIndex(matrix) like in Task_2_3;
        // create and use method FindFirstRowMaxIndex(matrix);
        // create and use method SwapColumns(matrix, FindDiagonalMaxIndex, FindFirstRowMaxIndex);
        SwapColumns(ref matrix, FindDiagonalMaxIndex, FindFirstRowMaxIndex);
        // end
    }

    public delegate int FindElementDelegate(int[,] matrix);

    public int FindFirstRowMaxIndex(int[,] matrix)
    {
        int maxIndex = 0;
        int maxValue = matrix[0, 0];

        for (int j = 1; j < matrix.GetLength(1); j++)
        {
            if (matrix[0, j] > maxValue)
            {
                maxValue = matrix[0, j];
                maxIndex = j;
            }
        }

        return maxIndex;
    }

    public void SwapColumns(ref int[,] matrix, FindElementDelegate findDiagonalMaxIndex, FindElementDelegate findFirstRowMaxIndex)
    {
        int[] C = new int[matrix.GetLength(0)]; 
        int diagonalIndex = findDiagonalMaxIndex(matrix); 
        int firstrowIndex = findFirstRowMaxIndex(matrix); 

        if (diagonalIndex == firstrowIndex) return; 

        for (int i = 0; i < matrix.GetLength(0); i++)
        {
            C[i] = matrix[i, diagonalIndex];
            matrix[i, diagonalIndex] = matrix[i, firstrowIndex];
            matrix[i, firstrowIndex] = C[i];
        }
    }

    public void Task_3_7(ref int[,] B, int[,] C)
    {
        // code here

        // create and use public delegate CountPositive(matrix, index);
        // use CountRowPositive(matrix, rowIndex) from Task_2_7
        // use CountColumnPositive(matrix, colIndex) from Task_2_7
        // create and use method InsertColumn(matrixB, CountRow, matrixC, CountColumn);

        // end
    }

    public void Task_3_10(ref int[,] matrix)
    {
        // code here
        // create and use public delegate FindIndex(matrix);
        // create and use method FindMaxBelowDiagonalIndex(matrix);
        // create and use method FindMinAboveDiagonalIndex(matrix);
        // use RemoveColumn(matrix, columnIndex) from Task_2_10
        // create and use method RemoveColumns(matrix, findMaxBelowDiagonalIndex, findMinAboveDiagonalIndex)
        FindIndex searchArea = default(FindIndex);
    
        RemoveColumns(ref matrix, FindMaxBelowDiagonalIndex, FindMinAboveDiagonalIndex);
        // end
    }

    public delegate int FindIndex(int[,] matrix);

    public int FindMaxBelowDiagonalIndex(int[,] matrix)
    {
        int maxIndex = 0;
        int maxValue = matrix[0, 0];

        for (int i = 0; i < matrix.GetLength(0); i++)
        {
            for (int j = 0; j <= i; j++)
            {
                if (matrix[i, j] > maxValue)
                {
                    maxValue = matrix[i, j];
                    maxIndex = j; 
                }
            }
        }

        return maxIndex;
    }

    public int FindMinAboveDiagonalIndex(int[,] matrix)
    {
        int minIndex = 0;
        int minValue = matrix[0, 1];
        int size = matrix.GetLength(0);

        for (int i = 0; i < size; i++)
        {
            for (int j = i + 1; j < size; j++) 
            {
                if (matrix[i, j] < minValue)
                {
                    minValue = matrix[i, j];
                    minIndex = j;
                }
            }
        }

        return minIndex;
    }

    public void RemoveColumns(ref int[,] matrix, FindIndex findMaxBelowDiagonalIndex, FindIndex findMinAboveDiagonalIndex)
    {
        int maxIndex = findMaxBelowDiagonalIndex(matrix);
        int minIndex = findMinAboveDiagonalIndex(matrix);
        
        if (maxIndex == minIndex)
        {
            RemoveColumn(ref matrix, maxIndex);
        }
        else if (maxIndex > minIndex)
        {
            RemoveColumn(ref matrix, maxIndex);
            RemoveColumn(ref matrix, minIndex);
        }
        else
        {
            RemoveColumn(ref matrix, minIndex);
            RemoveColumn(ref matrix, maxIndex);
        }
    }

    public void Task_3_13(ref int[,] matrix)
    {
        // code here

        // use public delegate FindElementDelegate(matrix) from Task_3_6
        // create and use method RemoveRows(matrix, findMax, findMin)

        // end
    }

    public void Task_3_22(int[,] matrix, out int[] rows, out int[] cols)
    {
        // code here
        // create and use public delegate GetNegativeArray(matrix);
        // use GetNegativeCountPerRow(matrix) from Task_2_22
        // use GetMaxNegativePerColumn(matrix) from Task_2_22
        // create and use method FindNegatives(matrix, searcherRows, searcherCols, out rows, out cols);
        rows = null;
        cols = null;

        FindNegatives(matrix, CountNegativeInRows, FindMaxNegativePerColumn, out rows, out cols);
        // end
    }

    public delegate int[] GetNegativeArray(int[,] matrix);

    public void FindNegatives(int[,] matrix, GetNegativeArray searcherRows, GetNegativeArray searcherCols, out int[] rows, out int[] cols)
    {
        rows = searcherRows(matrix);
        cols = searcherCols(matrix);
    }

    public void Task_3_27(int[,] A, int[,] B)
    {
        // code here

        // create and use public delegate ReplaceMaxElement(matrix, rowIndex, max);
        // use ReplaceMaxElementOdd(matrix) from Task_2_27
        // use ReplaceMaxElementEven(matrix) from Task_2_27
        // create and use method EvenOddRowsTransform(matrix, replaceMaxElementOdd, replaceMaxElementEven);

        // end
    }

    public void Task_3_28a(int[] first, int[] second, ref int answerFirst, ref int answerSecond)
    {
        // code here
        // create public delegate IsSequence(array, left, right);
        // create and use method FindIncreasingSequence(array, A, B); similar to FindSequence(array, A, B) in Task_2_28a
        // create and use method FindDecreasingSequence(array, A, B); similar to FindSequence(array, A, B) in Task_2_28a
        // create and use method DefineSequence(array, findIncreasingSequence, findDecreasingSequence);
        answerFirst = DefineSequence(first, findIncreasingSequence, findDecreasingSequence);
        answerSecond = DefineSequence(second, findIncreasingSequence, findDecreasingSequence);
        // end
    }

    public delegate bool IsSequence(int[] arr, int left, int right);

    public bool findIncreasingSequence(int[] arr, int a, int b)
    {
        int c = FindSequence(arr, a, b);

        return (c == 1);
    }

    public bool findDecreasingSequence(int[] arr, int a, int b)
    {
        int c = FindSequence(arr, a, b);

        return (c == -1);
    }

    public int DefineSequence(int[] arr, IsSequence findIncreasingSequence, IsSequence findDecreasingSequence)
    {
        if (findIncreasingSequence(arr, 0, arr.Length - 1))
            return 1;
        else if (findDecreasingSequence(arr, 0, arr.Length - 1))
            return -1;
        else
            return 0;
    }

    public void Task_3_28c(int[] first, int[] second, ref int[] answerFirstIncrease, ref int[] answerFirstDecrease, ref int[] answerSecondIncrease, ref int[] answerSecondDecrease)
    {
        // code here
        // create public delegate IsSequence(array, left, right);
        // use method FindIncreasingSequence(array, A, B); from Task_3_28a
        // use method FindDecreasingSequence(array, A, B); from Task_3_28a
        // create and use method FindLongestSequence(array, sequence);
        answerFirstIncrease = FindLongestSequence(first, findIncreasingSequence);
        answerFirstDecrease = FindLongestSequence(first, findDecreasingSequence);
       
        answerSecondIncrease = FindLongestSequence(second, findIncreasingSequence);
        answerSecondDecrease = FindLongestSequence(second, findDecreasingSequence);
        // end
    }

    public int[] FindLongestSequence(int[] arr, IsSequence sequence)
    {
        int maxLength = 0;
        int start = 0, end = 0;

        for (int i = 0; i < arr.Length; i++)
        {
            for (int j = i + 1; j < arr.Length; j++)
            {
                if (sequence(arr, i, j))
                {
                    int currentLength = j - i + 1;

                    if (currentLength > maxLength)
                    {
                        maxLength = currentLength;
                        start = i;
                        end = j;
                    }
                }
            }
        }

        return new int[] {start, end};
    }
    #endregion
    #region bonus part
    public double[,] Task_4(double[,] matrix, int index)
    {
        // MatrixConverter[] mc = new MatrixConverter[4]; - uncomment me

        // code here

        // create public delegate MatrixConverter(matrix);
        // create and use method ToUpperTriangular(matrix);
        // create and use method ToLowerTriangular(matrix);
        // create and use method ToLeftDiagonal(matrix); - start from the left top angle
        // create and use method ToRightDiagonal(matrix); - start from the right bottom angle

        // end

        return matrix;
    }
    #endregion
}
