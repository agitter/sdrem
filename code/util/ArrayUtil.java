/*
 * Copyright (c) 2009, Anthony Gitter
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without modification,
 * are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 * 3. Neither the name of Carnegie Mellon University nor the names of its contributors
 *    may be used to endorse or promote products derived from this software without
 *    specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS ``AS IS''
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 * DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
 * OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

package util;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;

public class ArrayUtil {
	
	/**
	 * Takes two matrices of the same dimension and adds every element in
	 * the second matrix to the first matrix, overwriting the original values.
	 * Returns a reference to the first
	 * matrix.  If one of the matrices is null, a reference to the other
	 * is returned.
	 * @param matrix1
	 * @param matrix2
	 * @return a reference to the first matrix (unless it is null)
	 */
	public static int[][] sumMatrices(int[][] matrix1, int[][] matrix2)
	{
		if(matrix1 == null && matrix2 == null)
		{
			throw new IllegalArgumentException("Both matrices are null");
		}
		
		if(matrix1 == null)
		{
			return matrix2;
		}
		
		if(matrix2 == null)
		{
			return matrix1;
		}
		
		if(matrix1.length != matrix2.length)
		{
			throw new IllegalArgumentException("Matrix dimensions are not the same");
		}
		
		for(int i = 0; i < matrix1.length; i++)
		{
			if(matrix1[i].length != matrix2[i].length)
			{
				throw new IllegalArgumentException("Matrix dimensions are not the same");
			}
		}
		
		// We have two equally sized valid matrices so sum the elements
		for(int i = 0; i < matrix1.length; i++)
		{
			for(int j = 0; j < matrix1[i].length; j++)
			{
				matrix1[i][j] += matrix2[i][j];
			}
		}
			
		return matrix1;
	}
	
	/**
	 * Takes two matrices of the same dimension and adds every element in
	 * the second matrix to the first matrix, overwriting the original values.
	 * Returns a reference to the first
	 * matrix.  If one of the matrices is null, a reference to the other
	 * is returned.
	 * @param matrix1
	 * @param matrix2
	 * @return a reference to the first matrix (unless it is null)
	 */
	public static Integer[][] sumMatrices(Integer[][] matrix1, Integer[][] matrix2)
	{
		if(matrix1 == null && matrix2 == null)
		{
			throw new IllegalArgumentException("Both matrices are null");
		}
		
		if(matrix1 == null)
		{
			return matrix2;
		}
		
		if(matrix2 == null)
		{
			return matrix1;
		}
		
		if(matrix1.length != matrix2.length)
		{
			throw new IllegalArgumentException("Matrix dimensions are not the same");
		}
		
		for(int i = 0; i < matrix1.length; i++)
		{
			if(matrix1[i].length != matrix2[i].length)
			{
				throw new IllegalArgumentException("Matrix dimensions are not the same");
			}
		}
		
		// We have two equally sized valid matrices so sum the elements
		for(int i = 0; i < matrix1.length; i++)
		{
			for(int j = 0; j < matrix1[i].length; j++)
			{
				matrix1[i][j] = new Integer(matrix1[i][j].intValue() + matrix2[i][j].intValue());
			}
		}
			
		return matrix1;
	}
	
	/**
	 * Takes two matrices of the same dimension and adds every element in
	 * the second matrix to the first matrix, overwriting the original values.
	 * Returns a reference to the first
	 * matrix.  If one of the matrices is null, a reference to the other
	 * is returned.
	 * @param matrix1
	 * @param matrix2
	 * @return a reference to the first matrix (unless it is null)
	 */
	public static double[][] sumMatrices(double[][] matrix1, double[][] matrix2)
	{
		if(matrix1 == null && matrix2 == null)
		{
			throw new IllegalArgumentException("Both matrices are null");
		}
		
		if(matrix1 == null)
		{
			return matrix2;
		}
		
		if(matrix2 == null)
		{
			return matrix1;
		}
		
		if(matrix1.length != matrix2.length)
		{
			throw new IllegalArgumentException("Matrix dimensions are not the same");
		}
		
		for(int i = 0; i < matrix1.length; i++)
		{
			if(matrix1[i].length != matrix2[i].length)
			{
				throw new IllegalArgumentException("Matrix dimensions are not the same");
			}
		}
		
		// We have two equally sized valid matrices so sum the elements
		for(int i = 0; i < matrix1.length; i++)
		{
			for(int j = 0; j < matrix1[i].length; j++)
			{
				matrix1[i][j] += matrix2[i][j];
			}
		}
			
		return matrix1;
	}
	
	/**
	 * Divide every element in the matrix by an integer.
	 * @param matrix
	 * @param constant
	 * @return A reference to the matrix
	 */
	public static int[][] divByConst(int[][] matrix, int constant)
	{
		for(int i = 0; i < matrix.length; i++)
		{
			for(int j = 0; j < matrix[i].length; j++)
			{
				matrix[i][j] /= constant;
			}
		}
		
		return matrix;
	}
	
	/**
	 * Divide every element in the matrix by a double.
	 * @param matrix
	 * @param constant
	 * @return A reference to the matrix
	 */
	public static double[][] divByConst(double[][] matrix, double constant)
	{
		for(int i = 0; i < matrix.length; i++)
		{
			for(int j = 0; j < matrix[i].length; j++)
			{
				matrix[i][j] /= constant;
			}
		}
		
		return matrix;
	}
	
	// TODO more error checking
	/**
	 * Assumes all rows are the same length and the matrix has at least
	 * one column.
	 * @param matrix
	 * @return a new double matrix with the same values as the Integer array
	 */
	public static double[][] integerToDouble(Integer[][] matrix)
	{
		double[][] dMatrix = new double[matrix.length][matrix[0].length];
		
		for(int i = 0; i < matrix.length; i++)
		{
			for(int j = 0; j < matrix[i].length; j++)
			{
				dMatrix[i][j] = matrix[i][j].doubleValue();
			}
		}
		
		return dMatrix;
	}
	
	/**
	 * Assumes all rows are the same length and the matrix has at least
	 * one column.
	 * @param matrix
	 * @return a new double matrix with the same values as the int array
	 */
	public static double[][] intToDouble(int[][] matrix)
	{
		double[][] dMatrix = new double[matrix.length][matrix[0].length];
		
		for(int i = 0; i < matrix.length; i++)
		{
			for(int j = 0; j < matrix[i].length; j++)
			{
				dMatrix[i][j] = matrix[i][j];
			}
		}
		
		return dMatrix;
	}
	
	/**
	 * Assumes all rows are the same length and the matrix has at least
	 * one column.
	 * @param matrix
	 * @return a new Integer matrix with the same values as the int array
	 */
	public static Integer[][] intToInteger(int[][] matrix)
	{
		Integer[][] integerMatrix = new Integer[matrix.length][matrix[0].length];
		
		for(int i = 0; i < matrix.length; i++)
		{
			for(int j = 0; j < matrix[i].length; j++)
			{
				integerMatrix[i][j] = new Integer(matrix[i][j]);
			}
		}
		
		return integerMatrix;
	}
	
	
	public static void printMatrix(int[][] matrix)
	{
		for(int i = 0; i < matrix.length; i++)
		{
			for(int j = 0; j < matrix[i].length; j++)
			{
				System.out.print(matrix[i][j] + "\t");
			}
			System.out.println();
		}
	}
	
	public static void printMatrix(double[][] matrix)
	{
		for(int i = 0; i < matrix.length; i++)
		{
			for(int j = 0; j < matrix[i].length; j++)
			{
				System.out.print(matrix[i][j] + "\t");
			}
			System.out.println();
		}
	}
	
	/**
	 * Return the median element
	 * @param values
	 * @return
	 */
	public static double median(Collection<Double> values)
	{
		if(values.size() == 0)
		{
			throw new IllegalArgumentException("Cannot take median of an empty collection");
		}
		
		Double[] sortedVals = new Double[0];
		sortedVals = values.toArray(sortedVals);
		Arrays.sort(sortedVals);
		
		if(sortedVals.length % 2 == 1) // odd
		{
			int middle = (int)Math.floor(sortedVals.length / 2.0);
			return sortedVals[middle];
		}
		else // even
		{
			int lower = (int)Math.floor((sortedVals.length-1) / 2.0);
			int upper = (int)Math.ceil((sortedVals.length-1) / 2.0);
			
			return (sortedVals[lower] + sortedVals[upper]) / 2;
		}
	}
	
	/**
	 * Return the maximum element
	 * @param values
	 * @return
	 */
	public static double max(Collection<Double> values)
	{
		if(values.size() == 0)
		{
			throw new IllegalArgumentException("Cannot take maximum of an empty collection");
		}
		
		// Default to the first value
		double max = values.iterator().next();
		
		for(double cur : values)
		{
			if(cur > max)
			{
				max = cur;
			}
		}
		
		return max;
	}
	
	/**
	 * Return the minimum element
	 * @param values
	 * @return
	 */
	public static double min(Collection<Double> values)
	{
		if(values.size() == 0)
		{
			throw new IllegalArgumentException("Cannot take minimum of an empty collection");
		}
		
		// Default to the first value
		double min = values.iterator().next();
		
		for(double cur : values)
		{
			if(cur < min)
			{
				min = cur;
			}
		}
		
		return min;
	}
	
	/**
	 * Return the average value
	 * @param values
	 * @return
	 */
	public static double average(Collection<Double> values)
	{
		if(values.size() == 0)
		{
			throw new IllegalArgumentException("Cannot take average of an empty collection");
		}
		
		double sum = 0;
		
		for(double cur : values)
		{
			sum += cur;
		}
		
		return sum / values.size();
	}
}
