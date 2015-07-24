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

package alg;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.zip.GZIPInputStream;

/**
 * A class that contains a String representation of the
 * edges and vertices on a path and the path weight.  Used
 * to compare paths stored to a file before creating Path objects.
 *
 */
public class StringPath implements Comparable<StringPath>
{
	/** A String representation of a path */
	private String path;
	private double weight;
	
	/**
	 * Construct the StringPath from the String written by a StringPath's toString
	 * @param pathAndWeight
	 */
	public StringPath(String pathAndWeight)
	{
		String[] parts = pathAndWeight.trim().split("\t");
		path = parts[0];
		weight = Double.parseDouble(parts[1]);
	}
	
	/**
	 * Construct the StringPath from the list of edges and vertices and the path weight.
	 * The path weight is generally given without including the target weight if
	 * the StringPath is to be written to a file when precomputing paths.
	 * @param edges
	 * @param vertices
	 * @param weight
	 */
	public StringPath(List<Edge> edges, List<Vertex> vertices, double weight)
	{
		this.weight = weight;
		
		// Build the String representation of the path
		StringBuffer buf = new StringBuffer();
		// Traverse all edges and vertices
		for(int e = 0; e < edges.size(); e++)
		{
			// The edges aren't oriented so we need to use the
			// vertex list to print the edge endpoints in the
			// correct order
			Edge pathEdge = edges.get(e);
			if(!(pathEdge.containsVertex(vertices.get(e)) &&
					pathEdge.containsVertex(vertices.get(e+1))))
			{
				throw new IllegalStateException("Edge " + pathEdge + " was expected to " + 
						"contain " + vertices.get(e) + " and " + vertices.get(e+1));
			}

			buf.append(vertices.get(e)).append(":");
			buf.append(pathEdge.getType()).append(":");
			buf.append(vertices.get(e+1)).append("|");
		}

		// Remove the extra "|"
		buf.deleteCharAt(buf.length()-1);
		path = buf.toString();
	}
	
	/**
	 * Construct a path using the edges, vertices, and path weight from a Path object
	 * @param path
	 */
	public StringPath(Path path)
	{
		this(Arrays.asList(path.getEdges()),
				Arrays.asList(path.getVertices()),
				path.maxWeight());
	}
	
	
	public String getPath()
	{
		return path;
	}
	
	public double getWeight()
	{
		return weight;
	}
	
	/**
	 * Parse the String representation of the path to obtain
	 * the name of the last vertex
	 * @return
	 */
	public String getTarget()
	{
		int ind = path.lastIndexOf(':');
		
		if(ind < 0)
		{
			throw new IllegalStateException("No target found: " + this);
		}
		
		return path.substring(ind + 1);
	}
	
	/**
	 * StringPaths are typically created before the target weight is known
	 * so update the weight by multiplying it by the target weight
	 * @param targetWeight
	 */
	public void addTargetWeight(double targetWeight)
	{
		weight *= targetWeight;
	}

	public String toString()
	{
		return path + "\t" + weight;
	}
	
	/**
	 * Compare this StringPath with another StringPath using
	 * the weight to make the comparion.  If the weights are
	 * equal, use String representation of the paths to
	 * break ties.
	 */
	public int compareTo(StringPath other)
	{
		int weightComp = Double.compare(weight, other.weight);
		if(weightComp == 0)
		{
			return path.compareTo(other.path);
		}
		else
		{
			return weightComp;
		}
	}
	
	/**
	 * Create a list of StringPath objects from a file that contains the toString
	 * of StringPath objects.  Assumes no header and an uncompressed file.
	 * @param pathFile
	 * @return
	 * @throws IOException
	 */
	public static ArrayList<StringPath> loadFile(String pathFile)
		throws IOException
	{
		return loadFile(pathFile, false);
	}
	
	/**
	 * Create a list of StringPath objects from a file that contains the toString
	 * of StringPath objects.  Assumes no header and allows the user to specifiy if
	 * the file is compressed.
	 * @return
	 * @throws IOException
	 */
	public static ArrayList<StringPath> loadFile(String pathFile, boolean compressed)
		throws IOException
	{
		BufferedReader reader;
		if(compressed)
		{
			reader = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(pathFile))));
		}
		else
		{
			reader = new BufferedReader(new FileReader (pathFile));
		}
		
		ArrayList<StringPath> paths = new ArrayList<StringPath>();
		
		String line;
		while((line = reader.readLine()) != null)
		{
			paths.add(new StringPath(line));
		}
		reader.close();
				
		return paths;
	}
}