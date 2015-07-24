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
import java.util.Collection;
import java.util.HashMap;

/**
 * A vertex in a graph that may take a weight.  The vertex may
 * take an additional weight that is applied only when it
 * functions as a target.
 *
 */

public class Vertex {

	private static final double DEFAULT_WEIGHT = 1;
	
	private double weight, targWeight;
	/** If the vertex id is [name]_[info], the name is used.  If there is
	 * no info, the id and name are the same */
	private String name;
	/** An id is a gene name appended with _[info] where info is
	 * used to differentiate this vertex from another with the same name */
	private String id;
	/** A non-negative, numerical id used for this vertex when it appears
	 * in a lp_solve linear program.  Assigned automatically and cannot
	 * be changed. */
	private int lpId;
	/** A counter used to assign each vertex a unique LP id */
	private static int nextLpId = 0;
	
	/** A map to track which vertex names are already associated with
	 * existing vertices */
	private static HashMap<String, Vertex> vertices = new HashMap<String, Vertex>();
	
	private Vertex(double weight, double targWeight, String id)
	{
		this.weight = weight;
		this.targWeight = targWeight;
		this.id = id;
		lpId = nextLpId++;
		
		if(id.contains("_"))
		{
			int index = id.indexOf('_');
			name = id.substring(0, index);
		}
		else
		{
			name = id;
		}
	}

	/**
	 * Create a new Vertex or return the existing Vertex with this id.
	 * The Vertex's weights need to be set separately.
	 * @param id The id will be saved in upper case
	 * @return
	 */
	public static Vertex createVertex(String id)
	{
		id = id.toUpperCase().trim();
		
		// If this vertex exists return a reference to it
		if(vertices.containsKey(id))
		{
			return vertices.get(id);
		}
		// Otherwise create the new Vertex with default weights
		else
		{
			Vertex v = new Vertex(DEFAULT_WEIGHT, DEFAULT_WEIGHT, id);
			vertices.put(id, v);
			return v;
		}
	}
	
	
	/**
	 * 
	 * @return a reference to the set of all vertices that have been created
	 */
	public static Collection<Vertex> getVertices()
	{
		return vertices.values();
	}
	
	public void setWeight(double w)
	{
		weight = w;
	}
	
	public void setTargetWeight(double w)
	{
		targWeight = w;
	}
	
	public double getWeight()
	{
		return weight;
	}
	
	public double getTargetWeight()
	{
		return targWeight;
	}
	
	public String getName()
	{
		return name;
	}
	
	/**
	 * Ignores the id of the vertices, if one is set
	 * @param otherVertex
	 * @return true if this vertex and otherVertex have the same name
	 */
	public boolean sameName(Vertex otherVertex)
	{
		return name.equals(otherVertex.getName());
	}
	
	public String getId()
	{
		return id;
	}
	
	public int getLpId()
	{
		return lpId;
	}
	
	/**
	 * 
	 * @return the number of vertex names (keys) in the static vertex map
	 */
	public static int uniqueVertices()
	{
		return vertices.keySet().size();
	}
	
	public String detailedToString()
	{
		return id + " : " + weight + " : " + targWeight;
	}
	
	public String toString()
	{
		return id;
	}
}
