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

import util.MapUtil;

/**
 * A vertex in a graph that may take a weight.  The vertex may
 * take an additional weight that is applied only when it
 * functions as a target.
 *
 */

public class Vertex implements Comparable<Vertex>{

	private static double defaultWeight = 1;
	private static double defaultTargWeight = 1;
	
	/** Track the largest vertex weight seen so far among all vertices.
	 * Does not change even if the vertex with this maximum weight has its
	 * weight reduced later. */
	private static double maxWeight = Double.NEGATIVE_INFINITY;
	/** Track the largest target weight seen so far among all vertices.
	 * Does not change even if the vertex with this maximum weight has its
	 * target weight reduced later. */
	private static double maxTargWeight = Double.NEGATIVE_INFINITY;
	
	private double weight, targWeight;
	/** The vertex id and name are now identical */
	private String name;
	/** The vertex id and name are now identical */
	private String id;
	/** A non-negative, numerical id used for this vertex when it appears
	 * in a lp_solve linear program.  Assigned automatically and cannot
	 * be changed. */
	private int lpId;
	/** A counter used to assign each vertex a unique LP id */
	private static int nextLpId = 0;
	
	/** If setNodePriors has been called, then the priors provided there
	 * are stored in this map and will be used to set the weight of all
	 * new nodes. */
	private static HashMap<String, Double> priors = null;
	
	/** A map to track which vertex names are already associated with
	 * existing vertices */
	private static HashMap<String, Vertex> vertices = new HashMap<String, Vertex>();
	
	private Vertex(double weight, double targWeight, String id)
	{
		setWeight(weight);
		setTargetWeight(targWeight);
		this.id = id;
		lpId = nextLpId++;
		// No longer store separate name and id
		name = id;
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
			Vertex v;
			// Use the specific prior provided if there is one for this node
			if(priors != null && priors.containsKey(id))
			{
				v = new Vertex(priors.get(id), defaultTargWeight, id);
			}
			else
			{
				v = new Vertex(defaultWeight, defaultTargWeight, id);
			}
			vertices.put(id, v);
			return v;
		}
	}
	
	// TODO Need to have this make any paths in the graph stale
	/**
	 * Updates the weights of all existing vertices, but does not affect their
	 * target weights.  Nodes with a specific weight listed will have that weight set,
	 * and all other nodes will use the new default weight.  All future vertices
	 * that are created will have their weights set in the same manner.  Should not be
	 * called after paths have been found.  If the priorFile is null or the empty
	 * string, all nodes will be assigned the default prior.
	 * @param priorFile a list of nodes and their prior (aka weight)
	 * @param defaultPrior this value is used for all nodes that do not have a specific
	 * weight listed in the prior file
	 */
	public static void setNodePriors(String priorFile, double defaultPrior)
	{
		defaultWeight = defaultPrior;
		
		if(priorFile == null || priorFile.equals(""))
		{
			// Reset all node-specific priors and use the default prior for all nodes
			priors = new HashMap<String, Double>();
			System.out.println("Using default node prior " + defaultPrior + " for all nodes");
		}
		else
		{
			HashMap<String, String> strPriors = MapUtil.loadMap(priorFile, 0, 1, false);
			System.out.println("Using " + strPriors.size() + " node-specific " + 
					"priors and default prior " + defaultPrior + " for all other nodes");
			
			// Convert the priors to doubles
			priors = new HashMap<String, Double>(strPriors.size());
			for(String id : strPriors.keySet())
			{
				priors.put(id.toUpperCase().trim(), Double.valueOf(strPriors.get(id)));
			}
		}
		
		// Iterate through all existing vertices and update their weights/priors
		for(Vertex v : getVertices())
		{
			if(priors.containsKey(v.id))
			{
				v.setWeight(priors.get(v.id));
			}
			else
			{
				v.setWeight(defaultWeight);
			}
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
	
	/**
	 * Obtain a reference to the Vertex with the specified id or return
	 * null if no such Vertex exists
	 * @param id
	 * @return the Vertex, if found, or null
	 */
	public static Vertex findVertex(String id)
	{
		id = id.toUpperCase().trim();
		
		// If this Vertex exists return a reference to it
		if(vertices.containsKey(id))
		{
			return vertices.get(id);
		}
		else
		{
			return null;
		}
	}
	
	/**
	 * @return The largest weight that has been assigned to a vertex
	 */
	public static double getMaxWeight()
	{
		return maxWeight;
	}
	
	/**
	 * @return The largest target weight that has been assigned to a vertex
	 */
	public static double getMaxTargetWeight()
	{
		return maxTargWeight;
	}
	
	/**
	 * Also updates the maximum weight seen among all vertices if necessary
	 * @param w
	 */
	public void setWeight(double w)
	{
		weight = w;
		
		if(w > maxWeight)
		{
			maxWeight = w;
		}
	}
	
	/**
	 * Also updates the maximum target weight seen among all vertices if necessary
	 * @param w
	 */
	public void setTargetWeight(double w)
	{
		targWeight = w;
		
		if(w > maxTargWeight)
		{
			maxTargWeight = w;
		}
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
	 * @return true if this vertex and otherVertex have the same name (case-sensitive)
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

	/**
	 * Compare based on node id
	 */
	public int compareTo(Vertex other)
	{
		return this.getId().compareTo(other.getId());
	}
}
