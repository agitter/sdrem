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

import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;

/**
 * Stores a pathway from SGD.  The vertices are small molecules(?)
 * along the pathway that are operated on by enzymes (genes) on the directed
 * edges.
 *
 */
public class SGDPathway {
	
	private HashMap<String, HashSet<PathEdge>> edges;
	private String name, id;
	/** A unique id is assigned to each edge so that if the metabolic pathway
	 * is converted to a gene pathway different occurrences of a single
	 * gene along different metabolic pathway edges can be differentiated */
	private static int edgeId = 0;
	
	public SGDPathway()
	{
		edges = new HashMap<String, HashSet<PathEdge>>();
	}
	
	public void addEdge(String source, String target, HashSet<String> edgeGenes)
	{
		PathEdge pe = new PathEdge(edgeGenes, target, edgeId++);
		
		HashSet<PathEdge> edgeSet;
		if(edges.containsKey(source))
		{
			edgeSet = edges.get(source);
		}
		else
		{
			edgeSet = new HashSet<PathEdge>();
		}
		
		edgeSet.add(pe);
		edges.put(source, edgeSet);
	}
	
	/**
	 * Convert this metabolic pathway to a gene pathway.
	 * @param pathCounter the id of this pathway, which ensures
	 * that if this pathway is merged its genes do not become
	 * associated with genes/edges in the other pathway
	 * @return
	 */
	public Graph toGenePathway(int pathCounter)
	{
		Graph genePath = new Graph();
		HashSet<Vertex> allVertices = new HashSet<Vertex>();
		
		// Iterate through all metabolic pathway edges, then find the edges
		// that come after each edge
		for(Map.Entry<String, HashSet<PathEdge>> entry : edges.entrySet())
		{
			for(PathEdge edge : entry.getValue())
			{
				if(edges.containsKey(edge.target))
				{
					// Find all the edges leaving this metabolite
					for(PathEdge nextEdge : edges.get(edge.target))
					{
						// Create a directed edge from all source genes
						// (genes on the first metabolic pathway edge)
						// to all target genes (unless the source and
						// target have the same name)
						for(String src : edge.genes)
						{
							for(String targ : nextEdge.genes)
							{
								if(!src.equalsIgnoreCase(targ))
								{
									Vertex srcVert = Vertex.createVertex(src + 
											"_id" + edge.id + "pa" + pathCounter);
									Vertex targVert = Vertex.createVertex(targ + 
											"_id" + nextEdge.id + "pa" + pathCounter);
									
									genePath.addEdge(new DirEdge(srcVert, targVert));
									allVertices.add(srcVert);
									allVertices.add(targVert);
								}
							}
						}
					}
				}
			}
		}
		
		genePath.addSources(allVertices);
		
		return genePath;
	}
	
	public int numEdges()
	{
		return edges.size();
	}
	
	public void setName(String newName)
	{
		name = newName;
	}
	
	public void setId(String newId)
	{
		id = newId;
	}
	
	public String toString()
	{
		return id + " :: " + edges;
	}

	/**
	 * Stores the genes along the edge as well as the molecule at the
	 * end of the directed edge.
	 *
	 */
	private class PathEdge
	{
		HashSet<String> genes;
		String target;
		int id;
		
		/**
		 * @param geneSet
		 * @param targ
		 */
		PathEdge(HashSet<String> geneSet, String targ, int id)
		{
			// Copy the genes into a new set
			genes = new HashSet<String>(geneSet);
			target = targ;
			this.id = id;
		}
		
		public String toString()
		{
			return genes + " : " + target;
		}
	}
}