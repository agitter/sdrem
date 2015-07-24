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

import java.util.ArrayList;
import java.util.List;

/**
 * A weighted connection between two vertices.  Edges can
 * be directed or undirected.  Edge directions are allowed
 * to change unless the edge's direction is fixed.
 */

public abstract class Edge {
	
	protected static final double DEFAULT_WEIGHT = 1;
	/** Track the largest edge weight seen so far among all types of edges */
	protected static double maxWeight = Double.NEGATIVE_INFINITY;
	
	/** The edge is not oriented. */
	public static final int UNORIENTED = 0;
	/** The edge is oriented from vertex 1 to vertex 2 */
	public static final int FORWARD = 1;
	/** The edge is oriented from vertex 2 to vertex 1 */
	public static final int BACKWARD = 2;
	
	protected double weight;
	protected Vertex v1, v2;
	/** The edge's orientation: UNORIENTED, FORWARD, or BACKWARD */
	protected int orientation;
	/** The orientation cannot be changed */
	protected boolean fixed;
	
	/**
	 * Track which paths wish to use the edge in the forward direction
	 * (from vertex 1 to vertex 2)
	 */
	protected ArrayList<Path> fwdPaths;
	
	/**
	 * When the satisfied path cache is enabled, the number of satisfied
	 * paths going through this edge in the direction it is currently
	 * oriented is stored.  When the cache is disabled, -1 is stored.
	 */
	protected int satisfiedPathCache;
	
	protected Edge(Vertex v1, Vertex v2, double weight, int orient, boolean fixed)
	{
		this.v1 = v1;
		this.v2 = v2;
		this.weight = weight;
		orientation = orient;
		this.fixed = fixed;
		
		fwdPaths = new ArrayList<Path>();
		satisfiedPathCache = -1;
		
		if(weight > maxWeight)
		{
			maxWeight = weight;
		}
	}
	
	/**
	 * @return The largest weight that has been assigned to an edge object
	 */
	public static double getMaxWeight()
	{
		return maxWeight;
	}
	
	public double getWeight()
	{
		return weight;
	}
	
	/**
	 * Undirected edges arbitrarily use vertex 1 as the source.
	 * This does not affect their ability to be used in either direction
	 * by a path.
	 * @return
	 */
	public Vertex getSource()
	{
		if(orientation == BACKWARD)
		{
			return v2;
		}
		else
		{
			return v1;
		}
	}
	
	/**
	 * Undirected edges arbitrarily use vertex 2 as the target.
	 * This does not affect their ability to be used in either direction
	 * by a path.
	 * @return
	 */
	public Vertex getTarget()
	{
		if(orientation == BACKWARD)
		{
			return v1;
		}
		else
		{
			return v2;
		}
	}
	
	/**
	 * Returns whether the edge is being used in the forward or backward
	 * direction if this vertex is desired as the source
	 * @param source
	 * @return
	 */
	public int findDirection(Vertex source)
	{
		if(v1 == source)
		{
			return FORWARD;
		}
		else if(v2 == source)
		{
			return BACKWARD;
		}
		else
		{
			throw new IllegalArgumentException(source.getName() + " is not a vertex of this edge "
					+ printVertices());
		}
		
	}
	
	/** UNORIENTED, FORWARD, or BACKWARD */
	public int getOrientation()
	{
		return orientation;
	}
	
	public boolean isOriented()
	{
		return (!(orientation == UNORIENTED));
	}
	
	public boolean isFixed()
	{
		return fixed;
	}
	
	public boolean containsVertex(Vertex v)
	{
		return (v1 == v) || (v2 == v);
	}
	
	public String printVertices()
	{
		return getSource().getName() + ":" + getTarget().getName();
	}

	public String toString()
	{
		return printVertices();
	}
	
	/**
	 * @return the sum of the max weights of all paths that use this edge in
	 * the forward direction
	 */
	public double fwdPathMaxWeights()
	{
		double sum = 0;
		for(Path p : fwdPaths)
		{
			sum += p.maxWeight();
		}
		
		return sum;
	}
	
	
	/**
	 * @return the sum of the active weights of all paths that use this edge in
	 * the forward direction
	 */
	public double fwdPathWeights()
	{
		double sum = 0;
		for(Path p : fwdPaths)
		{
			sum += p.weight();
		}
		
		return sum;
	}
	
	/**
	 * 
	 * @return the number of paths that use this edge in the forward direction
	 */
	public int numFwdPaths()
	{
		return fwdPaths.size();
	}
	
	/**
	 * @return the number of paths that use this edge in the forward direction
	 * and are satisfied (every edge along them is directed in the needed
	 * direction to connect the source and target)
	 */
	public int numSatisfiedFwdPaths()
	{
		int count = 0;
		for(Path p : fwdPaths)
		{
			if(p.isConnected())
			{
				count++;
			}
		}
		
		return count;
	}
	
	/**
	 * It is recommended that this method is called before updating
	 * edge use statistics for a large number of paths.
	 * After calling this method once, the number of satisfied paths
	 * through this edge that use it in the direction it is oriented
	 * is calculated and cached.  To maintain correctness,
	 * disableCachedSatisfiedPaths must be called before changing
	 * any edge orientation because the cached value will
	 * no longer be accurate once the orientation changes.
	 *
	 */
	public void enableCachedSatisfiedPaths()
	{
		// Ensure the number will be calculated by resetting the
		// cached path value
		satisfiedPathCache = -1;
		satisfiedPathCache = numConsistentSatisfiedPaths();
	}
	
	/**
	 * This method clears the satisfied path cache
	 * and must be called before changing any orientations
	 * in the network.
	 *
	 */
	public void disableCachedSatisfiedPaths()
	{
		satisfiedPathCache = -1;
	}
	
	/**
	 * Associates this edge with a path
	 * @param p the path
	 * @param source the vertex the path requires to be the source of the edge,
	 * i.e. the first vertex of the edge that will be encountered when following
	 * the path
	 */
	public abstract void assocPath(Path p, Vertex source);
	
	/**
	 * Removes this edge's association with the path
	 * @param p the path
	 */
	public abstract void removePath(Path p);
	
	/**
	 * @param p
	 * @return the set of paths that want to use this edge in a different
	 * direction than path p
	 */
	public abstract List<Path> getConflictingPaths(Path p);
	
	/**
	 * @return the number of path-path pairs that disagree on the orienation
	 * this edge should take
	 */
	public abstract int countConflicts();
	
	/**
	 * Calculates the new global score - old global score if this edge
	 * were to be flipped.  Returns 0 if the edge direction is fixed or
	 * the edge is unoriented.
	 * @return
	 */
	public abstract double computeFlipDelta();
	
	/**
	 * @return true if this edge is used in one or more paths.
	 */
	public abstract boolean isUsed();
	
	/**
	 * Considers all paths, not just paths that are connected from
	 * source to target.
	 * @return the number of paths that use this edge in the direction it
	 * is oriented, or the number of paths that use it in either direction
	 * if it is not oriented
	 */
	public abstract int numConsistentPaths();
	
	/**
	 * Considers only paths that are connected from
	 * source to target.
	 * @return the number of paths that use this edge in the direction it
	 * is oriented, or the number of paths that use it in either direction
	 * if it is not oriented
	 */
	public abstract int numConsistentSatisfiedPaths();
	
	/**
	 * Considers all paths, not just paths that are connected from
	 * source to target.
	 * @return the number of paths that use this edge in the direction 
	 * opposite the one in which it
	 * is oriented, or 0 if it is not oriented
	 */
	public abstract int numInconsistentPaths();
	
	/**
	 * Return "dir" or "undir"
	 * @return
	 */
	public abstract String getType();
}
