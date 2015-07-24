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
import java.util.HashSet;
import java.util.List;
import java.util.Set;


public class DirEdge extends Edge {

	public DirEdge(Vertex v1, Vertex v2)
	{
		this(v1, v2, DEFAULT_WEIGHT);
	}

	public DirEdge(Vertex v1, Vertex v2, double weight)
	{
		// Directed edges are always oriented in the forward direction
		// and always fixed
		super(v1, v2, weight, FORWARD, true);
	}

	/**
	 * Associates this edge with a path
	 * @param p the path
	 * @param source the vertex the path requires to be the source of the edge
	 * i.e. the first vertex of the edge that will be encountered when following
	 * the path
	 */
	@Override
	public void assocPath(Path p, Vertex source)
	{
		if(source == v1)
		{
			fwdPaths.add(p);
		}
		else
		{
			throw new IllegalArgumentException(source.getName() + " is not the source vertex of this directed edge");
		}
	}

	/**
	 * Removes this edge's association with the path
	 * @param p the path
	 */
	@Override
	public void removePath(Path p)
	{
		fwdPaths.remove(p);
	}

	/**
	 * Returns the set of paths that want to use this edge in a different
	 * direction that path p
	 * @param p
	 * @return
	 */
	@Override
	public List<Path> getConflictingPaths(Path p)
	{
		// Directed edges can't be used in the backward direction so
		// all paths using them use them in the same direction
		return new ArrayList<Path>();
	}

	/**
	 * Returns the number of path-path pairs that disagree on the orienation
	 * this edge should take
	 * @return
	 */
	@Override
	public int countConflicts()
	{
		return 0;
	}

	/**
	 * Calculates the new global score - old global score if this edge
	 * were to be flipped.  Returns 0 if the edge direction is fixed or
	 * the edge is unoriented.
	 * @return
	 */
	@Override
	public double computeFlipDelta()
	{
		return 0;
	}

	/**
	 * @return true if this edge is used in one or more paths.
	 */
	@Override
	public boolean isUsed()
	{
		return (fwdPaths.size() > 0);
	}
	
	/**
	 * Considers all paths, not just paths that are connected from
	 * source to target.
	 * @return the number of paths that use this edge in the direction it
	 * is oriented (which must be the forward direction)
	 */
	@Override
	public int numConsistentPaths()
	{
		return numFwdPaths();
	}
	
	/**
	 * Considers only paths that are connected from
	 * source to target.
	 * @return the number of paths that use this edge in the direction it
	 * is oriented (which must be the forward direction)
	 */
	@Override
	public int numConsistentSatisfiedPaths()
	{
		// If the satisfied path cache is not -1, the cache is enabled
		// so return the cached value instead of recalculating it
		if(satisfiedPathCache != -1)
		{
			return satisfiedPathCache;
		}
		
		return numSatisfiedFwdPaths();
	}
	
	/**
	 * Considers all paths, not just paths that are connected from
	 * source to target.
	 * @return 0 because no path can use this edge in the BACKWARD direction
	 */
	@Override
	public int numInconsistentPaths()
	{
		return 0;
	}
	
	/**
	 * Return "dir" to specify that this is a directed edge
	 */
	@Override
	public String getType()
	{
		return "dir";
	}
}
