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
import java.util.Observable;
import java.util.Observer;



public class UndirEdge extends Edge {

	/**
	 * Track which paths wish to use the edge in the backward direction
	 * (from vertex 2 to vertex 1)
	 */
	protected ArrayList<Path> backPaths;
	
	/**
	 * Tracks how many times this edge's flip method has been called.
	 * This may not correspond to the number of actual flips since
	 * the flip method will not flip a fixed or unoriented edge.
	 */
	protected int flipCount;
	
	/**
	 * An index used to map this edge to its index in the set
	 * of linear program variables.  The id is typically (but not always)
	 * based on the 0-indexed position of the edge in the list
	 * of conflict edges.  The id should be non-negative.
	 */
	protected int id;
	
	/** Used to notify the Graph object when the orientation changes */
	protected OrientationPublisher publisher = new OrientationPublisher();
	
	public UndirEdge(Vertex v1, Vertex v2)
	{
		this(v1, v2, DEFAULT_WEIGHT);
	}
	
	public UndirEdge(Vertex v1, Vertex v2, double weight)
	{
		this(v1, v2, weight, UNORIENTED, false);
	}
	
	public UndirEdge(Vertex v1, Vertex v2, double weight, int orient, boolean fixed)
	{
		super(v1, v2, weight, orient, fixed);

		backPaths = new ArrayList<Path>();
		flipCount = 0;
		id = -1;
	}
	
	/**
	 * Associates this edge with a path.  Undirected edges can be associated with
	 * paths that wish to use them in the forward direction and others that wish
	 * to use them in the backward direction.
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
		else if(source == v2)
		{
			backPaths.add(p);
		}
		else
		{
			throw new IllegalArgumentException(source.getName() + " is not a vertex of this edge");
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
		backPaths.remove(p);
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
		if(fwdPaths.contains(p))
		{
			return backPaths;
		}
		else if(backPaths.contains(p))
		{
			return fwdPaths;
		}
		else
		{
			return new ArrayList<Path>();
		}
	}
	
	/**
	 * @return the number of path-path pairs that disagree on the orienation
	 * this edge should take
	 */
	@Override
	public int countConflicts()
	{
		return fwdPaths.size() * backPaths.size();
	}
	
	/**
	 * If the edge is not fixed, changes the direction
	 * from FORWARD to BACKWARD or vice versa.  Does not affect
	 * unoriented edges.
	 *
	 */
	public void flip()
	{
		flipCount++;
		
		if(orientation == FORWARD)
		{
			setOrientation(BACKWARD);
		}
		else if(orientation == BACKWARD)
		{
			setOrientation(FORWARD);
		}
	}
	
	public void resetFlipCount()
	{
		flipCount = 0;
	}
	
	public int getFlipCount()
	{
		return flipCount;
	}
	
	/**
	 * If the edge is not fixed, randomly sets the orientation
	 * to FORWARD or BACKWARD.
	 *
	 */
	public void randOrient()
	{
		if(Math.random() < 0.5)
		{
			setOrientation(FORWARD);
		}
		else
		{
			setOrientation(BACKWARD);
		}
	}

	/**
	 * If the edge is not fixed, set it to be oriented forward, oriented
	 * backward, or unoriented.
	 * Notifies all observers that the orientation was set.
	 * @param orient
	 */
	public void setOrientation(int orient)
	{
		if(!fixed)
		{
			orientation = orient;
			publisher.publish(this);
		}
	}
	
	// TODO need to take action with the paths that wanted to use the edge in the
	// other direction?  Do we allow a path to be fixed if there are conflicts?
	// Maybe use publishers?  But would need a lot because of all the paths.
	/**
	 * Fix the current edge orientation if it is forward or backward.
	 * An edge cannot be fixed as unoriented.
	 *
	 */
	public void fix()
	{
		if(orientation == UNORIENTED)
		{
			throw new IllegalArgumentException("Cannot fix an edge's direction as UNORIENTED");
		}
		
		fixed = true;
	}
	
	public void unfix()
	{
		fixed = false;
	}
	
	/**
	 * Fixes an edge's orientation if it is not already fixed and there are no
	 * conflicts on this edge.  If the edge isn't used, defaults to being fixed
	 * in the forward orientation.
	 * @return true if the edge orientation is fixed now but wasn't before
	 */
	public boolean fixIfNoConflicts()
	{
		if(fixed)
		{
			return false;
		}
		if(backPaths.isEmpty())
		{
			setOrientation(FORWARD);
			fix();
			return true;
		}
		else if(fwdPaths.isEmpty())
		{
			setOrientation(BACKWARD);
			fix();
			return true;
		}
		else
		{
			return false;
		}
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
		if(fixed || orientation == UNORIENTED)
		{
			return 0;
		}
		
		// Determine which paths are using the edge in its current
		// direction and which would use it in its flipped direction
		ArrayList<Path> curPaths, flipPaths;
		if(orientation == FORWARD)
		{
			curPaths = fwdPaths;
			flipPaths = backPaths;
		}
		else
		{
			curPaths = backPaths;
			flipPaths = fwdPaths;
		}
		
		double flipWeights = 0, curWeights = 0;
		
		// Calculate the sum of the path weights for paths using the edge
		// in its current direction
		for(Path p : curPaths)
		{
			curWeights += p.weight();
		}
		
		// Calculate the sum of the path weights for paths using the edge
		// in its flipped direction
		flip();
		for(Path p : flipPaths)
		{
			flipWeights += p.weight();
		}
		flip();
		// We want the flip count to only track external calls to flip
		flipCount -= 2;
		
		return flipWeights - curWeights;
	}
	
	/**
	 * @return true if this edge is used in one or more paths.
	 */
	@Override
	public boolean isUsed()
	{
		return ((fwdPaths.size() + backPaths.size()) > 0);
	}
	
	// TODO only allow the id to be set once?
	/**
	 * User must ensure that the new id is not already in use
	 * by another edge, this is not checked.
	 * @param newId the non-negative value to use for the new id
	 */
	public void setId(int newId)
	{
		if(newId < 0)
		{
			throw new IllegalArgumentException("Edge ids must be non-negative: " + newId);
		}
		
		id = newId;
	}
	
	public int getId()
	{
		return id;
	}
	
	/**
	 * @return the sum of the max weights of all paths that use this edge in
	 * the backward direction
	 */
	public double backPathMaxWeights()
	{
		double sum = 0;
		for(Path p : backPaths)
		{
			sum += p.maxWeight();
		}
		
		return sum;
	}
	
	/**
	 * @return the sum of the active weights of all paths that use this edge in
	 * the backward direction
	 */
	public double backPathWeights()
	{
		double sum = 0;
		for(Path p : backPaths)
		{
			sum += p.weight();
		}
		
		return sum;
	}
	
	
	/**
	 * 
	 * @return the number of paths that use this edge in the back direction
	 */
	public int numBackPaths()
	{
		return backPaths.size();
	}
	
	/**
	 * 
	 * @return the number of paths that use this edge in the back direction
	 * and are satisfied (every edge along them is directed in the needed
	 * direction to connect the source and target)
	 */
	public int numSatisfiedBackPaths()
	{
		int count = 0;
		for(Path p : backPaths)
		{
			if(p.isConnected())
			{
				count++;
			}
		}
		
		return count;
	}
	
	
	/**
	 * Considers all paths, not just paths that are connected from
	 * source to target.
	 * @return the number of paths that use this edge in the direction it
	 * is oriented, or the number of paths that use it in either direction
	 * if it is not oriented
	 */
	@Override
	public int numConsistentPaths()
	{
		if(orientation == UNORIENTED)
		{
			return numFwdPaths() + numBackPaths();
		}
		else if(orientation == FORWARD)
		{
			return numFwdPaths();
		}
		else if(orientation == BACKWARD)
		{
			return numBackPaths();
		}
		else
		{
			throw new IllegalStateException(orientation + " is not a valid orientation"
					+ " for edge " + this);
		}
	}
	
	/**
	 * Considers only paths that are connected from
	 * source to target.
	 * @return the number of paths that use this edge in the direction it
	 * is oriented, or the number of paths that use it in either direction
	 * if it is not oriented
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
		
		// If the cache is not enabled, calculate the number
		// of satisfied paths that use this edge in a direction
		// consistent with its current orientation.
		if(orientation == UNORIENTED)
		{
			return numSatisfiedFwdPaths() + numSatisfiedBackPaths();
		}
		else if(orientation == FORWARD)
		{
			return numSatisfiedFwdPaths();
		}
		else if(orientation == BACKWARD)
		{
			return numSatisfiedBackPaths();
		}
		else
		{
			throw new IllegalStateException(orientation + " is not a valid orientation"
					+ " for edge " + this);
		}
	}
	
	/**
	 * Considers all paths, not just paths that are connected from
	 * source to target.
	 * @return the number of paths that use this edge in the direction 
	 * opposite the one in which it
	 * is oriented, or 0 if it is not oriented
	 */
	@Override
	public int numInconsistentPaths()
	{
		if(orientation == UNORIENTED)
		{
			return 0;
		}
		else if(orientation == FORWARD)
		{
			return numBackPaths();
		}
		else if(orientation == BACKWARD)
		{
			return numFwdPaths();
		}
		else
		{
			throw new IllegalStateException(orientation + " is not a valid orientation"
					+ " for edge " + this);
		}
	}
	
	/**
	 * Return "undir" to specify that this is a undirected edge
	 */
	@Override
	public String getType()
	{
		return "undir";
	}
	
	/**
	 * Any observers will be notified when this edge's orientation changes.
	 * @param o
	 */
	public void addOrientationObserver(Observer o)
	{
		publisher.addObserver(o);
	}
	
	/**
	 * Any observers will be notified when this edge's orientation changes.
	 * @param o
	 */
	public void removeOrientationObserver(Observer o)
	{
		publisher.deleteObserver(o);
	}
	
	protected class OrientationPublisher extends Observable
	{
		public void publish(Edge e)
		{
		      setChanged();
		      notifyObservers(e);
		}
	}
}
