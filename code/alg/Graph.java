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

import java.io.File;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.Collection;
import java.util.Collections;
import java.util.GregorianCalendar;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Observable;
import java.util.Observer;
import java.util.PriorityQueue;
import java.util.Set;
import java.util.Stack;
import java.util.concurrent.Callable;
import java.util.concurrent.ConcurrentHashMap;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.PriorityBlockingQueue;
import java.util.zip.GZIPOutputStream;

/**
 * A weighted mixed graph containing vertices of type V and
 * edges of type E.  Two vertices v1 and v2 may be connected by
 * up to three edges: v1 -> v2, v2 -> v1, and v1 -- v2.  The graph
 * contains sources and targets and can generate paths from
 * sources to targets.
 * <p>
 * A graph can also be used to represent a complex pathway,
 * such as KEGG pathways.
 */

public class Graph implements Observer{

	// Two different maps are needed because it's possible that v1 may
	// have two unique edges to v2.  The edge map associated with v1
	// cannot contain the key v2 twice
	/** A map to store edges whose direction cannot be changed */
	protected Map<Vertex, Map<Vertex, Edge>> dirMap;
	/** A map to store edges that can be directed and are initially undirected */
	protected Map<Vertex, Map<Vertex, Edge>> undirMap;

	/** The set of vertices that function as sources */
	protected ArrayList<Vertex> sources;
	/** The set of vertices that function as targets */
	protected Set<Vertex> targets;

	/** The set of directed edges */
	protected ArrayList<DirEdge> dirEdges;
	/** The set of undirected edges */
	protected ArrayList<UndirEdge> undirEdges;

	/** A map from vertices to degree */
	protected Map<Vertex, Integer> degreeCache;
	
	/** Track the number of times a search branch is aborted early and at what
	 * depth it is aborted */
	protected BigInteger[] aborts;

	public Graph()
	{
		dirMap = new HashMap<Vertex, Map<Vertex, Edge>>();
		undirMap = new HashMap<Vertex, Map<Vertex, Edge>>();

		sources = new ArrayList<Vertex>();
		targets = new HashSet<Vertex>();

		dirEdges = new ArrayList<DirEdge>();
		undirEdges = new ArrayList<UndirEdge>();

		degreeCache = new HashMap<Vertex, Integer>();
	}

	public void addSource(Vertex v)
	{
		sources.add(v);
	}

	public void addSources(Collection<Vertex> v)
	{
		sources.addAll(v);
	}

	public void clearSources()
	{
		sources.clear();
	}

	public void addTarget(Vertex v)
	{
		targets.add(v);
	}

	/**
	 * Add a single directed edge to the graph
	 * @param e
	 * @throws IllegalArgumentException if the edge already exists
	 */
	public void addEdge(DirEdge e)
	{
		addEdge(e, dirMap);
		dirEdges.add(e);
	}

	/**
	 * Add an undirected edge to the graph.
	 * Also adds the graph as an observer so that it is notified if the edge's
	 * orientation changes.
	 * @param e
	 * @throws IllegalArgumentException if the edge already exists
	 */
	public void addEdge(UndirEdge e)
	{
		addEdge(e, undirMap);
		undirEdges.add(e);
		e.addOrientationObserver(this);
	}


	/**
	 * Add edge e to the edge map.  Will add undirected edges to the map
	 * twice, from source to target and target to source.
	 * @param e
	 * @param edgeMap
	 */
	private void addEdge(Edge e, Map<Vertex, Map<Vertex, Edge>> edgeMap)
	{
		// Get the source and target of the edge
		// This source and target are not to be confused with the sources and
		// targets in the Graph
		Vertex source = e.getSource();
		Vertex target = e.getTarget();

		// No matter what the edge's orientation is, it will be used to connect
		// the source to the target
		if(!edgeMap.containsKey(source))
		{
			edgeMap.put(source, new HashMap<Vertex, Edge>());
		}

		Map<Vertex, Edge> vertMap = edgeMap.get(source);
		if(vertMap.containsKey(target))
		{
			throw new IllegalArgumentException("This type of edge from " + source.getName()
					+ " to " + target.getName() + " already exists");
		}

		vertMap.put(target, e);

		// Now check if the edge's orientation is UNORIENTED.  If so, add the edge
		// in the reverse direction as well
		if(e.getOrientation() == Edge.UNORIENTED)
		{
			if(!edgeMap.containsKey(target))
			{
				edgeMap.put(target, new HashMap<Vertex, Edge>());
			}

			vertMap = edgeMap.get(target);
			if(vertMap.containsKey(source))
			{
				throw new IllegalArgumentException("This type of edge from " + target.getName()
						+ " to " + source.getName() + " already exists");
			}

			vertMap.put(source, e);	
		}
	}

	/**
	 * Inserts all sources, targets, and edges from the other graph into this
	 * graph.  If an edges already exists in this graph, the original edge
	 * remains and the identical edge from the other graph is discarded.  Thus
	 * the weight of the original edge will be the weight after the merge.
	 * @param otherGraph
	 */
	public void mergeGraph(Graph otherGraph)
	{
		// Make sure we don't end up with two copies of a Vertex
		// in the list of sources
		for(Vertex otherSrc : otherGraph.sources)
		{
			if(!sources.contains(otherSrc))
			{
				sources.add(otherSrc);
			}
		}

		// The targets are a set
		for(Vertex otherTarg : otherGraph.targets)
		{
			targets.add(otherTarg);
		}

		// addEdge throws IllegalArgumentException when the edge already exists
		// so we need to catch and ignore these exceptions
		for(DirEdge otherDirEdge : otherGraph.dirEdges)
		{
			try
			{
				addEdge(otherDirEdge);
			}
			catch(IllegalArgumentException e)
			{System.err.println("Skipping edge " + otherDirEdge);}
		}
		for(UndirEdge otherUndirEdge : otherGraph.undirEdges)
		{
			try
			{
				addEdge(otherUndirEdge);
			}
			catch(IllegalArgumentException e)
			{System.err.println("Skipping edge " + otherUndirEdge);}
		}
	}

	/**
	 * Remove directed edge e
	 * @param e
	 */
	public void removeEdge(DirEdge e)
	{
		removeEdge(e, dirMap);
		dirEdges.remove(e);
	}

	/**
	 * Remove undirected edge e.
	 * Also removes this graph as an observer.
	 * @param e
	 */
	public void removeEdge(UndirEdge e)
	{
		removeEdge(e, undirMap);
		undirEdges.remove(e);
		e.removeOrientationObserver(this);
	}
	
	/**
	 * Remove edge e.
	 * Also removes this graph as an observer if the edge
	 * is undirected.
	 * @param e
	 */
	public void removeEdge(Edge e)
	{
		if(e instanceof DirEdge)
		{
			removeEdge((DirEdge) e);
		}
		else
		{
			removeEdge((UndirEdge) e);
		}
	}

	/**
	 * Remove edge e from the edge map if it is present in the edge map.
	 * @param e
	 * @param edgeMap
	 */
	private void removeEdge(Edge e, Map<Vertex, Map<Vertex, Edge>> edgeMap)
	{
		// Get the source and target
		// The edge's source and target are not related to the Graph's sources
		// and targets
		Vertex source = e.getSource();
		Vertex target = e.getTarget();

		// First check if the edge is used to connect the source and target
		if(edgeMap.containsKey(source))
		{
			Map<Vertex, Edge> vertMap = edgeMap.get(source);
			// Verify that the source is connected to the target with this particular edge
			if(vertMap.containsKey(target) && vertMap.get(target) == e)
			{
				// If it is, remove the edge
				vertMap.remove(target);
			}
		}

		// If the edge was unoriented when it was added, it may have been used to
		// connect the target to the source as well.  Check for this.
		if(edgeMap.containsKey(target))
		{
			Map<Vertex, Edge> vertMap = edgeMap.get(target);
			// We know the target has neighbors, now see if the source is one of them
			// and the edge is used to make the connection
			if(vertMap.containsKey(source) && vertMap.get(source) == e)
			{
				// If it is, remove the edge
				vertMap.remove(source);
			}
		}
	}

	/**
	 * Refreshes the set of edges by deleting any existing references
	 * to the edge and re-adding the edge.  This method is useful when an
	 * edge's orientation has changed.  Directed edges' orientations are fixed
	 * so they do not need to be updated.
	 * @param e
	 */
	public void updateEdge(DirEdge e)
	{
		removeEdge(e, dirMap);
		addEdge(e, dirMap);
	}

	/**
	 * Refreshes the set of edges by deleting any existing references
	 * to the edge and re-adding the edge.  This method is useful when an
	 * edge's orientation has changed.  Directed edges' orientations are fixed
	 * so they do not need to be updated.  Does not affect this graph's status
	 * as an observer for the edge.
	 * @param e
	 */
	public void updateEdge(UndirEdge e)
	{
		removeEdge(e, undirMap);
		addEdge(e, undirMap);
	}


	// TODO correctness testing
	/**
	 * Remove all edges in the graph that are not in a path from the source containing
	 * at least minDepth edges and no more than maxDepth edges.
	 * Paths contain unique edges but not necessarily
	 * unique vertices because two vertices can be connected by an undirected edge
	 * and a directed edge.
	 * Paths may use unoriented edges in either direction but will respect the orientation
	 * of an edge that has been assigned a direction even if the direction is not fixed
	 * (as long as the graph's edges have been properly updated after every change
	 * in orientation).
	 * By default, does not require the paths to end at a target.
	 * @param minDepth
	 * @param maxDepth
	 */
	public void pruneGraphNaive(int minDepth, int maxDepth)
	{
		pruneGraphNaive(minDepth, maxDepth, true);
	}


	/**
	 * Remove all edges in the graph that are not in a path from the source containing
	 * at least minDepth edges and no more than maxDepth edges.
	 * Paths contain unique edges but not necessarily
	 * unique vertices because two vertices can be connected by an undirected edge
	 * and a directed edge.
	 * Paths may use unoriented edges in either direction but will respect the orientation
	 * of an edge that has been assigned a direction even if the direction is not fixed
	 * (as long as the graph's edges have been properly updated after every change
	 * in orientation).
	 * @param minDepth
	 * @param maxDepth
	 * @param ignoreTargs if true, the set of target vertices will be ignored.  Rather
	 * the min and maxDepth will be used exclusively to find paths from a source to any
	 * vertex in the depth range.  If false, paths must be in the specified depth
	 * range and end at a target.
	 */
	public void pruneGraphNaive(int minDepth, int maxDepth, boolean ignoreTargs)
	{
		long start = System.currentTimeMillis();

		HashSet<Edge> edgesToRemove = new HashSet<Edge>();
		edgesToRemove.addAll(getAllEdges());

		// Iterate through all sources, searching for all paths from
		// that source to any target (or node at the specified depth)
		for(Vertex s : sources)
		{
			// No need to continue searching if all edges have been used in at
			// least one path already
			if(edgesToRemove.size() > 0)
			{
				Stack<Edge> eStack = new Stack<Edge>();
				Stack<Vertex> vStack = new Stack<Vertex>();
				vStack.push(s);
				pruneGraphNaiveHelper(edgesToRemove, eStack, vStack, 0, minDepth, maxDepth, ignoreTargs);
			}
		}

		System.out.println("Pruning " + edgesToRemove.size() + " edges");

		// After all paths have been found, see if any edges should be removed
		// because they were not on any of the paths.
		for(Edge e : edgesToRemove)
		{
			// TODO can this be done without instanceof?
			if(e instanceof UndirEdge)
			{
				removeEdge((UndirEdge) e);
			}
			else
			{
				removeEdge((DirEdge) e);
			}
		}

		long stop = System.currentTimeMillis();
		System.out.println("Time (ms) to perform DFS pruning: " +
				(stop - start) + "\n");
	}

	/**
	 * Helper method for pruneGraph
	 * @param edgesToRemove
	 * @param eStack
	 * @param vStack
	 * @param curDepth
	 * @param minDepth
	 * @param maxDepth
	 * @param ignoreTargs
	 */
	private void pruneGraphNaiveHelper(HashSet<Edge> edgesToRemove,
			Stack<Edge> eStack, Stack<Vertex> vStack, int curDepth, int minDepth, int maxDepth,
			boolean ignoreTargs)
	{
		if(edgesToRemove.size() > 0 && curDepth < maxDepth)
		{
			// Find all neighbors connected via directed edges
			// and then all neighbors connected via undirected edges
			pruneGraphNaiveHelper(edgesToRemove, eStack, vStack, curDepth, minDepth, maxDepth, dirMap, ignoreTargs);
			pruneGraphNaiveHelper(edgesToRemove, eStack, vStack, curDepth, minDepth, maxDepth, undirMap, ignoreTargs);
		}
	}

	/**
	 * Helper method for pruneGraphHelper
	 * @param edgesToRemove
	 * @param eStack
	 * @param vStack
	 * @param curDepth
	 * @param minDepth
	 * @param maxDepth
	 * @param edgeMap
	 * @param ignoreTargs
	 */
	private void pruneGraphNaiveHelper(HashSet<Edge> edgesToRemove,
			Stack<Edge> eStack, Stack<Vertex> vStack, int curDepth, int minDepth, int maxDepth,
			Map<Vertex, Map<Vertex, Edge>> edgeMap, boolean ignoreTargs)
	{
		if(edgesToRemove.size() > 0 && curDepth < maxDepth)
		{	
			curDepth++;			
			Vertex curVertex = vStack.peek();

			// First visit all of the neighbors connected with the type
			// of edge contained in edgeMap
			if(edgeMap.containsKey(curVertex))
			{
				Map<Vertex, Edge> neighbors = edgeMap.get(curVertex);

				// At each directed neighbor, push the edge and neighbor and recurse.
				// Also remove all edges from the set of edgesToRemove if a path
				// in the specified depth range is found
				for(Map.Entry<Vertex, Edge> pair : neighbors.entrySet())
				{
					Edge curEdge = pair.getValue();
					Vertex nextVertex = pair.getKey();

					// Paths cannot contain cycles so do not add this edge if
					// it leads us to a vertex we have already visited
					// in this path
					if(!vStack.contains(nextVertex))
					{
						eStack.push(curEdge);
						vStack.push(nextVertex);

						// We want a path to this vertex if:
						// 1) The path is in the depth range
						// 2) (optional) the path ends at a target
						// We already know the current depth is <= the max depth
						if((curDepth >= minDepth) &&
								(ignoreTargs || targets.contains(nextVertex)))
						{
							// The edges in this valid path should not be pruned
							// from the network
							edgesToRemove.removeAll(eStack);
						}

						pruneGraphNaiveHelper(edgesToRemove, eStack, vStack, curDepth,
								minDepth, maxDepth, ignoreTargs);

						eStack.pop();
						vStack.pop();
					}
				} // end loop through neighbors
			} // check vertex is in map
		} // check depth and edges remaining
	}


	/**
	 * Remove all edges in the graph that are more than maxDepth away from a source
	 * Paths contain unique edges but not necessarily
	 * unique vertices because two vertices can be connected by an undirected edge
	 * and a directed edge.
	 * Paths may use unoriented edges in either direction but will respect the orientation
	 * of an edge that has been assigned a direction even if the direction is not fixed
	 * (as long as the graph's edges have been properly updated after every change
	 * in orientation).
	 * Does not require the paths to end at a target.
	 * @param maxDepth
	 */
	public void pruneGraph(long maxDepth)
	{
		// Use BFS to get the distance of every vertex from its closest source
		HashMap<Vertex, Integer> vertexDepth = bfs();
		
		ArrayList<DirEdge> dirToRemove = new ArrayList<DirEdge>();
		for(DirEdge e : dirEdges)
		{
			// The source depth tells the number of edges that are needed to
			// get from a graph source to this edge's source. 
			// Adding 1 gives this edge's depth.
			if(((long) vertexDepth.get(e.getSource()) + 1) > maxDepth)
			{
				dirToRemove.add(e);
			}
		}
		for(DirEdge e : dirToRemove)
		{
			removeEdge(e);
		}
		
		ArrayList<UndirEdge> undirToRemove = new ArrayList<UndirEdge>();
		for(UndirEdge e : undirEdges)
		{
			// The min depth tells the number of edges that are needed to
			// get from a source to this edge.  Adding 1 gives this edge's
			// depth.  Use the source and target depth because
			// these edges are undirected.
			int srcDepth = vertexDepth.get(e.getSource());
			int targDepth = vertexDepth.get(e.getTarget());
			if((Math.min((long) srcDepth, (long) targDepth) + 1) > maxDepth)
			{
				undirToRemove.add(e);
			}
		}
		for(UndirEdge e : undirToRemove)
		{
			removeEdge(e);
		}
		
		System.out.println("BFS pruning with length " + maxDepth + 
				" removed "+ dirToRemove.size() + " directed edges" +
				" and " + undirToRemove.size() + " undirected edges");
	}

	// TODO update javadoc
	/**
	 * Paths contain unique edges but not necessarily
	 * unique vertices because two vertices can be connected by an undirected edge
	 * and a directed edge.
	 * Paths may use unoriented edges in either direction but will respect the orientation
	 * of an edge that has been assigned a direction even if the direction is not fixed
	 * (as long as the graph's edges have been properly updated after every change
	 * in orientation).
	 */
	public HashMap<Vertex, Integer> bfs()
	{
		long start = System.currentTimeMillis();

		HashMap<Vertex, Integer> vertexDepth = new HashMap<Vertex, Integer>();
		for(Vertex v : getVertices())
		{
			vertexDepth.put(v, Integer.MAX_VALUE);
		}

		// Iterate through all sources, initiating a BFS at each one
		for(Vertex s : sources)
		{
			HashSet<Vertex> visited = new HashSet<Vertex>();
			HashSet<Vertex> active = new HashSet<Vertex>();
			active.add(s);
			vertexDepth.put(s,0);

			bfsHelper(vertexDepth, active, visited);
		}

		long stop = System.currentTimeMillis();
		System.out.println("Time (ms) to perform BFS: " +
				(stop - start) + "\n");

//		for(Vertex v : vertexDepth.keySet())
//		{
//			System.out.println(v + "\t" + vertexDepth.get(v));
//		}

		return vertexDepth;
	}

	// TODO javadoc
	/**
	 * 
	 */
	private void bfsHelper(
			HashMap<Vertex, Integer> vertexDepth, HashSet<Vertex> active,
			HashSet<Vertex> visited)
	{
		if(active.size() > 0)
		{
			HashSet<Vertex> nextActive = new HashSet<Vertex>();

			for(Vertex curVertex : active)
			{
				// Should be the same for all active vertices at this level
				int curDepth = vertexDepth.get(curVertex);

				// First visit all of the neighbors connected with directed edges
				if(dirMap.containsKey(curVertex))
				{
					// At each directed neighbor, TODO document
					for(Vertex neighbor : dirMap.get(curVertex).keySet())
					{
						if(curDepth + 1 < vertexDepth.get(neighbor))
						{
							vertexDepth.put(neighbor, curDepth + 1);
						}

						// Only need to visit a vertex if it has not been visited previously
						if(!visited.contains(neighbor) && !active.contains(neighbor))
						{
							nextActive.add(neighbor);
						}
					} // end loop through neighbors
				} // check vertex is in map


				// Next visit all of the neighbors connected with undirected edges
				if(undirMap.containsKey(curVertex))
				{
					// At each directed neighbor, TODO document
					for(Vertex neighbor : undirMap.get(curVertex).keySet())
					{
						if(curDepth + 1 < vertexDepth.get(neighbor))
						{
							vertexDepth.put(neighbor, curDepth + 1);
						}

						// Only need to visit a vertex if it has not been visited previously
						if(!visited.contains(neighbor) && !active.contains(neighbor))
						{
							nextActive.add(neighbor);
						}
					} // end loop through neighbors
				} // check vertex is in map

				visited.add(curVertex);
			} // iterate through active vertices

			active = null;
			bfsHelper(vertexDepth, nextActive, visited);
		} // check if there are active vetices
	}



	/**
	 * Counts all possible paths from sources to targets.  Paths must use less than
	 * or equal to maxDepth edges.  Paths contain unique edges but not necessarily
	 * unique vertices because two vertices can be connected by an undirected edge
	 * and a directed edge.
	 * Paths may use unoriented edges in either direction but will respect the orientation
	 * of an edge that has been assigned a direction even if the direction is not fixed
	 * (as long as the graph's edges have been properly updated after every change
	 * in orientation).  Faster than finding all paths and then counting them.
	 * @param maxDepth
	 * @return
	 */
	public BigInteger countPaths(int maxDepth)
	{
		return countPaths(maxDepth, false);
	}


	/**
	 * Counts all possible paths from sources to targets or vertices at a specified
	 * depth.  Paths must use less than
	 * or equal to maxDepth edges.  Paths contain unique edges but not necessarily
	 * unique vertices because two vertices can be connected by an undirected edge
	 * and a directed edge.
	 * Paths may use unoriented edges in either direction but will respect the orientation
	 * of an edge that has been assigned a direction even if the direction is not fixed
	 * (as long as the graph's edges have been properly updated after every change
	 * in orientation).  Faster than finding all paths and then counting them.
	 * @param maxDepth
	 * @param ignoreTargs if true, the set of target vertices will be ignored.  Rather
	 * the maxDepth will be used exclusively to find paths from a source to any
	 * vertex at the maximum depth.
	 * @return
	 */
	public BigInteger countPaths(int maxDepth, boolean ignoreTargs)
	{
		BigInteger count = new BigInteger("0");

		// Iterate through all sources, searching for all paths from
		// that source to any target
		for(Vertex s : sources)
		{
			Stack<Vertex> vStack = new Stack<Vertex>();
			vStack.push(s);
			count = countPathHelper(count, vStack, maxDepth, ignoreTargs);
		}

		return count;
	}

	/**
	 * Helper method for countPaths
	 * @param count
	 * @param vStack
	 * @param depthLeft
	 * @param ignoreTargs
	 * @return BigInteger the updated count
	 */
	private BigInteger countPathHelper(BigInteger count,
			Stack<Vertex> vStack, int depthLeft, boolean ignoreTargs)
	{
		if(depthLeft > 0)
		{
			// Find all neighbors connected via directed edges
			// and then all neighbors connected via undirected edges
			count = countPathHelper(count, vStack, depthLeft, dirMap, ignoreTargs);
			count = countPathHelper(count, vStack, depthLeft, undirMap, ignoreTargs);
		}

		return count;
	}

	/**
	 * Helper method for countPathHelper
	 * @param count
	 * @param vStack
	 * @param depthLeft
	 * @param edgeMap
	 * @param ignoreTargs
	 * @return BigInteger the updated count
	 */
	private BigInteger countPathHelper(BigInteger count,
			Stack<Vertex> vStack, int depthLeft,
			Map<Vertex, Map<Vertex, Edge>> edgeMap, boolean ignoreTargs)
	{
		// depthLeft is the number of edges that can still be appended to the
		// path while respecting the max depth
		if(depthLeft > 0)
		{		
			Vertex curVertex = vStack.peek();

			// First visit all of the neighbors connected with the type
			// of edge contained in edgeMap
			if(edgeMap.containsKey(curVertex))
			{
				Map<Vertex, Edge> neighbors = edgeMap.get(curVertex);

				// At each directed neighbor, push the edge and neighbor and recurse.
				// Also add a new path if the neighbor is a target.
				for(Map.Entry<Vertex, Edge> pair : neighbors.entrySet())
				{
					Vertex nextVertex = pair.getKey();

					// Paths cannot contain cycles so do not add this edge if
					// it leads us to a vertex we have already visited
					// in this path
					if(!vStack.contains(nextVertex))
					{
						vStack.push(nextVertex);

						// We want a path to this vertex if either:
						// 1) we're ignoring target vertices and found a path
						//  of depth maxDepth
						// 2) we're using targets and found a target
						if((ignoreTargs && depthLeft == 1) ||
								(!ignoreTargs && targets.contains(nextVertex)))
						{
							// We found a new path so add to the counter
							count = count.add(BigInteger.ONE);
						}

						count = countPathHelper(count, vStack, depthLeft - 1, ignoreTargs);

						vStack.pop();
					}
				} // end loop through neighbors
			} // check vertex is in map
		} // check depth

		return count;
	}


	/**
	 * Find all possible paths from sources to targets.  Paths must use less than
	 * or equal to maxDepth edges.  Paths contain unique edges but not necessarily
	 * unique vertices because two vertices can be connected by an undirected edge
	 * and a directed edge.
	 * Paths may use unoriented edges in either direction but will respect the orientation
	 * of an edge that has been assigned a direction even if the direction is not fixed
	 * (as long as the graph's edges have been properly updated after every change
	 * in orientation).
	 * @param maxDepth
	 * @return
	 */
	public ArrayList<Path> findPaths(int maxDepth)
	{
		return findPaths(maxDepth, false);
	}


	/**
	 * Find all possible paths from sources to targets or vertices at a specified
	 * depth.  Paths must use less than
	 * or equal to maxDepth edges.  Paths contain unique edges but not necessarily
	 * unique vertices because two vertices can be connected by an undirected edge
	 * and a directed edge.
	 * Paths may use unoriented edges in either direction but will respect the orientation
	 * of an edge that has been assigned a direction even if the direction is not fixed
	 * (as long as the graph's edges have been properly updated after every change
	 * in orientation).
	 * @param maxDepth
	 * @param ignoreTargs if true, the set of target vertices will be ignored.  Rather
	 * the maxDepth will be used exclusively to find paths from a source to any
	 * vertex at the maximum depth.
	 * @return
	 */
	public ArrayList<Path> findPaths(int maxDepth, boolean ignoreTargs)
	{
		ArrayList<Path> paths = new ArrayList<Path>();

		// Iterate through all sources, searching for all paths from
		// that source to any target
		for(Vertex s : sources)
		{
			Stack<Edge> eStack = new Stack<Edge>();
			Stack<Vertex> vStack = new Stack<Vertex>();
			vStack.push(s);
			findPathHelper(paths, eStack, vStack, maxDepth, ignoreTargs);
		}

		System.err.println("Updating path statistics");
		long statStart = System.currentTimeMillis();

		// After all paths have been found and created, update their cached
		// statistics
		beginStatsUpdate();
		for(Path p : paths)
		{
			p.updateCachedStats();
		}
		endStatsUpdate();
		long statEnd = System.currentTimeMillis();
		System.err.println((statEnd - statStart) + " ms to update path stats");

		return paths;
	}

	/**
	 * Helper method for findPaths
	 * @param paths
	 * @param eStack
	 * @param vStack
	 * @param depthLeft
	 * @param ignoreTargs
	 */
	private void findPathHelper(ArrayList<Path> paths,
			Stack<Edge> eStack, Stack<Vertex> vStack, int depthLeft, boolean ignoreTargs)
	{
		if(depthLeft > 0)
		{
			// Find all neighbors connected via directed edges
			// and then all neighbors connected via undirected edges
			findPathHelper(paths, eStack, vStack, depthLeft, dirMap, ignoreTargs);
			findPathHelper(paths, eStack, vStack, depthLeft, undirMap, ignoreTargs);
		}
	}

	/**
	 * Helper method for findPathHelper
	 * @param paths
	 * @param eStack
	 * @param vStack
	 * @param depthLeft
	 * @param edgeMap
	 * @param ignoreTargs
	 */
	private void findPathHelper(ArrayList<Path> paths,
			Stack<Edge> eStack, Stack<Vertex> vStack, int depthLeft,
			Map<Vertex, Map<Vertex, Edge>> edgeMap, boolean ignoreTargs)
	{
		// depthLeft is the number of edges that can still be appended to the
		// path while respecting the max depth
		if(depthLeft > 0)
		{		
			Vertex curVertex = vStack.peek();

			// First visit all of the neighbors connected with the type
			// of edge contained in edgeMap
			if(edgeMap.containsKey(curVertex))
			{
				Map<Vertex, Edge> neighbors = edgeMap.get(curVertex);

				// At each directed neighbor, push the edge and neighbor and recurse.
				// Also add a new path if the neighbor is a target.
				for(Map.Entry<Vertex, Edge> pair : neighbors.entrySet())
				{
					Edge curEdge = pair.getValue();
					Vertex nextVertex = pair.getKey();

					// Paths cannot contain cycles so do not add this edge if
					// it leads us to a vertex we have already visited
					// in this path
					if(!vStack.contains(nextVertex))
					{
						eStack.push(curEdge);
						vStack.push(nextVertex);

						// We want a path to this vertex if either:
						// 1) we're ignoring target vertices and found a path
						//  of depth maxDepth
						// 2) we're using targets and found a target
						if((ignoreTargs && depthLeft == 1) ||
								(!ignoreTargs && targets.contains(nextVertex)))
						{
							// TODO here we could not add the path if it's weight
							// is below the threshold
							Path newPath = new Path(this, eStack, vStack);
							paths.add(newPath);

							// TODO temporarily output how many paths have been found
							// to monitor path finding on large instances
							if(paths.size() % 10000 == 0)
							{
								Calendar now = new GregorianCalendar();
								System.err.println("Paths found: " + paths.size() + " at "
										+ (now.get(Calendar.MONTH) + 1) + "/"
										+ now.get(Calendar.DATE) + "/"
										+ now.get(Calendar.YEAR) + " "
										+ now.get(Calendar.HOUR_OF_DAY) + ":"
										+ now.get(Calendar.MINUTE));
							}
						}

						findPathHelper(paths, eStack, vStack, depthLeft - 1, ignoreTargs);

						eStack.pop();
						vStack.pop();
					}
				} // end loop through neighbors
			} // check vertex is in map
		} // check depth
	}


	/**
	 * Find paths from sources to targets.  Paths must use less than
	 * or equal to maxDepth edges.  Paths contain unique edges but not necessarily
	 * unique vertices because two vertices can be connected by an undirected edge
	 * and a directed edge.
	 * Paths may use unoriented edges in either direction but will respect the orientation
	 * of an edge that has been assigned a direction even if the direction is not fixed
	 * (as long as the graph's edges have been properly updated after every change
	 * in orientation).  If more than maxPaths paths exist, only the top maxPaths paths
	 * will be returned as determined by the path weight.  Will not exhaustively search
	 * all branches if it can be determined that all potential paths down a certain
	 * path will have weights less than current set of top paths.
	 * @param maxDepth
	 * @param maxPaths
	 * @return
	 */
	public PriorityQueue<Path> findTopPaths(int maxDepth, int maxPaths)
	{
		return findTopPaths(maxDepth, maxPaths, false);
	}

	/**
	 * Find all possible paths from sources to targets or vertices at a specified
	 * depth.  Paths must use less than
	 * or equal to maxDepth edges.  Paths contain unique edges but not necessarily
	 * unique vertices because two vertices can be connected by an undirected edge
	 * and a directed edge.
	 * Paths may use unoriented edges in either direction but will respect the orientation
	 * of an edge that has been assigned a direction even if the direction is not fixed
	 * (as long as the graph's edges have been properly updated after every change
	 * in orientation).  If more than maxPaths paths exist, only the top maxPaths paths
	 * will be returned as determined by the path weight.  Will not exhaustively search
	 * all branches if it can be determined that all potential paths down a certain
	 * path will have weights less than current set of top paths.
	 * @param maxDepth
	 * @param maxPaths
	 * @param ignoreTargs if true, the set of target vertices will be ignored.  Rather
	 * the maxDepth will be used exclusively to find paths from a source to any
	 * vertex at the maximum depth.
	 * @return
	 */
	public PriorityQueue<Path> findTopPaths(int maxDepth, int maxPaths, boolean ignoreTargs)
	{
		// PathWeight is the current ranking metric that makes the most sense to use
		// at this point because some others are dependent on the other paths that have
		// been found already which means the set of paths returned would be highly
		// dependent on the order of network traversal.  Degree-based metrics
		// are less useful for sorting paths.
		PriorityQueue<Path> paths = new PriorityQueue<Path>(1001, Path.getComparator("PathWeight"));

		// Iterate through all sources, searching for all paths from
		// that source to any target
		for(Vertex s : sources)
		{
			// Reset the counts for every new source
			aborts = new BigInteger[maxDepth];
			for(int d = 0; d < maxDepth; d++)
			{
				aborts[d] = BigInteger.ZERO;
			}
			
			Stack<Edge> eStack = new Stack<Edge>();
			Stack<Vertex> vStack = new Stack<Vertex>();
			vStack.push(s);
			findTopPathHelper(paths, eStack, vStack, s.getWeight(), maxDepth, maxPaths, ignoreTargs);
		
			System.out.print("Aborts for source " + s);
			for(int d = maxDepth - 1; d >= 0; d--)
			{
				System.out.print("\td=" + (maxDepth - d) + ":" + aborts[d]);
			}
			System.out.println();
		}

		System.err.println("Updating path statistics");
		long statStart = System.currentTimeMillis();

		// After all paths have been found and created, update their cached
		// statistics
		beginStatsUpdate();
		for(Path p : paths)
		{
			p.updateCachedStats();
		}
		endStatsUpdate();
		long statEnd = System.currentTimeMillis();
		System.err.println((statEnd - statStart) + " ms to update path stats");

		return paths;
	}

	/**
	 * Helper method for findTopPaths
	 * @param paths
	 * @param eStack
	 * @param vStack
	 * @param stackWeight the product of all weights of the edges and vertices on the stacks
	 * @param depthLeft
	 * @param maxPaths
	 * @param ignoreTargs
	 */
	private void findTopPathHelper(PriorityQueue<Path> paths,
			Stack<Edge> eStack, Stack<Vertex> vStack, double stackWeight,
			int depthLeft, int maxPaths, boolean ignoreTargs)
	{
		if(depthLeft > 0)
		{
			// Find all neighbors connected via directed edges
			// and then all neighbors connected via undirected edges
			findTopPathHelper(paths, eStack, vStack, stackWeight, depthLeft, maxPaths, dirMap, ignoreTargs);
			findTopPathHelper(paths, eStack, vStack, stackWeight, depthLeft, maxPaths, undirMap, ignoreTargs);
		}
	}

	/**
	 * Helper method for findTopPathHelper
	 * @param paths
	 * @param eStack
	 * @param vStack
	 * @param stackWeight the product of all weights of the edges and vertices on the stacks
	 * @param depthLeft
	 * @param maxPaths
	 * @param edgeMap
	 * @param ignoreTargs
	 */
	private void findTopPathHelper(PriorityQueue<Path> paths,
			Stack<Edge> eStack, Stack<Vertex> vStack, double stackWeight,
			int depthLeft, int maxPaths,
			Map<Vertex, Map<Vertex, Edge>> edgeMap, boolean ignoreTargs)
	{
		// depthLeft is the number of edges that can still be appended to the
		// path while respecting the max depth
		if(depthLeft > 0)
		{		
			Vertex curVertex = vStack.peek();

			// First visit all of the neighbors connected with the type
			// of edge contained in edgeMap
			if(edgeMap.containsKey(curVertex))
			{
				Map<Vertex, Edge> neighbors = edgeMap.get(curVertex);

				// At each directed neighbor, push the edge and neighbor and recurse.
				// Also add a new path if the neighbor is a target.
				for(Map.Entry<Vertex, Edge> pair : neighbors.entrySet())
				{
					Edge curEdge = pair.getValue();
					Vertex nextVertex = pair.getKey();

					// Paths cannot contain cycles so do not add this edge if
					// it leads us to a vertex we have already visited
					// in this path
					if(!vStack.contains(nextVertex))
					{
						eStack.push(curEdge);
						vStack.push(nextVertex);
						
						// The new stack weight includes the weights of all edges
						// and vertices that were previously on the stack and
						// adds the weight of the edge and vertex most recently
						// placed on the stack.  Does not include any target weights.
						double newStackWeight = stackWeight * curEdge.getWeight() * nextVertex.getWeight();
						
						// We want a path to this vertex if either:
						// 1) we're ignoring target vertices and found a path
						//  of depth maxDepth
						// 2) we're using targets and found a target
						if((ignoreTargs && depthLeft == 1) ||
								(!ignoreTargs && targets.contains(nextVertex)))
						{
							Path newPath = new Path(this, eStack, vStack);
							
							// If we have fewer than the maximum number of paths,
							// add the path.  Otherwise check whether this path
							// ranks in the top maxPaths number of paths
							if(paths.size() < maxPaths)
							{
								paths.offer(newPath);
							}
							else
							{
								// The lowest-scoring top-ranked path
								Path lowPath = paths.peek();
								
								// Remove the old lowest-scoring path if
								// the new path scores higher
								if(paths.comparator().compare(newPath, lowPath) > 0)
								{
									paths.poll();
									paths.offer(newPath);
									
									// Remove references to this path so that the conflict
									// edges and edge use statistics do not consider it
									lowPath.deletePath();
								}
								else
								{
									// If the new path should not be added to the priority queue
									// remove reference to the new path so that the conflict
									// edges and edge use statistics do not consider it
									newPath.deletePath();
								}
							}
						}  // end conditional where a potential new top path is found

						// The recursive call should be made if the maximum number
						// of paths has not been found yet or if there is a possibility
						// that the branch in the search tree may yield higher weight
						// paths
						if(!abortBranch(maxPaths, paths, newStackWeight, depthLeft - 1))
						{
							findTopPathHelper(paths, eStack, vStack, newStackWeight, depthLeft - 1, maxPaths, ignoreTargs);
						}
						else
						{
							aborts[depthLeft-1] = aborts[depthLeft-1].add(BigInteger.ONE);
						}
						
						eStack.pop();
						vStack.pop();
					}
				} // end loop through neighbors
			} // check vertex is in map
		} // check depth
	}
	
	/**
	 * Determine whether findTopPaths should continue searching the current branch
	 * in the search tree.  The search will be aborted if it is guaranteed that
	 * no possible paths in the current branch will be added to the set of top paths.
	 * @param maxPaths the maximum number of top weight paths to be tracked
	 * @param paths the priority queue of top-ranked paths
	 * @param stackWeight the weight of all edges and vertices on the stacks in the current
	 * search tree branch
	 * @param depthLeft the number of edges that can still be added to the path so that
	 * it is not longer than the maximum allowable path length
	 * @return
	 */
	private boolean abortBranch(int maxPaths, PriorityQueue<Path> paths, double stackWeight, int depthLeft)
	{
		// Stop searching if the maximum depth has been reached
		if(depthLeft < 1)
		{
			return true;
		}
		
		// Continue searching if we have not found the maximum number of paths desired yet
		if(paths.size() < maxPaths)
		{
			return false;
		}
		
		// The path weight of the path from the set of top paths that has the lowest weight.
		// Abort the search if it is impossible for the current branch to find a path
		// with weight higher than this.
		double weightThreshold = paths.peek().getMetricValue("PathWeight");
		
		// The link includes adding another edge and vertex.  Optimal refers
		// to best possible case in the branch of the search tree.
		double optimalLinkWeight = Edge.getMaxWeight() * Vertex.getMaxWeight();
		// The best case weight of a path that could be found in this branch
		// of the search tree
		double optimalWithoutTarget = stackWeight;		
		
		if(optimalLinkWeight < 1)
		{
			// If the best case next edge + next vertex will decrease the path weight,
			// then the highest possible path weight will only add one more edge and vertex
			optimalWithoutTarget *= optimalLinkWeight;
		}
		else
		{
			// If the link weight improves the path weight, the highest possible path
			// weight could use the maximum number of edges
			optimalWithoutTarget *= Math.pow(optimalLinkWeight, depthLeft);
		}
		
		// A path will always end at a single target
		double optimalWeight = optimalWithoutTarget * Vertex.getMaxTargetWeight();
		
		// Abort the search if the best possible path we could find will not be placed
		// in the set of top-ranked paths
		return optimalWeight <= weightThreshold;
	}
	

	
	
	
	
	
	
	

	// Double check the following methods carefully
	// Variables for the multithreaded methods
	/** Read only copy of the sources */
	protected List<Vertex> sourcesMt;
	/** Read only copy of the targets */
	protected Set<Vertex> targetsMt;
	/** Read only copy of the edges*/
	protected Map<Vertex, Map<Vertex, Edge>> dirMapMt, undirMapMt;
	/** A map from sources to the number of branches that were terminated early */
	protected ConcurrentHashMap<Vertex, ArrayList<BigInteger>> abortsMt;
	/** A global heap that tracks the top M paths across all sources */
	protected PriorityBlockingQueue<Path> pathsMt;
	/** A lock for the findTopPathsMt method */
	protected Object findTopPathsLock = new Object();
	
	/** A map from targets to writers that will write the paths that end at that target */
	protected ConcurrentHashMap<Vertex, PrintWriter> writersMt;
	
	/**
	 * A class that initiates a DFS with a maximum number of top paths to be stored.
	 *
	 */
	public class TopPathsRun implements Callable
	{
		private Vertex runSource;
		private int runMaxDepth, runMaxPaths;
		private boolean runIgnoreTargs;
		
		public TopPathsRun(Vertex source, int maxDepth, int maxPaths, boolean ignoreTargs)
		{
			runSource = source;
			runMaxDepth = maxDepth;
			runMaxPaths = maxPaths;
			runIgnoreTargs = ignoreTargs;
		}
		
		/**
		 * Initiate the DFS from the specified source
		 * @return true if the DFS completed without an Exception
		 */
		public Boolean call() throws Exception
		{
			System.out.println("Starting DFS from " + runSource);
			
			// Initialize the counts for the new source
			ArrayList<BigInteger> srcAborts = new ArrayList<BigInteger>(runMaxDepth);
			for(int d = 0; d < runMaxDepth; d++)
			{
				srcAborts.add(BigInteger.ZERO);
			}
			abortsMt.put(runSource, srcAborts);
			
			Stack<Edge> eStack = new Stack<Edge>();
			Stack<Vertex> vStack = new Stack<Vertex>();
			vStack.push(runSource);
			findTopPathHelperMt(runSource, eStack, vStack, runSource.getWeight(), runMaxDepth, runMaxPaths, runIgnoreTargs);
			
			// Return true if the search terminated without an exception
			return true;
		}
	}
	
	/**
	 * Find paths from sources to targets.  Paths must use less than
	 * or equal to maxDepth edges.  Paths contain unique edges but not necessarily
	 * unique vertices because two vertices can be connected by an undirected edge
	 * and a directed edge.  Multithreaded version.
	 * Paths may use unoriented edges in either direction but will respect the orientation
	 * of an edge that has been assigned a direction even if the direction is not fixed
	 * (as long as the graph's edges have been properly updated after every change
	 * in orientation).  If more than maxPaths paths exist, only the top maxPaths paths
	 * will be returned as determined by the path weight.  Will not exhaustively search
	 * all branches if it can be determined that all potential paths down a certain
	 * path will have weights less than current set of top paths.
	 * @param maxDepth
	 * @param maxPaths
	 * @return
	 */
	public PriorityBlockingQueue<Path> findTopPathsMt(int maxDepth, int maxPaths)
	{
		return findTopPathsMt(maxDepth, maxPaths, false);
	}

	/**
	 * Find all possible paths from sources to targets or vertices at a specified
	 * depth.  Paths must use less than
	 * or equal to maxDepth edges.  Paths contain unique edges but not necessarily
	 * unique vertices because two vertices can be connected by an undirected edge
	 * and a directed edge.
	 * Paths may use unoriented edges in either direction but will respect the orientation
	 * of an edge that has been assigned a direction even if the direction is not fixed
	 * (as long as the graph's edges have been properly updated after every change
	 * in orientation).  If more than maxPaths paths exist, only the top maxPaths paths
	 * will be returned as determined by the path weight.  Will not exhaustively search
	 * all branches if it can be determined that all potential paths down a certain
	 * path will have weights less than current set of top paths.
	 * Multithreaded version.
	 * @param maxDepth
	 * @param maxPaths
	 * @param ignoreTargs if true, the set of target vertices will be ignored.  Rather
	 * the maxDepth will be used exclusively to find paths from a source to any
	 * vertex at the maximum depth.
	 * @return
	 */
	public PriorityBlockingQueue<Path> findTopPathsMt(int maxDepth, int maxPaths, boolean ignoreTargs)
	{
		try
		{
			// Each time this method is called, freeze the graph state so that multiple
			// threads can access it without locking in a read-only manner
			sourcesMt = Collections.unmodifiableList(sources);
			targetsMt = Collections.unmodifiableSet(targets);
			dirMapMt = Collections.unmodifiableMap(dirMap);
			undirMapMt = Collections.unmodifiableMap(undirMap);
			
			abortsMt = new ConcurrentHashMap<Vertex, ArrayList<BigInteger>>();

			// PathWeight is the current ranking metric that makes the most sense to use
			// at this point because some others are dependent on the other paths that have
			// been found already which means the set of paths returned would be highly
			// dependent on the order of network traversal.  Degree-based metrics
			// are less useful for sorting paths.
			pathsMt = new PriorityBlockingQueue<Path>(1001, Path.getComparator("PathWeight"));

			// Setup the DFS from each source
			Collection<Callable<Boolean>> dfsRuns = new ArrayList<Callable<Boolean>>(sourcesMt.size());
			for(Vertex s : sourcesMt)
			{
				dfsRuns.add(new TopPathsRun(s, maxDepth, maxPaths, ignoreTargs));
			}

			// Setup the ExcecutorService using one thread per processor available.
			// Because of lock contention when trying to modify the heap, it may not
			// help performance greatly to add more threads or even more processors past
			// a certain point.
			int numProcs = Runtime.getRuntime().availableProcessors();
			System.out.println("Number of processors available to the Java Virtual Machine: " + numProcs);
			ExecutorService executor = Executors.newFixedThreadPool(numProcs);

			// Submit the DFS runs and block until they complete
			List<Future<Boolean>> results = executor.invokeAll(dfsRuns);
			// Shutdown the threads in the thread pool
			executor.shutdown();

			// Verify that all of the runs completed successfully
			for(Future<Boolean> result : results)
			{
				if(!result.get())
				{
					throw new IllegalStateException("At least one DFS did not run to completion");
				}
			}

			// Display how many braches were terminated early
			for(Vertex s : sourcesMt)
			{
				System.out.print("Aborts for source " + s);
				for(int d = maxDepth - 1; d >= 0; d--)
				{
					System.out.print("\td=" + (maxDepth - d) + ":" + abortsMt.get(s).get(d));
				}
				System.out.println();
			}

			System.err.println("Updating path statistics");
			long statStart = System.currentTimeMillis();

			// After all paths have been found and created, update their cached
			// statistics
			beginStatsUpdate();
			for(Path p : pathsMt)
			{
				p.updateCachedStats();
			}
			endStatsUpdate();
			long statEnd = System.currentTimeMillis();
			System.err.println((statEnd - statStart) + " ms to update path stats");

			return pathsMt;
		}
		// TODO Reconsider how to handle these exceptions
		catch (InterruptedException e)
		{
			System.err.println("Multithreaded DFS was interrupted");
			throw new RuntimeException(e);
		}
		catch (ExecutionException e)
		{
			System.err.println("Multithreaded DFS exited with an exception");
			throw new RuntimeException(e);
		}
	}

	/**
	 * Helper method for findTopPathsMt.  Multithreaded version.
	 * @param source
	 * @param eStack
	 * @param vStack
	 * @param stackWeight the product of all weights of the edges and vertices on the stacks
	 * @param depthLeft
	 * @param maxPaths
	 * @param ignoreTargs
	 */
	private void findTopPathHelperMt(Vertex source,
			Stack<Edge> eStack, Stack<Vertex> vStack, double stackWeight,
			int depthLeft, int maxPaths, boolean ignoreTargs)
	{
		if(depthLeft > 0)
		{
			// Find all neighbors connected via directed edges
			// and then all neighbors connected via undirected edges
			findTopPathHelperMt(source, eStack, vStack, stackWeight, depthLeft, maxPaths, dirMapMt, ignoreTargs);
			findTopPathHelperMt(source, eStack, vStack, stackWeight, depthLeft, maxPaths, undirMapMt, ignoreTargs);
		}
	}

	/**
	 * Helper method for findTopPathHelperMt.  Multithreaded version.
	 * @param source
	 * @param eStack
	 * @param vStack
	 * @param stackWeight the product of all weights of the edges and vertices on the stacks
	 * @param depthLeft
	 * @param maxPaths
	 * @param edgeMap
	 * @param ignoreTargs
	 */
	private void findTopPathHelperMt(Vertex source,
			Stack<Edge> eStack, Stack<Vertex> vStack, double stackWeight,
			int depthLeft, int maxPaths,
			Map<Vertex, Map<Vertex, Edge>> edgeMap, boolean ignoreTargs)
	{
		// depthLeft is the number of edges that can still be appended to the
		// path while respecting the max depth
		if(depthLeft > 0)
		{		
			Vertex curVertex = vStack.peek();

			// First visit all of the neighbors connected with the type
			// of edge contained in edgeMap
			if(edgeMap.containsKey(curVertex))
			{
				Map<Vertex, Edge> neighbors = edgeMap.get(curVertex);

				// At each directed neighbor, push the edge and neighbor and recurse.
				// Also add a new path if the neighbor is a target and is in the
				// top-ranked paths.
				for(Map.Entry<Vertex, Edge> pair : neighbors.entrySet())
				{
					Edge curEdge = pair.getValue();
					Vertex nextVertex = pair.getKey();

					// Paths cannot contain cycles so do not add this edge if
					// it leads us to a vertex we have already visited
					// in this path
					if(!vStack.contains(nextVertex))
					{
						eStack.push(curEdge);
						vStack.push(nextVertex);
						
						// The new stack weight includes the weights of all edges
						// and vertices that were previously on the stack and
						// adds the weight of the edge and vertex most recently
						// placed on the stack.  Does not include any target weights.
						double newStackWeight = stackWeight * curEdge.getWeight() * nextVertex.getWeight();
						
						// We want a path to this vertex if either:
						// 1) we're ignoring target vertices and found a path
						//  of depth maxDepth
						// 2) we're using targets and found a target
						//
						// AND the path's weight places in in the set of top-ranked
						// paths we've found so far
						if((ignoreTargs && depthLeft == 1) ||
								(!ignoreTargs && targetsMt.contains(nextVertex)))
						{
							// Don't want multiple threads to try to insert a path
							// into a full heap at the same time.
							// This block must be the only code that writes
							// to the heap, although other code is allowed to
							// read from it without acquiring this lock,
							// understanding that its state may change
							// at any time.
							// Do not create or delete paths outside of this block either
							// becaues the Path class is not thread-safe.
							synchronized(findTopPathsLock)
							{
								Path newPath = new Path(this, eStack, vStack);

								// If we have fewer than the maximum number of paths,
								// add the path.  Otherwise check whether this path
								// ranks in the top maxPaths number of paths
								if(pathsMt.size() < maxPaths)
								{
									pathsMt.offer(newPath);
								}
								else
								{
									// The lowest-scoring top-ranked path
									Path lowPath = pathsMt.peek();

									// Remove the old lowest-scoring path if
									// the new path scores higher
									if(pathsMt.comparator().compare(newPath, lowPath) > 0)
									{
										// Add the new path first because otherwise there is
										// a possibility that another thread will be in abortBranchMt
										// after the removal and believe that the threshold
										// for aborting the search is higher than it actually is,
										// leading it to abort when it shouldn't
										// Offering before polling may mean that the threads
										// reading the heap state may see more paths than there actually
										// are, but because we already have maxPaths paths in the heap
										// at this point we know that the heap is full
										// whether its size returns maxPaths or maxPaths+1 to the reading
										// thread.
										pathsMt.offer(newPath);
										pathsMt.poll();

										// Remove references to this path so that the conflict
										// edges and edge use statistics do not consider it
										lowPath.deletePath();
									}
									else
									{
										// If the new path should not be added to the priority queue
										// remove reference to the new path so that the conflict
										// edges and edge use statistics do not consider it
										newPath.deletePath();
									}
								}
							} // end synchronized block
						}  // end conditional where a potential new top path is found

						// The recursive call should be made if the maximum number
						// of paths has not been found yet or if there is a possibility
						// that the branch in the search tree may yield higher weight
						// paths.
						// Because the heap may be modified while another thread executes this method,
						// the method may not abort sometimes when it could have.  It will never
						// abort when it should not have though.
						if(!abortBranchMt(maxPaths, newStackWeight, depthLeft - 1))
						{
							findTopPathHelperMt(source, eStack, vStack, newStackWeight, depthLeft - 1, maxPaths, ignoreTargs);
						}
						else
						{
							ArrayList<BigInteger> srcAborts = abortsMt.get(source);
							srcAborts.set(depthLeft-1, srcAborts.get(depthLeft-1).add(BigInteger.ONE));
						}
						
						eStack.pop();
						vStack.pop();
					}
				} // end loop through neighbors
			} // check vertex is in map
		} // check depth
	}
	
	/**
	 * Determine whether findTopPathsMt should continue searching the current branch
	 * in the search tree.  The search will be aborted if it is guaranteed that
	 * no possible paths in the current branch will be added to the set of top paths.
	 * Multithreaded version, and because the top paths heap may be modified while
	 * this method executes it may return false when it could have returned true.
	 * However, it is always okay to continue the search even if it could have been aborted
	 * so this is safe.
	 * @param maxPaths the maximum number of top weight paths to be tracked
	 * @param stackWeight the weight of all edges and vertices on the stacks in the current
	 * search tree branch
	 * @param depthLeft the number of edges that can still be added to the path so that
	 * it is not longer than the maximum allowable path length
	 * @return true if the branch can safely be terminated
	 */
	private boolean abortBranchMt(int maxPaths, double stackWeight, int depthLeft)
	{
		// Stop searching if the maximum depth has been reached
		if(depthLeft < 1)
		{
			return true;
		}
		
		// Continue searching if we have not found the maximum number of paths desired yet
		if(pathsMt.size() < maxPaths)
		{
			return false;
		}
		
		// The path weight of the path from the set of top paths that has the lowest weight.
		// Abort the search if it is impossible for the current branch to find a path
		// with weight higher than this.
		double weightThreshold = pathsMt.peek().getMetricValue("PathWeight");
		
		// The link includes adding another edge and vertex.  Optimal refers
		// to best possible case in the branch of the search tree.
		// TODO could cache these values along with the rest of the graph state
		double optimalLinkWeight = Edge.getMaxWeight() * Vertex.getMaxWeight();
		// The best case weight of a path that could be found in this branch
		// of the search tree
		double optimalWithoutTarget = stackWeight;		
		
		if(optimalLinkWeight < 1)
		{
			// If the best case next edge + next vertex will decrease the path weight,
			// then the highest possible path weight will only add one more edge and vertex
			optimalWithoutTarget *= optimalLinkWeight;
		}
		else
		{
			// If the link weight improves the path weight, the highest possible path
			// weight could use the maximum number of edges
			optimalWithoutTarget *= Math.pow(optimalLinkWeight, depthLeft);
		}
		
		// A path will always end at a single target
		double optimalWeight = optimalWithoutTarget * Vertex.getMaxTargetWeight();
		
		// Abort the search if the best possible path we could find will not be placed
		// in the set of top-ranked paths
		return optimalWeight <= weightThreshold;
	}
	
	
	/**
	 * A class that initiates a DFS from a singe source that writes
	 * all paths to a target.
	 *
	 */
	public class StorePathsRun implements Callable
	{
		private Vertex runSource;
		private int runMaxDepth;
		
		public StorePathsRun(Vertex source, int maxDepth)
		{
			runSource = source;
			runMaxDepth = maxDepth;
		}
		
		/**
		 * Initiate the DFS from the specified source
		 * @return true if the DFS completed without an Exception
		 */
		public Boolean call() throws Exception
		{
			long start = System.currentTimeMillis();
			System.out.println("Starting DFS from " + runSource);
			
			Stack<Edge> eStack = new Stack<Edge>();
			Stack<Vertex> vStack = new Stack<Vertex>();
			vStack.push(runSource);
			// The source weight is the initial weight of all vertices and edges on the stacks
			storePathHelperMt(runSource, eStack, vStack, runSource.getWeight(), runMaxDepth);
			
			long stop = System.currentTimeMillis();
			System.out.println("DFS from " + runSource + " took " + (stop - start) + "ms");
			
			// Return true if the search terminated without an exception
			return true;
		}
	}


	/**
	 * Writes all paths from sources to targets to files.  
	 * Paths to a specific target are written in their own file.  Paths must use less than
	 * or equal to maxDepth edges.  Paths contain unique edges but not necessarily
	 * unique vertices because two vertices can be connected by an undirected edge
	 * and a directed edge.  Multithreaded.
	 * Paths may use unoriented edges in either direction but will respect the orientation
	 * of an edge that has been assigned a direction even if the direction is not fixed
	 * (as long as the graph's edges have been properly updated after every change
	 * in orientation).
	 * @param outDir the directory in which the path files will be written
	 * @param maxDepth
	 */
	public void storePathsMt(String outDir, int maxDepth)
	{
		try
		{
			// Each time this method is called, freeze the graph state so that multiple
			// threads can access it without locking in a read-only manner
			sourcesMt = Collections.unmodifiableList(sources);
			targetsMt = Collections.unmodifiableSet(targets);
			dirMapMt = Collections.unmodifiableMap(dirMap);
			undirMapMt = Collections.unmodifiableMap(undirMap);

			if(!outDir.endsWith(File.separator))
			{
				outDir = outDir + File.separator;
			}
			
			// Create the writers for each target
			writersMt = new ConcurrentHashMap<Vertex, PrintWriter>();
			for(Vertex t : targetsMt)
			{
				// Write compressed files
				PrintWriter tWriter = new PrintWriter(new GZIPOutputStream(new FileOutputStream(outDir + t + ".txt.gz")));
				writersMt.put(t, tWriter);
			}

			// Setup the DFS from each source
			Collection<Callable<Boolean>> dfsRuns = new ArrayList<Callable<Boolean>>(sourcesMt.size());
			for(Vertex s : sourcesMt)
			{
				dfsRuns.add(new StorePathsRun(s, maxDepth));
			}

			// Setup the ExcecutorService using one thread per processor available.
			// If execution becomes IO-bound there may be little speedup beyond a
			// certain number of processors.
			int numProcs = Runtime.getRuntime().availableProcessors();
			System.out.println("Number of processors available to the Java Virtual Machine: " + numProcs);
			ExecutorService executor = Executors.newFixedThreadPool(numProcs);

			// Submit the DFS runs and block until they complete
			List<Future<Boolean>> results = executor.invokeAll(dfsRuns);
			// Shutdown the threads in the thread pool
			executor.shutdown();

			// Verify that all of the runs completed successfully
			for(Future<Boolean> result : results)
			{
				if(!result.get())
				{
					throw new IllegalStateException("At least one DFS did not run to completion");
				}
			}
			
			// Close all of the writers
			for(PrintWriter w : writersMt.values())
			{
				w.close();
			}
		}
		catch(IOException e)
		{
			e.printStackTrace();
		}
		catch(InterruptedException e)
		{
			System.err.println("Multithreaded DFS was interrupted");
			throw new RuntimeException(e);
		}
		catch(ExecutionException e)
		{
			System.err.println("Multithreaded DFS exited with an exception");
			throw new RuntimeException(e);
		}
	}

	/**
	 * Helper method for storePathsMt.  Multithreaded.
	 * @param source
	 * @param eStack
	 * @param vStack
	 * @param stackWeight the product of all weights of the edges and vertices on the stacks
	 * @param depthLeft
	 */
	private void storePathHelperMt(Vertex source,
			Stack<Edge> eStack, Stack<Vertex> vStack, double stackWeight,
			int depthLeft)
	{
		if(depthLeft > 0)
		{
			// Find all neighbors connected via directed edges
			// and then all neighbors connected via undirected edges
			storePathHelperMt(source, eStack, vStack, stackWeight, depthLeft, dirMapMt);
			storePathHelperMt(source, eStack, vStack, stackWeight, depthLeft, undirMapMt);
		}
	}

	/**
	 * Helper method for storePathHelperMt.  Multithreaded.
	 * @param source
	 * @param eStack
	 * @param vStack
	 * @param stackWeight the product of all weights of the edges and vertices on the stacks
	 * @param depthLeft
	 * @param edgeMap
	 */
	private void storePathHelperMt(Vertex source,
			Stack<Edge> eStack, Stack<Vertex> vStack, double stackWeight,
			int depthLeft, Map<Vertex, Map<Vertex, Edge>> edgeMap)
	{
		// depthLeft is the number of edges that can still be appended to the
		// path while respecting the max depth
		if(depthLeft > 0)
		{		
			Vertex curVertex = vStack.peek();

			// First visit all of the neighbors connected with the type
			// of edge contained in edgeMap
			if(edgeMap.containsKey(curVertex))
			{
				Map<Vertex, Edge> neighbors = edgeMap.get(curVertex);

				// At each directed neighbor, push the edge and neighbor and recurse.
				// Also add a new path if the neighbor is a target
				for(Map.Entry<Vertex, Edge> pair : neighbors.entrySet())
				{
					Edge curEdge = pair.getValue();
					Vertex nextVertex = pair.getKey();

					// Paths cannot contain cycles so do not add this edge if
					// it leads us to a vertex we have already visited
					// in this path
					if(!vStack.contains(nextVertex))
					{
						eStack.push(curEdge);
						vStack.push(nextVertex);
						
						// The new stack weight includes the weights of all edges
						// and vertices that were previously on the stack and
						// adds the weight of the edge and vertex most recently
						// placed on the stack.  Does not include any target weights.
						double newStackWeight = stackWeight * curEdge.getWeight() * nextVertex.getWeight();
						
						// Write the path to this vertex if we found a target
						if(targetsMt.contains(nextVertex))
						{
							// Calculate the weight of the path
							Vertex target = nextVertex;
//							double pathWeight = newStackWeight * target.getTargetWeight();
							// Do not want the stored path weight to include any target weights
							// because these can change for different SDREM iterations
							double pathWeight = newStackWeight;

							// Write the path.  Lock using the writer so that no two threads
							// try to write to the same writer simultaneously
							PrintWriter tWriter = writersMt.get(target);
							synchronized(tWriter)
							{
								tWriter.println(new StringPath(eStack, vStack, pathWeight));
							} // end synchronized block

						}  // end conditional where a new path is found

						storePathHelperMt(source, eStack, vStack, newStackWeight, depthLeft - 1);
			
						eStack.pop();
						vStack.pop();
					}
				} // end loop through neighbors
			} // check vertex is in map
		} // check depth
	}
	
	
	/**
	 * 
	 * @param p the test path
	 * @return true if the graph contains the vertices in p in the
	 * same order as in p.  Compares vertices by name.
	 */
	public boolean containsPath(Path p)
	{
		return (matchingVertices(p) == p.getNumVertices());
	}


	// TODO use a map of vertex names to vertices to speed up this algorithm
	/**
	 * Search the graph for the longest possible subsequence of a Path p.
	 * Compares vertices by vertex name, not id.  Considers all possible
	 * start points in p.  Will work for paths that contain the same vertex
	 * multiple times if the corresponding cycles exist in the graph.
	 * @param p
	 * @return the number of vertices in p that can be found on a path
	 * in this graph (0 through p.length)
	 */
	public int matchingVertices(Path p)
	{
		int maxMatches = 0;

		// Iterate through all possible starting points in p
		for(int start = 0; start < p.getNumEdges(); start++)
		{
			// Check all vertices with either directed or undirected edges
			int matches = Math.max(matchingVerticesHelper(p.getVertices(), start, dirMap),
					matchingVerticesHelper(p.getVertices(), start, undirMap));

			if(matches > maxMatches)
			{
				maxMatches = matches;
			}
		}

		return maxMatches;
	}


	/**
	 * See if the starting point in the path has any outgoing edges in the graph
	 * @param path the vertices in the path
	 * @param pathInd the index of the vertex we are trying to find next
	 * @param edgeMap
	 * @return the number of vertices in the longest matching subpath
	 */
	private int matchingVerticesHelper(Vertex[] path, int pathInd,
			Map<Vertex, Map<Vertex, Edge>> edgeMap)
	{
		if(pathInd < path.length)
		{
			int maxMatches = 0;

			// Because we are matching by name, it is not possible to directly
			// find the edges in the maps that start at a particular vertex.  We must
			// iterate through all of them (this can be slow)
			for(Map.Entry<Vertex, Map<Vertex, Edge>> entry : edgeMap.entrySet())
			{
				Vertex curVertex = entry.getKey();
				// See if this vertex matches the current vertex in the path by name
				if(curVertex.sameName(path[pathInd]))
				{
					// Guaranteed that the curVertex will match path[pathInd], but let
					// the other helper method do the actual increment to the number
					// of matches for cleaner recursion
					int matches = matchingVerticesHelper(path, pathInd, curVertex, dirMap);
					if(matches > maxMatches)
					{
						maxMatches = matches;
					}

					matches = matchingVerticesHelper(path, pathInd, curVertex, undirMap);
					if(matches > maxMatches)
					{
						maxMatches = matches;
					}
				}
			}

			return maxMatches;
		}
		else
		{
			// There are no more possible starting points
			return 0;
		}
	}


	/**
	 * Check if the vertex in the Graph matches the vertex in the path,
	 * and if so follow the outgoing edges
	 * @param path the vertices in the path
	 * @param pathInd the index of the vertex we are trying to find next
	 * @param graphVertex the Vertex we are presently comparing with the vertex at
	 * path[pathInd]
	 * @param the set of edges to consider (directed or undirected)
	 * @return
	 */
	private int matchingVerticesHelper(Vertex[] path, int pathInd, Vertex graphVertex,
			Map<Vertex, Map<Vertex, Edge>> edgeMap)
	{
		if(pathInd < path.length)
		{
			// See if the path vertex matches the graph vertex
			if(graphVertex.sameName(path[pathInd]))
			{
				// We must have at least 1 match because
				// graphVertex matches path[pathInd]
				int maxMatches = 1;

				// If there was a match, get the neighbors and recurse
				// after making sure the graphVertex has neighbors
				if(edgeMap.containsKey(graphVertex))
				{
					for(Vertex neighbor : edgeMap.get(graphVertex).keySet())
					{
						// From the neighbor, we need to check it's neighbors along the
						// directed and undirected edges
						// We add 1 to the number of matches to account for the match of
						// graphVertex and path[pathInd]
						int matches =  1 + Math.max(
								matchingVerticesHelper(path, pathInd + 1, neighbor, dirMap),
								matchingVerticesHelper(path, pathInd + 1, neighbor, undirMap));

						if(matches > maxMatches)
						{
							maxMatches = matches;
						}
					}
				}

				return maxMatches;
			}
			else
			{
				// If they did not match, return 0 because the longest
				// matching subpath starting at pathInd has length 0
				return 0;
			}
		}
		else
		{
			// Base case
			return 0;
		}
	}

	
	/**
	 * Create Path objects from the collection of StringPaths
	 * @param storedPaths
	 * @return
	 */
	public ArrayList<Path> loadStoredPaths(Collection<StringPath> storedPaths)
	{
		ArrayList<Path> paths = new ArrayList<Path>();

		// For each stored path, find references to the Edge and Vertex
		// objects it uses.  Then use them to create a Path object.
		for(StringPath storedPath : storedPaths)
		{
			ArrayList<Edge> edges = new ArrayList<Edge>();
			Stack<Vertex> vertices = new Stack<Vertex>();
			
			for(String storedEdge : storedPath.getPath().split("\\|"))
			{
				String[] edgeParts = storedEdge.split(":");
				
				// Lookup the vertices of this edge
				Vertex v1 = Vertex.findVertex(edgeParts[0]);
				if(v1 == null)
				{
					throw new IllegalStateException("Cannot find vertex " + edgeParts[0]);
				}
				Vertex v2 = Vertex.findVertex(edgeParts[2]);
				if(v2 == null)
				{
					throw new IllegalStateException("Cannot find vertex " + edgeParts[2]);
				}
				
				// Add the vertices to the list
				if(vertices.size() == 0)
				{
					vertices.push(v1);
					vertices.push(v2);
				}
				else
				{
					// Confirm the first vertex of this edge is the same as the second
					// vertex of the last edge
					if(!v1.equals(vertices.peek()))
					{
						throw new IllegalStateException(storedPath.getPath() +
								" has disconnected consecutive edges");
					}
					vertices.push(v2);
				}
				
				// Find the edge and add it to the list
				Map<Vertex, Map<Vertex, Edge>> edgeMap;
				if(edgeParts[1].equals("dir"))
				{
					edgeMap = dirMap;
				}
				else if(edgeParts[1].equals("undir"))
				{
					edgeMap = undirMap;
				}
				else
				{
					throw new IllegalStateException(edgeParts[1] +
						" is not a recognized edge type");
				}
				
				if(!edgeMap.containsKey(v1))
				{
					throw new IllegalStateException("Cannot find edge " + storedEdge);
				}
				Map<Vertex, Edge> neighbors = edgeMap.get(v1);
				
				if(!neighbors.containsKey(v2))
				{
					throw new IllegalStateException("Cannot find edge " + storedEdge);
				}
				edges.add(neighbors.get(v2));
			} // Finish processing an edge from the StringPath
			
			Path newPath = new Path(this, edges, vertices);
			paths.add(newPath);
			
			// Verify the network used to store the paths is the same network
			// used now
			if(Math.abs(newPath.maxWeight() - storedPath.getWeight()) > 1e-10)
			{
				throw new IllegalStateException("Path did not have expected weight: " + 
						newPath + "\t" + newPath.maxWeight() + "\t" + storedPath.getWeight());
			}
		} // Finish loading a StringPath

		// After all paths have been loaded, update their cached statistics
		beginStatsUpdate();
		for(Path p : paths)
		{
			p.updateCachedStats();
		}
		endStatsUpdate();

		return paths;
	}

	/**
	 * An edge observed by the graph will notify the graph when its
	 * orientation changes so that the graph can update its map.
	 */
	public void update(Observable o, Object arg)
	{
		if(arg instanceof UndirEdge)
		{
			updateEdge((UndirEdge) arg);
		}
	}

	public ArrayList<DirEdge> getDirEdges()
	{
		return dirEdges;
	}

	public ArrayList<UndirEdge> getUndirEdges()
	{
		return undirEdges;
	}

	/**
	 * Get the outgoing directed edges from vertex v.
	 * @param v
	 * @return The outgoing edges or an empty list if
	 * the vertex is not in the graph
	 */
	public ArrayList<Edge> getDirEdges(Vertex v)
	{
		if(dirMap.containsKey(v))
		{
			return new ArrayList<Edge>(dirMap.get(v).values());
		}
		else
		{
			return new ArrayList<Edge>();
		}
	}

	/**
	 * Get the undirected edges that connect vertex v to some
	 * other vertex (and only the outgoing undirected edges if
	 * the undirected edges have been assigned an orientation)
	 * @param v
	 * @return The connected edges or an empty list if
	 * the vertex is not in the graph
	 */
	public ArrayList<Edge> getUndirEdges(Vertex v)
	{
		if(undirMap.containsKey(v))
		{
			return new ArrayList<Edge>(undirMap.get(v).values());
		}
		else
		{
			return new ArrayList<Edge>();
		}
	}

	/**
	 * 
	 * @return a new ArrayList containing all directed and undirected edges
	 */
	public ArrayList<Edge> getAllEdges()
	{
		ArrayList<Edge> all = new ArrayList<Edge>();
		all.addAll(dirEdges);
		all.addAll(undirEdges);
		return all;
	}

	public ArrayList<Vertex> getSources()
	{
		return sources;
	}

	public Set<Vertex> getTargets()
	{
		return targets;
	}

	// TODO couldn't a vertex appear in the returned list more than once?
	/**
	 * 
	 * @return a list of all vertices that have at least one
	 * outgoing edge or undirected edge (i.e. all keys in the edge maps)
	 */
	public ArrayList<Vertex> getInteriorVertices()
	{
		ArrayList<Vertex> vertices = new ArrayList<Vertex>();
		vertices.addAll(undirMap.keySet());
		vertices.addAll(dirMap.keySet());
		return vertices;		
	}

	// TODO couldn't a vertex appear in the returned list more than once?
	/**
	 * 
	 * @return a list of all vertices, obtained by taking both vertices
	 * from each edge
	 */
	public ArrayList<Vertex> getVertices()
	{
		ArrayList<Vertex> vertices = new ArrayList<Vertex>();
		for(Edge e : dirEdges)
		{
			vertices.add(e.getSource());
			vertices.add(e.getTarget());
		}
		for(Edge e : undirEdges)
		{
			vertices.add(e.getSource());
			vertices.add(e.getTarget());
		}
		return vertices;
	}


	/**
	 * Calculates the degree (in and out) for a vertex
	 * by iterating through all edges.  Does not use
	 * the cached degree.  The same neighbor connected by multiple edges
	 * (protein-protein and protein-DNA) will be counted multiple times.
	 * @param v
	 * @param byName if true, consider the degree of any
	 * vertices with the same name as this vertex instead
	 * of only those that are this exact Vertex object
	 * @return
	 */
	public int getDegree(Vertex v, boolean byName)
	{
		return getDegree(v, byName, false);
	}



	/**
	 * Calculates the degree (in and out) for a vertex
	 * by iterating through all edges.  Uses the cached
	 * value if one exists.  Do not rely on the cache value
	 * if edges have been added, oriented, or removed since
	 * the last time getDegree was called.
	 * The same neighbor connected by multiple edges
	 * (protein-protein and protein-DNA) will be counted multiple times.
	 * @param v
	 * @param byName if true, consider the degree of any
	 * vertices with the same name as this vertex instead
	 * of only those that are this exact Vertex object
	 * @param useCache if true, try to get the degree from
	 * the degree cache and only calculate it if the Vertex
	 * does not exist in the cache
	 * @return
	 */
	public int getDegree(Vertex v, boolean byName, boolean useCache)
	{
		if(useCache && degreeCache.containsKey(v))
		{
			return degreeCache.get(v).intValue();
		}

		int degree = 0;
		// Inefficient to allocate space for all edges, but 
		// this shouldn't be performance-critical
		for(Edge e : getAllEdges())
		{
			Vertex source = e.getSource();
			Vertex target = e.getTarget();

			if(byName)
			{
				if(v.sameName(source) || v.sameName(target))
				{
					degree++;
				}
			}
			else
			{
				if(v.equals(source) || v.equals(target))
				{
					degree++;
				}
			}
		}		

		degreeCache.put(v, Integer.valueOf(degree));

		return degree;
	}

	/**
	 * Calculates the degree (in and out) for a vertex
	 * by iterating through all edges.  Matches by name
	 * instead of id since only the name is given
	 * as a parameter.  Does not use a cache.
	 * The same neighbor connected by multiple edges
	 * (protein-protein and protein-DNA) will be counted multiple times.
	 * @param name
	 * @return
	 */
	public int getDegree(String name)
	{
		int degree = 0;
		// Inefficient to allocate space for all edges, but 
		// this shouldn't be performance-critical
		for(Edge e : getAllEdges())
		{
			Vertex source = e.getSource();
			Vertex target = e.getTarget();

			if(source.getName().equalsIgnoreCase(name) ||
					target.getName().equalsIgnoreCase(name))
			{
				degree++;
			}

		}		
		return degree;
	}

	// TODO should be called automatically
	/**
	 * Must be called after adding, removing, or orienting edges.
	 *
	 */
	public void clearDegreeCache()
	{
		degreeCache.clear();
	}


	/**
	 * Have each edge cache its number of satisfied forward and backward
	 * paths
	 *
	 */
	public void beginStatsUpdate()
	{
		for(Edge e : getAllEdges())
		{
			e.enableCachedSatisfiedPaths();
		}
	}

	// TODO combine with clearing vertex cache?
	/**
	 * Clear the edge cache for satisfied forward and backward paths
	 *
	 */
	public void endStatsUpdate()
	{
		for(Edge e : getAllEdges())
		{
			e.disableCachedSatisfiedPaths();
		}
	}


	/**
	 * Writes the undirected degree of each vertex.  Doesn't close the writer.
	 * @param writer
	 */
	public void writeUndirDegree(PrintWriter writer)
	{
		writer.println("Vertex\tDegree");
		for(Map.Entry<Vertex, Map<Vertex, Edge>> mapping : undirMap.entrySet())
		{
			writer.println(mapping.getKey() + "\t" + mapping.getValue().size());
		}
	}


	/**
	 * Writes a list of all directed and undirected edges.  Doesn't close the writer.
	 * @param writer
	 */
	public void writeEdges(PrintWriter writer)
	{
		for(DirEdge d : dirEdges)
		{
			writer.println(d.getSource() + "\t-->\t" + d.getTarget());
		}

		for(UndirEdge d : undirEdges)
		{
			writer.println(d.getSource() + "\t---\t" + d.getTarget());
		}
	}
	
	/**
	 * Writes a list of all directed and undirected edges.  Automatically creates
	 * a file writer
	 * @param filename
	 * @throws IOException
	 */
	public void writeEdges(String filename) throws IOException
	{
		PrintWriter writer = new PrintWriter(new FileWriter(filename));
		writeEdges(writer);
		writer.close();
	}

	/**
	 * Only prints undirected edges
	 */
	public void printEdges(String file) throws IOException
	{
		PrintWriter writer = new PrintWriter(new FileWriter(file));
		for(Vertex key : undirMap.keySet())
		{
			for(Vertex target : undirMap.get(key).keySet())
			{
				writer.println(key + "\t" + target + "\t" + undirMap.get(key).get(target).getOrientation());
			}
		}
		writer.close();
	}
}
