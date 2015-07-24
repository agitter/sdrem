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

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.Collection;
import java.util.GregorianCalendar;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Observable;
import java.util.Observer;
import java.util.Set;
import java.util.Stack;


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
							Path newPath = new Path(this, eStack, vStack);
							paths.add(newPath);


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
	 * 
	 * @param p the test path
	 * @return true if the graph contains the vertices in p in the
	 * same order as in p.  Compares vertices by name.
	 */
	public boolean containsPath(Path p)
	{
		return (matchingVertices(p) == p.getNumVertices());
	}


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
