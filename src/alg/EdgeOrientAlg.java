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
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Random;
import java.util.Scanner;
import java.util.Set;

import util.MapUtil;

/**
 * Runs the edge orientation algorithm of choice.
 *
 */
public class EdgeOrientAlg {

	private Graph graph;
	private ArrayList<Path> paths;
	/** A list of all edges whose orientation is in conflict. That
	 * is, at least one path wishes to use them in each direction. */
	private ArrayList<UndirEdge> conflictEdges;
	/** Stores the best found orientations for the conflict edges */
	private int[] savedOrientations;

	private static final int MAX_DEPTH = 5;
	/** The default number of times to restart the randomized orientation
	 * followed by local search in randPlusSearchSln */
	private static final int RAND_PLUS_SEARCH_RESTARTS = 10;

	/**
	 * @param g the Graph must have already been initialized with the
	 * edges and vertices
	 */
	public EdgeOrientAlg(Graph g)
	{
		graph = g;
	}


	/**
	 * Randomly orients the conflict edges with randomOrient and then
	 * performs local search with localSearchSln.  Repeats this
	 * process the default number of times.  Returns the global
	 * score of the best found configuration and sets the conflict edges
	 * to have this orientation.
	 * 
	 * @return
	 */
	public double randPlusSearchSln()
	{
		return randPlusSearchSln(RAND_PLUS_SEARCH_RESTARTS);
	}

	/**
	 * Randomly orients the conflict edges with randomOrient and then
	 * performs local search with localSearchSln.  Repeats this
	 * process the default number of times.  Returns the global
	 * score of the best found configuration and sets the conflict edges
	 * to have this orientation.
	 * 
	 * @return
	 */
	public double randPlusSearchSln(int iterations)
	{
		double bestGlobal = Double.NEGATIVE_INFINITY;

		for(int rep = 0; rep < iterations; rep++)
		{
			randomOrient();
			localSearchSln();

			// If a new best configuration is found, save it
			double newGlobal = globalScore();
			if(newGlobal > bestGlobal)
			{
				bestGlobal = newGlobal;
				saveConflictOrientations();
			}	
		}

		// Reload the best conflict edge orientations in case the caller
		// wishes to access it
		loadConflictOrientations();
		graphStateChanged();

		System.out.println("\nBest random + edge flip local search after " + 
				iterations + " iterations: " + bestGlobal);
		System.out.println("Max possible: " + maxGlobalScore() + "\n");

		return bestGlobal;		
	}


	/**
	 * Iteratively flips the edge
	 * that yields the larges global score increase until a local
	 * maximum is reached.  If any of the conflict edges has not been 
	 * oriented (it's orientation is Edge.UNORIENTED) all edges will
	 * be randomly oriented first.
	 * 
	 * Changes the state of the stored Graph.
	 * 
	 * 
	 * @return the global score
	 */
	public double localSearchSln()
	{
		// Randomly orient all conflict edges if any edges have not
		// been assigned an orientation
		// This automatically finds paths and the conflict edges if necessary
		if(!conflictEdgesOriented())
		{
			System.out.println("Conflict edges not oriented so performing " +
			"random orientation");
			randomOrient();
		}

		double global = globalScore();
		int counter = 0;

		// Don't try to flip edges if there are no conflict edges
		if(conflictEdges.size() > 0)
		{
			long start = System.currentTimeMillis();
			System.out.println("Beginning edge flip local search");
			double oldGlobal = Double.NEGATIVE_INFINITY;

			// Iterate as long as flipping an edge will yield an improvement
			while(oldGlobal < global)
			{
				counter++;
				oldGlobal = global;

				// Find the edge whose flip will yield the greatest
				// increase in global score
				UndirEdge bestEdge = conflictEdges.get(0);
				double bestEdgeDelta = Double.NEGATIVE_INFINITY;
				for(UndirEdge e : conflictEdges)
				{
					double curDelta = e.computeFlipDelta();
					if(bestEdgeDelta < curDelta)
					{
						bestEdge = e;
						bestEdgeDelta = curDelta;
					}
				}

				// Only perform the flip if it will give a positive change
				// in the global objective function
				if(bestEdgeDelta > 0)
				{
					bestEdge.flip();
					global += bestEdgeDelta;
				}

//				System.out.println("At iteration " + counter + " global score" +
//				" changed from " + oldGlobal + " to " + global);
			} // end local search loop

			long stop = System.currentTimeMillis();
			System.out.println("Finished local search");
			System.out.println("Time (ms): " + (stop - start) + "\n");
			
			graphStateChanged();

		} // end check for conflictEdges.size() > 0

		global = globalScore();
		System.out.println("\nEdge flip local search: " + global);
		System.out.println("Max possible: " + maxGlobalScore() + "\n");

		return global;
	}

	/**
	 * Assign edge orientations greedily by choosing the edge
	 * direction that supports the paths with the highest
	 * total weights.  Defaults to FORWARD if there is a tie.
	 * @return
	 */
	public double greedyEdgeSln()
	{
		// Find the paths and conflict edges if needed
		if(paths == null || conflictEdges == null)
		{
			findConflicts();
		}

		for(UndirEdge e : conflictEdges)
		{
			if(e.fwdPathMaxWeights() >= e.backPathMaxWeights())
			{
				e.setOrientation(Edge.FORWARD);
			}
			else
			{
				e.setOrientation(Edge.BACKWARD);
			}
		}
		
		graphStateChanged();

		// Report the results
		double global = globalScore();
		System.out.println("Greedy edge solution: " + global);
		System.out.println("Max possible: " + maxGlobalScore() + "\n");

		return global;
	}

	/**
	 * Randomly orients all conflict edges (undirected edges that
	 * have one or more paths that wish to use them in each direction).
	 * Path and conflict edges will be identified if they have not
	 * already been
	 *
	 */
	public void randomOrient()
	{
		if(paths == null || conflictEdges == null)
		{
			findConflicts();
		}

		// Randomly orient the edges that had conflicts
		for(UndirEdge e : conflictEdges)
		{
			e.randOrient();
			e.resetFlipCount();
		}
		
		graphStateChanged();
	}


	/**
	 * Stores all conflict edges in the graph and sets the class's
	 * conflictEdges variable.  Also fixes all edges that
	 * do not have conflicts.  Returns the number of edges that
	 * have conflicts.
	 *
	 */
	public int findConflicts()
	{
		if(paths == null)
		{
			findPaths();
		}

		long start = System.currentTimeMillis();
		System.out.println("Fixing edges without conflicts and finding conflict edges");
		ArrayList<UndirEdge> edges = graph.getUndirEdges();
		conflictEdges = new ArrayList<UndirEdge>();

		int fixCount = 0, usedCount = 0;
		for(UndirEdge e: edges)
		{
			if(e.isUsed())
			{
				usedCount++;
			}

			// fixIfNoConflicts() returns true iff the edge was not
			// already fixed but is now fixed because there were no conflicts.
			// We also want to count edges that were already fixed.
			if(e.isFixed() || e.fixIfNoConflicts())
			{
				fixCount++;
			}
			else
			{
				conflictEdges.add(e);
			}
		}
		long stop = System.currentTimeMillis();
		System.out.println(usedCount + " of " + edges.size() + " edges were used in at least one path");
		System.out.println(fixCount + " of " + edges.size() + " edges did not have conflicts");
		System.out.println("Time (ms): " + (stop - start) + "\n");

		return conflictEdges.size();
	}

	public ArrayList<UndirEdge> getConflictEdges()
	{
		if(conflictEdges == null)
		{
			findConflicts();
		}

		return conflictEdges;
	}

	/**
	 * 
	 * @return true iff all conflict edges have been oriented as either
	 * FORWARD or BACKWARD
	 */
	private boolean conflictEdgesOriented()
	{
		for(UndirEdge e: conflictEdges)
		{
			if(e.getOrientation() == Edge.UNORIENTED)
			{
				return false;
			}
		}

		return true;
	}


	/**
	 * Use the graph to find paths containing up to the default
	 * max depth number of edges
	 * from sources to targets.  The paths are stored by the
	 * edge orientation algorithm and by default not visible
	 * to the caller.  Sets the class's paths variable.
	 * @return the number of paths
	 *
	 */
	public int findPaths()
	{
		return findPaths(MAX_DEPTH);
	}


	/**
	 * Use the graph to find paths containing up to depth edges
	 * from sources to targets.  The paths are stored by the
	 * edge orientation algorithm and by default not visible
	 * to the caller.  Sets the class's paths variable.
	 * @param depth
	 * @retun the number of paths
	 */
	public int findPaths(int depth)
	{
		long start = System.currentTimeMillis();
		System.out.println("Finding paths up to depth " + depth);

		paths = graph.findPaths(depth);

		long stop = System.currentTimeMillis();
		System.out.println("Found " + paths.size() + " paths using depth " + depth);
		System.out.println("Time (ms): " + (stop - start) + "\n");

		return paths.size();
	}
	
	/**
	 * Return a reference to the list of paths.  Finds all
	 * paths first if they have not been found already.  Will return
	 * all paths, satisfied and unsatisfied.
	 * @return
	 */
	public ArrayList<Path> getPaths()
	{
		if(paths == null)
		{
			findPaths();
		}
		
		return paths;
	}
	
	
	/**
	 * Return a new list containing only satisfied (connected) paths.  Finds all
	 * paths first if they have not been found already.
	 * @return
	 */
	public ArrayList<Path> getSatisfiedPaths()
	{
		ArrayList<Path> satisfied = new ArrayList<Path>();
		
		for(Path p : getPaths())
		{
			if(p.isConnected())
			{
				satisfied.add(p);
			}
		}
		
		return satisfied;
	}
	
	/**
	 * Must be called after the state of the graph
	 * has been changed to maintain consistency
	 * among data structures.  Updates
	 * path edge use counts and clears the degree cache.
	 *
	 */
	public void graphStateChanged()
	{
		graph.clearDegreeCache();
		updateEdgeUses();
	}
	
	/**
	 * Update the edge use statistics for all paths.
	 * Should be called after any changes in edge
	 * orientation.
	 *
	 */
	public void updateEdgeUses()
	{
		if(paths == null)
		{
			findPaths();
		}
		
		graph.beginStatsUpdate();
		for(Path p : paths)
		{
			p.updateEdgeUses();
		}
		graph.endStatsUpdate();
	}
	
	/**
	 * Restore all unfixed edges to the unoriented state.  All edges
	 * that aren't conflict edges were fixed if findConflicts was called.
	 * Does nothing if paths have not been found yet.
	 */
	public void clearOrientations()
	{
		if(paths != null)
		{
			for(Path p : paths)
			{
				p.clearOrientations();
			}		
			
			graphStateChanged();
		}
	}

	/**
	 * Calculates the global score for the network orientation,
	 * which is the sum of the weights of all satisfied paths.
	 * If the paths in the graph have not alread been found,
	 * they will be searched for first.
	 * @return the global network score
	 */
	public double globalScore()
	{
		if(paths == null)
		{
			findPaths();
		}

		double global = 0;
		for(Path p: paths)
		{
			global += p.weight();
		}

		return global;
	}

	/**
	 * Calculates the maximum possible score for the network orientation,
	 * which is the score that would be obtained if there were no conflicts
	 * and every path could be optimally oriented.  This score is rarely
	 * feasible.
	 * @return
	 */
	public double maxGlobalScore()
	{
		if(paths == null)
		{
			findPaths();
		}

		double max = 0;
		for(Path p: paths)
		{
			max += p.maxWeight();
		}

		return max;
	}

	/**
	 * Stores the orientations of the conflict edges.  Calling this
	 * method will overwrite the previous saved orientations
	 * stored in the savedOrientations member variable.
	 * @return A reference to the saved orienations.
	 */
	public int[] saveConflictOrientations()
	{
		if(conflictEdges == null)
		{
			throw new IllegalStateException("Cannot save the conflict edges' orientations " +
			"because the conflict edges have not yet been identified");
		}

		// Create a new array to hold the orientations
		savedOrientations = new int[conflictEdges.size()];

		// Save each edges' orientation
		for(int e = 0; e < savedOrientations.length; e++)
		{
			savedOrientations[e] = conflictEdges.get(e).getOrientation();
		}

		return savedOrientations;
	}
	
	/**
	 * Stores the orientations of the conflict edges to a file.
	 */
	public void writeConflictOrientations(String file) throws IOException
	{
		PrintWriter writer = new PrintWriter(new FileWriter(file));
		writeConflictOrientations(writer);
		writer.close();
	}
	
	
	/**
	 * Stores the orientations of the conflict edges to a file.
	 */
	public void writeConflictOrientations(PrintWriter writer)
	{
		if(conflictEdges == null)
		{
			throw new IllegalStateException("Cannot save the conflict edges' orientations " +
			"because the conflict edges have not yet been identified");
		}

		// Create a new buffer to hold the orientations
		StringBuffer buf = new StringBuffer();

		// Save each edges' orientation
		for(int e = 0; e < conflictEdges.size(); e++)
		{
			buf.append(conflictEdges.get(e).getOrientation()).append(" ");
		}
		
		writer.println(buf.toString());
	}


	/**
	 * Load the orientations stored in the savedOrientations member
	 * variable into the conflictEdges.
	 *
	 */
	public void loadConflictOrientations()
	{
		if(savedOrientations == null)
		{
			throw new IllegalStateException("Cannot load the stored orientation " +
			"because no orientation has been saved");
		}

		loadConflictOrientations(savedOrientations);
	}


	/**
	 * Load the orientations stored in the savedOrientations member
	 * variable into the conflictEdges.
	 *
	 */
	public void loadConflictOrientations(int[] orientations)
	{
		if(orientations.length != conflictEdges.size())
		{
			throw new IllegalStateException("There are a different number of " +
			"orientations and conflict edges");
		}

		// Load each edges' orientation
		for(int e = 0; e < orientations.length; e++)
		{
			conflictEdges.get(e).setOrientation(orientations[e]);
		}
		
		graphStateChanged();
	}
	
	
	/**
	 * Load conflict orientations from file.  Assumes the orientations are
	 * in the space-separated format generated by writeConflictOrientations
	 * @param file
	 * @throws IOException
	 */
	public void loadConflictOrientations(String file) throws IOException
	{
		BufferedReader reader = new BufferedReader(new FileReader(file));
		String line = reader.readLine();
		String[] strParts = line.split(" ");
		int[] orient = new int[strParts.length];
		
		for(int o = 0; o < orient.length; o++)
		{
			orient[o] = Integer.parseInt(strParts[o]);
		}
		
		loadConflictOrientations(orient);
	}


	/**
	 * @param orientation1
	 * @param orientation2
	 * @return the number of edge orientation disagreements between the two
	 * orientations
	 */
	public int compareOrientations(int[] orientation1, int[] orientation2)
	{
		if(orientation1.length != orientation2.length)
		{
			throw new IllegalStateException("There are a different number of " +
			"edges in the orientations");
		}

		int dif = 0;
		// Count the differences
		for(int e = 0; e < orientation1.length; e++)
		{
			if(orientation1[e] != orientation2[e])
			{
				dif++;
			}
		}

		return dif;
	}


	/**
	 * Uses the writer to output statistics about the graph.  Currently
	 * only prints the degree of each vertex using undirected edges.
	 * @param writer
	 */
	public void writeStats(PrintWriter writer)
	{
		writer.println("Undirected degree");
		graph.writeUndirDegree(writer);
	}

	/**
	 * Print the edges in each path along with the length and max path weight
	 * @param writer
	 */
	public void writePathWeights(PrintWriter writer)
	{
		for(Path p: paths)
		{
			writer.println(p + "\t" + p.getNumEdges() + "\t" + p.maxWeight());
		}
	}

	
	public void printSortedPaths(String file, String metric, boolean connectedOnly)
		throws IOException
	{
		PrintWriter writer = new PrintWriter(new FileWriter(file));
		
		ArrayList<Path> toSort;
		if(connectedOnly)
		{
			toSort = getSatisfiedPaths();
		}
		else
		{
			toSort = getPaths();
		}
		
		Collections.sort(toSort, Path.getComparator(metric));
		Collections.reverse(toSort); // Order largest to smallest
			
		for(Path p : toSort)
		{
			writer.println(p + "\t" + p.isConnected() + "\t" + p.getMetricValue(metric));
		}
		writer.close();
	}
	
	
	/**
	 * Write all edges that are on a source-target path.  It also
	 * specifies if the edge was originally undirected or directed (protein-protein or
	 * protein-DNA), if it is presently oriented, and its weight.  Finds paths
	 * if they have not already been found.  Only uses satisfied paths.
	 * @param outFile
	 */
	public void writePathEdges(String outFile)
	{
		if(paths == null)
		{
			findPaths();
		}
		
		try
		{
			HashSet<Edge> pathEdges = new HashSet<Edge>();
			
			// Iterate through all satisfied paths and add the edges on the path
			for(Path curPath : getSatisfiedPaths())
			{
				for(Edge curEdge : curPath.getEdges())
				{
					pathEdges.add(curEdge);
				}
			}
			
			PrintWriter writer = new PrintWriter(new FileWriter(outFile));
			writer.println("Source\tType\tTarget\tOriented\tWeight");
			
			for(DirEdge dEdge : graph.getDirEdges())
			{
				if(pathEdges.contains(dEdge))
				{
					String src = dEdge.getSource().getName();
					String targ = dEdge.getTarget().getName();
					
					writer.println(src + "\tpd\t" + targ + "\t" + dEdge.isOriented() + "\t" + dEdge.getWeight());
				}
			}
			
			for(UndirEdge uEdge : graph.getUndirEdges())
			{
				if(pathEdges.contains(uEdge))
				{
					String src = uEdge.getSource().getName();
					String targ = uEdge.getTarget().getName();
					
					writer.println(src + "\tpp\t" + targ + "\t" + uEdge.isOriented() + "\t" + uEdge.getWeight());
				}
			}
			writer.close();
		}
		catch(IOException e)
		{
			e.printStackTrace();
		}
	}
}
