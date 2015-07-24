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
import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.Collection;
import java.util.Collections;
import java.util.Date;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.PriorityQueue;
import java.util.Random;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.zip.GZIPInputStream;
import java.util.zip.GZIPOutputStream;

import util.MapUtil;
import util.StrOps;

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

	/** Used for brute force solution **/
	int configsChecked;

	private static final int MAX_DEPTH = 5;
	/** The default number of times to restart the randomized rounding
	 * portion of the minSatSln algorithm */
	private static final int BERTSIMAS_ROUND_RESTARTS = 100;
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



	//	TODO make a version that runs it multiple times?
	/**
	 * Orients all edges that are not in conflict, then randomly assigns
	 * directions to all edges in conflict.  Changes the state of the
	 * stored Graph.  Does not apply local search after the random orientation.
	 * 
	 * 
	 * @param writer A summary of the random orientations will be written to
	 * this writer.  Not used presently.
	 * @param prefix A prefix that will be prepended to every line of output
	 * (for instance, the number of edges, sources, targets, etc.)  Not used presently.
	 * @return the global score
	 */
	public double randomSln(PrintWriter writer, String prefix) throws IOException
	{
		// TODO decide what to do about output
		// This String will be prepended to every line of output
		//prefix += "\t" + paths.size() + "\t" + usedCount + "\t" + fixCount + "\t" + edges.size();

		// Randomly orient the edges
		// This will automatically find paths and conflict edges if needed
		randomOrient();
		graphStateChanged();

		double global = globalScore();
		System.out.println("Score from random orientation : "+ global);
		System.out.println("Max possible: " + maxGlobalScore() + "\n");

		return global;
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

		// TODO check input
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


	// TODO make a version that runs it multiple times?
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


			// Check if any edges were flipped more than once
//			int numMultiFlip = 0;
//			for(UndirEdge e : conflictEdges)
//			{
//			if(e.getFlipCount() > 1)
//			{
//			numMultiFlip++;
//			System.out.println("Edge " + e + " flipped " +
//			e.getFlipCount() + " times");
//			}
//			}
		} // end check for conflictEdges.size() > 0

		global = globalScore();
		System.out.println("\nEdge flip local search: " + global);
		System.out.println("Max possible: " + maxGlobalScore() + "\n");

		return global;
	}

	
	/**
	 * Uses toulbar2 (http://carlit.toulouse.inra.fr/cgi-bin/awki.cgi/ToolBarIntro),
	 * a weighted CSP solver, to oriented all conflict edges.
	 * 
	 * Presently the interface is very awkward.  The method
	 * first writes the edge orientation as a WCSP, then waits for the user
	 * to notify it that toulbar2 found a solution, and then sets the
	 * edges' orientations.
	 * 
	 * Changes the state of the stored Graph.
	 * 
	 * @param toulbar2Instance the file where the WCSP instance will be written
	 * @param toulbar2Sol where to expect the solution from toulbar2 when it finishes
	 * @return the global score
	 * @throws IOException
	 */
	public double toulbar2Sln(String toulbar2Instance, String toulbar2Sol) throws IOException
	{
		genToulbar2(toulbar2Instance);
		return scoreToulbar2(toulbar2Sol);
	}
	
	
	/**
	 * Generates an input file for the toulbar2 program
	 * @param toulbar2Instance the name of file to be written that will be
	 * used as input for toulbar2
	 * @throws IOException
	 */
	public void genToulbar2(String toulbar2Instance) throws IOException
	{
		// Make sure that the paths and conflict edges
		// have already been identified
		if(paths == null || conflictEdges == null)
		{
			findConflicts();
		}

		long start = System.currentTimeMillis();
		PrintWriter writer = new PrintWriter(new FileWriter(toulbar2Instance));

		// Identify those paths that use conflict edges.  The other paths
		// will be satisfied no matter what and do not need to be incluced
		// in the linear program
		ArrayList<Path> conflictPaths = new ArrayList<Path>(paths.size());
		for(Path p : paths)
		{
			if(p.hasConflicts())
			{
				conflictPaths.add(p);
			}
		}

		int numCp = conflictPaths.size();
		int numCe = conflictEdges.size();
		int numCv = numCp + numCe;

		// Print the XCSP 2.1 format headers
		writer.println("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
		writer.println("<instance>");
		writer.println("<presentation name=\"EdgeOrientation\" " +
				"format=\"XCSP 2.1\" type=\"WCSP\"/>\n");
		
		// The only variables are binary so only one domain is needed
		writer.println("<domains nbDomains=\"1\">");
		writer.println("<domain name=\"D0\" nbValues=\"2\">0..1</domain>");
		writer.println("</domains>\n");
		
		
		// Create a variable for each conflict edge.  Set the 
		// edge ids at the same time so that Paths can
		// easily formulate the linear program constraints
		// with their references to the edges they consist of
		writer.println("<variables nbVariables=\""+ numCe + "\">");
		for(int e = 0; e < conflictEdges.size(); e++)
		{
			conflictEdges.get(e).setId(e);
			writer.println("<variable name=\"E" + e + "\" domain=\"D0\"/>");
		}
		writer.println("</variables>\n");
		

		// Add relation tags for each conflict path
		writer.println("<relations nbRelations=\"" + numCp + "\">");
		for(int p = 0; p < numCp; p++)
		{
			conflictPaths.get(p).writeToulbar2Relation(writer, p);

			if((p + 1) % 10000 == 0)
			{
				System.out.println("Added relation for " + (p + 1) +
						" of " + numCp + " conflict paths");
			}
		}
		System.out.println("Added relations for all conflict paths");
		writer.println("</relations>\n");
		
		// Use (1000 * numCp) + 1 as a proxy to infinity because
		// this cost will never be achieved
		int maxCost = ((1000 * numCp) + 1);
		writer.println("<constraints nbConstraints=\"" + numCp +
			"\" maximalCost=\"" + maxCost + "\">");
		// Add constraint tags for each conflict path
		for(int p = 0; p < numCp; p++)
		{
			conflictPaths.get(p).writeToulbar2Constraint(writer, p);

			if((p + 1) % 10000 == 0)
			{
				System.out.println("Added constraint for " + (p + 1) +
						" of " + numCp + " conflict paths");
			}
		}
		System.out.println("Added constraints for all conflict paths");

		writer.println("</constraints>");
		writer.println("</instance>");
		writer.close();

		long stop = System.currentTimeMillis();
		System.out.println("Time (ms) to create the CSP: " +
				(stop - start) + "\n");
	}
	
	// TODO javadoc
	public double scoreToulbar2PlusSearch(String toulbar2Sol) throws IOException
	{
		scoreToulbar2(toulbar2Sol);
		return localSearchSln();
	}

	// TODO javadoc
	public double scoreToulbar2(String toulbar2Sol) throws IOException
	{
		// Make sure that the paths and conflict edges
		// have already been identified
		if(paths == null || conflictEdges == null)
		{
			findConflicts();
		}
		
		int numCe = conflictEdges.size();		
		
		// Read the solution file
		BufferedReader reader = new BufferedReader(new FileReader(toulbar2Sol));
		
		// The solution should be printed on one line
		String solution = reader.readLine().trim();
		String solParts[] = solution.split(" ");
		reader.close();
		
		// Make sure the solution gives the expected number of edge orientations
		if(solParts.length != numCe)
		{
			throw new IllegalArgumentException("Expected " + numCe + " edge orientation " +
					"assignments but found " + solParts.length + " in the file " +
					toulbar2Sol);
		}
		
		System.out.println("Read " + solParts.length + " edge orientations from " +
				toulbar2Sol);
		
		// Make the edge assignemnts.  1 corresponds to forward.
		for(int e = 0; e < numCe; e++)
		{
			if(Integer.parseInt(solParts[e]) == 0)
			{
				conflictEdges.get(e).setOrientation(Edge.BACKWARD);
			}
			else if(Integer.parseInt(solParts[e]) == 1)
			{
				conflictEdges.get(e).setOrientation(Edge.FORWARD);
			}
			else
			{
				throw new IllegalArgumentException(solParts[e] + " is not a valid " +
						"edge orientation");
			}
		}
		
		graphStateChanged();
		
		double global = globalScore();

		// Report the results
		System.out.println("\nToulbar2 WCSP approximation: " + global);
		System.out.println("Max possible: " + maxGlobalScore() + "\n");

		return global;
	}
	
	/**
	 * Creates an instance of a MATLAB linear program that is a relaxation
	 * of the edge orientation problem.  The constraint matrix can
	 * be generated as a full or sparse MATLAB matrix (sparse is recommended).
	 */
	public void genMatlabLp() throws IOException
	{
		boolean sparse = true;

		// Make sure that the paths and conflict edges
		// have already been identified
		if(paths == null || conflictEdges == null)
		{
			findConflicts();
		}

		long start = System.currentTimeMillis();
		System.out.println("Creating the linear program");

		// Set ids based on the array list index so that Paths can
		// easily formulate the linear program constraints
		// with their references to the edges they consist of
		for(int e = 0; e < conflictEdges.size(); e++)
		{
			conflictEdges.get(e).setId(e);
		}

		// Identify those paths that use conflict edges.  The other paths
		// will be satisfied no matter what and do not need to be incluced
		// in the linear program
		ArrayList<Path> conflictPaths = new ArrayList<Path>(paths.size());
		for(Path p : paths)
		{
			if(p.hasConflicts())
			{
				conflictPaths.add(p);
			}
		}

		int numCp = conflictPaths.size();
		int numCe = conflictEdges.size();
		int numCv = numCp + numCe;

		// Construct the objective function.  Path are weighted with their maximum
		// possible weight and variables are not weighted.
		// Path variables appear in order before edge variables
		StringBuffer buf = new StringBuffer();

		// Iterate through all conflict paths to get their maximum weight
		for(Path p : conflictPaths)
		{
			buf.append(p.maxWeight());
			buf.append("\n");
		}

		// Iterate through all edge variables to add 0 as their
		// coefficient in the objective function
		for(int e = 1; e < numCe; e++)
		{
			buf.append("0\n");
		}
		// Started at index 1 so add the final 0 without a trailing tab
		buf.append("0");

		// Write the objective function, the f matrix in MATLAB
		PrintWriter writer = new PrintWriter(new FileWriter(
		"../testData/simulated/matlab_f_matrix.txt"));
		writer.println(buf.toString());
		writer.close();

		// Create the file for the A MATLAB matrix, the constraint matrix
		PrintWriter writerA;
		if(sparse)
		{
			writerA = new PrintWriter(new FileWriter(
			"../testData/simulated/matlab_A_matrix_sparse.txt"));
		}
		else
		{
			writerA = new PrintWriter(new FileWriter(
			"../testData/simulated/matlab_A_matrix.txt"));
		}

		// Create the file for the B MATLAB matrix, the constraint matrix
		// right hand side
		PrintWriter writerB = new PrintWriter(new FileWriter(
		"../testData/simulated/matlab_B_matrix.txt"));

		// Add contraints for each conflict path
		int constraintRow = 1;
		for(int p = 0; p < numCp; p++)
		{
			constraintRow = conflictPaths.get(p).writeMatlabConstraints(writerA, writerB,
					p, numCp, numCe, constraintRow, sparse);

			if((p + 1) % 1000 == 0)
			{
				System.out.println("Added " + constraintRow + " constraints for "
						+ (p + 1) + " of " + numCp + " conflict paths");
			}
		}
		System.out.println("Added " + (constraintRow - 1) +
		" constraints for all conflict paths");

		writerA.close();
		writerB.close();

		// Create the upper bound and lower bound MATLAB matrices
		PrintWriter writerUb = new PrintWriter(new FileWriter(
		"../testData/simulated/matlab_ub_matrix.txt"));
		PrintWriter writerLb = new PrintWriter(new FileWriter(
		"../testData/simulated/matlab_lb_matrix.txt"));

		StringBuffer bufferUb = new StringBuffer();
		StringBuffer bufferLb = new StringBuffer();

		// Set the upper and lower bounds on all path and edge variables
		// They all must be <= 1 because this is a relaxation of
		// a 0-1 integer program.
		for(int v = 1; v < numCv; v++)
		{
			bufferUb.append("1\n");
			bufferLb.append("0\n");
		}
		bufferUb.append("1");
		bufferLb.append("0");

		writerUb.println(bufferUb.toString());
		writerLb.println(bufferLb.toString());
		writerUb.close();
		writerLb.close();

		long stop = System.currentTimeMillis();
		System.out.println("Time (ms) to write the linear program: " +
				(stop - start) + "\n");
	}

	/**
	 * Exhaustively examine every possible conflict edge orientation
	 * configuration and return the score of the optimal
	 * configuration.  This is only tractable for graphs with 20 or so
	 * conflict edges.
	 * @return
	 */
	public double bruteForceSln()
	{
		// Make sure that the paths and conflict edges
		// have already been identified
		if(paths == null || conflictEdges == null)
		{
			findConflicts();
		}

		if(conflictEdges.size() < 1 || conflictEdges.size() > 25)
		{
			throw new IllegalArgumentException("Brute force algorithm " +
					"is not recommended for " + conflictEdges.size() + " conflicts");
		}

		configsChecked = 0;
		long start = System.currentTimeMillis();

		// Initially orient all edges in the forward direction
		for(UndirEdge e : conflictEdges)
		{
			e.setOrientation(Edge.FORWARD);
		}

		// Now recursively check all possible edge configurations.
		// Store the score and acutal orientation when a new best
		// orientation is found.
		bruteForceHelper(0, Double.NEGATIVE_INFINITY);
		long stop = System.currentTimeMillis();

		// Report the results.  The saved orientation must be loaded first.
		loadConflictOrientations();
		graphStateChanged();
		double global = globalScore();
		System.out.println("Brute force search took (ms): " + (stop - start));
		System.out.println("Brute force configurations checked: " + configsChecked);
		System.out.println("Brute force solution: " + global);
		System.out.println("Max possible: " + maxGlobalScore() + "\n");

		return global;
	}


	private double bruteForceHelper(int edgeInd, double bestScore)
	{
		// If this isn't the last edge in the list, recurse
		if(edgeInd < (conflictEdges.size() - 1))
		{
			bestScore = bruteForceHelper(edgeInd + 1, bestScore);

			// After returning from the recursion, flip the edge and then
			// make the recursive call again
			conflictEdges.get(edgeInd).flip();
			bestScore = bruteForceHelper(edgeInd + 1, bestScore);

			// Return the edge to its original orientation
			conflictEdges.get(edgeInd).flip();

			return bestScore;
		}
		else if(edgeInd == (conflictEdges.size() - 1))
		{
			// The last edge in the list calculates the score for the current
			// orienation and the orientation which it is flipped.
			// If either yields a new high score, the score and
			// configuration are stored
			double curScore = globalScore();
//			for(UndirEdge e : conflictEdges)
//			{
//			System.out.print(e.getOrientation() + " ");
//			}
//			System.out.println();
			if(curScore > bestScore)
			{
				bestScore = curScore;
				saveConflictOrientations();
			}
			configsChecked++;

			conflictEdges.get(edgeInd).flip();

			curScore = globalScore();
//			for(UndirEdge e : conflictEdges)
//			{
//			System.out.print(e.getOrientation() + " ");
//			}
//			System.out.println();
			if(curScore > bestScore)
			{
				bestScore = curScore;
				saveConflictOrientations();
			}
			configsChecked++;

			// Return the edge to its original orientation
			conflictEdges.get(edgeInd).flip();

			return bestScore;
		}
		else
		{
			throw new IllegalArgumentException("There are " + conflictEdges.size() +
					" conflict edges so " + edgeInd + " is not a valid index");
		}
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

		// TODO make these global variables?
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
	 * Use the graph to find paths containing up to depth edges
	 * from sources to targets.  The paths are stored by the
	 * edge orientation algorithm and by default not visible
	 * to the caller.  Sets the class's paths variable.
	 * If more than maxPaths paths exist, only the top maxPaths paths
	 * will be returned as determined by the path weight.
	 * Runs a multithreaded DFS unless all paths are being enumerated.
	 * @param depth
	 * @param maxPaths less than or equal to 0, calls findPaths without a limit on the number
	 * of paths
	 * @retun the number of paths stored
	 */
	public int findTopPaths(int depth, int maxPaths)
	{
		if(maxPaths <= 0)
		{
			return findPaths(depth);
		}
		
		long start = System.currentTimeMillis();
		System.out.println("Finding up to " + maxPaths + " paths up to depth "
				+ depth + " sorting by path weight");

		// Convert from PriorityBlockingQueue to ArrayList
		paths = new ArrayList<Path>(graph.findTopPathsMt(depth, maxPaths));

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
	
	
	// TODO this is not safe if paths were not fully enumerated.
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
	 * Finds paths at depth 2 through 6 and reports how many paths
	 * were found and how long the search took.  Does not affect
	 * the stored paths of the EdgeOrientAlg.
	 *
	 */
	public void findPathsTest()
	{
		long start = System.currentTimeMillis();

		for(int depth = 2; depth <= 6; depth++)
		{
			ArrayList<Path> localPaths = graph.findPaths(depth);
			System.out.println("Found " + localPaths.size() + " paths using depth " + depth);
			System.out.println("Total elapsed time (ms): " + (System.currentTimeMillis() - start));
		}
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
	
	/**
	 * Print the edges in each path along with the length and max satisfied edge uses
	 * @param writer
	 */
	public void writeMaxSatEdgeUses(PrintWriter writer)
	{
		for(Path p: paths)
		{
			writer.println(p + "\t" + p.getNumEdges() + "\t" + p.maxSatEdgeUses());
		}
	}

	
	public void printPaths(String file) throws IOException
	{
		PrintWriter writer = new PrintWriter(new FileWriter(file));
		for(Path p : paths)
		{
			writer.println(p + "\t" + p.isConnected() + "\t" + p.maxWeight());
		}
		writer.close();
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
	
	/**
	 * Write a compressed file containing the StringPath representation of
	 * all satisfied paths
	 * @param outFile
	 */
	public void writeSatisfiedStringPaths(String outFile)
	{
		try
		{
			PrintWriter writer = new PrintWriter(new GZIPOutputStream(new FileOutputStream(outFile)));
			
			for(Path p : getSatisfiedPaths())
			{
				writer.println(new StringPath(p));
			}
			
			writer.close();
		}
		catch(IOException e)
		{
			e.printStackTrace();
		}
	}



	/**
	 * Store all paths to all targets 
	 * @param sourceFile
	 * @param targFile any target weights will not be included in the stored path weights
	 * @param edgeFile
	 * @param priorFile node priors
	 * @param pathLength
	 * @param defaultNodePrior
	 * @param runName use to name the log file and output directory
	 * @param baseDir Parent directory of the directory where the paths
	 * will be written.  A separate file is written for each target.
	 * The log file will be created here as well.
	 */
	public static void storePaths(String sourceFile, String targFile,
			String edgeFile, String priorFile, int pathLength, double defaultNodePrior,
			String runName, String baseDir)
	{
		try
		{
			if(!baseDir.endsWith("/"))
			{
				baseDir = baseDir + "/";
			}
			
			PrintStream stdout = System.out;
			FileOutputStream fos = new FileOutputStream(baseDir + runName + ".out");
			PrintStream ps = new PrintStream(fos);
			System.setOut(ps);
			
			// Print starting time
			long start = System.currentTimeMillis();
			Calendar calendar = Calendar.getInstance();
			Date now = calendar.getTime();
			System.out.println(now);
			
			// Create the directory if it doesn't exist
			String outDir = baseDir + runName + "/";
			File outFile = new File(outDir);
			if(!outFile.exists())
			{
				outFile.mkdir();
			}
			
			// Print configuration
			System.out.println("\nNetwork configuration:");
			System.out.println("Source file: " + sourceFile);
			System.out.println("Target file: " + targFile);
			System.out.println("Edge file: " + edgeFile);
			System.out.println("Node prior file: " + priorFile);
			System.out.println("Path length: " + pathLength);
			System.out.println("Default node prior: " + defaultNodePrior + "\n");
			
			// Set the node priors
			Vertex.setNodePriors(priorFile, defaultNodePrior);
			
			// Create the graph and save the paths
			Graph g = new Graph();
			DataLoader.readEdgesFromEda(g, edgeFile);
			DataLoader.readSources(g, sourceFile);
			DataLoader.readTargets(g, targFile);
			g.storePathsMt(outDir, pathLength);
			
			// Print elapsed time
			calendar = Calendar.getInstance();
			now = calendar.getTime();
			System.out.println("\n" + now);
			long stop = System.currentTimeMillis();
			System.out.println("Total time (s): " + ((stop - start)/1000));
			ps.close();
			// Reset the standard output stream
			System.setOut(stdout);
		}
		catch(IOException e)
		{
			e.printStackTrace();
		}
	}
	
	
	/**
	 * Take the stored paths to all targets written by storePaths and keep only
	 * the top-ranked subset of those paths.  Multithreaded.
	 * @param numPaths the maximum number of high-scoring paths to keep for each target
	 * @param storedPathsDir
	 * @param filteredPathsDir
	 */
	public static void filterStoredPathsMt(int numPaths, String storedPathsDir,
			String filteredPathsDir)
	{
		try
		{
			// Print starting time
			long start = System.currentTimeMillis();
			Calendar calendar = Calendar.getInstance();
			Date now = calendar.getTime();
			System.out.println(now);
			
			System.out.println("\nStored paths: " + storedPathsDir);
			System.out.println("Filtered paths: " + filteredPathsDir);
			System.out.println("Paths to keep: " + numPaths);

			if(!filteredPathsDir.endsWith(File.separator))
			{
				filteredPathsDir = filteredPathsDir + File.separator;
			}
			StrOps.createDirectory(filteredPathsDir);

			// Obtain the set of original target files
			File storedPaths = new File(storedPathsDir);
			File[] origFiles = storedPaths.listFiles(new FilenameFilter()
			{
				public boolean accept(File dir, String name)
				{
					return name.endsWith(".txt.gz");
				}
			});

			// Setup a task for each target
			Collection<Callable<Boolean>> filterRuns = new ArrayList<Callable<Boolean>>(origFiles.length);
			for(File f : origFiles)
			{
				filterRuns.add(new FilterStoredPathsRun(f.getPath(),
						filteredPathsDir + f.getName(),
						numPaths));
			}

			// Setup the ExcecutorService using one thread per processor available.
			// If execution becomes IO-bound there may be little speedup beyond a
			// certain number of processors.
			int numProcs = Runtime.getRuntime().availableProcessors();
			System.out.println("Number of processors available to the Java Virtual Machine: " + numProcs + "\n");
			ExecutorService executor = Executors.newFixedThreadPool(numProcs);

			// Submit the runs and block until they complete
			List<Future<Boolean>> results = executor.invokeAll(filterRuns);
			// Shutdown the threads in the thread pool
			executor.shutdown();

			// Verify that all of the runs completed successfully
			for(Future<Boolean> result : results)
			{
				if(!result.get())
				{
					throw new IllegalStateException("At least one path filtering task completed abnormally");
				}
			}
			
			// Print elapsed time
			calendar = Calendar.getInstance();
			now = calendar.getTime();
			System.out.println("\n" + now);
			long stop = System.currentTimeMillis();
			System.out.println("Total time (s): " + ((stop - start)/1000));
		}
		catch(IOException e)
		{
			e.printStackTrace();
		}
		catch(InterruptedException e)
		{
			System.err.println("Path filtering was interrupted");
			throw new RuntimeException(e);
		}
		catch(ExecutionException e)
		{
			System.err.println("Path filtering exited with an exception");
			throw new RuntimeException(e);
		}
	}
	
	/**
	 * A class that initiates filtering of one target's stored paths, writing only
	 * the top-ranked paths to a new file.
	 */
	public static class FilterStoredPathsRun implements Callable
	{
		String inFile, outFile;
		int numPaths;
		
		/**
		 * @param inFile the file with the full set of paths
		 * @param outFile the file with the filtered paths
		 * @param numPaths maximum number of paths to keep
		 */
		public FilterStoredPathsRun(String inFile, String outFile, int numPaths)
		{
			this.inFile = inFile;
			this.outFile = outFile;
			this.numPaths = numPaths;
		}
		
		public Boolean call() throws Exception
		{
			// Store the top-ranked paths in a heap
			PriorityQueue<StringPath> topPaths = new PriorityQueue<StringPath>(numPaths);
			
			// Read through the full list of paths storing up to numPaths
			BufferedReader reader = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(inFile))));
			String line;
			int count = 0;
			while((line = reader.readLine()) != null)
			{
				StringPath newPath = new StringPath(line);
				count++;
				
				// Add paths regardless of weight until we hit the desired number
				if(topPaths.size() < numPaths)
				{
					topPaths.add(newPath);
				}
				else
				{
					StringPath lowPath = topPaths.peek();
					// If the new path has a higher weight replace the previous
					// lowest weight path in the heap
					if(lowPath.compareTo(newPath) < 0)
					{
						topPaths.poll();
						topPaths.add(newPath);
					}
				}
			}
			reader.close();
			
			// Write the top-ranked paths
			PrintWriter writer = new PrintWriter(new GZIPOutputStream(new FileOutputStream(outFile)));
			for(StringPath p : topPaths)
			{
				writer.println(p);
			}
			writer.close();
			
			System.out.println((new File(inFile)).getName() + " stored paths:\t" + count);
			
			// Return true if the filtering terminated without an exception
			return true;
		}
	}
		
	/**
	 * Load top-ranking paths from the set of paths that have been precomputed
	 * for each target.  The paths are stored by the
	 * edge orientation algorithm and by default not visible
	 * to the caller.  Sets the class's paths variable.
	 * If more than maxPaths paths have been stored, only the top maxPaths paths
	 * will be returned as determined by the path weight.
	 * Adds the target weight to the stored path weight.
	 * @param maxPaths the number of top-ranked stored paths to load, or -1 if all
	 * paths should be loaded
	 * @param storedPathsDir the directory where the stored paths are located
	 * @return the number of paths stored
	 * @throws IOException
	 */
	public int loadStoredTopPaths(int maxPaths, String storedPathsDir)
		throws IOException
	{

		long start = System.currentTimeMillis();

		if(maxPaths < 0)
		{
			System.out.println("Loading all paths from " + storedPathsDir);
		}
		else
		{
			System.out.println("Loading up to " + maxPaths + " paths from " +
					storedPathsDir);
		}

		if(!storedPathsDir.endsWith(File.separator))
		{
			storedPathsDir = storedPathsDir + File.separator;
		}

		// Before loading the paths, verify that all the necessary files are there
		// Track the target Vertex so that the StringPaths can have their weights
		// updated to include the target weight
		HashMap<File, Vertex> targetFiles = new HashMap<File, Vertex>();
		for(Vertex t : graph.getTargets())
		{
			File newFile = new File(storedPathsDir + t + ".txt.gz");
			if(!newFile.exists())
			{
				throw new IllegalStateException("No stored paths for " + t);
			}
			targetFiles.put(newFile, t);
		}

		// Create a heap that finds the top-scoring maxPaths paths
		// Do not create the Path objects until the top-scoring paths have been found
		int stored = 0;
		PriorityQueue<StringPath> topPaths = new PriorityQueue<StringPath>();
		for(File targFile : targetFiles.keySet())
		{
			BufferedReader reader = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(targFile))));

			String line;
			while((line = reader.readLine()) != null)
			{
				// Create the StringPath and then update its weight with the
				// target weight
				StringPath newPath = new StringPath(line);
				Vertex target = targetFiles.get(targFile);
				newPath.addTargetWeight(target.getTargetWeight());
				stored++;

				// Add paths regardless of weight until we hit the desired number
				// or if no bound was placed on the number of paths to load
				if(maxPaths < 0 || topPaths.size() < maxPaths)
				{
					topPaths.add(newPath);
				}
				else
				{
					StringPath lowPath = topPaths.peek();
					// If the new path has a higher weight replace the previous
					// lowest weight path in the heap
					if(lowPath.compareTo(newPath) < 0)
					{
						topPaths.poll();
						topPaths.add(newPath);
					}
				}
			}
			reader.close();
		}

		paths = graph.loadStoredPaths(topPaths);

		long stop = System.currentTimeMillis();
		System.out.println("Loaded " + paths.size() + " out of " +
				stored + " stored paths for the " + targetFiles.size() + " targets");
		System.out.println("Time (ms): " + (stop - start) + "\n");

		return paths.size();

	}
}
