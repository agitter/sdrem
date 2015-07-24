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
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;

import util.ArrayUtil;

/**
 * Methods for evaluating the results of an edge orientation
 * algorithm.
 *
 */
public class OrientEval {

	private EdgeOrientAlg edgeOrient;
	/** The physical interaction network */
	private Graph network;
	/** The merged gold standard pathways */
	private Graph goldStandard;
	
	/** Maps genes to their aliases and description */
	private static HashMap<String, String> geneInfo;
	/** Maps genes to their aliases and description */
	private static HashMap<String, String> geneSyn;
	private static String geneInfoFile = "C:/Users/agitter/Documents/HostVirus/edgeOrientation/testData/sgd/orf_info.txt";
	
	
	/** Maps protein pairs to Strings of experimental evidence for
	 * their interaction */
	private static HashMap<String, String> ppiEvidence;
	private static String ppiEvidFile = "C:/Users/agitter/Documents/HostVirus/pnm/factor/data/PPI/BioGrid_unique_ints_w_evidence.txt";
	
	public OrientEval(Graph network, Graph goldStandard, EdgeOrientAlg edgeOrient)
	{
		this.network = network;
		this.goldStandard = goldStandard;
		this.edgeOrient = edgeOrient;
	}
	
	public static void main(String[] b)
	{
		try {
			System.out.println(geneInfo("YKR035W-A"));
			System.out.println(geneInfo("YOL090W"));
			System.out.println(edgeEvidence("YDR097C", "YOL090W"));
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	
	/**
	 * Validate the paths in the network (before or after orientation) finding
	 * for each path the longest matching subpath in the gold standard pathways.
	 * Then find the average path weights for all longest matching subpaths of
	 * the same length.  The max weight of the path is used, not the active weight.
	 * @param connectedOnly only consider network paths that are connected.  Set this
	 * to false to ignore the effects of the edge orientation.
	 * @return a structure with two arrays.  One array counts the number of times
	 * a path of a given length (num vertices) has a longest subpath in the
	 * gold standard graph of a particular size. Dimensions are: <br>
	 * [# vertices in path][# vertices in longest subpath] <br>
	 * The other array gives the average weight of all paths
	 * that are of a certain length and have a longest subpath in the gold
	 * standard of a certain length. Dimensions are: <br>
	 * [# vertices in path][# vertices in longest subpath]
	 */
	public DoubleArray validatePaths(boolean connectedOnly)
	{
		if(connectedOnly)
		{
			System.out.println("Validating connected network paths");
		}
		else
		{
			System.out.println("Validating all network paths");
		}
		
		
		// Commented out code is the older, slower way that requires
		// enumerating all gold standard subpaths
		
		// First find the longest path in the network so we can
		// create a 2d array that counts the number of times
		// a path of a given length (num vertices) has a longest subpath in the
		// gold standard graph of a particular size
		int maxLength = 0;
		ArrayList<Path> netPaths = edgeOrient.getPaths();
		for(Path netPath : netPaths)
		{
			if(netPath.getNumVertices() > maxLength)
			{
				maxLength = netPath.getNumVertices();
			}
		}
		
		// Dimensions are [# vertices in path][# vertices in longest subpath]
		int[][] maxSubpath = new int[maxLength + 1][maxLength + 1];
		double[][] avgWeight = new double[maxLength + 1][maxLength + 1];
		
		System.out.println(netPaths.size() + " network paths.  " + 
				"Max number of vertices: " + maxLength);
		
		// Could optimize by separating by length and start vertex
		// Now store all of the subpaths up to max length
//		ArrayList<Path> gsSubpaths = new ArrayList<Path>();
//		for(int l = 1; l < maxLength; l++)
//		{
//			gsSubpaths.addAll(goldStandard.findPaths(l, true));
//		}
		
		int pathsChecked = 0, matchingPaths = 0;
		for(Path netPath : netPaths)
		{
			if(!connectedOnly || (connectedOnly && netPath.isConnected()))
			{
				pathsChecked++;
				
				// Iterate through all gold standard paths
				// to see if any of them has the same vertices as
				// this network path.  Stop when a match is found (there
				// could be multiple matches but we don't care)
//				Iterator<Path> gsIter = gsSubpaths.iterator();
//				boolean matches = false;
//				while(!matches && gsIter.hasNext())
//				{
//					matches = netPath.pathMatches(gsIter.next());
//				}
//				
//				if(matches)
//				{
//					matchingPaths++;
//				}
								
				int matches = goldStandard.matchingVertices(netPath);
				maxSubpath[netPath.getNumVertices()][matches]++;

				avgWeight[netPath.getNumVertices()][matches] += netPath.maxWeight();
			}
		}
		
		matchingPaths = maxSubpath[maxLength][maxLength];
		if(connectedOnly)
		{
			System.out.println(matchingPaths + " of " + pathsChecked + " connected network " +
				"paths were found in the gold standard");
		}
		else
		{
			System.out.println(matchingPaths + " of " + pathsChecked + " network " +
			"paths were found in the gold standard");
		}
		
		
		// Calculate the average weights from the count and sum of the weihts
		for(int i = 0; i < maxSubpath.length; i++)
		{
			for(int j = 0; j < maxSubpath[i].length; j++)
			{
				if(maxSubpath[i][j] > 1)
				{
					avgWeight[i][j] /= maxSubpath[i][j];
				}
			}
		}
		
		ArrayUtil.printMatrix(maxSubpath);
		System.out.println("---------------------------");
		ArrayUtil.printMatrix(avgWeight);
		
		DoubleArray result = new DoubleArray(maxSubpath, avgWeight);
		
		return result;
	}
	
	
	/**
	 * Validate the paths in the network (before or after orientation) finding
	 * for each path the longest matching subpath in the gold standard pathways.
	 * Then for all those paths with the desired length and longest subpath,
	 * see how many of them rank in the top X of sorted paths.  The sorted paths only
	 * include paths of the same length as the desired length.  Paths are sorted
	 * by various criteria described below.  Many of the sorting criteria is not
	 * sufficient to unambiguously sort the paths so some results may be misleading.
	 * @param pathLength paths with this number of vertices will be considered
	 * @param longestSubpathLength paths with the maximum longest subpath >= this
	 * value will be considered
	 * @param connectedOnly only consider network paths that are connected.  Edge
	 * orientation will still affect the statistics the paths are ranked by
	 * (e.g. vertex degree)
	 * @return a 2d array where the rows gives the rank threshold and the columns
	 * give the sort criteria.  The rows are:
	 * <br>0: Top X, where X is the number of paths meeting the length and subpath criteria
	 * <br>1: Top 1%
	 * <br>2: Top 5%
	 * <br>3: Top 10%
	 * <br>4: Top 25%
	 * <br>5: Top 50%
	 * <br>6: Top 100
	 * <br>7: Top 1000
	 * <p>
	 * Except for the first column, the columns vary the sort criteria:
	 * <br>0: Number of paths in the rank threshold
	 * <br>1: Path weight (ignoring whether the path is connected)
	 * <br>2: Max edge weight
	 * <br>3: Average edge weight
	 * <br>4: Min edge weight
	 * <br>5: Max edge use
	 * <br>6: Average edge use
	 * <br>7: Min edge use
	 * <br>8: Max edge use (satisfied paths only)
	 * <br>9: Average edge use (satisfied paths only)
	 * <br>10: Min edge use (satisfied paths only)
	 * <br>11: Max vertex degree
	 * <br>12: Average vertex degree
	 * <br>13: Min vertex degree
	 */
	public int[][] validateRankedPaths(int pathLength, int longestSubpathLength, boolean connectedOnly)
	{
		if(connectedOnly)
		{
			System.out.println("Validating connected network paths by rank");
		}
		else
		{
			System.out.println("Validating all network paths by rank");
		}
		System.out.println("Paths of length " + pathLength + " with longest subpaths" +
				" with at least " + longestSubpathLength + " vertices");
		
		// First find the longest path in the network so we can
		// create a 2d array that counts the number of times
		// a path of a given length (num vertices) has a longest subpath in the
		// gold standard graph of a particular size
		int maxLength = 0;
		ArrayList<Path> netPaths = edgeOrient.getPaths();
		for(Path netPath : netPaths)
		{
			if(netPath.getNumVertices() > maxLength)
			{
				maxLength = netPath.getNumVertices();
			}
		}
		
		System.out.println(netPaths.size() + " network paths.  " + 
				"Max number of vertices: " + maxLength);
		
		// Stores those paths of the desired length that meet the longest
		// subpath criteria
		ArrayList<Path> pathsOfInterest = new ArrayList<Path>();
		ArrayList<Path> allPaths = new ArrayList<Path>();
		for(Path netPath : netPaths)
		{
			if(!connectedOnly || (connectedOnly && netPath.isConnected()))
			{
				// Only rank paths that have the desired path length
				int length = netPath.getNumVertices();
				if(length == pathLength)
				{
					allPaths.add(netPath);
					
					int matches = goldStandard.matchingVertices(netPath);
					
					if(matches >= longestSubpathLength)
					{
						pathsOfInterest.add(netPath);
					}
				}
			}
		}
		
		System.out.println(pathsOfInterest.size() + " paths of interest out of " +
				allPaths.size());
		
		// Set the threshold sizes.  Round when finding the top n%
		int numThresholds = 8;
		int[] thresholds = new int[numThresholds];
		thresholds[0] = pathsOfInterest.size();
		thresholds[1] = (int) Math.round(allPaths.size() * 0.01); // Top 1%
		thresholds[2] = (int) Math.round(allPaths.size() * 0.05); // Top 5%
		thresholds[3] = (int) Math.round(allPaths.size() * 0.1); // Top 10%
		thresholds[4] = (int) Math.round(allPaths.size() * 0.25); // Top 25%
		thresholds[5] = (int) Math.round(allPaths.size() * 0.50); // Top 50%
		thresholds[6] = 100;
		thresholds[7] = 1000;
		
		// Set the sorting criteria
		int numCriteria = 13;
		String[] sortCriteria = new String[numCriteria];
		sortCriteria[0] = "MaxWeight";
		sortCriteria[1] = "MaxEdgeWeight";
		sortCriteria[2] = "AvgEdgeWeight";
		sortCriteria[3] = "MinEdgeWeight";
		sortCriteria[4] = "MaxUses";
		sortCriteria[5] = "AvgUses";
		sortCriteria[6] = "MinUses";
		sortCriteria[7] = "MaxSatUses";
		sortCriteria[8] = "AvgSatUses";
		sortCriteria[9] = "MinSatUses";
		sortCriteria[10] = "MaxDegree";
		sortCriteria[11] = "AvgDegree";
		sortCriteria[12] = "MinDegree";
		
		// Create the result array and store the thresholds
		// Make sure the desired threshold doesn't exceed the number of available
		// paths
		int[][] results = new int[numThresholds][numCriteria+1];
		for(int t = 0; t < numThresholds; t++)
		{
			thresholds[t] = Math.min(thresholds[t], allPaths.size());
			results[t][0] = thresholds[t];
		}
		
		// Populate the results for each sorting criteria / threshold combination
		for(int c = 0; c < numCriteria; c++)
		{
			// Sort the paths
			ArrayList<Path> sortedPaths = new ArrayList<Path>(allPaths);
			Collections.sort(sortedPaths, Path.getComparator(sortCriteria[c]));
			Collections.reverse(sortedPaths); // Order largest to smallest
			
			// For each threshold, take the top N sorted paths and save them
			// in a set
			for(int t = 0; t < numThresholds; t++)
			{
				HashSet<Path> topPaths = new HashSet<Path>();
				for(int p = 0; p < thresholds[t]; p++)
				{
					topPaths.add(sortedPaths.get(p));
				}
				
				if(topPaths.size() != thresholds[t])
				{
					throw new IllegalStateException("Not all top ranked paths are unique");
				}
				
				// Now see how many paths of interest are in the top ranked paths
				int hits = 0;
				for(int i = 0; i < pathsOfInterest.size(); i++)
				{
					if(topPaths.contains(pathsOfInterest.get(i)))
					{
						hits++;
					}
				}
				
				// Use c+1 because the first column is the number of paths
				// for each threshold
				results[t][c+1] = hits;
			}
		}
		
		return results;
	}
	
	/**
	 * Write out the information needed to plot a ROC curve, including
	 * TP, FN, FP, TN, TPR, FPR.  Also include the metric value (i.e. the path weight).
	 * @param pathLength the number of vertices a path must have to be included
	 * @param longestSubpathLength the number of consecutive vertices required for a partial match
	 * @param totalPaths the total number of paths with this length (all undirected paths in the graph)
	 * @param totalMatches the total number of those paths that are partial matches
	 * @param rankingCriteria e.g. PathWeight
	 * @param outFile
	 * @throws IOException
	 */
	public void writeROC(int pathLength, int longestSubpathLength, 
		int totalPaths, int totalMatches, String rankingCriteria, String outFile)
		throws IOException
	{
		System.out.println("Calculating ROC for paths of length " + pathLength +
				" vertices with longest subpaths" +
				" that match at least " + longestSubpathLength + " consecutive " +
				"gold standard vertices");
		
		// Only evaluate satisfied paths
		ArrayList<Path> predictedPaths = new ArrayList<Path>();
		HashSet<Path> matchPaths = new HashSet<Path>();
		
		for(Path netPath : edgeOrient.getSatisfiedPaths())
		{
			// Only interested in paths that have the desired number of vertices
			int length = netPath.getNumVertices();
			if(length == pathLength)
			{
				predictedPaths.add(netPath);
				
				// If the desired gold standard matching criteria is met, add
				// the path to the set of matches
				int matches = goldStandard.matchingVertices(netPath);

				if(matches >= longestSubpathLength)
				{
					matchPaths.add(netPath);
				}
			}
		}
		
		System.out.println(predictedPaths.size() + " satisifed network paths with "
				+ pathLength + " vertices");
		System.out.println(matchPaths.size() + " satisifed network paths with "
				+ pathLength + " vertices match gold standard paths");
		System.out.println(totalPaths + " possible paths with "
				+ pathLength + " vertices");
		
		// Sort the paths using the desired criteria
		ArrayList<Path> sortedPaths = new ArrayList<Path>(predictedPaths);
		Collections.sort(sortedPaths, Path.getComparator(rankingCriteria));
		Collections.reverse(sortedPaths); // Order largest to smallest
		
		PrintWriter writer = new PrintWriter(new FileWriter(outFile));
		writer.println("#Pred\tTP\tFN\tFP\tTN\tTPR\tFPR\t" + rankingCriteria);
		
		// Initialize
		int preds = 0;
		int tp = 0; // No predictions are made
		int fn = totalMatches; // All matching paths are not predicted to be matches
		int fp = 0; // No predictions mare made
		int tn = totalPaths - totalMatches; // Correctly predicted these to not be matches
		double tpr = ((double) tp) / (tp + fn);
		double fpr = ((double) fp) / (fp + tn);
		
		writer.println(preds + "\t" + tp + "\t" + fn + "\t" + fp + "\t" + tn + "\t"
				+ tpr + "\t" + fpr + "\t-1");
		
		// Iterate through all possible thresholds and output the TP, FN, FP, TN
		for(Path prediction : sortedPaths)
		{
			preds++;
			
			if(matchPaths.contains(prediction))
			{
				tp++;
				fn--;
			}
			else
			{
				fp++;
				tn--;
			}
			
			tpr = ((double) tp) / (tp + fn);
			fpr = ((double) fp) / (fp + tn);
			
			writer.println(preds + "\t" + tp + "\t" + fn + "\t" + fp + "\t" + tn + "\t"
					+ tpr + "\t" + fpr + "\t" + prediction.getMetricValue(rankingCriteria));
		}
		writer.close();
	}
	
	
	/**
	 * Count the number of paths that have the desired length and the subset
	 * of those that are partial gold standard matches.
	 * @param pathLength the number of vertices a path must have in order to
	 * be considered
	 * @param longestSubpathLength the number of consecutive vertices in a gold
	 * standard path that must match in order to say a path is a partial match
	 * @return an int array {# paths with desired length, # partial matches}
	 */
	public int[] countPathsOfInterest(int pathLength, int longestSubpathLength)
	{
		System.out.println("Counting paths of length " + pathLength +
				" vertices with longest subpaths" +
				" that match at least " + longestSubpathLength + " consecutive " +
				"gold standard vertices");
		
		// Count satisfied and unsatisfied paths
		ArrayList<Path> predictedPaths = new ArrayList<Path>();
		HashSet<Path> matchPaths = new HashSet<Path>();
		
		for(Path netPath : edgeOrient.getPaths())
		{
			// Only interested in paths that have the desired number of vertices
			int length = netPath.getNumVertices();
			if(length == pathLength)
			{
				predictedPaths.add(netPath);
				
				// If the desired gold standard matching criteria is met, add
				// the path to the set of matches
				int matches = goldStandard.matchingVertices(netPath);

				if(matches >= longestSubpathLength)
				{
					matchPaths.add(netPath);
				}
			}
		}
		
		System.out.println(predictedPaths.size() + " network paths with "
				+ pathLength + " vertices");
		System.out.println(matchPaths.size() + " network paths with "
				+ pathLength + " vertices match gold standard paths");
		
		int[] counts = {predictedPaths.size(), matchPaths.size()};
		return counts;
	}
	
	/**
	 * Determine how many gold standard subpaths are present in
	 * the physical interaction network.  This is useful in determining
	 * if the gold standard paths can potentially be recovered if
	 * different sources and targets were used or if the required
	 * edges aren't present.  Respects the orientation of any all
	 * edges in the network.
	 * @param subpathLength the max length of the linear subpaths
	 * (in edges) that will be searched for in the network
	 * @return an 2d array whose ith row gives the number of subpaths
	 * of length i+1 are present in the network and the total number
	 * of subpaths of that length
	 */
	public int[][] gsPathsInNetwork(int subpathLength)
	{
		// [i][0] is the number of subpaths of length i+1 in the network
		// [i][1] is the total number of subpaths of length i+1
		int[][] subpathStats = new int[subpathLength][2];
		
		// Iterate through all subpath lengths up the to max
		for(int len = 0; len < subpathLength; len++)
		{
			ArrayList<Path> gsSubpaths = goldStandard.findPaths(len+1, true);
			subpathStats[len][1] = gsSubpaths.size();
			
			for(Path gsPath : gsSubpaths)
			{
				if(network.containsPath(gsPath))
				{
					subpathStats[len][0]++;
					System.out.println("Matched: " + gsPath);
				}
			}
		}
		
		ArrayUtil.printMatrix(subpathStats);
		
		return subpathStats;
	}
	
	
	
	/**
	 * After orienting the network, calculate statistics for the
	 * conflict edges including the number that are present in the
	 * gold standard, present and oriented correctly, etc.
	 * Also calculate average weights in the interaction network,
	 * the average number of network paths the edges appear in,
	 * and the degree of their vertices.
	 * @return a 2d array of statistics indexed [i][j].  i takes values
	 * 0-3 where
	 * <br> 0 means the edge is not in the gold standard in either direction
	 * <br> 1 means the edge is only in the gold standard in
	 * the direction it was oriented
	 * <br> 2 means the edge is not in the gold standard in the direction
	 * it was oriented but is in the other direction
	 * <br> 3 means the edge is in the gold standard in both directions
	 * <p>
	 * The columns (j values) give different statistics
	 * <br> 0 is the count of conflict edges in this category
	 * <br> 1 is the average weight of the edges in this category
	 * <br> 2 is the average number of paths that use the edge
	 * in the direction it is oriented
	 * <br> 3 is the average number of paths that use the edge
	 * in the direction opposite the one in which it is oriented
	 * <br> 4 is the average degree of the source in the network
	 * <br> 5 is the average degree of the target in the network
	 */
	public double[][] validateConflictEdges()
	{
		double[][] stats = new double[4][6];
		
		ArrayList<UndirEdge> conflicts = edgeOrient.getConflictEdges();
		System.out.println(conflicts.size() + " conflict edges");
		
		for(UndirEdge conflict : conflicts)
		{
			// Create temporary paths in order to probe the gold
			// standard pathways
			ArrayList<Edge> edgeList = new ArrayList<Edge>(1);
			edgeList.add(conflict);
			
			ArrayList<Vertex> vertList = new ArrayList<Vertex>(2);
			vertList.add(conflict.getSource());
			vertList.add(conflict.getTarget());
			
			Path path = new Path(network, edgeList, vertList);
			
			ArrayList<Vertex> revVertList = new ArrayList<Vertex>(2);
			revVertList.add(conflict.getTarget());
			revVertList.add(conflict.getSource());
			
			Path revPath = new Path(network, edgeList, revVertList);
			
			if(!path.isConnected() || revPath.isConnected())
			{
				throw new IllegalStateException("Temporary path for " + conflict +
						" improperly constructed");
			}
			
			// See if the gold standard pathways contain the edge in either direction
			boolean containsPath = goldStandard.containsPath(path);
			boolean containsRevPath = goldStandard.containsPath(revPath);
			
			// Store the index into the array based the results
			int ind = 0;
			// Redundant case, but helps readability
			if(!containsPath && !containsRevPath)
			{
				// The vertices of the conflict edge are not linked in
				// the gold standard pathways
				ind = 0;
			}
			else if(containsPath && !containsRevPath)
			{
				// The conflict edge was oriented "correctly"
				ind = 1;
			}
			else if(!containsPath && containsRevPath)
			{
				// The conflict edge was oriented "incorrectly"
				ind = 2;
			}
			else if(containsPath && containsRevPath)
			{
				// The conflict edge could have been oriented in
				// either direction according to the gold standard
				// pathways.  The gold standard could contain this edge
				// as an undirected edge or there could be different
				// directed edges between these vertices in each direction
				// (perhaps in different individual gold standard paths)
				ind = 3;
			}
			
			// Add this edge's count and weight to the statistics
			stats[ind][0]++;
			stats[ind][1] += conflict.getWeight();
			
			// Add the number of paths that use it in the direction
			// it is oriented and the number that don't
			stats[ind][2] += conflict.numConsistentPaths();
			stats[ind][3] += conflict.numInconsistentPaths();
			
			// Add the degree of the source and target
			stats[ind][4] += network.getDegree(conflict.getSource(), false);
			stats[ind][5] += network.getDegree(conflict.getTarget(), false);
			
			
			// Remove the paths from the edge's association list
			path.deletePath();
			revPath.deletePath();
		}
		
		// Calculate average weight and path use
		for(int i = 0; i < 4; i++)
		{
			if(stats[i][0] > 0)
			{
				stats[i][1] /= stats[i][0];
				stats[i][2] /= stats[i][0];
				stats[i][3] /= stats[i][0];
				stats[i][4] /= stats[i][0];
				stats[i][5] /= stats[i][0];
			}
		}
		
		ArrayUtil.printMatrix(stats);
		
		return stats;
	}
	
	
	
	
	
	
	// TODO should this calculation discard the target's weight when calculating
	// the path weights?
	/**
	 * For each target in the graph, calculate the sum of all
	 * path weights terminating at this target and the percentage
	 * of the global weight that the target's paths contribute.
	 *
	 * @param writer
	 */
	public void writeTargetSupport(PrintWriter writer)
	{
		Set<Vertex> targets = network.getTargets();

		// Create a map from targets to indices for the arrays that will store
		// the path weights
		HashMap<Vertex, Integer> indMap = new HashMap<Vertex, Integer>(targets.size());

		// Store the target's names, the indices, and the initial weights
		double[] actualWeight = new double[targets.size()];
		double[] maxWeight = new double[targets.size()];
		String[] names = new String[targets.size()];
		int[] numPaths = new int[targets.size()];
		int[] connectedPaths = new int[targets.size()];
		int index = 0;
		for(Vertex target : targets)
		{
			actualWeight[index] = 0;
			maxWeight[index] = 0;
			names[index] = target.getName();
			numPaths[index] = 0;
			connectedPaths[index] = 0;

			indMap.put(target, Integer.valueOf(index));

			index++;
		}

		// Iterate through all paths and update the weight arrays and the number
		// of satisfied paths array for that target.  Also store the global
		// weight and max weight
		double global = 0;
		double maxGlobal = 0;
		for(Path p : edgeOrient.getPaths())
		{
			int targIndex = indMap.get(p.getTarget()).intValue();
			double pathMax = p.maxWeight();

			numPaths[targIndex]++;

			// See if the path is connected from source to target
			if(p.isConnected())
			{
				connectedPaths[targIndex]++;

				// TODO this assumes path weight when not connected is 0.  Is that OK?
				// If it's connected the path's actual weight is the max weight
				// If not, the path weight is 0
				actualWeight[targIndex] += pathMax;
				global += pathMax;
			}

			// Add the max weight
			maxWeight[targIndex] += p.maxWeight();
			maxGlobal += pathMax;
		}

		// Print the compiled statistics
		writer.println("Name\tConnected Paths\tTotal Paths\tCP/TP\tCP/Global Paths\t" +
				"Connected Path Weight\tMax Total Path Weight\tCPW/MTPW\t" +
		"CPW/Global Weight\tMTPW/Max Global Weight");
		for(int i = 0; i < targets.size(); i++)
		{
			writer.println(names[i] + "\t" + connectedPaths[i] + "\t" +
					numPaths[i] + "\t" + (((double)connectedPaths[i])/numPaths[i]) + "\t" +
					(((double)connectedPaths[i])/edgeOrient.getPaths().size()) + "\t" + actualWeight[i] + "\t" +
					maxWeight[i] + "\t" + (actualWeight[i]/maxWeight[i]) + "\t" +
					(actualWeight[i]/global) + "\t" + (maxWeight[i]/maxGlobal));
		}
	}
	
	/**
	 * Write the number of paths using each conflict edge,
	 * the direction they use the edge in,
	 * and the sum of their weights.  If the edge has been
	 * oriented, will also print the orientation and the
	 * weight of the paths for the chosen orientation
	 * @param writer
	 */
	public void writeEdgeUse(PrintWriter writer)
	{
		writer.println("Sum of all path weights: " + edgeOrient.maxGlobalScore());
		writer.println("Edge\tOrientation\tOriented path max weight\t" +
				"Oriented path active weight\t" +
				"Fwd paths\tBack paths\tFwd path weight\t" +
				"Back path weight\tMin path max weight\t" +
		"Fwd path active weight\tBack path active weight");
		for(UndirEdge e : edgeOrient.getConflictEdges())
		{
			double fwdMaxWeight = e.fwdPathMaxWeights();
			double backMaxWeight = e.backPathMaxWeights();
			double fwdWeight = e.fwdPathWeights();
			double backWeight = e.backPathWeights();

			String orient = " \t ";
			if(e.getOrientation() == Edge.UNORIENTED)
			{
				orient = "unoriented\t" + (fwdMaxWeight + backMaxWeight) + "\t" +
				(fwdWeight + backWeight);
			}
			else if(e.getOrientation()== Edge.FORWARD)
			{
				orient = "forward\t" + fwdMaxWeight + "\t" + fwdWeight;
			}
			else if(e.getOrientation()== Edge.BACKWARD)
			{
				orient = "backward\t" + backMaxWeight + "\t" + backWeight;
			}

			writer.println(e + "\t" + orient + "\t" + e.numFwdPaths() + "\t" + e.numBackPaths() + "\t" +
					fwdMaxWeight + "\t" + backMaxWeight + "\t" + Math.min(fwdMaxWeight, backMaxWeight) +
					"\t" + fwdWeight + "\t" + backWeight);
		}
	}


	// TODO merge with writeTargetSupport?
	/**
	 * Writes the number of paths that run through each vertex,
	 * their weight, and their active weight after the orientation
	 * @param writer
	 */
	public void writeVertexUse(PrintWriter writer)
	{
		ArrayList<Vertex> vertices = new ArrayList<Vertex>(Vertex.getVertices());

		// Create a map from vertices to indices for the arrays that will store
		// the path weights
		HashMap<Vertex, Integer> indMap = new HashMap<Vertex, Integer>(vertices.size());

		// Store the vertice's names, the indices, and the initial weights
		double[] actualWeight = new double[vertices.size()];
		double[] maxWeight = new double[vertices.size()];
		String[] names = new String[vertices.size()];
		int[] numPaths = new int[vertices.size()];
		int[] connectedPaths = new int[vertices.size()];
		int index = 0;
		for(Vertex v : vertices)
		{
			actualWeight[index] = 0;
			maxWeight[index] = 0;
			names[index] = v.getName();
			numPaths[index] = 0;
			connectedPaths[index] = 0;

			indMap.put(v, Integer.valueOf(index));

			index++;
		}

		// Iterate through all paths and update the weight arrays and the number
		// of satisfied paths array for that target.  Also store the global
		// weight and max weight
		double global = 0;
		double maxGlobal = 0;
		for(Path p : edgeOrient.getPaths())
		{
			// Iterate through all vertices in the path
			for(Vertex v : p.getVertices())
			{
				int targIndex = indMap.get(v).intValue();
				double pathMax = p.maxWeight();

				numPaths[targIndex]++;

				// See if the path is connected from source to target
				if(p.isConnected())
				{
					connectedPaths[targIndex]++;

					// TODO this assumes path weight when not connected is 0.  Is that OK?
					// If it's connected the path's actual weight is the max weight
					// If not, the path weight is 0
					actualWeight[targIndex] += pathMax;
					global += pathMax;
				}

				// Add the max weight
				maxWeight[targIndex] += p.maxWeight();
				maxGlobal += pathMax;
			}
		}

		// Print the compiled statistics
		writer.println("Name\tConnected Paths\tTotal Paths\tCP/TP\tCP/Global Paths\t" +
				"Connected Path Weight\tMax Total Path Weight\tCPW/MTPW\t" +
		"CPW/Global Weight\tMTPW/Max Global Weight");
		for(int i = 0; i < vertices.size(); i++)
		{
			// Only output vertices that are used in paths
			if(numPaths[i] > 0)
			{
				writer.println(names[i] + "\t" + connectedPaths[i] + "\t" +
						numPaths[i] + "\t" + (((double)connectedPaths[i])/numPaths[i]) + "\t" +
						(((double)connectedPaths[i])/edgeOrient.getPaths().size()) + "\t" + actualWeight[i] + "\t" +
						maxWeight[i] + "\t" + (actualWeight[i]/maxWeight[i]) + "\t" +
						(actualWeight[i]/global) + "\t" + (maxWeight[i]/maxGlobal));
			}
		}
	}
	
	
	/**
	 * Writes the path details to a file for the top paths
	 * according to some criteria.
	 * The network must have been oriented already before running.  The criteria
	 * and respective files are:
	 * <br>Top 20 paths (4 or 5 edges) ranked by path weight that are completely in
	 * gold standard: in_gs_top_20_pw.out
	 * <br>Top 20 paths (5 edges) ranked by path weight that have between 2 and 4 consecutive edges
	 * in the gold standard: partial_gs_top_20_pw.out
	 * <br>Top 20 paths (5 edges) ranked by path weight that have no edges in the gold standard:
	 * no_gs_top_20_pw.out
	 * @param outDir the directory where the output files will be written
	 */
	public void topPathDetails(String outDir)
	{
		// Length is in vertices here
		// Should probably be at least 4, 3 starts to not make sense and will
		// not include any partial paths
		int pathLength = 6;
		
		// Stores those paths of the desired length into different lists
		// depending on how many consecutive vertices are in a gold
		// standard path
		ArrayList<Path> allPaths = new ArrayList<Path>();
		ArrayList<Path> inPaths = new ArrayList<Path>();
		ArrayList<Path> partialPaths = new ArrayList<Path>();
		ArrayList<Path> noMatchPaths = new ArrayList<Path>();
		ArrayList<Path> ignorePaths = new ArrayList<Path>();
		for(Path netPath : edgeOrient.getPaths())
		{
			// Only consider paths that are connected after orientation
			if(netPath.isConnected())
			{
				// Includes paths that aren't the correct length
				allPaths.add(netPath);
				
				// Only rank paths that have the desired path length
				int length = netPath.getNumVertices();
				if(length == pathLength)
				{				
					int matches = goldStandard.matchingVertices(netPath);
					
					// Paths that are completely in the gold standard
					if(matches == pathLength)
					{
						inPaths.add(netPath);
					}
					// Paths with at least 2 consecutive edges in the gold standard
					// that are not completely in the gold standard
					else if(matches < pathLength && matches > 2)
					{
						partialPaths.add(netPath);
					}
					// Paths that have no edges in the gold standard
					else if(matches <= 1)
					{
						noMatchPaths.add(netPath);
					}
					// Not interested if a single edge matches
					else if(matches == 2)
					{
						ignorePaths.add(netPath);
					}
				}
				// Also consider slightly shorter paths that are still
				// completely in the gold standard
				else if(length == pathLength - 1)
				{				
					int matches = goldStandard.matchingVertices(netPath);
					
					// Shorter paths that are completely in the gold standard
					if(matches == pathLength - 1)
					{
						inPaths.add(netPath);
					}
				}
			}
		}
		
		System.out.println("All paths: " + allPaths.size());
		System.out.println("In paths: " + inPaths.size());
		System.out.println("Partial paths: " + partialPaths.size());
		System.out.println("No match paths: " + noMatchPaths.size());
		System.out.println("Ignored paths: " + ignorePaths.size());
		
		try
		{
			// Top paths by path weight that are completely in the gold standard
			PrintWriter writer = new PrintWriter(new FileWriter(outDir + "in_gs_top_20_pw.out"));
			
			Collections.sort(inPaths, Path.getComparator("PathWeight"));
			Collections.reverse(inPaths); // Order largest to smallest
			int limit = Math.min(inPaths.size(), 20);
			for(int i = 0; i < limit; i++)
			{
				writePathDetails(inPaths.get(i), writer);
//				System.out.println(inPaths.get(i) + "\t" + inPaths.get(i).weight());
			}
			writer.close();
			
			
			
			// Top paths by path weight that are partially in the gold standard
			writer = new PrintWriter(new FileWriter(outDir + "partial_gs_top_20_pw.out"));
			
			Collections.sort(partialPaths, Path.getComparator("PathWeight"));
			Collections.reverse(partialPaths); // Order largest to smallest
			limit = Math.min(partialPaths.size(), 20);
			for(int i = 0; i < limit; i++)
			{
				writePathDetails(partialPaths.get(i), writer);
			}
			writer.close();
			
			
			
			// Top paths by path weight that have no edges in the gold standard
			writer = new PrintWriter(new FileWriter(outDir + "no_gs_top_20_pw.out"));
			
			Collections.sort(noMatchPaths, Path.getComparator("PathWeight"));
			Collections.reverse(noMatchPaths); // Order largest to smallest
			limit = Math.min(noMatchPaths.size(), 20);
			for(int i = 0; i < limit; i++)
			{
				writePathDetails(noMatchPaths.get(i), writer);
			}
			writer.close();
		}
		catch(IOException e)
		{
			e.printStackTrace();
		}
	}
		
	
	public void testPathDetails()
	{
		
		try {
			PrintWriter writer = new PrintWriter(new FileWriter("C:/Users/agitter/Documents/HostVirus/edgeOrientation/testData/text.out"));
			
			Path p = edgeOrient.getPaths().get(0);
			writePathDetails(p, writer);
			
			p = edgeOrient.getPaths().get(10);
			writePathDetails(p, writer);
			
			p = edgeOrient.getPaths().get(100);
			writePathDetails(p, writer);
			
			p = edgeOrient.getPaths().get(1000);
			writePathDetails(p, writer);
			
			writer.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	
	}
	
	
	/**
	 * Write gene info for every vertex in the path and the
	 * experimental evidence supporting every edge
	 * @param p The path of interest
	 * @param writer
	 */
	public void writePathDetails(Path pa, PrintWriter writer)
	{
		try
		{
			Vertex[] vertices = pa.getVertices();
			for(int p = 0; p < vertices.length - 1; p++)
			{
				writer.print(vertices[p].getName() + "\t");
				writer.print(geneInfo(vertices[p].getName()) + "\t");
				writer.print("<" + synonym(vertices[p].getName()) + "," + synonym(vertices[p+1].getName()) + ">\t");

				if(goldStandard.containsPath(createTmpPath(vertices[p], vertices[p+1])))
				{
					writer.print("Yes\t");
				}
				else
				{
					writer.print("No\t");
				}
				
				// Also check of the edge appears in the gold standard in the other direction
				// This will be useful for checking if an edge is undirected or if
				// the edge is present but the orientation is wrong.
				if(goldStandard.containsPath(createTmpPath(vertices[p+1], vertices[p])))
				{
					writer.print("Yes\t");
				}
				else
				{
					writer.print("No\t");
				}
				
				writer.println(edgeEvidence(vertices[p].getName(), vertices[p+1].getName()));
			}
			writer.print(vertices[vertices.length-1].getName() + "\t");
			writer.println(geneInfo(vertices[vertices.length-1].getName()) + "\n\n");
		}
		catch(IOException e)
		{
			e.printStackTrace();
		}
	}
	
	// TODO comment
	public void evalNoMatch(int pathLength, String directedEdges, String outFile)
		throws IOException
	{
		String rankingMetric = "PathWeight";
		
		// TODO remove the array lists that are used?
		// Stores those paths of the desired length into different lists
		// depending on how many consecutive vertices are in a gold
		// standard path
		ArrayList<Path> inPaths = new ArrayList<Path>();
		ArrayList<Path> partialPaths = new ArrayList<Path>();
		ArrayList<Path> noMatchPaths = new ArrayList<Path>();
		ArrayList<Path> ignorePaths = new ArrayList<Path>();

		// Only use satisfied paths
		for(Path netPath : edgeOrient.getSatisfiedPaths())
		{// Only rank paths that have the desired path length
			int length = netPath.getNumVertices();
			if(length == pathLength)
			{				
				int matches = goldStandard.matchingVertices(netPath);

				// Paths that are completely in the gold standard
				if(matches == pathLength)
				{
					inPaths.add(netPath);
				}
				// Paths with at least 2 consecutive edges in the gold standard
				// that are not completely in the gold standard
				else if(matches < pathLength && matches > 2)
				{
					partialPaths.add(netPath);
				}
				// Paths that have no edges in the gold standard
				else if(matches <= 1)
				{
					noMatchPaths.add(netPath);
				}
				// Not interested if a single edge matches
				else if(matches == 2)
				{
					ignorePaths.add(netPath);
				}
			}
		}

		System.out.println("In paths (all edges match): " + inPaths.size());
		System.out.println("Partial paths (> 1, < all edges match): " + partialPaths.size());
		System.out.println("No match paths (0 edges match): " + noMatchPaths.size());
		System.out.println("Ignored paths (1 edge matches): " + ignorePaths.size());

		System.out.println("\nLoading edges from " + directedEdges);
		
		// Load the known edges that are not in the gold standard
		HashMap<Vertex, HashSet<Vertex>> knownEdges = new HashMap<Vertex, HashSet<Vertex>>();
		BufferedReader reader = new BufferedReader(new FileReader(directedEdges));
		String line;
		while((line = reader.readLine()) != null)
		{
			String[] parts = line.split("\t");
			if(parts.length != 3 || !parts[1].equals("->"))
			{
				throw new IOException("Improperly formatted file");
			}
			
			Vertex source = Vertex.createVertex(parts[0]);
			Vertex target = Vertex.createVertex(parts[2]);
			
			HashSet<Vertex> targSet;
			if(knownEdges.containsKey(source))
			{
				targSet = knownEdges.get(source);
				targSet.add(target);
				knownEdges.put(source, targSet);
			}
			else
			{
				targSet = new HashSet<Vertex>();
				targSet.add(target);
				knownEdges.put(source, targSet);
			}
			
			
			System.out.println("Loaded edge " + synonym(parts[0]) +
					" -> " + synonym(parts[2]));
		}
		reader.close();

		// Top paths by path weight that have no edges in the gold standard
		PrintWriter writer = new PrintWriter(new FileWriter(outFile));
		writer.println("Edge matches\t" + rankingMetric + "\tPath");

		Collections.sort(noMatchPaths, Path.getComparator(rankingMetric));
		Collections.reverse(noMatchPaths); // Order largest to smallest

		// Loop through all paths and count how many of the known edges
		// (in the correct direction) they contain
		for(Path p : noMatchPaths)
		{
			int hits = 0;
			
			Vertex[] vertices = p.getVertices();
			for(int v = 0; v < vertices.length - 1; v++)
			{
				// See if the current vertex and its directed
				// edge to the next vertex are in the collection
				// of known edges
				if(knownEdges.containsKey(vertices[v]) &&
						knownEdges.get(vertices[v]).contains(vertices[v+1]))
				{
					hits++;
				}
			}
			
			writer.print(hits + "\t" + p.getMetricValue(rankingMetric));
			
			for(Vertex v : vertices)
			{
				writer.print("\t" + synonym(v.getName()));
			}
			writer.println();
		}
		
		writer.close();
	}
	
	/**
	 * Creates a temporary Path from a source Vertex to a target Vertex
	 * @param src
	 * @param targ
	 * @return
	 */
	private Path createTmpPath(Vertex src, Vertex targ)
	{
		ArrayList<Edge> edgeList = new ArrayList<Edge>(1);
		edgeList.add(new DirEdge(src, targ));
		
		ArrayList<Vertex> vertList = new ArrayList<Vertex>(2);
		vertList.add(src);
		vertList.add(targ);
		
		return new Path(network, edgeList, vertList);
	}
	
	
	/**
	 * Returns a String that contains tab separated list of the gene's
	 * standard name, aliases, and SGD description
	 * @param gene
	 * @return
	 * @throws IOException
	 */
	private static String geneInfo(String orfName) throws IOException
	{
		if(geneInfo == null)
		{
			loadGeneInfo();
		}
		
		// Lookup the gene info
		if(geneInfo.containsKey(orfName))
		{
			return geneInfo.get(orfName);
		}
		else
		{
			return "-\t-\t-";
		}
	}
	
	
	/**
	 * Returns a String that contains the standard gene name if one exists
	 * or the ORF name if one does not exist
	 * @param gene
	 * @return
	 * @throws IOException
	 */
	private static String synonym(String orfName) throws IOException
	{
		if(geneSyn == null)
		{
			loadGeneInfo();
		}
		
		// Lookup the gene info
		if(geneSyn.containsKey(orfName))
		{
			return geneSyn.get(orfName);
		}
		else
		{
			return orfName;
		}
	}
	
	private static void loadGeneInfo() throws IOException
	{
		geneInfo = new HashMap<String, String>();
		geneSyn = new HashMap<String, String>();
		
		BufferedReader reader = new BufferedReader(new FileReader(geneInfoFile));
		
		String line;
		while((line = reader.readLine()) != null)
		{
			int index = line.indexOf('\t');
			String orf = line.substring(0, index);
			geneInfo.put(orf, line.substring(index + 1));
			
			int index2 = line.indexOf('\t', index + 1);
			String syn = line.substring(index + 1, index2);
			if(syn.equals("-"))
			{
				syn = orf;
			}
			geneSyn.put(orf, syn);
		}
		
		reader.close();
	}
	
	/**
	 * Returns a String that contains a list of types of experiments
	 * that support an undirected PPI between these proteins
	 * @param p1
	 * @param p2
	 * @return
	 * @throws IOException
	 */
	private static String edgeEvidence(String p1, String p2) throws IOException
	{
		if(ppiEvidence == null)
		{
			loadPpiEvidence();
		}
		
		String interaction = sortInt(p1, p2);
		
		// Lookup the evidence String
		if(ppiEvidence.containsKey(interaction))
		{
			return ppiEvidence.get(interaction);
		}
		else
		{
			return "-";
		}
	}
	
	/**
	 * Load evidence of PPI from a file to be used with edgeEvidence
	 * @throws IOException
	 */
	private static void loadPpiEvidence() throws IOException
	{
		ppiEvidence = new HashMap<String, String>();
		
		BufferedReader reader = new BufferedReader(new FileReader(ppiEvidFile));
		
		// Read the header, which gives the interaction types
		String line = reader.readLine();
		line = line.replace("# ", "");
		String[] header = line.split("\t");
		
		while((line = reader.readLine()) != null)
		{
			String[] parts = line.split("\t");
			String interaction = sortInt(parts[0], parts[1]);
			
			// Loop through all interaction types, building up a String of evidence
			StringBuffer evBuf = new StringBuffer();
			for(int p = 2; p < parts.length; p++)
			{
				// Only include evidence types that aren't absent
				if(!parts[p].equals("0"))
				{
					evBuf.append(header[p]).append(" (");
					evBuf.append(parts[p]).append("), ");
				}
			}
			
			// If there is no evidence, write "-"
			if(evBuf.length() == 0)
			{
				evBuf.append('-');
			}
			// Otherwise remove the trailing ", "
			else
			{
				evBuf.delete(evBuf.length() - 2, evBuf.length() - 1);
			}
			
			// Store the evidence in the map
			ppiEvidence.put(interaction, evBuf.toString());
		}
		
		reader.close();
	}
	
	// Copied from pnm project (preprocess/PPI.java)
	/**
	 * Sorts protein 1 and protein 2 so that the lexicographically lesser
	 * one is ordered first.  Returns the pair as a single tab-delimited String.
	 * @param p1
	 * @param p2
	 * @return
	 */
	private static String sortInt(String p1, String p2)
	{
		if(p1.compareToIgnoreCase(p2) <= 0)
		{
			return p1 + "\t" + p2;
		}
		else
		{
			return p2 + "\t" + p1;
		}
	}
	
	class DoubleArray
	{
		public int[][] iArray;
		public double[][] dArray;
		
		public DoubleArray(int[][] i, double[][] d)
		{
			iArray = i;
			dArray = d;
		}
	}
}
