package alg;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collection;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.Set;

import util.ArrayUtil;
import util.Entrez;
import util.MapUtil;
import util.StatUtil;
import util.StrOps;

/**
 * This class provides a variety of methods for generating
 * activity scores for targets.
 *
 */
public class TargetScores {

	public static final String SYN_FILE = "../testData/binding/SGD_standardToOrf.txt";
	/** Map from standard names to ORF names */
	public static HashMap<String, String> synMap = null;
	public static final String REV_SYN_FILE = "../testData/binding/SGD_orfToStandard.txt";
	/** Map from ORF names to standard names */
	public static HashMap<String, String> revSynMap = null;

	/**
	 * Uses the top ranked satisfied paths only by some criteria to obtain a
	 * distribution of metrics (path weight, edge use, etc.) for the real
	 * targets and the specified number of random targets.  Does not
	 * presently use this distribution to generate a score, but prints
	 * the distributions to a file and returns the targets' percentiles
	 * within the distribution.  Use a default path length of 5.
	 * @param runs
	 * @param pathThreshold the number of top ranked paths to be considered
	 * when calculating scores.  Use -1 to indicate the number should be chosen automatically
	 * (theshold will be 5 * (number of real targets + number of random targets))
	 * @param randTargsRatio the ratio of random targets to real targets (i.e. 2 gives
	 * twice as many random targets as real targets).  Must be > 0.
	 * @param rankingMetric one of the ranking metrics than can be used
	 * to create a Path comparator
	 * @param edgeFile
	 * @param srcFile
	 * @param targFile
	 * @param runName a String that is appended after "score" in the filename
	 * to identify runs
	 * @param outDir the output directory, the filename will be generated automatically.
	 * If the empty string "", then no output will be written
	 * @return a map from target names to their percentile in the random distribution
	 */
	public static HashMap<String, Double> rankedPathScore(int runs, int pathThreshold, int randTargsRatio,
				String rankingMetric, String edgeFile, String srcFile, String targFile,
				String runName, String outDir)
	{
		return rankedPathScore(runs, pathThreshold, 1, rankingMetric, edgeFile, srcFile,
				targFile, runName, 5, outDir);
	}
	
	/**
	 * Uses the top ranked satisfied paths only by some criteria to obtain a
	 * distribution of metrics (path weight, edge use, etc.) for the real
	 * targets and an equal number of random targets.  Does not
	 * presently use this distribution to generate a score, but prints
	 * the distributions to a file and returns the targets' percentiles
	 * within the distribution.  Defaults to enumerating all paths.
	 * @param runs
	 * @param pathThreshold
	 * @param rankingMetric one of the ranking metrics than can be used
	 * to create a Path comparator
	 * @param edgeFile
	 * @param srcFile
	 * @param targFile
	 * @param runName a String that is appended after "score" in the filename
	 * to identify runs
	 * @param pathLength the max path length
	 * @param outDir the output directory, the filename will be generated automatically.
	 * If the empty string "", then no output will be written
	 * @return a map from target names to their percentile in the random distribution
	 */
	public static HashMap<String, Double> rankedPathScore(int runs, int pathThreshold,
			String rankingMetric, String edgeFile, String srcFile, String targFile,
			String runName, int pathLength, String outDir)
	{
		return rankedPathScore(runs, pathThreshold, 1, rankingMetric, edgeFile, srcFile,
				targFile, runName, pathLength, outDir);
	}
	
	/**
	 * Uses the top ranked satisfied paths only by some criteria to obtain a
	 * distribution of metrics (path weight, edge use, etc.) for the real
	 * targets and an equal number of random targets.  Does not
	 * presently use this distribution to generate a score, but prints
	 * the distributions to a file and returns the targets' percentiles
	 * within the distribution.  Defaults to enumerating all paths.
	 * @param runs
	 * @param pathThreshold
	 * @param rankingMetric one of the ranking metrics than can be used
	 * to create a Path comparator
	 * @param edgeFile
	 * @param srcFile
	 * @param targFile
	 * @param runName a String that is appended after "score" in the filename
	 * to identify runs
	 * @param pathLength the max path length
	 * @param outDir the output directory, the filename will be generated automatically.
	 * If the empty string "", then no output will be written
	 * @return a map from target names to their percentile in the random distribution
	 */
	public static HashMap<String, Double> rankedPathScore(int runs, int pathThreshold,
			float randTargsRatio,
			String rankingMetric, String edgeFile, String srcFile, String targFile,
			String runName, int pathLength, String outDir)
	{
		return rankedPathScore(runs, pathThreshold, randTargsRatio, -1,
				rankingMetric, edgeFile, srcFile,
				targFile, runName, pathLength, "", outDir);
	}


	/**
	 * Uses the top ranked satisfied paths only by some criteria to obtain a
	 * distribution of metrics (path weight, edge use, etc.) for the real
	 * targets and the specified number of random targets.  Does not
	 * presently use this distribution to generate a score, but prints
	 * the distributions to a file and returns the targets' percentiles
	 * within the distribution
	 * @param runs
	 * @param pathThreshold the number of top ranked paths to be considered
	 * when calculating scores.  Use -1 to indicate the number should be chosen automatically
	 * (theshold will be 5 * (number of real targets + number of random targets))
	 * @param randTargsRatio the ratio of random targets to real targets (i.e. 2 gives
	 * twice as many random targets as real targets).  Must be > 0.
	 * @param maxPaths the maximum number of paths to enumerate.  -1 for all paths.
	 * @param rankingMetric one of the ranking metrics than can be used
	 * to create a Path comparator
	 * @param edgeFile
	 * @param srcFile
	 * @param targFile
	 * @param runName a String that is appended after "score" in the filename
	 * to identify runs
	 * @param pathLength the max path length
	 * @param storedPathsDir the name of a directory that contains precomputed paths
	 * @param outDir the output directory, the filename will be generated automatically.
	 * If the empty string "", then no output will be written
	 * @return a map from target names to their percentile in the random distribution
	 */
	public static HashMap<String, Double> rankedPathScore(int runs, int pathThreshold,
			float randTargsRatio, int maxPaths,
			String rankingMetric, String edgeFile, String srcFile, String targFile,
			String runName, int pathLength, String storedPathsDir, String outDir)
	{
		if(randTargsRatio <= 0)
		{
			throw new IllegalArgumentException("Must have random"
					+ " targets.  " + randTargsRatio + " is not valid.");
		}

		System.out.println("Runs: " + runs);

		long start = System.currentTimeMillis();

		try
		{
			HashMap<Vertex, Integer> realMap = DataLoader.createTargetMap(targFile);
			int numRealTargs = realMap.keySet().size();
			int numRandTargs = Math.round(numRealTargs * randTargsRatio);
			System.out.println("Real targets: " + numRealTargs + "\tRandom targets: " + numRandTargs);

			double[] targScores = new double[numRealTargs];
			double[] realDist = new double[numRealTargs * runs];
			double[] randDist = new double[numRandTargs * runs];
			for(int i = 0; i < numRealTargs; i++)
			{
				targScores[i] = 0;
			}

			// Check if the path threshold should be calculated automatically
			if(pathThreshold < 0)
			{
				pathThreshold = 5 * (numRealTargs + numRandTargs);
				System.out.println("Automatically setting path threshold");
			}
			System.out.println("Number of top paths used in calculations: " + pathThreshold);


			int distInd = 0;
			for(int r = 0; r < runs; r++)
			{
				System.out.println("Run " + r);

				// Reset all state in between runs.
				// It is inconvenient to remove targets and paths and track
				// state correctly
				Graph g = new Graph();

				DataLoader.readEdgesFromEda(g, edgeFile);
				DataLoader.readSources(g, srcFile);
				// Load new random targets every iteration
				// Use the targets with stored paths, if available, as the
				// set of possible random targets
				HashMap<Vertex, Integer> randMap;
				if(storedPathsDir == null || storedPathsDir.equals(""))
				{
					randMap =
						DataLoader.readRealRandTargets(g, targFile, randTargsRatio);
				}
				else
				{
					randMap =
						DataLoader.readRealRandStoredTargets(g, targFile, storedPathsDir, randTargsRatio);
				}

				EdgeOrientAlg orient = new EdgeOrientAlg(g);

				// Find new paths only if a directory with previously stored paths
				// was not provided
				if(storedPathsDir == null || storedPathsDir.equals(""))
				{
					orient.findTopPaths(pathLength, maxPaths);
				}
				else
				{
					orient.loadStoredTopPaths(maxPaths, storedPathsDir);
				}
				orient.findConflicts();				
				orient.randPlusSearchSln(1);

				// Create an array to store the scores for the real targets
				// and another for the random targets
				double[] realScores = new double[numRealTargs];
				for(int i = 0; i < numRealTargs; i++)
				{
					realScores[i] = 0;
				}

				double[] randScores = new double[numRandTargs];
				for(int i = 0; i < numRandTargs; i++)
				{
					randScores[i] = 0;
				}

				// Sort the satisfied paths by the desired metric
				ArrayList<Path> sortedPaths = new ArrayList<Path>(orient.getSatisfiedPaths());
				Collections.sort(sortedPaths, Path.getComparator(rankingMetric));
				Collections.reverse(sortedPaths); // Order largest to smallest

				// For the top thresh paths, add the score to the appropriate array
				int thresh = Math.min(pathThreshold, sortedPaths.size());
				for(int p = 0; p < thresh; p++)
				{
					Path curPath = sortedPaths.get(p);
					Vertex curTarg = curPath.getTarget();

					// See if it's a real or random
					double[] curScores;
					int scoreInd;
					if(realMap.keySet().contains(curTarg))
					{
						curScores = realScores;
						scoreInd = realMap.get(curTarg);
					}
					else if(randMap.keySet().contains(curTarg))
					{
						curScores = randScores;
						scoreInd = randMap.get(curTarg);
					}
					else
					{
						throw new IllegalStateException("Target " + curTarg + " is not real" +
						" or random");
					}

					// TODO use average for some metrics?
					// Add this path's cached value for this metric to the total score
					curScores[scoreInd] += curPath.getMetricValue(rankingMetric);
				} // end iterating through top paths

				// Copy/merge the real and random scores
				for(int i = 0; i < numRealTargs; i++)
				{
					targScores[i] += realScores[i];
					realDist[i + numRealTargs * r] = realScores[i];
					distInd++;
				}

				for(int i = 0; i < numRandTargs; i++)
				{
					randDist[i + numRandTargs * r] = randScores[i];
				}

			} // done with runs


			if((runs * numRealTargs) != distInd)
			{
				throw new IllegalStateException("Should have " + (runs * numRealTargs) +
						" data points in " +
						"the distribution.  Have " + distInd + ".");
			}


			// Take the average for the real targets
			for(int i = 0; i < numRealTargs; i++)
			{
				targScores[i] = targScores[i] / runs;
			}

			// Calculate the percentage of random target scores the average
			// real target scores are greater than
			int[] numGreater  = new int[numRealTargs];
			for(int i = 0; i < numRealTargs; i++)
			{
				numGreater[i] = 0;
			}
			for(int r = 0; r < randDist.length; r++)
			{
				for(int t = 0; t < numRealTargs; t++)
				{
					if(targScores[t] > randDist[r])
					{
						numGreater[t]++;
					}
				}
			}

			double[] percentGreater = new double[numRealTargs];
			for(int t = 0; t < numRealTargs; t++)
			{
				percentGreater[t] = numGreater[t] / ((double) randDist.length);
			}


			// Print the real scores and the distributions if a output
			// directory was specified
			if(!outDir.equals(""))
			{
				String out = outDir + "scores" + runName + "_r" + randTargsRatio + 
				"_" + rankingMetric + "_" 
				+ runs + "_" + pathThreshold + ".txt";
				PrintWriter writer = new PrintWriter(new FileWriter(out));
				writer.println("Target\tAverage score\t# Rand greater than\t% Rand greater than");
				for(Vertex targ : realMap.keySet())
				{
					int ind = realMap.get(targ);
					writer.println(targ.getName() + "\t" + targScores[ind] +
							"\t" + numGreater[ind] + "\t" + percentGreater[ind]);
				}

				writer.println("\n\nReal dist.\tRand dist.");

				int minDistLen = Math.min(realDist.length, randDist.length);
				for(int d = 0; d < minDistLen; d++)
				{
					writer.println(realDist[d] + "\t" + randDist[d]);
				}

				// Now print the rest of the longer distribution if needed
				if(realDist.length > randDist.length)
				{
					for(int d = minDistLen; d < realDist.length; d++)
					{
						writer.println(realDist[d]);
					}
				}
				else if(realDist.length < randDist.length)
				{
					for(int d = minDistLen; d < randDist.length; d++)
					{
						writer.println("\t" + randDist[d]);
					}
				}

				writer.close();
			}

			// Create a map that contains the targets' percentiles in
			// the random distribution
			HashMap<String, Double> percentiles = new HashMap<String, Double>();
			for(Vertex targ : realMap.keySet())
			{
				int ind = realMap.get(targ);
				percentiles.put(targ.getName(), percentGreater[ind]);
			}

			long time = (System.currentTimeMillis() - start) / 1000;
			System.out.println("Time to generate scores (s): " + time);

			return percentiles;
		}
		catch(Exception e)
		{
			e.printStackTrace();
			System.exit(1);
			return null;
		}
	}

	/**
	 * Uses the top ranked satisfied paths only by some criteria to score all nodes in the
	 * graph.  Prints statistics for each node, including the scores and if
	 * they are a source/target.  The statistics are averaged over the specified
	 * number of runs.  Uses a default path length of 5.  Enumerates all paths.
	 * @param runs
	 * @param pathThreshold the number of top ranked paths to be considered
	 * when calculating scores.  Use -1 to indicate the number should be chosen automatically
	 * (theshold will be 5 * number of real targets)
	 * @param rankingMetric one of the ranking metrics than can be used
	 * to create a Path comparator.  Will be used to select the top paths.
	 * @param edgeFile
	 * @param srcFile
	 * @param targFile
	 * @param runName a String that is appended after "nodeScore" in the filename
	 * to identify runs
	 * @param outDir the output directory, the filename will be generated automatically.
	 * If the empty string "", then no output will be written
	 * @return a NodeScoreInfo object that summarizes the output
	 */
	public static NodeScoreInfo rankedPathNodeScore(int runs, int pathThreshold,
			String rankingMetric, String edgeFile, String srcFile, String targFile,
			String runName, String outDir)
	{
		return rankedPathNodeScore(runs, pathThreshold, rankingMetric, edgeFile,
				srcFile, targFile, runName, 5, outDir);
	}
	
	/**
	 * Uses the top ranked satisfied paths only by some criteria to score all nodes in the
	 * graph.  Prints statistics for each node, including the scores and if
	 * they are a source/target.  The statistics are averaged over the specified
	 * number of runs.  Enumerates all paths.
	 * @param runs
	 * @param pathThreshold the number of top ranked paths to be considered
	 * when calculating scores.  Use -1 to indicate the number should be chosen automatically
	 * (theshold will be 5 * number of real targets)
	 * @param rankingMetric one of the ranking metrics than can be used
	 * to create a Path comparator.  Will be used to select the top paths.
	 * @param edgeFile
	 * @param srcFile
	 * @param targFile
	 * @param runName a String that is appended after "nodeScore" in the filename
	 * to identify runs
	 * @param pathLength
	 * @param outDir the output directory, the filename will be generated automatically.
	 * If the empty string "", then no output will be written
	 * @return a NodeScoreInfo object that summarizes the output
	 */
	public static NodeScoreInfo rankedPathNodeScore(int runs, int pathThreshold,
			String rankingMetric, String edgeFile, String srcFile, String targFile,
			String runName, int pathLength, String outDir)
	{
		return rankedPathNodeScore(runs, pathThreshold, rankingMetric, edgeFile,
				srcFile, targFile, runName, pathLength, -1, "", outDir);
	}
	
	/**
	 * Uses the top ranked satisfied paths only by some criteria to score all nodes in the
	 * graph.  Prints statistics for each node, including the scores and if
	 * they are a source/target.  The statistics are averaged over the specified
	 * number of runs.  Assumes that node priors have already been set externally.
	 * @param runs
	 * @param pathThreshold the number of top ranked paths to be considered
	 * when calculating scores.  Use -1 to indicate the number should be chosen automatically
	 * (theshold will be 5 * number of real targets)
	 * @param rankingMetric one of the ranking metrics than can be used
	 * to create a Path comparator.  Will be used to select the top paths.
	 * @param edgeFile
	 * @param srcFile
	 * @param targFile
	 * @param runName a String that is appended after "nodeScore" in the filename
	 * to identify runs
	 * @param pathLength the max path length
	 * @param maxPaths the number of paths to enumerate during the DFS.  Use -1 to enumerate
	 * all paths.
	 * @param storedPathsDir the name of a directory that contains precomputed paths
	 * @param outDir the output directory, the filename will be generated automatically.
	 * If the empty string "", then no output will be written
	 * @return a NodeScoreInfo object that summarizes the output
	 */
	public static NodeScoreInfo rankedPathNodeScore(int runs, int pathThreshold,
			String rankingMetric, String edgeFile, String srcFile, String targFile,
			String runName, int pathLength, int maxPaths, String storedPathsDir,
			String outDir)
	{
		if(revSynMap == null)
		{
			// Don't want the code to crash if the hard-coded file isn't found
			File mapFile = new File(REV_SYN_FILE);
			if(mapFile.exists())
			{
				revSynMap = MapUtil.loadMap(REV_SYN_FILE, 0, 1, false);
			}
			else
			{
				revSynMap = new HashMap<String, String>();
			}
		}

		System.out.println("Runs: " + runs);

		long start = System.currentTimeMillis();

		try
		{
			HashMap<Vertex, Integer> realMap = DataLoader.createTargetMap(targFile);
			int numRealTargs = realMap.keySet().size();

			// Check if the path threshold should be calculated automatically
			if(pathThreshold < 0)
			{
				pathThreshold = 5 * numRealTargs;
				System.out.println("Automatically setting path threshold");
			}
			System.out.println("Number of top paths used in calculations: " + pathThreshold);

			Graph g = new Graph();
			DataLoader.readEdgesFromEda(g, edgeFile);
			DataLoader.readSources(g, srcFile);
			DataLoader.readTargets(g, targFile);

			EdgeOrientAlg orient = new EdgeOrientAlg(g);
			// Find new paths only if a directory with previously stored paths
			// was not provided
			if(storedPathsDir == null || storedPathsDir.equals(""))
			{
				orient.findTopPaths(pathLength, maxPaths);
			}
			else
			{
				orient.loadStoredTopPaths(maxPaths, storedPathsDir);
			}
			orient.findConflicts();		

			// Determine the number of nodes in the network.   This isn't exactly equal
			// to the number of Vertex objects because some of the vertices may have
			// different ids but the same name.
			HashSet<String> vertices = new HashSet<String>();
			for(Vertex v : Vertex.getVertices())
			{
				vertices.add(v.getName());
			}
			int numNodes = vertices.size();

			// Create a map from node names to indices and store all names
			HashMap<String, Integer> nodeMap = new HashMap<String, Integer>();
			ArrayList<String> tmpNames = new ArrayList<String>(vertices);
			String[] names = new String[tmpNames.size()];
			names = tmpNames.toArray(names);
			for(int n = 0; n < numNodes; n++)
			{
				nodeMap.put(names[n], n);
			}

			// Create two sets of arrays, one to track the score in the current run
			// and another to track the final score (the average of all run scores)
			double[][] curScores;
			double[][] finalScores = new double[numNodes][4];

			// Track the names of the source and target vertices
			HashSet<String> sourceNames = new HashSet<String>();
			for(Vertex v : g.getSources())
			{
				sourceNames.add(v.getName());
			}
			HashSet<String> targetNames = new HashSet<String>();
			for(Vertex v : g.getTargets())
			{
				targetNames.add(v.getName());
			}

			// Begin the runs
			for(int r = 0; r < runs; r++)
			{
				// Clear the current scores with every run
				curScores = new double[numNodes][4];

				System.out.println("Run " + r);

				// The set of sources and targets does not change, thus it is
				// not necessary to find paths at every iteration
				orient.randPlusSearchSln(1);

				// Sort the satisfied paths by the desired metric
				ArrayList<Path> sortedPaths = new ArrayList<Path>(orient.getSatisfiedPaths());
				Collections.sort(sortedPaths, Path.getComparator(rankingMetric));
				Collections.reverse(sortedPaths); // Order largest to smallest

				// For the top thresh paths, add the score to the appropriate array
				int thresh = Math.min(pathThreshold, sortedPaths.size());
				double topWeight = 0;
				for(int p = 0; p < thresh; p++)
				{
					Path curPath = sortedPaths.get(p);

					// Get all nodes along the path
					Vertex[] curNodes = curPath.getVertices();

					for(Vertex curNode : curNodes)
					{
						// Add this path's count and weight for the nodes in it
						// Note that this is not dependent on the ranking metric,
						// because we want a weighted estimate of the fraction
						// of top paths through a node.
						int ind = nodeMap.get(curNode.getName());
						curScores[ind][0] += 1;
						curScores[ind][1] += curPath.getMetricValue("PathWeight");

						// The scores involving all paths must also include these
						// top ranked paths
						curScores[ind][2] += 1;
						curScores[ind][3] += curPath.getMetricValue("PathWeight");
					}

					topWeight += curPath.getMetricValue("PathWeight");
				} // end iterating through top paths

				// Now iterate through the remaining paths
				double totalWeight = topWeight;
				for(int p = thresh; p < sortedPaths.size(); p++)
				{
					Path curPath = sortedPaths.get(p);

					// Get all nodes along the path
					Vertex[] curNodes = curPath.getVertices();

					for(Vertex curNode : curNodes)
					{
						// Add this path's count and weight for the nodes in it
						int ind = nodeMap.get(curNode.getName());
						curScores[ind][2] += 1;
						curScores[ind][3] += curPath.getMetricValue("PathWeight");
					}

					totalWeight += curPath.getMetricValue("PathWeight");
				} // end iterating through top paths


				// Sanity check
				if(Math.abs(totalWeight - orient.globalScore()) > (0.001 * totalWeight))
				{
					throw new IllegalStateException("Total path weights do not match: " + 
							totalWeight + "\t" + orient.globalScore());
				}

				// Now calculate the scores for this run and aggregate the value
				// for the final score
				for(int n = 0; n < numNodes; n++)
				{
					// Another sanity check
					if(curScores[n][0] > thresh)
					{
						throw new IllegalStateException("More paths than threshold: " + 
								curScores[n][0] + "\t" + thresh);
					}

					// Divide by the number of top paths
					curScores[n][0] /= thresh;
					// Divide by the weight of the top paths
					curScores[n][1] /= topWeight;
					// Divide by the number of all satisfied paths
					curScores[n][2] /= sortedPaths.size();
					// Divide by the weight of the all satisfied paths
					curScores[n][3] /= totalWeight;

					for(int i = 0; i < 4; i++)
					{
						finalScores[n][i] += curScores[n][i];
					}
				}

			} // done with runs

			// Take the average over the number of runs to obtain the final scores
			finalScores = ArrayUtil.divByConst(finalScores, runs);

			NodeScoreInfo info = null;

			// Print the real scores and the distributions if a output
			// directory was specified
			if(!outDir.equals(""))
			{
				String out = outDir + "nodeScores" + runName + 
				"_" + rankingMetric + "_" 
				+ runs + "_" + pathThreshold + ".txt";
				PrintWriter writer = new PrintWriter(new FileWriter(out));
				writer.println("ORF name\tStandard name\tSource\tTarget\tDegree\t" + 
						"% top " + pathThreshold + " paths through node\t" +
						"Weighted % top " + pathThreshold + " paths through node\t" +
				"% all paths through node\tWeighted % all paths through node");

				HashMap<String, String> infoMap = new HashMap<String, String>();

				for(int n = 0; n < numNodes; n++)
				{
					StringBuffer buf = new StringBuffer(names[n]);
					buf.append("\t");

					// If this node isn't in the synonym map there is most likely a problem
					// (such as the targets file containing blank entries)
					if(!revSynMap.containsKey(names[n]))
					{
						// No longer want to throw an exception, because some of the binding
						// data has atypical targets
//						throw new IllegalStateException("Node " + names[n] + " not present in" +
//						" the synonyms map");
						buf.append(names[n]);
					}
					else
					{
						buf.append(revSynMap.get(names[n]));
					}

					buf.append("\t");

					if(sourceNames.contains(names[n]))
					{
						buf.append("Y").append("\t");
					}
					else
					{
						buf.append("N").append("\t");
					}

					if(targetNames.contains(names[n]))
					{
						buf.append("Y").append("\t");
					}
					else
					{
						buf.append("N").append("\t");
					}

					// The edges in the graph don't change from run to run
					// so it's okay to use the last run's graph.
					buf.append(g.getDegree(names[n]));

					for(int i = 0; i < 4; i++)
					{
						buf.append("\t").append(finalScores[n][i]);
					}

					infoMap.put(names[n], buf.toString());

					writer.println(buf.toString());
				}

				info = new NodeScoreInfo(nodeMap, finalScores, infoMap);

				writer.close();
			}

			long time = (System.currentTimeMillis() - start) / 1000;
			System.out.println("Time to generate scores (s): " + time);

			return info;
		}
		catch(Exception e)
		{
			e.printStackTrace();
			System.exit(1);
			return null;
		}
	}
	
	
	// TODO document
	// Based on rankedPathNodeScore
	/**
	 * @param runs
	 * @param edgeFile
	 * @param srcFile
	 * @param targFile
	 * @param nodePriorsFile
	 * @param runName a String that is appended after "nodeScore" in the filename
	 * to identify runs
	 * @param pathLength the max path length
	 * @param maxPaths the maximum number of paths to store when finding paths
	 * in the network.  Use -1 to allow all paths.
	 * @param outDir the output directory, the filename will be generated automatically.
	 * @return a NodeScoreInfo object that summarizes the output
	 */
	public static void writeNodeScores(int runs,
			String edgeFile, String srcFile, String targFile,
			String nodePriorsFile, String runName, int pathLength, int maxPaths, String outDir)
	{
		System.out.println("Runs: " + runs);		

		try
		{
			Graph g = new Graph();
			DataLoader.readEdgesFromEda(g, edgeFile);
			DataLoader.readSources(g, srcFile);
			DataLoader.readTargets(g, targFile);
			Vertex.setNodePriors(nodePriorsFile, 0.5);

			// Only need to find paths once, then just generate node scores using
			// new random restarts of local search
			EdgeOrientAlg orient = new EdgeOrientAlg(g);
			orient.findTopPaths(pathLength, maxPaths);
			orient.findConflicts();	

//			int[] pathThresholds = {1, 2, 3, 4, 100, -1};
			int[] pathThresholds = {100, 1000, 5000, 10000, 50000, 100000, 500000, 1000000, -1};
			int numThresholds = pathThresholds.length;
			
			// Determine the number of nodes in the network.   This isn't exactly equal
			// to the number of Vertex objects because some of the vertices may have
			// different ids but the same name.
			HashSet<String> vertices = new HashSet<String>();
			for(Vertex v : Vertex.getVertices())
			{
				vertices.add(v.getName());
			}
			int numNodes = vertices.size();

			// Create a map from node names to indices and store all names
			HashMap<String, Integer> nodeMap = new HashMap<String, Integer>();
			ArrayList<String> tmpNames = new ArrayList<String>(vertices);
			String[] names = new String[tmpNames.size()];
			names = tmpNames.toArray(names);
			for(int n = 0; n < numNodes; n++)
			{
				nodeMap.put(names[n], n);
			}

			// Create arrays to track the scores in the current run
			double[][] scores;
			double[][] weightedScores;
			double[] pathWeightSum;

			// Track the names of the source and target vertices
			HashSet<String> sourceNames = new HashSet<String>();
			for(Vertex v : g.getSources())
			{
				sourceNames.add(v.getName());
			}
			HashSet<String> targetNames = new HashSet<String>();
			for(Vertex v : g.getTargets())
			{
				targetNames.add(v.getName());
			}

			// Begin the runs.  Each new run corresponds to a new starting point
			// and local search
			for(int r = 1; r <= runs; r++)
			{
				long start = System.currentTimeMillis();

				// Clear the current scores with every run
				scores = new double[numNodes][numThresholds];
				weightedScores = new double[numNodes][numThresholds];
				pathWeightSum = new double[numThresholds];

				System.out.println("Run " + r);

				orient.randPlusSearchSln(1);

				// Sort the satisfied paths by the desired metric
				ArrayList<Path> sortedPaths = new ArrayList<Path>(orient.getSatisfiedPaths());
				Collections.sort(sortedPaths, Path.getComparator("PathWeight"));
				Collections.reverse(sortedPaths); // Order largest to smallest

				// Make sure we never set a threshold that is greater than
				// the nubmer of satisfied paths
				for(int t = 0; t < numThresholds; t++)
				{
					pathThresholds[t] = Math.min(pathThresholds[t], sortedPaths.size());
				}
				// The last threshold is all paths
				pathThresholds[pathThresholds.length-1] = sortedPaths.size();

				for(int p = 0; p < sortedPaths.size(); p++)
				{
					Path curPath = sortedPaths.get(p);

					// Get all nodes along the path
					Vertex[] curNodes = curPath.getVertices();

					for(Vertex curNode : curNodes)
					{
						// Add this path's count and weight for the nodes in it
						int ind = nodeMap.get(curNode.getName());
						for(int t = 0; t < numThresholds; t++)
						{
							if(p < pathThresholds[t])
							{
								scores[ind][t] += 1;
								weightedScores[ind][t] += curPath.getMetricValue("PathWeight");
							}
						}
					}

					for(int t = 0; t < numThresholds; t++)
					{
						if(p < pathThresholds[t])
						{
							pathWeightSum[t] += curPath.getMetricValue("PathWeight");
						}
					}
				} // end iterating through paths


				// Sanity check
				double totalWeight = pathWeightSum[numThresholds-1];
				if(Math.abs(totalWeight - orient.globalScore()) > (0.001 * totalWeight))
				{
					throw new IllegalStateException("Total path weights do not match: " + 
							totalWeight + "\t" + orient.globalScore());
				}

				// Now calculate the scores for this run and aggregate the value
				// for the final score
				for(int n = 0; n < numNodes; n++)
				{
					for(int t = 0; t < numThresholds; t++)
					{
						// Another sanity check
						if(scores[n][t] > pathThresholds[t])
						{
							throw new IllegalStateException("More paths than threshold: " + 
									scores[n][t] + "\t" + pathThresholds[t]);
						}

						// Divide by the number of top paths
						scores[n][t] /= pathThresholds[t];
						// Divide by the weight of the top paths
						weightedScores[n][t] /= pathWeightSum[t];
					}
				}


				// Print the real scores
				String out = outDir + "nodeScores_" + runName + 
				"_PathWeight_k" + pathLength + "_run" + r + "of"
				+ runs + ".txt";
				PrintWriter writer = new PrintWriter(new FileWriter(out));
				writer.print("Gene\tSource\tTarget\tDegree");

				for(int t = 0; t < numThresholds; t++)
				{
					writer.print("\tTop " + pathThresholds[t]);
				}
				for(int t = 0; t < numThresholds; t++)
				{
					writer.print("\tTop " + pathThresholds[t] + " weighted");
				}
				writer.println();

				for(int n = 0; n < numNodes; n++)
				{
					StringBuffer buf = new StringBuffer(names[n]);
					buf.append("\t");

					if(sourceNames.contains(names[n]))
					{
						buf.append("Y\t");
					}
					else
					{
						buf.append("N\t");
					}

					if(targetNames.contains(names[n]))
					{
						buf.append("Y\t");
					}
					else
					{
						buf.append("N\t");
					}

					// The edges in the graph don't change from run to run
					// so it's okay to use the last run's graph.
					buf.append(g.getDegree(names[n]));

					for(int t = 0; t < numThresholds; t++)
					{
						buf.append("\t").append(scores[n][t]);
					}
					for(int t = 0; t < numThresholds; t++)
					{
						buf.append("\t").append(weightedScores[n][t]);
					}

					writer.println(buf.toString());
				}

				// Lastly write out the total path weight
				StringBuffer buf = new StringBuffer("Total weight");
				buf.append("\t-\t-\t-");
				for(int t = 0; t < numThresholds; t++)
				{
					buf.append("\t-");
				}
				for(int t = 0; t < numThresholds; t++)
				{
					buf.append("\t").append(pathWeightSum[t]);
				}
				writer.println(buf.toString());
				writer.close();

				long time = (System.currentTimeMillis() - start) / 1000;
				System.out.println("Time to generate scores (s): " + time);
			} // done with runs

		}
		catch(Exception e)
		{
			e.printStackTrace();
			System.exit(1);
		}
	}


	/**
	 * Orients the network using local search and calculates the ratio of total
	 * satisfied path weight attributed to a target over all weight divided by
	 * the number of targets.  Does not
	 * presently use this distribution to generate a score, but prints
	 * the distributions to a file.
	 * @param runs
	 * @param edgeFile
	 * @param srcFile
	 * @param targFile
	 * @runName will be used in the automatic filename creation
	 * @param outDir the output directory, the filename will be generated automatically
	 */
	public static void weightRatioScore(int runs,
			String edgeFile, String srcFile, String targFile,
			String runName, String outDir)
	{
		weightRatioScore(runs, 1, edgeFile, srcFile, targFile, runName, outDir);
	}

	// Outdated, does not allow float randTargsRatio, maxPaths, or storedPathsDir
	/**
	 * Orients the network using local search and calculates the ratio of total
	 * satisfied path weight attributed to a target over all weight divided by
	 * the number of targets.  Uses a number of random targets equal to 
	 * (num real targets * randTargsRatio).  Does not
	 * presently use this distribution to generate a score, but prints
	 * the distributions to a file.
	 * @param runs
	 * @param randTargsRatio the ratio of random targets to real targets (i.e. 2 gives
	 * twice as many random targets as real targets).  Must be > 0.
	 * @param edgeFile
	 * @param srcFile
	 * @param targFile
	 * @runName will be used in the automatic filename creation
	 * @param outDir the output directory, the filename will be generated automatically
	 */
	public static void weightRatioScore(int runs, int randTargsRatio,
			String edgeFile, String srcFile, String targFile,
			String runName, String outDir)
	{
		if(randTargsRatio < 1)
		{
			throw new IllegalArgumentException("Must have at least as many random"
					+ " targets as real targets.  " + randTargsRatio + " is not valid.");
		}		

		System.out.println("Runs: " + runs);

		long start = System.currentTimeMillis();

		try
		{
			HashMap<Vertex, Integer> realMap = DataLoader.createTargetMap(targFile);
			int numRealTargs = realMap.keySet().size();
			int numRandTargs = numRealTargs * randTargsRatio;

			double[] targScores = new double[numRealTargs];
			double[] realDist = new double[numRealTargs * runs];
			double[] randDist = new double[numRandTargs * runs];
			for(int i = 0; i < numRealTargs; i++)
			{
				targScores[i] = 0;
			}

			int distInd = 0;
			for(int r = 0; r < runs; r++)
			{
				System.out.println("Run " + r);

				// Reset all state in between runs.
				// It is inconvenient to remove targets and paths and track
				// state correctly
				Graph g = new Graph();

				DataLoader.readEdgesFromEda(g, edgeFile);
				DataLoader.readSources(g, srcFile);
				// Load new random targets every iteration
				HashMap<Vertex, Integer> randMap = DataLoader.readRealRandTargets(g, targFile, randTargsRatio);

				EdgeOrientAlg orient = new EdgeOrientAlg(g);

				orient.findPaths(5);
				orient.findConflicts();				
				orient.randPlusSearchSln(1);

				// Create an array to store the scores for the real targets
				// and another for the random targets
				double[] realScores = new double[numRealTargs];
				for(int i = 0; i < numRealTargs; i++)
				{
					realScores[i] = 0;
				}

				double[] randScores = new double[numRandTargs];
				for(int i = 0; i < numRandTargs; i++)
				{
					randScores[i] = 0;
				}

				// For all paths, add the weight to the appropriate array
				for(Path curPath : orient.getPaths())
				{
					Vertex curTarg = curPath.getTarget();

					// See if it's a real or random
					double[] curScores;
					int scoreInd;
					if(realMap.keySet().contains(curTarg))
					{
						curScores = realScores;
						scoreInd = realMap.get(curTarg);
					}
					else if(randMap.keySet().contains(curTarg))
					{
						curScores = randScores;
						scoreInd = randMap.get(curTarg);
					}
					else
					{
						throw new IllegalStateException("Target " + curTarg + " is not real" +
						" or random");
					}

					// Add this path's weight to the total target weight.  Important to
					// use weight instead of cached weight because not all of
					// the paths will be satisfied
					curScores[scoreInd] += curPath.weight();
				} // end iterating through top paths

				// Store the total weight of all satisfied paths
				double totalWeight = orient.globalScore();

				// Calculate the ratio and divide by the number of targets (real and random)
				for(int i = 0; i < numRealTargs; i++)
				{
					// Compute the true desired metric:
					// weight for this target / (total weight / number of targets)
					realScores[i] = (realScores[i] * (numRealTargs + numRandTargs)) / totalWeight;

					targScores[i] += realScores[i];
					realDist[i + numRealTargs * r] = realScores[i];

					distInd++;
				}

				for(int i = 0; i < numRandTargs; i++)
				{
					// Compute the true desired metric:
					// weight for this target / (total weight / number of targets)
					randScores[i] = (randScores[i] * (numRealTargs + numRandTargs)) / totalWeight;
					randDist[i + numRandTargs * r] = randScores[i];
				}

			} // done with runs


			if((runs * numRealTargs) != distInd)
			{
				throw new IllegalStateException("Should have " + (runs * numRealTargs) +
						" data points in " +
						"the distribution.  Have " + distInd + ".");
			}


			// Take the average for the real targets
			for(int i = 0; i < numRealTargs; i++)
			{
				targScores[i] = targScores[i] / runs;
			}

			// Calculate the percentage of random target scores the average
			// real target scores are greater than
			int[] numGreater  = new int[numRealTargs];
			for(int i = 0; i < numRealTargs; i++)
			{
				numGreater[i] = 0;
			}
			for(int r = 0; r < randDist.length; r++)
			{
				for(int t = 0; t < numRealTargs; t++)
				{
					if(targScores[t] > randDist[r])
					{
						numGreater[t]++;
					}
				}
			}

			double[] percentGreater = new double[numRealTargs];
			for(int t = 0; t < numRealTargs; t++)
			{
				percentGreater[t] = numGreater[t] / ((double) randDist.length);
			}


			// Print the real scores and the distributions
			String out = outDir + "scores" + runName + "_r" + randTargsRatio + 
			"_WeightedRatio_" + runs + ".txt";
			PrintWriter writer = new PrintWriter(new FileWriter(out));
			writer.println("Target\tAverage score\t# Rand greater than\t% Rand greater than");
			for(Vertex targ : realMap.keySet())
			{
				int ind = realMap.get(targ);
				writer.println(targ.getName() + "\t" + targScores[ind] +
						"\t" + numGreater[ind] + "\t" + percentGreater[ind]);
			}

			writer.println("\n\nReal dist.\tRand dist.");

			int minDistLen = Math.min(realDist.length, randDist.length);
			for(int d = 0; d < minDistLen; d++)
			{
				writer.println(realDist[d] + "\t" + randDist[d]);
			}

			// Now print the rest of the longer distribution if needed
			if(realDist.length > randDist.length)
			{
				// This should never be the case since randTargsRatio is a positive int
				for(int d = minDistLen; d < realDist.length; d++)
				{
					writer.println(realDist[d]);
				}
			}
			else if(realDist.length < randDist.length)
			{
				for(int d = minDistLen; d < randDist.length; d++)
				{
					writer.println("\t" + randDist[d]);
				}
			}

			writer.close();

			long time = (System.currentTimeMillis() - start) / 1000;
			System.out.println("Time to generate scores (s): " + time);
		}
		catch(Exception e)
		{
			e.printStackTrace();
			System.exit(1);
		}
	}




//	TODO automatic naming
	/**
	 * Adjusts the prior belief on target activities in the TF-gene binding data
	 * grid according to the support of the targets in the network, as determined
	 * by rankedPathScore.  Does not use node scores
	 * @param percentiles a map of target names in ORF form to percentile scores
	 * (such as those from rankedPathScore)
	 * @param threshold a threshold on the percentile scores.  Anything above the
	 * threshold has its prior increased, anything below decreases
	 * @param oldPrior the binding data file
	 * @param newPrior
	 */
	public static void adjustPrior(HashMap<String, Double> percentiles,
			double threshold, String oldPrior, String newPrior)
	{
		adjustPrior(percentiles, threshold, null, 0, oldPrior, newPrior);
	}


	/**
	 * Adjusts the prior belief on target activities in the TF-gene binding data
	 * grid according to the support of the targets in the network, as determined
	 * by rankedPathScore.  Also incorporates node scores from rankedPathNodeScore
	 * so that if a node has a high score and is not already a source or target
	 * its binding prior can increase.  Uses a default minimum prior of 0, which will not
	 * affect the prior adjustments.
	 * @param targetScores a map of target names in ORF form to percentile scores
	 * (such as those from rankedPathScore)
	 * @param targetThresh a threshold on the target scores.  Anything above the
	 * threshold has its prior increased, anything below decreases.
	 * @param nodeScores a map of node names in ORF form to (weighted) percentages of
	 * paths (all or top only) the node appears in.  This map should not contain
	 * any sources or targets.  Will be ignored if null.
	 * @param nodeThresh a threshold on the node scores.  Any node surpassing this
	 * threshold will have its binding prior increase if it is a TF.
	 * @param oldPrior the binding data file
	 * @param newPrior
	 */
	public static void adjustPrior(HashMap<String, Double> targetScores,
			double targetThresh, HashMap<String, Double> nodeScores, double nodeThresh,
			String oldPrior, String newPrior)
	{
		adjustPrior(targetScores, targetThresh, nodeScores, nodeThresh, oldPrior, newPrior, 0);
	}

	// TODO automatic naming
	/**
	 * Adjusts the prior belief on target activities in the TF-gene binding data
	 * grid according to the support of the targets in the network, as determined
	 * by rankedPathScore.  Also incorporates node scores from rankedPathNodeScore
	 * so that if a node has a high score and is not already a source or target
	 * its binding prior can increase.
	 * @param targetScores a map of target names in ORF form to percentile scores
	 * (such as those from rankedPathScore)
	 * @param targetThresh a threshold on the target scores.  Anything above the
	 * threshold has its prior increased, anything below decreases.
	 * @param nodeScores a map of node names in ORF form to (weighted) percentages of
	 * paths (all or top only) the node appears in.  This map should not contain
	 * any sources or targets.  Will be ignored if null.
	 * @param nodeThresh a threshold on the node scores.  Any node surpassing this
	 * threshold will have its binding prior increase if it is a TF.
	 * @param oldPrior the binding data file
	 * @param newPrior
	 * @param minPrior no prior will be changed to be less than this minimum value
	 */
	public static void adjustPrior(HashMap<String, Double> targetScores,
			double targetThresh, HashMap<String, Double> nodeScores, double nodeThresh,
			String oldPrior, String newPrior, double minPrior)
	{
		try
		{
			if(synMap == null)
			{
				// Don't want the code to crash if the hard-coded file isn't found
				File mapFile = new File(SYN_FILE);
				if(mapFile.exists())
				{
					synMap = MapUtil.loadMap(SYN_FILE, 1, 0, false);
				}
				else
				{
					synMap = new HashMap<String, String>();
				}				
			}

			BufferedReader reader = new BufferedReader(new FileReader(oldPrior));

			// The first row contains the TF names (except the first column)
			String line = reader.readLine();
			String tfs[] = line.split("\t");

			// Translate the TFs to their ORF name if possible
			for(int t = 1; t < tfs.length; t++)
			{
				if(synMap.containsKey(tfs[t]))
				{
					String tf = synMap.get(tfs[t]);
					if(tf.contains("|"))
					{
						// Use the first name in the map because the synonyms file
						// lists the best match first
						tf = tf.substring(0, tf.indexOf('|'));
					}

					tfs[t] = tf;
				}
			}

			// Create arrays to determine how to modify the bound genes' prior
			// for each TF
			int incrTargCount = 0;
			int decrTargCount = 0;
			boolean[] change = new boolean[tfs.length];
			boolean[] increase = new boolean[tfs.length];
			for(int t = 0; t < tfs.length; t++)
			{
				if(targetScores.containsKey(tfs[t]))
				{
					change[t] = true;

					if(targetScores.get(tfs[t]) >= targetThresh)
					{
						increase[t] = true;
						incrTargCount++;
					}
					else
					{
						increase[t] = false;
						decrTargCount++;
					}
				}
				else
				{
					change[t] = false;
					increase[t] = false;
				}

			}
			System.out.println("Increase prior of " + incrTargCount +
			" TFs because of target scores");
			System.out.println("Decrease prior of " + decrTargCount +
			" TFs because of target scores");


			// Adjust the above boolean arrays if there are node scores
			int incrNodeCount = 0;
			if(nodeScores != null)
			{
				for(int t = 0; t < tfs.length; t++)
				{
					if(nodeScores.containsKey(tfs[t]) &&
							nodeScores.get(tfs[t]) >= nodeThresh)
					{
						change[t] = true;
						increase[t] = true;
						incrNodeCount++;
					}
				}
			}
			System.out.println("Increase prior of " + incrNodeCount
					+ " TFs because of node scores");

			// Write out the header again and then process each line individually
			// using the change and increase arrays
			PrintWriter writer = new PrintWriter(new FileWriter(newPrior));
			writer.println(line);

			while((line = reader.readLine()) != null)
			{
				String[] parts = line.split("\t");
				if(parts.length != tfs.length)
				{
					throw new IllegalStateException("Row beginning with " + parts[0] + 
							" does not contain " + tfs.length + " columns");
				}

				StringBuffer buf = new StringBuffer(parts[0]);

				// Start after the gene name, which has already been written
				for(int t = 1; t < tfs.length; t++)
				{
					buf.append("\t");

					if(change[t])
					{
						double val = Double.parseDouble(parts[t]);
						// Don't require the value to be exactly 0.0
						if(Math.abs(val) < 0.000001)
						{
							// Write 0 again if the TF doesn't bind this gene
							buf.append(parts[t]);
						}
						// Nonzero values need to be adjusted
						else
						{
							if(increase[t])
							{
								// The new value is the average of val and 1
								buf.append((1 + val) / 2.0);
							}
							else
							{
								// The new value is the average of val and 0
								// unless it is less than the minimum value
								buf.append(Math.max(val / 2.0, minPrior));
							}
						}
					}
					else
					{
						// Write the original value if this TF wasn't a target or didn't
						// have a high node score
						buf.append(parts[t]);
					}
				} // done processing a line

				writer.println(buf.toString());
			} // done reading file

			writer.close();


		}
		catch(IOException e)
		{
			e.printStackTrace();
		}
	}

	/**
	 * Takes a NodeScoreInfo object from rankedPathNodeScore and returns a map from
	 * nodes that are not sources or targets to their value in the data grid in the
	 * NodeScoreInfo object, which contains the (weighted) % of top/all paths the
	 * node appears in.
	 * @param info
	 * @param column the column in the data grid to use.  See rankedPathNodeScore
	 * @return
	 */
	public static HashMap<String, Double> findNodesOfInterest(NodeScoreInfo info, int column)
	{
		HashMap<String, Double> map = new HashMap<String, Double>();

		// Loop through all nodes
		for(String node : info.infoMap.keySet())
		{
			// Find out if this node is a source or target
			String[] parts = info.infoMap.get(node).split("\t");
			if(parts[2].equals("N") && parts[3].equals("N"))
			{
				// The node is neither a source nor a target			
				// Get the data value (the %) for this node from the data grid
				double val = info.finalScores[info.nodeMap.get(node)][column];

				map.put(node, val);
			}
		}

		return map;
	}


	// Replaced by analyzeTorTargets, but kept in case results from the
	// old TOR runs need to be reproduced
	/**
	 * Takes a file with a list of targets determined to be important by the
	 * edge orientation and writes a file that contains their role in a
	 * rapamycin screen and the Zaman2008 review.  Also writes the targets
	 * that are in the TORC1 pathway but were not found.
	 * @param targets a file containing proteins detected by the edge orientation + 
	 * DREM model.  Can be targets or targets + internal nodes.  Should be
	 * ORF names.
	 * @param screen a screen that identified rapamycin-related genes
	 * @param torRelated a list of genes known to be involved in the TORC1 pathway.
	 * Uses standard names.
	 * @param out
	 */
	public static void analyzeTorTargetsSimple(String targets, String screen, 
			String torRelated, String out)
	{
		try
		{
			if(synMap == null)
			{
				synMap = MapUtil.loadMap(SYN_FILE, 1, 0, false);
			}

			if(revSynMap == null)
			{
				revSynMap = MapUtil.loadMap(REV_SYN_FILE, 0, 1, false);
			}

			HashMap<String, String> screenMap = MapUtil.loadMap(screen, 0, 2, true);

			HashSet<String> targList = MapUtil.loadSet(targets, false);
			HashSet<String> torc1SetStandard = MapUtil.loadSet(torRelated);

			// Translate from standard to ORF
			HashSet<String> torc1Set = new HashSet<String>();
			for(String gene : torc1SetStandard)
			{
				gene = gene.toUpperCase();
				if(!synMap.containsKey(gene))
				{
					throw new IllegalStateException(gene + " cannot be translated");
				}

				torc1Set.add(synMap.get(gene));
			}


			PrintWriter writer = new PrintWriter(new FileWriter(out));
			writer.println("Gene\tORF name\tRapamycin screen\tTORC1 pathway");

			// Loop through all the targets found
			for(String targ : targList)
			{
				StringBuffer buf = new StringBuffer();
				targ = targ.toUpperCase();
				String std = targ;
				if(revSynMap.containsKey(targ))
				{
					std = revSynMap.get(targ);
				}

				buf.append(std).append("\t").append(targ).append("\t");

				// Print the value in the rapamycin screen or N/A
				if(screenMap.containsKey(targ))
				{
					buf.append(screenMap.get(targ));
				}
				else
				{
					buf.append("N/A");
				}
				buf.append("\t");


				// Print if this target is in the TORC1 pathway set
				if(torc1Set.contains(targ))
				{
					buf.append("Y");
					torc1Set.remove(targ);
				}
				else
				{
					buf.append("N");
				}

				writer.println(buf);
			}

			// Now write the TORC1 nodes there were missed by the model
			writer.println("\n\nTORC1 pathway not found");
			for(String torc1 : torc1Set)
			{
				String std = torc1;
				if(revSynMap.containsKey(torc1))
				{
					std = revSynMap.get(torc1);
				}
				writer.println(std);
			}

			writer.close();
		}
		catch(IOException e)
		{
			e.printStackTrace();
		}

	}


	/**
	 * Examines the set of proteins that have been identified as important by the
	 * edge orientation and writes a file that contains their role in a
	 * citric acid screen and signaling databases.  Also writes the targets
	 * that are in the TOR pathway but were not found.  The important targets include
	 * the predefined sources, the DREM-identified targets, and other non-source non-target
	 * nodes that have a large enough fraction of top ranked paths through them.
	 * @param sources the predefined sources.  Should be
	 * ORF names.
	 * @param targets a file containing proteins detected by the edge orientation + 
	 * DREM model.  Should be ORF names.
	 * @param nodeScoresDir the directory containing the node scores file
	 * @param nodescoresPrefix the first part of the node scores file from edge orientation.
	 * One file in the nodeScoresDir should match this prefix.
	 * @param nodeScoreThresh a percent threshold.  Nodes with a % top paths >= this
	 * threshold are deemed important.
	 * @param zamanMembers a list of genes known to be involved in the TOR pathway.
	 * Uses standard names.
	 * @param sgdMembers a list of genes from SGD that are annotated with the GO term
	 * TOR signaling cascade
	 * @param xieScreen genes and their sensitivity to rapamycin.
	 * @param hillenmeyerScreen Hillenmeyer et al. rapamycin screen
	 * @param chanScreen Chan et al. rapamycin screen
	 * @param rapmycinPhen rapamycin phenotype annotation in SGD
	 * @param bindingData the binding data input to DREM, which gives the set of TFs considered
	 * @param out
	 */
	public static void analyzeTorTargets(String sources, String targets, 
			String nodeScoresDir, final String nodeScoresPrefix,
			double nodeScoreThresh, String zamanMembers, String sgdMembers,
			String xieScreen, String hillenmeyerScreen, String chanScreen,
			String rapamycinPhen,
			String bindingData, String out)
	{
		try
		{
			HashSet<String> srcList = MapUtil.loadSet(sources, false);
			// Load the targets but ignore their target weights
			Set<String> targList = MapUtil.loadMap(targets, 0, 1, false).keySet();

			// Find the node scores
			File nsDir = new File(nodeScoresDir);
			File[] nsArray = nsDir.listFiles(new FilenameFilter()
			{
				public boolean accept(File dir, String name)
				{
					return (name.matches(nodeScoresPrefix + ".*"));
				}
			});
			if(nsArray.length != 1)
			{
				throw new IOException("Could not find unique node score file");
			}


			// Load the node score of all proteins, including sources and targets
			// Use the % of top paths for the node score
			HashMap<String, String> nodeScoreMap = MapUtil.loadMap(nsArray[0].getCanonicalPath(), 0, 5, true);

			// Determine which proteins have node scores >= the threshold
			HashSet<String> importantNodes = new HashSet<String>();
			for(String node : nodeScoreMap.keySet())
			{
				double score = Double.parseDouble(nodeScoreMap.get(node));
				if(score >= nodeScoreThresh)
				{
					importantNodes.add(node);
				}
			}

			// Sources, targets, and proteins with high node scores are all
			// considered to be important in the response
			HashSet<String> important = new HashSet<String>();
			important.addAll(srcList);
			important.addAll(targList);
			important.addAll(importantNodes);

			HashSet<String> zamanStandard = MapUtil.loadSet(zamanMembers);
			// Translate from standard to ORF
			HashSet<String> zamanSet = new HashSet<String>();
			for(String gene : zamanStandard)
			{
				zamanSet.add(MapUtil.getOrf(gene.toUpperCase().trim()));
			}
			
			// Already has ORF names
			HashSet<String> sgdSet = MapUtil.loadSet(sgdMembers);

			// Col 0 is the ORF name
			// Col 2 is the rapamycin sensitivity
			HashMap<String, String> xieMap = MapUtil.loadMap(xieScreen, 0, 2, true);

			// Col 0 is the ORF name
			// Col 2 is the lowest p-value across all replicates
			// Some genes appear in more than one line so use a multimap
			HashMap<String, HashSet<String>> hillenmeyerMap = MapUtil.loadMultiMap(hillenmeyerScreen, 0, 1, true);
			
			// Col 0 is the ORF name
			// Col 3 is the rapamycin sensitivity
			HashMap<String, String> chanMap = MapUtil.loadMap(chanScreen, 0, 3, true);

			// Col 0 is the ORF name
			// Names may appear in more than one line
			Set<String> rapPhenSet = MapUtil.loadMap(rapamycinPhen, 0, 4, false).keySet();
			
			// Load the set of TFs
			BufferedReader reader = new BufferedReader(new FileReader(bindingData));
			String[] tfsStandard = reader.readLine().split("\t");
			reader.close();
			HashSet<String> tfSet = new HashSet<String>();
			// Skip the first column, which is the word "Gene"
			for(int t = 1; t < tfsStandard.length; t++)
			{
				tfSet.add(MapUtil.getOrf(tfsStandard[t]));
			}
			
			PrintWriter writer = new PrintWriter(new FileWriter(out));
			writer.println("Gene\tORF name\tRole\tTF\tTOR pathway (Zaman)" +
					"\tTOR pathway (SGD)\tXie screen\tHillenmeyer screen\t" +
					"Chan screen\tRapamycin phenotype (SGD)\tAt least one match");

			// Loop through all the important nodes found
			for(String protein : important)
			{
				boolean foundMatch = false;

				StringBuffer buf = new StringBuffer();
				String std = MapUtil.getStandard(protein);

				buf.append(std).append("\t").append(protein).append("\t");

				// Print if this protein is a source and/or target
				if(srcList.contains(protein) && targList.contains(protein))
				{
					buf.append("Src,Targ");
				}
				else if(srcList.contains(protein))
				{
					buf.append("Src");
				}
				else if(targList.contains(protein))
				{
					buf.append("Targ");
				}
				else
				{
					buf.append("Other");
				}
				buf.append("\t");
				
				// Print if this protein is a TF
				// All targets must be, but some other proteins may be as well
				if(tfSet.contains(protein))
				{
					buf.append("Y");
				}
				else
				{
					buf.append("N");
				}
				buf.append("\t");
				
				// Print if this protein is in the TOR pathway in the Zaman paper
				if(zamanSet.contains(protein))
				{
					foundMatch = true;
					buf.append("Y");
					zamanSet.remove(protein);
				}
				else
				{
					buf.append("N");
				}
				buf.append("\t");

				// Print if this protein is in the TOR pathway in SGD
				if(sgdSet.contains(protein))
				{
					foundMatch = true;
					buf.append("Y");
					sgdSet.remove(protein);
				}
				else
				{
					buf.append("N");
				}
				buf.append("\t");
				
				// Print the level of sensitivity in the Xie screen
				if(xieMap.containsKey(protein))
				{
					String level = xieMap.get(protein);
					buf.append(level);
					foundMatch = true;
				}
				else
				{
					buf.append("N/A");
				}
				buf.append("\t");
				
				// Print the p-value of this protein in the rapamycin screen
				// It will contribute as a "match" if it has a p-value <= 0.0001
				if(hillenmeyerMap.containsKey(protein))
				{
					double bestPval = Double.MAX_VALUE;
					HashSet<String> pvals = hillenmeyerMap.get(protein);
					for(String pvalStr : pvals)
					{
						double nextPval = Double.parseDouble(pvalStr);
						if(nextPval < bestPval)
						{
							bestPval = nextPval;
						}
					}
					buf.append(bestPval);
					if(bestPval <= 0.0001)
					{
						foundMatch = true;
					}
				}
				else
				{
					buf.append("N/A");
				}
				buf.append("\t");
				
				// Print the level of sensitivity in the Chan screen
				if(chanMap.containsKey(protein))
				{
					String level = chanMap.get(protein);
					buf.append(level);
					foundMatch = true;
				}
				else
				{
					buf.append("N/A");
				}
				buf.append("\t");
				
				// Print if the gene is annotated with a rapamycin sensitivity
				// phenotype in SGD
				if(rapPhenSet.contains(protein))
				{
					foundMatch = true;
					buf.append("Y");
				}
				else
				{
					buf.append("N");
				}
				buf.append("\t");

				if(foundMatch)
				{
					buf.append("Y");
				}
				else
				{
					buf.append("N");
				}
				buf.append("\t");

				writer.println(buf);
			}

			// Now write the nodes there were missed by the model
			writer.println("\n\nZaman not found");
			writer.println("Gene\tIs TF?\tXie screen\tHillenmeyer screen\tChan screen\tRapamycin phenotype (SGD)");
			for(String tor : zamanSet)
			{
				writer.print(MapUtil.getStandard(tor) + "\t");
				if(tfSet.contains(tor))
				{
					writer.print("Y\t");
				}
				else
				{
					writer.print("N\t");
				}
				
				// Print the level of sensitivity in the Xie screen
				if(xieMap.containsKey(tor))
				{
					String level = xieMap.get(tor);
					writer.print(level);
				}
				else
				{
					writer.print("N/A");
				}
				writer.print("\t");
				
				// Print the p-value of this protein in the rapamycin screen
				if(hillenmeyerMap.containsKey(tor))
				{
					double bestPval = Double.MAX_VALUE;
					HashSet<String> pvals = hillenmeyerMap.get(tor);
					for(String pvalStr : pvals)
					{
						double nextPval = Double.parseDouble(pvalStr);
						if(nextPval < bestPval)
						{
							bestPval = nextPval;
						}
					}
					
					writer.print(bestPval);
				}
				else
				{
					writer.print("N/A");
				}
				writer.print("\t");
				
				// Print the level of sensitivity in the Chan screen
				if(chanMap.containsKey(tor))
				{
					String level = chanMap.get(tor);
					writer.print(level);
				}
				else
				{
					writer.print("N/A");
				}
				writer.print("\t");
				
				// Print if the gene is annotated with a rapamycin sensitivity
				// phenotype in SGD
				if(rapPhenSet.contains(tor))
				{
					writer.println("Y");
				}
				else
				{
					writer.println("N");
				}
			}

			writer.println("\n\nSGD not found");
			writer.println("Gene\tIs TF?\tXie screen\tHillenmeyer screen\tChan screen\tRapaymcin phenotype (SGD)");
			for(String tor : sgdSet)
			{
				writer.print(MapUtil.getStandard(tor) + "\t");
				if(tfSet.contains(tor))
				{
					writer.print("Y\t");
				}
				else
				{
					writer.print("N\t");
				}
				
				// Print the level of sensitivity in the Xie screen
				if(xieMap.containsKey(tor))
				{
					String level = xieMap.get(tor);
					writer.print(level);
				}
				else
				{
					writer.print("N/A");
				}
				writer.print("\t");
				
				// Print the p-value of this protein in the rapamycin screen
				if(hillenmeyerMap.containsKey(tor))
				{
					double bestPval = Double.MAX_VALUE;
					HashSet<String> pvals = hillenmeyerMap.get(tor);
					for(String pvalStr : pvals)
					{
						double nextPval = Double.parseDouble(pvalStr);
						if(nextPval < bestPval)
						{
							bestPval = nextPval;
						}
					}
					
					writer.print(bestPval);
				}
				else
				{
					writer.print("N/A");
				}
				writer.print("\t");
				
				// Print the level of sensitivity in the Chan screen
				if(chanMap.containsKey(tor))
				{
					String level = chanMap.get(tor);
					writer.print(level);
				}
				else
				{
					writer.print("N/A");
				}
				writer.print("\t");
				
				// Print if the gene is annotated with a rapamycin sensitivity
				// phenotype in SGD
				if(rapPhenSet.contains(tor))
				{
					writer.println("Y");
				}
				else
				{
					writer.println("N");
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
	 * Examines the set of proteins that have been identified as important by the
	 * edge orientation and writes a file that contains their role in a
	 * citric acid screen and signaling databases.  Also writes the targets
	 * that are in the HOG pathway but were not found.  The important targets include
	 * the predefined sources, the DREM-identified targets, and other non-source non-target
	 * nodes that have a large enough fraction of top ranked paths through them.
	 * @param sources the predefined sources.  Should be
	 * ORF names.
	 * @param targets a file containing proteins detected by the edge orientation + 
	 * DREM model.  Should be ORF names.
	 * @param nodeScoresDir the directory containing the node scores file
	 * @param nodescoresPrefix the first part of the node scores file from edge orientation.
	 * One file in the nodeScoresDir should match this prefix.
	 * @param nodeScoreThresh a percent threshold.  Nodes with a % top paths >= this
	 * threshold are deemed important.
	 * @param citricAcidScreen a screen that identified citric acid stress-related genes,
	 * @param keggMembers a list of genes known to be involved in the HOG pathway
	 * based on the KEGG database.  Uses standard names.
	 * @param sciSigMembers a list of genes known to be involved in the HOG pathway
	 * based on the Science Signaling database.  Uses standard names.
	 * @param capaldiHogTfs TFs mentioned by Capaldi 2008 that are activated by Hog1
	 * @param sorbitolScreen genes and their lowest p-value over multiple replicates
	 * of a sorbitol screen
	 * @param krantzMembers members of the HOG pathway as depicted in Figure 1 of Krantz2009
	 * @param sgdOsmotic genes annotated with osmotic stress resistance on SGD
	 * @param out
	 */
	public static void analyzeHogTargets(String sources, String targets, 
			String nodeScoresDir, final String nodeScoresPrefix,
			double nodeScoreThresh, String citricAcidScreen, 
			String keggMembers, String sciSigMembers, String capaldiHogTfs, 
			String sorbitolScreen, String krantzMembers, String sgdOsmotic, String out)
	{
		// Load the sources and predictions
		HashMap<String, Set<String>> predictions = loadSDREM(sources,
				targets, nodeScoresDir, nodeScoresPrefix, nodeScoreThresh);

		// Call this version of analyzeHogTargets to load the gold standards
		// and do the comparison
		analyzeHogTargets(predictions.get("sources"),
				predictions.get("targets"),
				predictions.get("internal"),
				citricAcidScreen, 
				keggMembers,
				sciSigMembers,
				capaldiHogTfs, 
				sorbitolScreen,
				krantzMembers,
				sgdOsmotic,
				out);
	}

	/**
	 * Examines the set of proteins that have been identified as important by 
	 * any computational method and writes a file that contains their role in a
	 * citric acid screen and signaling databases.  Also writes the proteins
	 * that are in the HOG pathway but were not found.
	 * @param sources the set of sources (ORF names)
	 * @param targets the set of targets (ORF names)
	 * @param internalNodes the set of internal nodes included in the predicted
	 * network that are neither sources nor targets (ORF names)
	 * @param citricAcidScreen a screen that identified citric acid stress-related genes,
	 * @param keggMembers a list of genes known to be involved in the HOG pathway
	 * based on the KEGG database.  Uses standard names.
	 * @param sciSigMembers a list of genes known to be involved in the HOG pathway
	 * based on the Science Signaling database.  Uses standard names.
	 * @param capaldiHogTfs TFs mentioned by Capaldi 2008 that are activated by Hog1
	 * @param sorbitolScreen genes and their lowest p-value over multiple replicates
	 * of a sorbitol screen
	 * @param krantzMembers members of the HOG pathway as depicted in Figure 1 of Krantz2009
	 * @param sgdOsmotic genes annotated with osmotic stress resistance on SGD
	 * @param out
	 */
	public static void analyzeHogTargets(Set<String> sources, Set<String> targets,
			Set<String> internalNodes, String citricAcidScreen, 
			String keggMembers, String sciSigMembers, String capaldiHogTfs, 
			String sorbitolScreen, String krantzMembers, String sgdOsmotic, String out)
	{
		try
		{
			// Sources, targets, and proteins with high node scores are all
			// considered to be important in the response
			HashSet<String> important = new HashSet<String>();
			important.addAll(sources);
			important.addAll(targets);
			important.addAll(internalNodes);

			HashSet<String> screenStandard = MapUtil.loadSet(citricAcidScreen);
			// Translate from standard to ORF
			HashSet<String> screenSet = new HashSet<String>();
			for(String gene : screenStandard)
			{
				screenSet.add(MapUtil.getOrf(gene.toUpperCase().trim()));
			}

			HashSet<String> keggStandard = MapUtil.loadSet(keggMembers);
			// Translate from standard to ORF
			HashSet<String> keggSet = new HashSet<String>();
			for(String gene : keggStandard)
			{
				keggSet.add(MapUtil.getOrf(gene.toUpperCase().trim()));
			}

			HashSet<String> sciSigStandard = MapUtil.loadSet(sciSigMembers);
			// Translate from standard to ORF
			HashSet<String> sciSigSet = new HashSet<String>();
			for(String gene : sciSigStandard)
			{
				sciSigSet.add(MapUtil.getOrf(gene.toUpperCase().trim()));
			}

			HashSet<String> capaldiStandard = MapUtil.loadSet(capaldiHogTfs);
			// Translate from standard to ORF
			HashSet<String> capaldiSet = new HashSet<String>();
			for(String gene : capaldiStandard)
			{
				capaldiSet.add(MapUtil.getOrf(gene.toUpperCase().trim()));
			}

			// Even though each gene should only be listed once, this isn't
			// actually the case
			HashMap<String, HashSet<String>> sorbitolMap = MapUtil.loadMultiMap(sorbitolScreen, 0, 1, false);

			HashSet<String> krantzStandard = MapUtil.loadSet(krantzMembers);
			// Translate from standard to ORF
			HashSet<String> krantzSet = new HashSet<String>();
			for(String gene : krantzStandard)
			{
				krantzSet.add(MapUtil.getOrf(gene.toUpperCase().trim()));
			}
			
			HashSet<String> sgdOsmoticSet = MapUtil.loadSet(sgdOsmotic);


			PrintWriter writer = new PrintWriter(new FileWriter(out));
			writer.println("Gene\tORF name\tRole\tCitric acid screen\tHOG pathway (KEGG)" +
					"\tHOG pathway (Sci. Sig.)\tHog1-activated TFs (Capaldi)" +
					"\tSorbitol screen\tHOG pathway (Krantz)\tOsmotic stress (SGD)\t" +
					"At least one match");

			// Loop through all the important nodes found
			for(String protein : important)
			{
				boolean foundMatch = false;

				StringBuffer buf = new StringBuffer();
				String std = MapUtil.getStandard(protein);

				buf.append(std).append("\t").append(protein).append("\t");

				// Print if this protein is a source and/or target
				if(sources.contains(protein) && targets.contains(protein))
				{
					buf.append("Src,Targ");
				}
				else if(sources.contains(protein))
				{
					buf.append("Src");
				}
				else if(targets.contains(protein))
				{
					buf.append("Targ");
				}
				else
				{
					buf.append("Other");
				}
				buf.append("\t");

				// Print if this protein is in the citric acid screen
				if(screenSet.contains(protein))
				{
					foundMatch = true;
					buf.append("Y");
					screenSet.remove(protein);
				}
				else
				{
					buf.append("N");
				}
				buf.append("\t");


				// Print if this protein is in the KEGG HOG pathway
				if(keggSet.contains(protein))
				{
					foundMatch = true;
					buf.append("Y");
					keggSet.remove(protein);
				}
				else
				{
					buf.append("N");
				}
				buf.append("\t");


				// Print if this protein is in the Science Signaling HOG pathway
				if(sciSigSet.contains(protein))
				{
					foundMatch = true;
					buf.append("Y");
					sciSigSet.remove(protein);
				}
				else
				{
					buf.append("N");
				}
				buf.append("\t");


				// Print if this protein is a Hog1-activated TF
				if(capaldiSet.contains(protein))
				{
					foundMatch = true;
					buf.append("Y");
					capaldiSet.remove(protein);
				}
				else
				{
					buf.append("N");
				}
				buf.append("\t");

				// Print the p-value of this protein in the sorbitol screen
				// It will contribute as a "match" if it has a p-value <= 0.0001
				if(sorbitolMap.containsKey(protein))
				{
					double bestPval = Double.MAX_VALUE;
					HashSet<String> pvals = sorbitolMap.get(protein);
					for(String pvalStr : pvals)
					{
						double nextPval = Double.parseDouble(pvalStr);
						if(nextPval < bestPval)
						{
							bestPval = nextPval;
						}
					}
					buf.append(bestPval);
					if(bestPval <= 0.0001)
					{
						foundMatch = true;
					}
				}
				else
				{
					buf.append("N/A");
				}
				buf.append("\t");

				// Print if this protein is a HOG pathway member in Figure 1 of Krantz2009
				if(krantzSet.contains(protein))
				{
					foundMatch = true;
					buf.append("Y");
					krantzSet.remove(protein);
				}
				else
				{
					buf.append("N");
				}
				buf.append("\t");
				
				
				// Print if thsi protein is a gene that is annotated with hyperosmotic
				// stress resistance in SGD
				if(sgdOsmoticSet.contains(protein))
				{
					foundMatch = true;
					buf.append("Y");
					sgdOsmoticSet.remove(protein);
				}
				else
				{
					buf.append("N");
				}
				buf.append("\t");
				

				if(foundMatch)
				{
					buf.append("Y");
				}
				else
				{
					buf.append("N");
				}
				buf.append("\t");

				writer.println(buf);
			}

			// Now write the nodes there were missed by the model
			writer.println("\n\nHOG (KEGG) not found");
			for(String hog : keggSet)
			{
				writer.println(MapUtil.getStandard(hog));
			}

			writer.println("\n\nHOG (Sci. Sig.) not found");
			for(String hog : sciSigSet)
			{
				writer.println(MapUtil.getStandard(hog));
			}

			writer.close();
		}
		catch(IOException e)
		{
			e.printStackTrace();
		}
	}

	
	/**
	 * Examines the set of proteins that have been identified as important by 
	 * ResponseNet and writes a file that contains their role in a
	 * citric acid screen and signaling databases.  Also writes the proteins
	 * that are in the HOG pathway but were not found.
	 * @param edgeFile the table from the ResponseNet website that lists all edges
	 * between proteins and genes.  Used to determine the sources, targets, and internal
	 * nodes
	 * @param citricAcidScreen a screen that identified citric acid stress-related genes,
	 * @param keggMembers a list of genes known to be involved in the HOG pathway
	 * based on the KEGG database.  Uses standard names.
	 * @param sciSigMembers a list of genes known to be involved in the HOG pathway
	 * based on the Science Signaling database.  Uses standard names.
	 * @param capaldiHogTfs TFs mentioned by Capaldi 2008 that are activated by Hog1
	 * @param sorbitolScreen genes and their lowest p-value over multiple replicates
	 * of a sorbitol screen
	 * @param krantzMembers members of the HOG pathway as depicted in Figure 1 of Krantz2009
	 * @param sgdOsmotic genes annotated with osmotic stress resistance on SGD
	 * @param out
	 */
	public static void analyzeHogResponseNet(String edgeFile, String citricAcidScreen, 
			String keggMembers, String sciSigMembers, String capaldiHogTfs, 
			String sorbitolScreen, String krantzMembers, String sgdOsmotic, String out)
	{
		HashMap<String, HashSet<String>> predictions = loadResponseNet(edgeFile);
		analyzeHogTargets(predictions.get("sources"),
				predictions.get("targets"),
				predictions.get("internal"),
				citricAcidScreen, 
				keggMembers,
				sciSigMembers,
				capaldiHogTfs, 
				sorbitolScreen,
				krantzMembers,
				sgdOsmotic,
				out);
	}
	
	/**
	 * Examines the set of proteins that have been identified as important by 
	 * Physical Network Models and writes a file that contains their role in a
	 * citric acid screen and signaling databases.  Also writes the proteins
	 * that are in the HOG pathway but were not found.
	 * @param edgeFile the PNM output file that gives the directions of all edges
	 * that are on an active path.  Used to determine the sources and internal
	 * nodes
	 * @param potentialSources this file lists all proteins that were identified as
	 * knockout nodes and could potentially be sources
	 * @param citricAcidScreen a screen that identified citric acid stress-related genes,
	 * @param keggMembers a list of genes known to be involved in the HOG pathway
	 * based on the KEGG database.  Uses standard names.
	 * @param sciSigMembers a list of genes known to be involved in the HOG pathway
	 * based on the Science Signaling database.  Uses standard names.
	 * @param capaldiHogTfs TFs mentioned by Capaldi 2008 that are activated by Hog1
	 * @param sorbitolScreen genes and their lowest p-value over multiple replicates
	 * of a sorbitol screen
	 * @param krantzMembers members of the HOG pathway as depicted in Figure 1 of Krantz2009
	 * @param sgdOsmotic genes annotated with osmotic stress resistance on SGD
	 * @param out
	 */
	public static void analyzeHogPNM(String edgeFile, String potentialSources, 
			String citricAcidScreen, String keggMembers, String sciSigMembers, String capaldiHogTfs, 
			String sorbitolScreen, String krantzMembers, String sgdOsmotic, String out)
	{
		HashMap<String, HashSet<String>> predictions = loadPNM(edgeFile, potentialSources);
		analyzeHogTargets(predictions.get("sources"),
				predictions.get("targets"),
				predictions.get("internal"),
				citricAcidScreen, 
				keggMembers,
				sciSigMembers,
				capaldiHogTfs, 
				sorbitolScreen,
				krantzMembers,
				sgdOsmotic,
				out);
		
		System.out.println("PNM selected " + predictions.get("sources").size() + " sources");
	}
	
	/**
	 * Loads the set of proteins that have been identified as important by the
	 * edge orientation and writes a file that contains their role in a
	 * citric acid screen and signaling databases.  The important targets include
	 * the predefined sources, the DREM-identified targets, and other non-source non-target
	 * nodes that have a large enough fraction of top ranked paths through them.
	 * @param sources the predefined sources.  Should be ORF names.
	 * @param targets a file containing proteins detected by the edge orientation + 
	 * DREM model.  Should be ORF names.
	 * @param nodeScoresDir the directory containing the node scores file
	 * @param nodescoresPrefix the first part of the node scores file from edge orientation.
	 * One file in the nodeScoresDir should match this prefix.
	 * @param nodeScoreThresh a percent threshold.  Nodes with a % top paths >= this
	 * threshold are deemed important.
	 */
	public static HashMap<String, Set<String>> loadSDREM(String sources, String targets, 
			String nodeScoresDir, final String nodeScoresPrefix,
			double nodeScoreThresh)
	{
		try
		{
			HashMap<String, Set<String>> nodesMap = new HashMap<String, Set<String>>();
			
			HashSet<String> srcList = MapUtil.loadSet(sources, false);
			// Load the targets but ignore their target weights
			Set<String> targList = MapUtil.loadMap(targets, 0, 1, false).keySet();

			// Find the node scores
			File nsDir = new File(nodeScoresDir);
			File[] nsArray = nsDir.listFiles(new FilenameFilter()
			{
				public boolean accept(File dir, String name)
				{
					return (name.matches(nodeScoresPrefix + ".*"));
				}
			});
			if(nsArray.length != 1)
			{
				throw new IOException("Could not find unique node score file");
			}


			// Load the node score of all proteins, including sources and targets
			// Use the % of top paths for the node score
			HashMap<String, String> nodeScoreMap = MapUtil.loadMap(nsArray[0].getCanonicalPath(), 0, 5, true);
			System.out.println("Loaded node scores from " + nsArray[0].getCanonicalPath());
			
			// Determine which proteins have node scores >= the threshold
			HashSet<String> importantNodes = new HashSet<String>();
			for(String node : nodeScoreMap.keySet())
			{
				double score = Double.parseDouble(nodeScoreMap.get(node));
				if(score >= nodeScoreThresh)
				{
					importantNodes.add(node);
				}
			}
			
			importantNodes.removeAll(srcList);
			importantNodes.removeAll(targList);
			
			nodesMap.put("sources", srcList);
			nodesMap.put("targets", targList);
			nodesMap.put("internal", importantNodes);
			// Don't know anything about the genes in this method
//			nodesMap.put("genes", new HashSet<String>());
			
			return nodesMap;
		}
		catch(IOException e)
		{
			e.printStackTrace();
			return null;
		}
	}
	
	/**
	 * Load the predictions made by ResponseNet.  Stores all proteins and genes as
	 * ORF names.  Originally written to support the output from the ResponseNet
	 * server around 5/2011, then updated to support the slightly different files
	 * returned by the sever on 6/2012.
	 * @param edgeFile
	 * @return
	 */
	public static HashMap<String, HashSet<String>> loadResponseNet(String edgeFile)
	{
		try
		{
			HashMap<String, HashSet<String>> nodesMap = new HashMap<String, HashSet<String>>();
			HashSet<String> sources = new HashSet<String>();
			HashSet<String> targets = new HashSet<String>();
			HashSet<String> internal = new HashSet<String>();
			HashSet<String> genes = new HashSet<String>();
			
			BufferedReader reader = new BufferedReader(new FileReader(edgeFile));
			// Skip the header
			String line = reader.readLine();
			
			// Read all edges and determine which category the two endpoints belong to
			while((line = reader.readLine()) != null)
			{
				String[] parts = line.split("\t");
				
				String edgeType = parts[3];
				String prot1 = parts[4];
				String prot2 = parts[5];
				
				if(edgeType.equals("SToSources"))
				{
					sources.add(MapUtil.getOrf(prot2));
				}
				else if(edgeType.equals("PPI"))
				{
					internal.add(MapUtil.getOrf(prot1));
					internal.add(MapUtil.getOrf(prot2));
				}
				else if(edgeType.equals("TFToGene"))
				{
					targets.add(MapUtil.getOrf(prot1));
				}
				// Not sure why the server output now repeats this edge type
				else if(edgeType.equals("TargetToT") || edgeType.equals("TargetToT  TargetToT"))
				{
					if(prot1.endsWith("_Gene"))
					{
						prot1 = prot1.substring(0, prot1.lastIndexOf("_Gene"));
					}
					else if(prot1.endsWith("_g"))
					{
						prot1 = prot1.substring(0, prot1.lastIndexOf("_g"));
					}
					genes.add(MapUtil.getOrf(prot1));
				}
				else
				{
					throw new IllegalArgumentException("Unrecognzied edge type: " + edgeType);
				}
			}
			
			// Don't want sources or targets to be included in the set of internal nodes
			internal.removeAll(sources);
			internal.removeAll(targets);
			
			nodesMap.put("sources", sources);
			nodesMap.put("targets", targets);
			nodesMap.put("internal", internal);
			nodesMap.put("genes", genes);
			
			return nodesMap;
		}
		catch(IOException e)
		{
			e.printStackTrace();
			return null;
		}
	}
	
	
	/**
	 * Load the predictions made by Physical Network Models.  Stores all proteins
	 * as ORF names.  The set of targets is always empty because PNM doesn't use the
	 * concept of targets.  The set of sources returned is the intersection between
	 * the nodes in the PNM network and the sources in the sourceFile
	 * @param edgeFile
	 * @param sourceFile PNM does not output which nodes are sources (knocked out nodes)
	 * so this needs to be provided.
	 * @return
	 */
	public static HashMap<String, HashSet<String>> loadPNM(String edgeFile, String sourceFile)
	{
		try
		{
			HashSet<String> potentialSources = MapUtil.loadSet(sourceFile, false);
			
			HashMap<String, HashSet<String>> nodesMap = new HashMap<String, HashSet<String>>();
			HashSet<String> sources = new HashSet<String>();
			HashSet<String> targets = new HashSet<String>();
			HashSet<String> internal = new HashSet<String>();
			// Could load the genes, but don't need them at the present time
//			HashSet<String> genes = new HashSet<String>();
			
			BufferedReader reader = new BufferedReader(new FileReader(edgeFile));
			// Skip the header
			String line = reader.readLine();
			
			// Read all edges and store the endpoints
			while((line = reader.readLine()) != null)
			{
				String[] parts = line.split(" ");

				// Don't care about the type of edge or the direction
				String prot1 = parts[0];
				String prot2 = parts[2];
				
				internal.add(MapUtil.getOrf(prot1));
				internal.add(MapUtil.getOrf(prot2));
			}
			
			sources.addAll(MapUtil.intersection(internal, potentialSources));
			
			// Don't want sources or targets to be included in the set of internal nodes
			internal.removeAll(sources);
			internal.removeAll(targets);

			nodesMap.put("sources", sources);
			nodesMap.put("targets", targets);
			nodesMap.put("internal", internal);
//			nodesMap.put("genes", genes);
			
			return nodesMap;
		}
		catch(IOException e)
		{
			e.printStackTrace();
			return null;
		}
	}
	
	/**
	 * Load post-processed predictions from GeneReg (targets only).
	 * Saves all targets that regulate at least a certain number of genes.
	 * @param regFile 
	 * @param threshold the gene target threshold 
	 * @return
	 */
	public static HashMap<String, HashSet<String>> loadGeneReg(String regFile, int threshold)
	{	
		HashMap<String, HashSet<String>> nodesMap = new HashMap<String, HashSet<String>>();
		HashSet<String> targets = new HashSet<String>();

		HashMap<String, String> targInts = MapUtil.loadMap(regFile, 1, 2, true);
		for(String tf : targInts.keySet())
		{
			int interactions = Integer.parseInt(targInts.get(tf));
			if(interactions >= threshold)
			{
				targets.add(tf);
			}
		}

		nodesMap.put("sources", new HashSet<String>());
		nodesMap.put("targets", targets);
		nodesMap.put("internal", new HashSet<String>());

		return nodesMap;
	}

	/**
	 * Examines the set of proteins that have been identified as important by the
	 * edge orientation and writes a file that contains their p-values in a sorbitol
	 * screen and other information.  Excludes any nodes that have been identified as
	 * targets by DREM, the sources nodes, and any protein that is previously known to
	 * be HOG-related according to KEGG, Science Signaling, or Krantz2009.  For
	 * each internal node identified, uses the set of satisfied paths from an orientation
	 * to write the set of targets for paths that involve the node.
	 * @param sources the predefined sources.  Should be
	 * ORF names.
	 * @param targets a file containing proteins detected by the edge orientation + 
	 * DREM model.  Should be ORF names.
	 * @param nodeScoresDir the directory containing the node scores file
	 * @param nodescoresPrefix the first part of the node scores file from edge orientation.
	 * One file in the nodeScoresDir should match this prefix.
	 * @param nodeScoreThresh a percent threshold.  Nodes with a % top paths >= this
	 * threshold are deemed important.
	 * @param keggMembers a list of genes known to be involved in the HOG pathway
	 * based on the KEGG database.  Uses standard names.
	 * @param sciSigMembers a list of genes known to be involved in the HOG pathway
	 * based on the Science Signaling database.  Uses standard names.
	 * @param sorbitolScreen genes and their lowest p-value over multiple replicates
	 * of a sorbitol screen
	 * @param krantzMembers members of the HOG pathway as depicted in Figure 1 of Krantz2009
	 * @param tfList a list of proteins denoted as transcription factors per Harbison2004 (standard names)
	 * @param satisfiedPaths the set of satisfied paths after a network orientation
	 * that will be used to identify targets of the paths that contain the nodes of interest
	 * @param out
	 */
	public static void hogInternalNodes(String sources, String targets, 
			String nodeScoresDir, final String nodeScoresPrefix,
			double nodeScoreThresh,
			String keggMembers, String sciSigMembers, 
			String sorbitolScreen, String krantzMembers,
			String tfList,
			ArrayList<Path> satisfiedPaths, String out)
	{
		try
		{
			HashSet<String> srcList = MapUtil.loadSet(sources, false);
			// Load the targets but ignore their target weights
			Set<String> targList = MapUtil.loadMap(targets, 0, 1, false).keySet();

			// Find the node scores
			File nsDir = new File(nodeScoresDir);
			File[] nsArray = nsDir.listFiles(new FilenameFilter()
			{
				public boolean accept(File dir, String name)
				{
					return (name.matches(nodeScoresPrefix + ".*"));
				}
			});
			if(nsArray.length != 1)
			{
				throw new IOException("Could not find unique node score file");
			}


			// Load the node score of all proteins, including sources and targets
			// Use the % of top paths for the node score
			HashMap<String, String> nodeScoreMapTop = MapUtil.loadMap(nsArray[0].getCanonicalPath(), 0, 5, true);
			// Also so the degree for each node and the % of all paths through the node
			HashMap<String, String> nodeScoreMapAll = MapUtil.loadMap(nsArray[0].getCanonicalPath(), 0, 7, true);
			HashMap<String, String> nodeDegreeMap = MapUtil.loadMap(nsArray[0].getCanonicalPath(), 0, 4, true);

			// Determine which proteins have node scores >= the threshold
			HashSet<String> importantNodes = new HashSet<String>();
			for(String node : nodeScoreMapTop.keySet())
			{
				double score = Double.parseDouble(nodeScoreMapTop.get(node));
				if(score >= nodeScoreThresh)
				{
					importantNodes.add(node);
				}
			}

			HashSet<String> keggStandard = MapUtil.loadSet(keggMembers);
			// Translate from standard to ORF
			HashSet<String> keggSet = new HashSet<String>();
			for(String gene : keggStandard)
			{
				keggSet.add(MapUtil.getOrf(gene.toUpperCase().trim()));
			}

			HashSet<String> sciSigStandard = MapUtil.loadSet(sciSigMembers);
			// Translate from standard to ORF
			HashSet<String> sciSigSet = new HashSet<String>();
			for(String gene : sciSigStandard)
			{
				sciSigSet.add(MapUtil.getOrf(gene.toUpperCase().trim()));
			}

			HashMap<String, String> sorbitolMap = MapUtil.loadMap(sorbitolScreen, 0, 1, false);

			HashSet<String> krantzStandard = MapUtil.loadSet(krantzMembers);
			// Translate from standard to ORF
			HashSet<String> krantzSet = new HashSet<String>();
			for(String gene : krantzStandard)
			{
				krantzSet.add(MapUtil.getOrf(gene.toUpperCase().trim()));
			}
			
			HashSet<String> knownHog = new HashSet<String>();
			knownHog.addAll(keggSet);
			knownHog.addAll(sciSigSet);
			knownHog.addAll(krantzSet);

			
			HashSet<String> tfsStandard = MapUtil.loadSet(tfList);
			// Translate from standard to ORF
			HashSet<String> tfSet = new HashSet<String>();
			for(String gene : tfsStandard)
			{
				tfSet.add(MapUtil.getOrf(gene.toUpperCase().trim()));
			}

			PrintWriter writer = new PrintWriter(new FileWriter(out));
			writer.println("Gene\tORF name\tTF\t% top paths\t% all paths" +
			"\tDegree\tSorbitol screen\tConnected targets");

			// Loop through all the important nodes found
			for(String protein : importantNodes)
			{
				// Only interested in nodes that are not sources,
				// targets, or known HOG members
				if(!(srcList.contains(protein) || targList.contains(protein) ||
						knownHog.contains(protein)))
				{
					StringBuffer buf = new StringBuffer();
					String std = MapUtil.getStandard(protein);

					buf.append(std).append("\t").append(protein).append("\t");
					
					if(tfSet.contains(protein))
					{
						buf.append("Y");
					}
					else
					{
						buf.append("N");
					}
					buf.append("\t");

					// Write the % of paths through the node
					buf.append(nodeScoreMapTop.get(protein)).append("\t");
					buf.append(nodeScoreMapAll.get(protein)).append("\t");

					// Write the degree
					buf.append(nodeDegreeMap.get(protein)).append("\t");

					// Print the p-value of this protein in the sorbitol screen
					if(sorbitolMap.containsKey(protein))
					{
						String pval = sorbitolMap.get(protein);
						buf.append(pval);
					}
					else
					{
						buf.append("N/A");
					}
					buf.append("\t");

					// Find and print all targets connnected by paths through this node
					ArrayList<String> connectedTargs = findTargets(protein, satisfiedPaths);
					ArrayList<String> connectedTargsOrf = new ArrayList<String>();
					for(String ct : connectedTargs)
					{
						connectedTargsOrf.add(MapUtil.getStandard(ct));
					}
					Collections.sort(connectedTargsOrf);
					for(String cto : connectedTargsOrf)
					{
						buf.append(cto).append(":");
					}
					buf.deleteCharAt(buf.length()-1);

					writer.println(buf);
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
	 * Compare SDREM output by analyzing how many predicted targets and internal
	 * nodes are common.  The first file in the list is taken as the reference set
	 * of predicitons and all other predictions are compared to it.
	 * @param modelNames A name for each predicted model
	 * @param sourceFile Assume the same sources are used for all predictions
	 * @param targetFiles A list of predcted targets
	 * @param nodeScoresDirs A list of directories that contain the node scores
	 * @param nodeScoresPrefix Assume a common node score filename, but this could
	 * be relaxed later
	 * @param nodeScoreThreshs A list of node score thresholds to use
	 * @param summaryOut A summary of the overlaps
	 * @param fullOut The full set of predictions for each method and their overlaps
	 */
	public static void compareSDREMPredictions(String[] modelNames, 
			String sourceFile, String[] targetFiles, 
			String[] nodeScoresDirs, String nodeScoresPrefix,
			double[] nodeScoreThreshs, String summaryOut, String fullOut)
	{
		if(modelNames.length != targetFiles.length ||
				modelNames.length != nodeScoresDirs.length ||
				modelNames.length != nodeScoreThreshs.length)
		{
			throw new IllegalArgumentException("All input arrays must have the same length: " + modelNames.length);
		}
		
		try
		{
			// Load all of the predictions
			ArrayList<HashMap<String, Set<String>>> predictions = new ArrayList<HashMap<String, Set<String>>>();
			// Also track all proteins that are predicted in any method
			HashSet<String> allNodes = new HashSet<String>();
			
			for(int p = 0; p < targetFiles.length; p++)
			{
				HashMap<String, Set<String>> newPreds = loadSDREM(sourceFile, targetFiles[p],
						nodeScoresDirs[p], nodeScoresPrefix, nodeScoreThreshs[p]);
				predictions.add(newPreds);
				
				allNodes.addAll(newPreds.get("sources"));
				allNodes.addAll(newPreds.get("targets"));
				allNodes.addAll(newPreds.get("internal"));
			}
			
			// Loop through all proteins that have been predicted and print which models
			// they were included int
			PrintWriter fullWriter = new PrintWriter(new FileWriter(fullOut));
			fullWriter.print("Gene\tORF name");
			for(String name : modelNames)
			{
				fullWriter.print("\t" + name);
			}
			fullWriter.println();
			
			for(String gene : allNodes)
			{
				StringBuffer lineOut = new StringBuffer(MapUtil.getStandard(gene));
				lineOut.append("\t").append(gene);
				for(int p = 0; p < targetFiles.length; p++)
				{
					lineOut.append("\t");
					
					// Allow the possibility that a gene is in more than one of
					// the following categories
					HashMap<String, Set<String>> curModel = predictions.get(p);
					if(curModel.get("sources").contains(gene))
					{
						lineOut.append("S");
					}
					if(curModel.get("targets").contains(gene))
					{
						lineOut.append("T");
					}
					if(curModel.get("internal").contains(gene))
					{
						lineOut.append("I");
					}
				}
				fullWriter.println(lineOut);
			}
			fullWriter.close();
			
			// TODO separate this into TFs and non-TFs?
			// Compute the overlaps for the targets and internal nodes of the reference
			// model and all other models
			PrintWriter summaryWriter = new PrintWriter(new FileWriter(summaryOut));
			String refName = modelNames[0];
			HashSet<String> refPreds = new HashSet<String>();
			refPreds.addAll(predictions.get(0).get("targets"));
			refPreds.addAll(predictions.get(0).get("internal"));
			
			// Skip the first model, the refernce model
			for(int p = 1; p < targetFiles.length; p++)
			{
				summaryWriter.println("\t" + refName + "\t" + modelNames[p]);
				
				HashSet<String> curPreds = new HashSet<String>();
				curPreds.addAll(predictions.get(p).get("targets"));
				curPreds.addAll(predictions.get(p).get("internal"));
				
				summaryWriter.println("# pred\t" + refPreds.size() + "\t" + curPreds.size());
				
				// Compute the number of predictions in one model that are also found
				// in the other
				int overlap = MapUtil.intersection(refPreds, curPreds).size();
				summaryWriter.println("Overlap\t" + overlap + "\t" + overlap);
				summaryWriter.println("% overlap\t" + ((double)overlap/refPreds.size()) + "\t" + ((double)overlap/curPreds.size()));
				summaryWriter.println();
			}
			summaryWriter.close();
		}
		catch(IOException e)
		{
			e.printStackTrace();
		}
	}
	
	/**
	 * Compare SDREM output by analyzing how many predicted targets and internal
	 * nodes are common.  Pairwise comparisons are performed for all models.
	 * Does not consider any external lists of nodes to be used for validation
	 * @param modelNames A name for each predicted model
	 * @param sourceFiles The sets of sources given to SDREM in each run
	 * @param targetFiles A list of predcted targets
	 * @param nodeScoresDirs A list of directories that contain the node scores, the
	 * first of which will also be used to obtain the node degree
	 * @param nodeScoresPrefix Assume a common node score filename, but this could
	 * be relaxed later
	 * @param nodeScoreThreshs A list of node score thresholds to use
	 * @param summaryOut A summary of the overlaps
	 * @param fullOut The full set of predictions for each method and their overlaps
	 */
	public static void compareSDREMPredictionsHuman(String[] modelNames, 
			String[] sourceFiles, String[] targetFiles, 
			String[] nodeScoresDirs, String nodeScoresPrefix,
			double[] nodeScoreThreshs, String summaryOut, String fullOut)
	{
		
	}
	
	/**
	 * Compare SDREM output by analyzing how many predicted targets and internal
	 * nodes are common.  Pairwise comparisons are performed for all models.
	 * @param modelNames A name for each predicted model
	 * @param sourceFiles The sets of sources given to SDREM in each run
	 * @param targetFiles A list of predcted targets
	 * @param nodeScoresDirs A list of directories that contain the node scores, the
	 * first of which will also be used to obtain the node degree
	 * @param nodeScoresPrefix Assume a common node score filename, but this could
	 * be relaxed later
	 * @param nodeScoreThreshs A list of node score thresholds to use
	 * @param validationNames  A list of the validation sets being use.  May be null.
	 * @param validationFiles A list of the files containing the named validation sets.
	 * May be null, but if validationNames or validationFiles is null, the other must
	 * be as well.  The validation files must not have headers.
	 * @param summaryOut A summary of the overlaps
	 * @param fullOut The full set of predictions for each method and their overlaps
	 */
	public static void compareSDREMPredictionsHuman(String[] modelNames, 
			String[] sourceFiles, String[] targetFiles, 
			String[] nodeScoresDirs, String nodeScoresPrefix,
			double[] nodeScoreThreshs, 
			String[] validationNames, String[] validationFiles,
			String summaryOut, String fullOut)
	{
		if(modelNames.length != sourceFiles.length ||
				modelNames.length != targetFiles.length ||
				modelNames.length != nodeScoresDirs.length ||
				modelNames.length != nodeScoreThreshs.length)
		{
			throw new IllegalArgumentException("All input arrays must have the same length: " + modelNames.length);
		}
		
		// Make sure that if either parameter is given, they both are
		if(validationNames == null ^ validationFiles == null)
		{
			throw new IllegalArgumentException("Must provide both validation names and files");
		}
		boolean haveValidation = (validationNames != null);
		if(haveValidation && (validationNames.length != validationFiles.length))
		{
			throw new IllegalArgumentException("Must provide the same number of validation names and files");
		}
		
		try
		{
			// Load all of the predictions
			ArrayList<HashMap<String, Set<String>>> predictions = new ArrayList<HashMap<String, Set<String>>>();
			// Also track all proteins that are predicted in any method
			HashSet<String> allNodes = new HashSet<String>();
			
			for(int p = 0; p < targetFiles.length; p++)
			{
				HashMap<String, Set<String>> newPreds = loadSDREM(sourceFiles[p], targetFiles[p],
						nodeScoresDirs[p], nodeScoresPrefix, nodeScoreThreshs[p]);
				predictions.add(newPreds);
				
				allNodes.addAll(newPreds.get("sources"));
				allNodes.addAll(newPreds.get("targets"));
				allNodes.addAll(newPreds.get("internal"));
			}
			
			// Load the degree of all the nodes
			HashMap<String, String> degreeMap = loadDegreeMap(nodeScoresDirs, nodeScoresPrefix);
			
			// Load the validation sets if specified
			ArrayList<HashMap<String, HashSet<String>>> validationSets = null;
			if(haveValidation)
			{
				validationSets = new ArrayList<HashMap<String, HashSet<String>>>();
				for(int v = 0; v < validationNames.length; v++)
				{
					validationSets.add(MapUtil.loadSets(validationFiles[v], false));
				}
			}
			
			// Loop through all proteins that have been predicted and print which models
			// they were included in
			PrintWriter fullWriter = new PrintWriter(new FileWriter(fullOut));
			fullWriter.print("Gene\tGene id\tDegree");
			for(String name : modelNames)
			{
				fullWriter.print("\t" + name);
			}
			fullWriter.print("\t# models");
			
			if(haveValidation)
			{
				for(int v = 0; v < validationNames.length; v++)
				{
					fullWriter.print("\t" + validationNames[v]);
				}
			}
			fullWriter.println();
			
			for(String geneId : allNodes)
			{				
				String gene = Entrez.getSymbol(geneId);
				if(gene == null)
				{
					System.err.println("Could not lookup gene id " + geneId);
					gene = "LOC" + geneId;
				}
				StringBuffer lineOut = new StringBuffer(gene);
				lineOut.append("\t").append(geneId);
				
				// Write the degree of the gene
				if(!degreeMap.containsKey(geneId))
				{
					throw new IllegalStateException("Cannot find degree of " + geneId);
				}
				lineOut.append("\t").append(degreeMap.get(geneId));
				
				int numModels = 0;
				for(int p = 0; p < targetFiles.length; p++)
				{
					lineOut.append("\t");
					
					// Allow the possibility that a gene is in more than one of
					// the following categories
					boolean inModel = false;
					HashMap<String, Set<String>> curModel = predictions.get(p);
					if(curModel.get("sources").contains(geneId))
					{
						lineOut.append("S");
						inModel = true;
					}
					if(curModel.get("targets").contains(geneId))
					{
						lineOut.append("T");
						inModel = true;
					}
					if(curModel.get("internal").contains(geneId))
					{
						lineOut.append("I");
						inModel = true;
					}
					
					// Have to account for the possibility that a node is both
					// a source and a target
					if(inModel)
					{
						numModels++;
					}
				}
				lineOut.append("\t").append(numModels);
				
				if(haveValidation)
				{
					for(int v = 0; v < validationNames.length; v++)
					{
						// TODO Could check if this type of validation has more than one possible
						// sets and print a number for multiple validation sets and Y/N
						// otherwise
						int validHits = 0;
						// The current type of validation is a map from set names (e.g. screens)
						// to set members (e.g. hits in those screens)
						HashMap<String, HashSet<String>> curSet = validationSets.get(v);
						for(String validSet : curSet.keySet())
						{
							if(curSet.get(validSet).contains(geneId))
							{
								validHits++;
							}
						}
						lineOut.append("\t").append(validHits);
					}
				}
				
				fullWriter.println(lineOut);
			}
			fullWriter.close();

			// Compute the overlaps for the sources, targets, and internal nodes of
			// all pairs of models
			PrintWriter summaryWriter = new PrintWriter(new FileWriter(summaryOut));
			for(int r = 0; r < targetFiles.length - 1; r++)
			{
				// Call the first model being compared the reference
				String refName = modelNames[r];
				Set<String> refSources = predictions.get(r).get("sources");
				// Predictions are the internal nodes and targets
				HashSet<String> refPreds = new HashSet<String>();
				refPreds.addAll(predictions.get(r).get("targets"));
				refPreds.addAll(predictions.get(r).get("internal"));
				// Sources and predictions
				HashSet<String> refAll = new HashSet<String>(refPreds);
				refAll.addAll(refSources);

				// Skip the reference model
				for(int p = r + 1; p < targetFiles.length; p++)
				{
					summaryWriter.println("\t" + refName + "\t" + modelNames[p]);

					Set<String> curSources = predictions.get(p).get("sources");

					HashSet<String> curPreds = new HashSet<String>();
					curPreds.addAll(predictions.get(p).get("targets"));
					curPreds.addAll(predictions.get(p).get("internal"));

					HashSet<String> curAll = new HashSet<String>(curPreds);
					curAll.addAll(curSources);

					summaryWriter.println("# sources\t" + refSources.size() + "\t" + curSources.size());
					// Compute the number of sources in the reference model that are also found
					// in the other
					int sourceOvr = MapUtil.intersection(refSources, curSources).size();
					summaryWriter.println("Source overlap\t" + sourceOvr + "\t" + sourceOvr);
					summaryWriter.println("% source overlap\t" + ((double)sourceOvr/refSources.size()) + "\t" + ((double)sourceOvr/curSources.size()));

					summaryWriter.println("# targets\t" + predictions.get(r).get("targets").size() + "\t" + predictions.get(p).get("targets").size());
					summaryWriter.println("# internal\t" + predictions.get(r).get("internal").size() + "\t" + predictions.get(p).get("internal").size());
					summaryWriter.println("# predictions\t" + refPreds.size() + "\t" + curPreds.size());
					// Compute the number of predictions in the reference model that are also found
					// in the other
					int overlap = MapUtil.intersection(refPreds, curPreds).size();
					summaryWriter.println("Prediction overlap\t" + overlap + "\t" + overlap);
					summaryWriter.println("% prediction overlap\t" + ((double)overlap/refPreds.size()) + "\t" + ((double)overlap/curPreds.size()));

					summaryWriter.println("# model members\t" + refAll.size() + "\t" + curAll.size());
					// Compute the number of predictions in the reference model that are also found
					// in the other
					int allOvr = MapUtil.intersection(refAll, curAll).size();
					summaryWriter.println("Model overlap\t" + allOvr + "\t" + allOvr);
					summaryWriter.println("% model overlap\t" + ((double)allOvr/refAll.size()) + "\t" + ((double)allOvr/curAll.size()));

					summaryWriter.println();
				}
			}
			summaryWriter.close();
		}
		catch(IOException e)
		{
			e.printStackTrace();
		}
	}
	
	public static HashMap<String, String> loadDegreeMap(String[] nodeScoresDirs, final String nodeScoresPrefix)
		throws IOException
	{
		HashMap<String, String> degreeMap = new HashMap<String, String>();
		
		for(int d = 0; d < nodeScoresDirs.length; d++)
		{
			// Load the degrees from model d
			File nsDir = new File(nodeScoresDirs[d]);
			File[] nsArray = nsDir.listFiles(new FilenameFilter()
			{
				public boolean accept(File dir, String name)
				{
					return (name.matches(nodeScoresPrefix + ".*"));
				}
			});
			// Load the degrees from the node scores file
			HashMap<String, String> curDegrees = MapUtil.loadMap(nsArray[0].getCanonicalPath(), 0, 4, true);
			
			// Make sure that the degrees agree with previously loaded degrees during the merge
			for(String geneId : curDegrees.keySet())
			{
				if(degreeMap.containsKey(geneId))
				{
					if(!degreeMap.get(geneId).equals(curDegrees.get(geneId)))
					{
						throw new IllegalStateException("Degrees for " + geneId + " disagree");
					}
				}
				else
				{
					degreeMap.put(geneId, curDegrees.get(geneId));
				}
			}
			// Debug
//			System.out.println("Loaded " + (d+1) + " models and " + degreeMap.size() + " node degrees");
		}
		
		return degreeMap;
	}

	/**
	 * Analyzes iterations 1 through lastItr to see how the sets of targets change.
	 * The targets files should use ORF names.
	 * @param lastItr
	 * @param targWeightSuffix ".model.activities" or ".targets"
	 * @param targDir
	 * @param allTargsFile
	 * @param newTargsFile
	 * @param droppedTargsFile
	 */
	public static void analyzeTargetChanges(int lastItr, String targWeightSuffix,
			String targDir, String allTargsFile, String newTargsFile, 
			String droppedTargsFile)
	{
		try
		{
			if(!targDir.endsWith("/"))
			{
				targDir = targDir + "/";
			}

			PrintWriter allWriter = new PrintWriter(new FileWriter(allTargsFile));
			PrintWriter newWriter = new PrintWriter(new FileWriter(newTargsFile));
			PrintWriter droppedWriter = new PrintWriter(new FileWriter(droppedTargsFile));

			// pastTargs uses standard names
			HashSet<String> pastTargs = new HashSet<String>();
			// lastTargs uses ORF names
			Set<String> lastTargs = new HashSet<String>();
			for(int i = 1; i <= lastItr; i++)
			{
				// Use the map loader to read the targets file, then
				// disregard the target weights
				String curTargFile = targDir + i + targWeightSuffix;
				Set<String> curTargs = MapUtil.loadMap(curTargFile, 0, 1, false).keySet();


				// The lists to be printed will use standard names
				ArrayList<String> allTargs = new ArrayList<String>();
				ArrayList<String> newTargs = new ArrayList<String>();

				for(String curTarg : curTargs)
				{
					allTargs.add(MapUtil.getStandard(curTarg));
					if(!lastTargs.contains(curTarg))
					{
						newTargs.add(MapUtil.getStandard(curTarg));
					}
				}
				Collections.sort(allTargs);
				Collections.sort(newTargs);


				allWriter.print("Iteration " + i);
				for(int t = 0; t < allTargs.size(); t++)
				{
					allWriter.print("\t" + allTargs.get(t));
				}
				allWriter.println();


				newWriter.print("Iteration " + i);
				for(int t = 0; t < newTargs.size(); t++)
				{
					newWriter.print("\t" + newTargs.get(t));
				}
				newWriter.println();


				for(int t = 0; t < newTargs.size(); t++)
				{
					if(pastTargs.contains(newTargs.get(t)))
					{
						newWriter.print("\tPrev targ");
					}
					else
					{
						newWriter.print("\t-");
					}
				}
				newWriter.println();


				ArrayList<String> droppedTargs = new ArrayList<String>();
				for(String lastTarg : lastTargs)
				{
					if(!curTargs.contains(lastTarg))
					{
						droppedTargs.add(MapUtil.getStandard(lastTarg));
					}
				}
				Collections.sort(droppedTargs);

				droppedWriter.print("Iteration " + i);
				for(int t = 0; t < droppedTargs.size(); t++)
				{
					droppedWriter.print("\t" + droppedTargs.get(t));
				}
				droppedWriter.println();


				lastTargs = curTargs;

				// Don't add the current targets to the set of
				// all past targets until after checking if the current
				// new targets were previous targets
				for(String curTarg : curTargs)
				{
					pastTargs.add(MapUtil.getStandard(curTarg));
				}
			}

			allWriter.close();
			newWriter.close();
			droppedWriter.close();
		}
		catch(IOException e)
		{
			e.printStackTrace();
		}

	}

	/**
	 * Find all target nodes in the satisfied paths that go through a particular node
	 * @param node
	 * @param satisifedPaths
	 * @return
	 */
	public static ArrayList<String> findTargets(String node, ArrayList<Path> satisfiedPaths)
	{
		HashSet<String> targets = new HashSet<String>();

		for(Path p : satisfiedPaths)
		{
			for(Vertex v : p.getVertices())
			{
				if(node.equalsIgnoreCase(v.getName()))
				{
					targets.add(p.getTarget().getName());
				}
			}
		}

		ArrayList<String> targList = new ArrayList<String>();
		targList.addAll(targets);
		Collections.sort(targList);
		return targList;
	}
	
	/**
	 * Create a file that groups all paths that go through the node according
	 * to which target they end at.  Assumes yeast genes.
	 * @param node paths that do not contain this path will be ignored
	 * @param satisfiedPaths the set of paths to consider
	 * @param outFile
	 * @throws IOException
	 */
	public static void groupPathsByTarget(String node,
			ArrayList<Path> satisfiedPaths,
			String outFile) throws IOException
	{
		Vertex keyVert = Vertex.findVertex(MapUtil.getOrf(node));
		if(keyVert == null)
		{
			throw new IllegalStateException("Cannot find vertex " + node);
		}
		
		HashMap<Vertex, HashSet<Path>> pathMap = new HashMap<Vertex, HashSet<Path>>();
		for(Path p : satisfiedPaths)
		{
			if(p.containsVertex(keyVert))
			{
				Vertex targ = p.getTarget();
				HashSet<Path> paths;
				if(pathMap.containsKey(targ))
				{
					paths = pathMap.get(targ);
				}
				else
				{
					paths = new HashSet<Path>();
				}
				paths.add(p);
				pathMap.put(targ, paths);
			}
		}
		
		System.out.println("Found " + pathMap.size() + " targets downstream of " + node);
		
		PrintWriter writer = new PrintWriter(new FileWriter(outFile));
		for(Vertex target : pathMap.keySet())
		{
			writer.println(MapUtil.getStandard(target.getName()));
			for(Path p : pathMap.get(target))
			{
				for(Vertex v : p.getVertices())
				{
					writer.print(MapUtil.getStandard(v.getName()) + "\t");
				}
				writer.println();
			}
			writer.println();
		}
		writer.close();
	}

	
	/**
	 * Compare the three HOG gold standards.  Determine how many members
	 * are common/unique to the gold standards.
	 * @param kegg KEGG HOG members
	 * @param sciSig Science Signaling database HOG members
	 * @param literature HOG members from reviews and other literature
	 * @param summaryOut a summary of how many members are shared/unique
	 * @param unionOut a list of all gold standard members
	 */
	public static void hogGoldStandardOverlap(String kegg, String sciSig,
			String literature, String summaryOut, String unionOut)
	{
		try
		{
			// Load all gold standards
			HashSet<String> keggStandard = MapUtil.loadSet(kegg);
			// Translate from standard to ORF
			HashSet<String> keggSet = new HashSet<String>();
			for(String gene : keggStandard)
			{
				keggSet.add(MapUtil.getOrf(gene.toUpperCase().trim()));
			}
			
			HashSet<String> sciSigStandard = MapUtil.loadSet(sciSig);
			// Translate from standard to ORF
			HashSet<String> sciSigSet = new HashSet<String>();
			for(String gene : sciSigStandard)
			{
				sciSigSet.add(MapUtil.getOrf(gene.toUpperCase().trim()));
			}
			
			HashSet<String> literatureStandard = MapUtil.loadSet(literature);
			// Translate from standard to ORF
			HashSet<String> literatureSet = new HashSet<String>();
			for(String gene : literatureStandard)
			{
				literatureSet.add(MapUtil.getOrf(gene.toUpperCase().trim()));
			}
			
			HashSet<String> allMembers = new HashSet<String>();
			allMembers.addAll(keggSet);
			allMembers.addAll(sciSigSet);
			allMembers.addAll(literatureSet);
			
			PrintWriter unionWriter = new PrintWriter(new FileWriter(unionOut));
			unionWriter.println("Standard\tORF\tKEGG\tScience Signaling\tLiterature");
			
			// Use binary encoding to map gold standard members to their correct bins
			int[] matches = new int[8]; 
			for(String gene : allMembers)
			{
				int index = 0;
				
				unionWriter.print(MapUtil.getStandard(gene) + "\t" + gene + "\t");
				if(keggSet.contains(gene))
				{
					unionWriter.print("X");
					index += 1;
				}
				
				unionWriter.print("\t");
				if(sciSigSet.contains(gene))
				{
					unionWriter.print("X");
					index += 2;
				}
				
				unionWriter.print("\t");
				if(literatureSet.contains(gene))
				{
					unionWriter.print("X");
					index += 4;
				}
				unionWriter.println();
				
				matches[index]++;
			}
			
			unionWriter.close();
			
			PrintWriter summaryWriter = new PrintWriter(new FileWriter(summaryOut));
			String[] categories = {"None", "KEGG only", "Science Signaling only",
					"KEGG + Science Signaling", "Literature only", "KEGG + Literature",
					"Science Signaling + Literature", "KEGG + Science Signagling + Literature"};
			for(int i = 0; i < 8; i++)
			{
				summaryWriter.println(matches[i] + ": " + categories[i]);
			}
			summaryWriter.close();
		}
		catch(IOException e)
		{
			e.printStackTrace();
		}
	}
	
	/**
	 * Calculate how many gold standard sources each protein appears in, and annotate
	 * which are TFs
	 * @param gsMembersFile standard gene names
	 * @param tfList standard gene names
	 * @param outFile
	 */
	public static void goldStandardOverlapDetailed(String gsMembersFile, String tfList, String outFile)
	{
		try
		{
			// Load all of the sets of gold standard members
			HashMap<String, HashSet<String>> gsMap = MapUtil.loadSets(gsMembersFile, false);
			ArrayList<String> gsOrder = new ArrayList<String>(gsMap.keySet());
			Collections.sort(gsOrder);
			HashSet<String> gsMembers = new HashSet<String>();
			
			HashSet<String> tfs = MapUtil.loadSet(tfList, false);
			
			PrintWriter writer = new PrintWriter(new FileWriter(outFile));
			writer.print("Standard\tORF\tTF");
			for(String gs : gsOrder)
			{
				writer.print("\t" + gs);
				
				gsMembers.addAll(gsMap.get(gs));
			}
			writer.println("\tCount");
			
			int count;
			for(String protein : gsMembers)
			{
				count = 0;
				
				writer.print(protein + "\t" + MapUtil.getOrf(protein));
				
				if(tfs.contains(protein.toUpperCase()))
				{
					writer.print("\tY");
				}
				else
				{
					writer.print("\tN");
				}
				
				for(String gs : gsOrder)
				{
					writer.print("\t");
					if(gsMap.get(gs).contains(protein))
					{
						count++;
						writer.print("X");
					}
				}
				
				writer.println("\t" + count);
			}
			
			writer.close();
		}
		catch(IOException e)
		{
			e.printStackTrace();
		}
	}
	
	/**
	 * Calculates p-values for the overlap of the HOG predictions and the
	 * gold standards.  Separate p-values are calculated for targets identified
	 * by DREM and the other signaling proteins.  This version assumes no nodes
	 * were deleted (held out) from the interaction network).  The interaction
	 * network genes and gold standard genes are used to define the set of all genes.
	 * @param sources the predefined sources.  Should be
	 * ORF names.
	 * @param targets a file containing proteins detected by the edge orientation + 
	 * DREM model.  Should be ORF names.
	 * @param nodeScoresDir the directory containing the node scores file
	 * @param nodescoresPrefix the first part of the node scores file from edge orientation.
	 * One file in the nodeScoresDir should match this prefix.
	 * @param nodeScoreThresh a percent threshold.  Nodes with a % top paths >= this
	 * threshold are deemed important.
	 * @param keggMembers a list of genes known to be involved in the HOG pathway
	 * based on the KEGG database.  Uses standard names.
	 * @param sciSigMembers a list of genes known to be involved in the HOG pathway
	 * based on the Science Signaling database.  Uses standard names.
	 * @param krantzMembers members of the HOG pathway as depicted in Figure 1 of Krantz2009
	 * @param tfList all proteins that are considered TFs (i.e. potential targets).  Standard names.
	 * Hot1 and Sko1 are hard-coded to be in this list.  Hog1 is removed from this list.
	 * @param physicalNetwork the PPI and protein-DNA network is used to build the list of all genes
	 * @param out
	 */
	public static void hogPredictionSignificance(String sources, String targets, 
			String nodeScoresDir, final String nodeScoresPrefix,
			double nodeScoreThresh,
			String keggMembers, String sciSigMembers,
			String krantzMembers, String tfList, String physicalNetwork,
			String out)
	{
		hogPredictionSignificance(sources, targets, nodeScoresDir, nodeScoresPrefix,
				nodeScoreThresh, null, keggMembers, sciSigMembers, krantzMembers,
				tfList, physicalNetwork, out);
	}
	
	/**
	 * Calculates p-values for the overlap of the HOG predictions and the
	 * gold standards.  Separate p-values are calculated for targets identified
	 * by DREM and the other signaling proteins.  The interaction
	 * network genes and gold standard genes are used to define the set of all genes.
	 * @param sources the predefined sources.  Should be
	 * ORF names.
	 * @param targets a file containing proteins detected by the edge orientation + 
	 * DREM model.  Should be ORF names.
	 * @param nodeScoresDir the directory containing the node scores file
	 * @param nodescoresPrefix the first part of the node scores file from edge orientation.
	 * One file in the nodeScoresDir should match this prefix.
	 * @param nodeScoreThresh a percent threshold.  Nodes with a % top paths >= this
	 * threshold are deemed important.
	 * @param deletedNodesFile a file listing the nodes that were deleted and held out
	 * of the interaction network.  These will be removed from the gold standard and
	 * physical network during
	 * evaluation because it is not possible for the algorithm to recover them.
	 * Should be specified as ORF names.
	 * @param keggMembers a list of genes known to be involved in the HOG pathway
	 * based on the KEGG database.  Uses standard names.
	 * @param sciSigMembers a list of genes known to be involved in the HOG pathway
	 * based on the Science Signaling database.  Uses standard names.
	 * @param krantzMembers members of the HOG pathway as depicted in Figure 1 of Krantz2009
	 * @param tfList all proteins that are considered TFs (i.e. potential targets).  Standard names.
	 * Hot1 and Sko1 are hard-coded to be in this list.  Hog1 is removed from this list.
	 * @param physicalNetwork the PPI and protein-DNA network is used to build the list of all genes
	 * @param out
	 */
	public static void hogPredictionSignificance(String sources, String targets, 
			String nodeScoresDir, final String nodeScoresPrefix,
			double nodeScoreThresh, String deletedNodesFile,
			String keggMembers, String sciSigMembers,
			String krantzMembers, String tfList, String physicalNetwork,
			String out)
	{
		try
		{
			// Load the targets but ignore their target weights
			Set<String> targList = MapUtil.loadMap(targets, 0, 1, false).keySet();
			HashSet<String> srcList = MapUtil.loadSet(sources, false);

			// Find the node scores
			File nsDir = new File(nodeScoresDir);
			File[] nsArray = nsDir.listFiles(new FilenameFilter()
			{
				public boolean accept(File dir, String name)
				{
					return (name.matches(nodeScoresPrefix + ".*"));
				}
			});
			if(nsArray.length != 1)
			{
				throw new IOException("Could not find unique node score file");
			}


			// Load the node score of all proteins, including sources and targets
			// Use the % of top paths for the node score
			HashMap<String, String> nodeScoreMap = MapUtil.loadMap(nsArray[0].getCanonicalPath(), 0, 5, true);

			// Determine which proteins have node scores >= the threshold
			HashSet<String> importantNodes = new HashSet<String>();
			for(String node : nodeScoreMap.keySet())
			{
				double score = Double.parseDouble(nodeScoreMap.get(node));
				if(score >= nodeScoreThresh)
				{
					importantNodes.add(node);
				}
			}
			// Want to exclude the sources and targets from the list
			importantNodes.removeAll(srcList);
			importantNodes.removeAll(targList);

			// Use the other overloaded version of the method to perform the comparison
			hogPredictionSignificance(srcList, targList, importantNodes, deletedNodesFile,
					keggMembers, sciSigMembers, krantzMembers, tfList, physicalNetwork,
					out);
		}
		catch(IOException e)
		{
			e.printStackTrace();
		}
	}
	
	/**
	 * Calculates p-values for the overlap of the HOG predictions and the
	 * gold standards.  Separate p-values are calculated for targets (TF)
	 * and the other signaling proteins.  The interaction
	 * network genes and gold standard genes are used to define the set of all genes.
	 * @param sources the predefined sources (ORF names)
	 * @param targets a file containing proteins detected by edge orientation + 
	 * DREM or another computational method (ORF names)
	 * @param internalNodes the set of internal nodes included in the predicted
	 * network that are neither sources nor targets (ORF names)
	 * @param nodeScoresDir the directory containing the node scores file
	 * @param nodescoresPrefix the first part of the node scores file from edge orientation.
	 * One file in the nodeScoresDir should match this prefix.
	 * @param nodeScoreThresh a percent threshold.  Nodes with a % top paths >= this
	 * threshold are deemed important.
	 * @param deletedNodesFile a file listing the nodes that were deleted and held out
	 * of the interaction network.  These will be removed from the gold standard during
	 * evaluation because it is not possible for the algorithm to recover them.
	 * Should be specified as ORF names.
	 * @param keggMembers a list of genes known to be involved in the HOG pathway
	 * based on the KEGG database.  Uses standard names.
	 * @param sciSigMembers a list of genes known to be involved in the HOG pathway
	 * based on the Science Signaling database.  Uses standard names.
	 * @param krantzMembers members of the HOG pathway as depicted in Figure 1 of Krantz2009.  Uses standard names.
	 * @param tfList all proteins that are considered TFs (i.e. potential targets).  Standard names.
	 * Hot1 and Sko1 are hard-coded to be in this list.  Hog1 is removed from this list.
	 * @param physicalNetwork the PPI and protein-DNA network is used to build the list of all genes
	 * @param out
	 */
	public static void hogPredictionSignificance(Set<String> sources,
			Set<String> targets, Set<String> internalNodes,
			String deletedNodesFile,
			String keggMembers, String sciSigMembers,
			String krantzMembers, String tfList, String physicalNetwork,
			String out)
	{
		try
		{
			// Want to exclude the sources and targets from the list of internal nodes
			internalNodes.removeAll(sources);
			internalNodes.removeAll(targets);			

			// Load the list of TFs
			HashSet<String> allTfsStandard = MapUtil.loadSet(tfList, false);
			HashSet<String> allTfs = new HashSet<String>();
			for(String tf : allTfsStandard)
			{
				allTfs.add(MapUtil.getOrf(tf));
			}
			// Manually add the two TFs with HOG-specific binding data
			allTfs.add("YMR172W"); // Hot1
			allTfs.add("YNL167C"); // Sko1
			// Remove Hog1, which is a MAPK but has binding data
			allTfs.remove("YLR113W");
			
			// Load the set of all genes from the network
			// This will be complemented by the set of all gold standard genes
			Set<String> edges1 = MapUtil.loadMap(physicalNetwork, 0, 1, false).keySet();
			Set<String> edges2 = MapUtil.loadMap(physicalNetwork, 2, 1, false).keySet();
			HashSet<String> allGenes = MapUtil.union(edges1, edges2);
//			System.out.println(edges1.size());
//			System.out.println(edges2.size());

			// Load the gold standard
			HashSet<String> keggStandard = MapUtil.loadSet(keggMembers);
			// Translate from standard to ORF
			HashSet<String> goldStandard = new HashSet<String>();
			for(String gene : keggStandard)
			{
				goldStandard.add(MapUtil.getOrf(gene.toUpperCase().trim()));
			}

			HashSet<String> sciSigStandard = MapUtil.loadSet(sciSigMembers);
			// Translate from standard to ORF
			for(String gene : sciSigStandard)
			{
				goldStandard.add(MapUtil.getOrf(gene.toUpperCase().trim()));
			}

			HashSet<String> krantzStandard = MapUtil.loadSet(krantzMembers);
			// Translate from standard to ORF
			for(String gene : krantzStandard)
			{
				goldStandard.add(MapUtil.getOrf(gene.toUpperCase().trim()));
			}
			
			// Load the deleted nodes and remove them from the network
			// and the gold standard
			HashSet<String> deletedNodes;
			if(deletedNodesFile != null && !deletedNodesFile.equals(""))
			{				
				deletedNodes= MapUtil.loadSet(deletedNodesFile, false);
				goldStandard.removeAll(deletedNodes);
				allGenes.removeAll(deletedNodes);
			}
			else
			{
				deletedNodes = new HashSet<String>();
			}
			
			// Exclude the known sources from the gold standard since they were
			// given as input to the algorithm
			goldStandard.removeAll(sources);
			
			// Add all gold standard members to the set of all genes in case
			// some of them were not in the interaction netowrk
			allGenes.addAll(goldStandard);
			
			HashSet<String> gsTfs = MapUtil.intersection(allTfs, goldStandard);
			HashSet<String> gsSignaling = MapUtil.subtract(goldStandard, gsTfs);
			HashSet<String> gsNoTargs = MapUtil.subtract(goldStandard, targets);
			
			PrintWriter writer = new PrintWriter(new FileWriter(out));

			if(deletedNodes.size() > 0)
			{
				writer.println("Deleted nodes (exluded from overlaps): " + deletedNodes.size());
			}
			
			writer.println("Sources (excluded from overlaps): " + sources.size());
			writer.println("Total TFs: " + allTfs.size());
			writer.println("Total genes: " + allGenes.size() + "\n");
			
			writer.println("Only consider targets as predicted TFs");
			writer.println("Targets: " + targets.size());
			writer.println("TFs in gold standard: " + gsTfs.size());
			HashSet<String> overlap = MapUtil.intersection(targets, gsTfs);
			writer.println("Overlap: " + overlap.size());
			double pval = StatUtil.overlapSignificance(overlap, targets, gsTfs, allTfs);
			writer.println("pval: " + pval);
//			for(String g : goldStandardTfs)
//			{
//				System.out.println(MapUtil.getStandard(g));
//			}
//			System.out.println();
			
			writer.println("Signaling predictions: " + internalNodes.size());
			writer.println("Gold standard members (minus targets): " + gsNoTargs.size());
			overlap = MapUtil.intersection(internalNodes, gsNoTargs);
			writer.println("Overlap: " + overlap.size());
			pval = StatUtil.overlapSignificance(overlap, internalNodes, gsNoTargs, allGenes);
			writer.println("pval: " + pval);
			
			
			writer.println("\nConsider all predictions together");
			HashSet<String> allPred = MapUtil.union(targets, internalNodes);
			HashSet<String> predTfs = MapUtil.intersection(allPred, allTfs);
			HashSet<String> predNonTfs = MapUtil.subtract(allPred, predTfs);
			
			writer.println("Predicted TFs: " + predTfs.size());
			writer.println("TFs in gold standard: " + gsTfs.size());
			overlap = MapUtil.intersection(predTfs, gsTfs);
			writer.println("Overlap: " + overlap.size());
			pval = StatUtil.overlapSignificance(overlap, predTfs, gsTfs, allTfs);
			writer.println("pval: " + pval);
			
			writer.println("Predicted non-TFs: " + predNonTfs.size());
			writer.println("Non-TFs in gold standard: " + gsSignaling.size());
			overlap = MapUtil.intersection(predNonTfs, gsSignaling);
			writer.println("Overlap: " + overlap.size());
			pval = StatUtil.overlapSignificance(overlap, predNonTfs, gsSignaling, allGenes);
			writer.println("pval: " + pval);
			
			writer.close();
		}
		catch(IOException e)
		{
			e.printStackTrace();
		}
	}
	
	/**
	 * Calculates p-values for the overlap of the HOG predictions and the
	 * gold standards.  Separate p-values are calculated for targets identified
	 * by ResponseNet and the other signaling proteins.
	 * @param edgeFile the table from the ResponseNet website that lists all edges
	 * between proteins and genes.  Used to determine the sources, targets, and internal
	 * nodes
	 * @param keggMembers a list of genes known to be involved in the HOG pathway
	 * based on the KEGG database.  Uses standard names.
	 * @param sciSigMembers a list of genes known to be involved in the HOG pathway
	 * based on the Science Signaling database.  Uses standard names.
	 * @param krantzMembers members of the HOG pathway as depicted in Figure 1 of Krantz2009
	 * @param tfList all proteins that are considered TFs (i.e. potential targets).  Standard names.
	 * @param physicalNetwork the PPI and protein-DNA network is used to obtain the list of all genes
	 * @param out
	 */
	public static void responseNetHogSignificance(String edgeFile,
			String keggMembers, String sciSigMembers,
			String krantzMembers, String tfList, String physicalNetwork,
			String out)
	{
		HashMap<String, HashSet<String>> predictions = loadResponseNet(edgeFile);

		// Use the other overloaded version of the method to perform the comparison
		// Assume no nodes were held out from the interaction network
		hogPredictionSignificance(predictions.get("sources"),
				predictions.get("targets"),
				predictions.get("internal"),
				null,
				keggMembers, sciSigMembers, krantzMembers, tfList, physicalNetwork,
				out);
	}
	
	/**
	 * Calculates p-values for the overlap of the PNM HOG predictions and the
	 * gold standards
	 * @param edgeFile the PNM output file that gives the directions of all edges
	 * that are on an active path.  Used to determine the sources and internal
	 * nodes
	 * @param potentialSources this file lists all proteins that were identified as
	 * knockout nodes and could potentially be sources
	 * @param keggMembers a list of genes known to be involved in the HOG pathway
	 * based on the KEGG database.  Uses standard names.
	 * @param sciSigMembers a list of genes known to be involved in the HOG pathway
	 * based on the Science Signaling database.  Uses standard names.
	 * @param krantzMembers members of the HOG pathway as depicted in Figure 1 of Krantz2009
	 * @param tfList all proteins that are considered TFs (i.e. potential targets).  Standard names.
	 * @param physicalNetwork the PPI and protein-DNA network is used to obtain the list of all genes
	 * @param out
	 */
	public static void pnmHogSignificance(String edgeFile, String potentialSources,
			String keggMembers, String sciSigMembers,
			String krantzMembers, String tfList, String physicalNetwork,
			String out)
	{
		HashMap<String, HashSet<String>> predictions = loadPNM(edgeFile, potentialSources);

		// Use the other overloaded version of the method to perform the comparison
		// Assume no nodes were held out from the interaction network
		hogPredictionSignificance(predictions.get("sources"),
				predictions.get("targets"),
				predictions.get("internal"),
				null,
				keggMembers, sciSigMembers, krantzMembers, tfList, physicalNetwork,
				out);
	}
	
	/**
	 * Calculates p-values for the overlap of the GeneReg HOG TF predictions and the
	 * gold standards.  Uses a threshold to only evaluate the most invovled regulators
	 * as determined by the number of genes they regulate.
	 * @param regFile the post-processed GeneReg output that tells how many genes
	 * each TF regulates
	 * @param threshold every TF that regulates at least this many genes is included
	 * in the GeneReg predictions
	 * @param keggMembers a list of genes known to be involved in the HOG pathway
	 * based on the KEGG database.  Uses standard names.
	 * @param sciSigMembers a list of genes known to be involved in the HOG pathway
	 * based on the Science Signaling database.  Uses standard names.
	 * @param krantzMembers members of the HOG pathway as depicted in Figure 1 of Krantz2009
	 * @param tfList all proteins that are considered TFs (i.e. potential targets).  Standard names.
	 * @param physicalNetwork the PPI and protein-DNA network is used to obtain the list of all genes
	 * @param out
	 */
	public static void geneRegSignificance(String regFile, int threshold,
			String keggMembers, String sciSigMembers,
			String krantzMembers, String tfList, String physicalNetwork,
			String out)
	{
		HashMap<String, HashSet<String>> predictions = loadGeneReg(regFile, threshold);

		// Use the other overloaded version of the method to perform the comparison
		// Assume no nodes were held out from the interaction network
		hogPredictionSignificance(predictions.get("sources"),
				predictions.get("targets"),
				predictions.get("internal"),
				null,
				keggMembers, sciSigMembers, krantzMembers, tfList, physicalNetwork,
				out);
	}
	
	/**
	 * Computes a p-value for the overlap between TOR predictions and all
	 * TOR screens and known pathway models
	 * @param sources the predefined sources.  Should be
	 * ORF names.
	 * @param targets a file containing proteins detected by the edge orientation + 
	 * DREM model.  Should be ORF names.
	 * @param nodeScoresDir the directory containing the node scores file
	 * @param nodescoresPrefix the first part of the node scores file from edge orientation.
	 * One file in the nodeScoresDir should match this prefix.
	 * @param nodeScoreThresh a percent threshold.  Nodes with a % top paths >= this
	 * threshold are deemed important.
	 * @param zamanMembers a list of genes known to be involved in the TOR pathway.
	 * Uses standard names
	 * @param sgdMembers a list of genes from SGD that are annotated with the GO term
	 * TOR signaling cascade
	 * @param xieScreen genes and their sensitivity to rapamycin.
	 * @param hillenmeyerScreen Hillenmeyer et al. rapamycin screen
	 * @param chanScreen Chan et al. rapamycin screen
	 * @param rapmycinPhen rapamycin phenotype annotation in SGD
	 * @param physicalNetwork the PPI and protein-DNA network is used to obtain the list of all genes
	 * @param out
	 */
	public static void torPredictionSignificance(String sources, String targets, 
			String nodeScoresDir, final String nodeScoresPrefix,
			double nodeScoreThresh, String zamanMembers, String sgdMembers,
			String xieScreen, String hillenmeyerScreen, String chanScreen,
			String rapamycinPhen, String physicalNetwork, String out)
	{
		try
		{
			// The physical network contributes to the set of all genes
			Set<String> edges1 = MapUtil.loadMap(physicalNetwork, 0, 1, false).keySet();
			Set<String> edges2 = MapUtil.loadMap(physicalNetwork, 2, 1, false).keySet();
			HashSet<String> allGenes = MapUtil.union(edges1, edges2);
			
			HashSet<String> srcList = MapUtil.loadSet(sources, false);
			allGenes.addAll(srcList);
			
			// Load the targets but ignore their target weights
			Set<String> targList = MapUtil.loadMap(targets, 0, 1, false).keySet();

			// Find the node scores
			File nsDir = new File(nodeScoresDir);
			File[] nsArray = nsDir.listFiles(new FilenameFilter()
			{
				public boolean accept(File dir, String name)
				{
					return (name.matches(nodeScoresPrefix + ".*"));
				}
			});
			if(nsArray.length != 1)
			{
				throw new IOException("Could not find unique node score file");
			}


			// Load the node score of all proteins, including sources and targets
			// Use the % of top paths for the node score
			HashMap<String, String> nodeScoreMap = MapUtil.loadMap(nsArray[0].getCanonicalPath(), 0, 5, true);

			// Determine which proteins have node scores >= the threshold
			HashSet<String> importantNodes = new HashSet<String>();
			for(String node : nodeScoreMap.keySet())
			{
				double score = Double.parseDouble(nodeScoreMap.get(node));
				if(score >= nodeScoreThresh)
				{
					importantNodes.add(node);
				}
			}

			// One set to contain the predicted TFs and signaling nodes
			// Specifically exclude the sources
			HashSet<String> predictions = new HashSet<String>();
			predictions.addAll(targList);
			predictions.addAll(importantNodes);
			predictions.removeAll(srcList);
			allGenes.addAll(predictions);

			// Load all of the known TOR evidence
			HashSet<String> known = new HashSet<String>();
			int upperBound = 0;

			// Col 0 is the ORF name
			// Col 2 is the rapamycin sensitivity
			HashMap<String, String> xieMap = MapUtil.loadMap(xieScreen, 0, 2, true);
			allGenes.addAll(xieMap.keySet());
			known.addAll(xieMap.keySet());
			System.out.println("Add Xie screen: " + known.size());
			upperBound += xieMap.keySet().size();
			
			// Col 0 is the ORF name
			// Col 2 is the lowest p-value across all replicates
			HashMap<String, HashSet<String>> hillenmeyerMap = MapUtil.loadMultiMap(hillenmeyerScreen, 0, 1, true);
			// Also want only the set of significant genes
			HashSet<String> hillenmeyerSig = new HashSet<String>();
			for(String gene : hillenmeyerMap.keySet())
			{
				// Some genes appear more than once, so add them if any time they appear
				// they have a low p-value
				HashSet<String> pvals = hillenmeyerMap.get(gene);
				for(String pval : pvals)
				{
					if(Double.parseDouble(pval) <= 0.0001)
					{
						hillenmeyerSig.add(gene);
					}
				}
			}
			allGenes.addAll(hillenmeyerMap.keySet());
			known.addAll(hillenmeyerSig);
			System.out.println("Add Hillenmeyer screen: " + known.size());
			upperBound += hillenmeyerSig.size();
			
			// Col 0 is the ORF name
			// Col 3 is the rapamycin sensitivity
			HashMap<String, String> chanMap = MapUtil.loadMap(chanScreen, 0, 3, true);
			allGenes.addAll(chanMap.keySet());
			known.addAll(chanMap.keySet());
			System.out.println("Add Chan screen: " + known.size());
			upperBound += chanMap.keySet().size();

			// Col 0 is the ORF name
			// Names may appear in more than one line
			Set<String> rapPhenSet = MapUtil.loadMap(rapamycinPhen, 0, 4, false).keySet();
			allGenes.addAll(rapPhenSet);
			known.addAll(rapPhenSet);
			System.out.println("Add SGD rapamycin phenotype: " + known.size());
			upperBound += rapPhenSet.size();
			
			// Exclude the sources
			known.removeAll(srcList);
			System.out.println("Remove sources: " + known.size());
			
			// First use only the screen data and ignore the Zaman and SGD models
			PrintWriter writer = new PrintWriter(new FileWriter(out));
			writer.println("Ignoring Zaman and SGD TOR models");
			
			HashSet<String> overlap = MapUtil.intersection(predictions, known);
			writer.println("Predictions: " + predictions.size());
			writer.println("TOR genes: " + known.size());
			writer.println("Overlap: " + overlap.size());
			writer.println("All genes: " + allGenes.size());
			
			double pval = StatUtil.overlapSignificance(overlap, predictions, known, allGenes);
			writer.println("p-value: " + pval);
			
			// Now retry, but include the Zaman and SGD models			
			HashSet<String> zamanStandard = MapUtil.loadSet(zamanMembers);
			// Translate from standard to ORF
			HashSet<String> zamanSet = new HashSet<String>();
			for(String gene : zamanStandard)
			{
				zamanSet.add(MapUtil.getOrf(gene.toUpperCase().trim()));
			}
			allGenes.addAll(zamanSet);
			known.addAll(zamanSet);
			System.out.println("Add Zaman model: " + known.size());
			upperBound += zamanSet.size();
			
			// Already has ORF names
			HashSet<String> sgdSet = MapUtil.loadSet(sgdMembers);
			allGenes.addAll(sgdSet);
			known.addAll(sgdSet);
			System.out.println("Add SGD model: " + known.size());
			upperBound += sgdSet.size();
			
			// Exclude the sources
			known.removeAll(srcList);
			System.out.println("Remove sources: " + known.size());
			
			// The number of genes that would be in the gold standard if there was
			// no agreement among any of the sources
			System.out.println("If no overlap, known TOR genes size: " + upperBound);
			
			writer.println("\n\nIncluding Zaman and SGD TOR models");
			
			overlap = MapUtil.intersection(predictions, known);
			writer.println("Predictions: " + predictions.size());
			writer.println("TOR genes: " + known.size());
			writer.println("Overlap: " + overlap.size());
			writer.println("All genes: " + allGenes.size());
			
			pval = StatUtil.overlapSignificance(overlap, predictions, known, allGenes);
			writer.println("p-value: " + pval);
			
			writer.close();
		}
		catch(IOException e)
		{
			e.printStackTrace();
		}

	}
	
	/**
	 * For each target, write how many of the top-ranked paths end at it,
	 * how many of all satisfied paths used in the run end at it, and various target
	 * scores.
	 * @param targFile the targets and their weights provided by DREM
	 * @param targScoresFile the network-derived target scores
	 * @param priorsFile the TF activity priors
	 * @param pathFile the satisifed paths
	 * @param topPaths
	 * @param outFile
	 */
	public static void targetFrequency(String targFile, String targScoresFile,
			String priorsFile, String pathFile, int topPaths, String outFile)
	{
		try
		{
			// Load the targets and their target weights in the network that were
			// derived from DREM
			HashMap<String, String> targWeights = MapUtil.loadMap(targFile, 0, 1, false);
			Set<String> targets = targWeights.keySet();
			
			// Load the network-derived target scores
			BufferedReader reader = new BufferedReader(new FileReader(targScoresFile));
			HashMap<String, String> targNetScore = new HashMap<String, String>();
			
			// Skip the header;
			reader.readLine();
			String line = reader.readLine();
			do
			{
				String[] parts = line.split("\t");
				targNetScore.put(parts[0], parts[3]);
				
				line = reader.readLine();
			} while(line != null && !line.trim().equals(""));
			reader.close();
			
			if(!MapUtil.sameElements(targets, targNetScore.keySet()))
			{
				throw new IllegalStateException("Target weights file and target network " +
						"scores file do not have the same targets");
			}
			
			// Load the TF activity priors
			HashMap<String, Double> targPriors = DREMInterface.readPrior(priorsFile);
			
			HashMap<String, Integer> allPathFreq = new HashMap<String, Integer>();
			HashMap<String, Integer> topPathFreq = new HashMap<String, Integer>();
			for(String target : targets)
			{
				allPathFreq.put(target, 0);
				topPathFreq.put(target, 0);
			}
			
			// Load the StringPaths and sort in descending order
			ArrayList<StringPath> paths = StringPath.loadFile(pathFile, true);
			Collections.sort(paths);
			Collections.reverse(paths);
			
			// Count how many paths end at each target
			for(int p = 0; p < paths.size(); p++)
			{
				// Find which target the path ends at
				String target = paths.get(p).getTarget();
				
				if(!targets.contains(target))
				{
					throw new IllegalStateException("Unrecognized target: " + target);
				}
				
				allPathFreq.put(target, allPathFreq.get(target) + 1);
				
				if(p < topPaths)
				{
					topPathFreq.put(target, topPathFreq.get(target) + 1);
				}
			}
			
			PrintWriter writer = new PrintWriter(new FileWriter(outFile));
			writer.println("Gene symbol\tGene id\tActivity prior\tTarget weight\t" +
					"Target network score\tTop " + 
					topPaths + " paths\tAll satisfied paths");
			for(String target : targets)
			{
				StringBuffer buf = new StringBuffer();
				String sym = Entrez.getSymbol(target);

				buf.append(sym).append("\t").append(target).append("\t");
				buf.append(targPriors.get(target)).append("\t");
				buf.append(targWeights.get(target)).append("\t");
				buf.append(targNetScore.get(target)).append("\t");
				buf.append(topPathFreq.get(target)).append("\t");
				buf.append(allPathFreq.get(target));
				writer.println(buf);
			}
			
			writer.close();
		}
		catch(IOException e)
		{
			e.printStackTrace();
		}
	}
	
	/**
	 * Determine which gold standard proteins lie on source-target paths
	 * @param srcFile
	 * @param targFile
	 * @param edgeFile
	 * @param pathLength
	 * @param goldStandardFile ORF names, no header
	 * @param outFile
	 */
	public static void goldStandardOnPaths(String srcFile, String targFile,
			String edgeFile, int pathLength, String goldStandardFile,
			String outFile)
	{
		try
		{
			// Load the graph and find paths
			Graph g = new Graph();
			DataLoader.readEdgesFromEda(g, edgeFile);
			DataLoader.readSources(g, srcFile);
			DataLoader.readTargets(g, targFile);
	
			EdgeOrientAlg orient = new EdgeOrientAlg(g);
			int numPaths = orient.findPaths(pathLength);
			
			// Construct the set of all vertices on the source-target paths
			// Consider all paths without regard to whether they could all
			// be simultaneously satisfied
			HashSet<String> allVerts = new HashSet<String>();
			for(Path p : orient.getPaths())
			{
				for(Vertex v : p.getVertices())
				{
					allVerts.add(v.getName());
				}
			}
			System.out.println(numPaths + " paths contain " + allVerts.size() + " unique vertices");
			
			// Load the gold standard and determine which gold standard members are on a path
			HashSet<String> goldStandard = MapUtil.loadSet(goldStandardFile, false);
			
			PrintWriter writer = new PrintWriter(new FileWriter(outFile));
			writer.println("Standard\tORF\tOn path");
			int onPath = 0;
			for(String protein : goldStandard)
			{
				writer.print(MapUtil.getStandard(protein) + "\t" + protein + "\t");
				if(allVerts.contains(protein))
				{
					onPath++;
					writer.println("Y");
				}
				else
				{
					writer.println("N");
				}
			}
			System.out.println(onPath + " of " + goldStandard.size() + " gold standard "
					+ "proteins are on a path");
			writer.close();
		}
		catch(IOException e)
		{
			e.printStackTrace();
		}
	}
	

	// This version does not take the number of top-ranked paths as a parameter,
	// it is fixed to 1000
	/**
	 * Compare different methods of predicting genetic screens.
	 * Consider sources and nodes with high node scores (which may include targets)
	 * @param edgeFile
	 * @param srcFile
	 * @param targFile
	 * @param nodePriorsFile
	 * @param defaultPrior
	 * @param nodeScoresFile
	 * @param nodeScoreThresh
	 * @param satPathsFile all of the satisifed paths from the SDREM output
	 * @param screenFile a list of screen hits and screen sources that is used to
	 * calculate AUC and compare ranking criteria
	 * @param secondScreenFile an optional second list of screen hits and screen sources
	 * that will be added to the full output file for single/pair hit predictions.  Set
	 * as null or "" to ignore.
	 * @param fullOutFile the predicted deletion scores for all genes of interest
	 * @param summaryOutFile a summary of how the ranking metrics perform at predicting
	 * the known screen hits
	 * @param pairsOutFile the predicted genetic interactions for all pairs of genes of
	 * interest.  Set as null or "" to ignore.
	 */
	public static void predictGeneticScreens(String edgeFile, String srcFile, String targFile,
			String nodePriorsFile, double defaultPrior,
			String nodeScoresFile,
			double nodeScoreThresh,
			String satPathsFile,
			String screenFile,
			String secondScreenFile,
			String fullOutFile,
			String summaryOutFile,
			String pairsOutFile)
	{
		try
		{
			boolean usePairs = (pairsOutFile != null && !pairsOutFile.equals(""));
			boolean twoScreens = (secondScreenFile != null && !secondScreenFile.equals(""));
			
			// Load the satisfied paths
			Vertex.setNodePriors(nodePriorsFile, defaultPrior);
			Graph graph = new Graph();
			DataLoader.readEdgesFromEda(graph, edgeFile);
			DataLoader.readSources(graph, srcFile);
			DataLoader.readTargets(graph, targFile);
			ArrayList<Path> satPaths = graph.loadStoredPaths(StringPath.loadFile(satPathsFile, true));
			System.out.println("Loaded " + satPaths.size() + " satisfied paths");

			// Store the 1000 top-ranked satisfied paths
			Collections.sort(satPaths, Path.getComparator("PathWeight"));
			Collections.reverse(satPaths); // Order largest to smallest
			int thresh = Math.min(1000, satPaths.size());
			
			ArrayList<Path> topPaths = new ArrayList<Path>(thresh);
			for(int t = 0; t < thresh; t++)
			{
				topPaths.add(satPaths.get(t));
			}
			
			HashSet<String> sources = MapUtil.loadSet(srcFile, false);
			HashSet<String> sourcesAndInternal = new HashSet<String>(sources);

			// Load the node score of all proteins, including sources and targets
			// Use the % of top paths for the node score
			HashMap<String, String> nodeScoreMapStr = MapUtil.loadMap(nodeScoresFile, 0, 5, true);
			HashMap<String, Double> nodeScoreMap = new HashMap<String, Double>(nodeScoreMapStr.size());
			
			// Determine which proteins have node scores >= the threshold
			for(String node : nodeScoreMapStr.keySet())
			{
				double score = Double.parseDouble(nodeScoreMapStr.get(node));
				nodeScoreMap.put(node, score);
				if(score >= nodeScoreThresh)
				{
					sourcesAndInternal.add(node);
				}
			}
			System.out.println("Scoring " + sourcesAndInternal.size() + " sources and high scoring nodes");
			
			// Also load the degree of all proteins
			HashMap<String, String> degreeMapStr = MapUtil.loadMap(nodeScoresFile, 0, 4, true);
			HashMap<String, Integer> degreeMap = new HashMap<String, Integer>(degreeMapStr.size());
			for(String node : degreeMapStr.keySet())
			{
				degreeMap.put(node, Integer.parseInt(degreeMapStr.get(node)));
			}
			
			// Convert the protein names to Vertex objects
			ArrayList<Vertex> nodesToScore = new ArrayList<Vertex>();
			for(String nodeId : sourcesAndInternal)
			{
				Vertex node = Vertex.findVertex(nodeId);
				if(node == null)
				{
					throw new IllegalStateException("Cannot find vertex " + nodeId + " in the graph");
				}
				nodesToScore.add(node);
			}
			// Sort the nodes so that the output file is identical across runs
			// with the same settings
			Collections.sort(nodesToScore);
			
			// Calculate the predicted genetic screen scores using the different metrics
			// Name the variables using all/top-ranked satisfied paths, source-target pairs
			// or target connectivity only (single source mode), unweighted/weighted,
			// and average/min (only applicable for weighted)
			// First the four unweighted metrics
			HashMap<String, Double> allStUnw = reachabilityScore(satPaths,
					nodesToScore,
					false,
					false,
					false);
			HashMap<String, Double> topStUnw = reachabilityScore(topPaths,
					nodesToScore,
					false,
					false,
					false);
			HashMap<String, Double> allTargUnw = reachabilityScore(satPaths,
					nodesToScore,
					true,
					false,
					false);
			HashMap<String, Double> topTargUnw = reachabilityScore(topPaths,
					nodesToScore,
					true,
					false,
					false);
			
			// Then the four weighted metrics that aggregate using the average
			HashMap<String, Double> allStAvg = reachabilityScore(satPaths,
					nodesToScore,
					false,
					true,
					true);
			HashMap<String, Double> topStAvg = reachabilityScore(topPaths,
					nodesToScore,
					false,
					true,
					true);
			HashMap<String, Double> allTargAvg = reachabilityScore(satPaths,
					nodesToScore,
					true,
					true,
					true);
			HashMap<String, Double> topTargAvg = reachabilityScore(topPaths,
					nodesToScore,
					true,
					true,
					true);
			
			// Then the four weighted metrics that aggregate using the min
			HashMap<String, Double> allStMin = reachabilityScore(satPaths,
					nodesToScore,
					false,
					true,
					false);
			HashMap<String, Double> topStMin = reachabilityScore(topPaths,
					nodesToScore,
					false,
					true,
					false);
			HashMap<String, Double> allTargMin = reachabilityScore(satPaths,
					nodesToScore,
					true,
					true,
					false);
			HashMap<String, Double> topTargMin = reachabilityScore(topPaths,
					nodesToScore,
					true,
					true,
					false);
			
			// Control for the node score and degree being used for tie-breaking
			HashMap<String, Double> nodeScoreOnly = new HashMap<String, Double>();
			for(Vertex node : nodesToScore)
			{
				// 1 means that 100% of the targets are still connected after the
				// gene is removed (therefore need to check node score to rank)
				nodeScoreOnly.put(node.toString(), 1.0);
			}
			
			// Create an upper bound scoring metric that ranks all screen hits highest
			HashMap<String, HashSet<String>> screenHitsFull = MapUtil.loadMultiMap(screenFile, 0, 1, false);
			Set<String> screenHits = screenHitsFull.keySet();
			// Store the number of screens each hit appears in
			HashMap<String, Integer> screenHitsCounts = new HashMap<String, Integer>();
			// Default to 0
			for(Vertex node : nodesToScore)
			{
				String nodeId = node.toString();
				screenHitsCounts.put(nodeId, 0);
			}
			for(String gene : screenHits)
			{
				screenHitsCounts.put(gene, screenHitsFull.get(gene).size());
			}
			// Use the number of screens a hit appears in to rank genes that appear
			// in more screens higher by giving them lower scores
			HashSet<String> hitsInNodesToScore = new HashSet<String>();
			HashMap<String, Double> upperBound = new HashMap<String, Double>();
			for(Vertex node : nodesToScore)
			{
				String nodeId = node.toString();
				if(screenHits.contains(nodeId))
				{
					// Presence in multiple screens gives a lower (higher ranking)
					// score but reserve the worst possible score (1) for genes
					// that are not screen hits
					upperBound.put(nodeId, 1/(1.0 + screenHitsCounts.get(nodeId)));
					hitsInNodesToScore.add(nodeId);
				}
				else
				{
					// 1 means that all targets are still connected after this
					// gene is removed, this is the worst possible score (least
					// severe "phenotype")
					upperBound.put(nodeId, 1.0);
				}
			}
			
			// Load the second set of screen hits if one was provided
			HashMap<String, HashSet<String>> screenHitsFull2 = null;
			HashMap<String, Integer> screenHitsCounts2 = null;
			if(twoScreens)
			{
				screenHitsFull2 = MapUtil.loadMultiMap(secondScreenFile, 0, 1, false);
				// Also store the number of screens each hit appears in
				screenHitsCounts2 = new HashMap<String, Integer>();
				// Default to 0
				for(Vertex node : nodesToScore)
				{
					String nodeId = node.toString();
					screenHitsCounts2.put(nodeId, 0);
				}
				for(String gene : screenHitsFull2.keySet())
				{
					screenHitsCounts2.put(gene, screenHitsFull2.get(gene).size());
				}
			}
						
			
			ArrayList<HashMap<String, Double>> scoresList = new ArrayList<HashMap<String, Double>>();
			scoresList.add(allStUnw);
			scoresList.add(topStUnw);
			scoresList.add(allTargUnw);
			scoresList.add(topTargUnw);
			scoresList.add(allStAvg);
			scoresList.add(topStAvg);
			scoresList.add(allTargAvg);
			scoresList.add(topTargAvg);
			scoresList.add(allStMin);
			scoresList.add(topStMin);
			scoresList.add(allTargMin);
			scoresList.add(topTargMin);
			scoresList.add(nodeScoreOnly);
			scoresList.add(upperBound);
			String[] metricNamesArr = {"allStUnw",
					"topStUnw",
					"allTargUnw",
					"topTargUnw",
					"allStAvg",
					"topStAvg",
					"allTargAvg",
					"topTargAvg",
					"allStMin",
					"topStMin",
					"allTargMin",
					"topTargMin",
					"nodeScore",
					"upperBound"};
			ArrayList<String> metricNames = new ArrayList<String>(Arrays.asList(metricNamesArr));
			// Rank the nodes for each metric
			ArrayList<HashMap<String, Integer>> ranksList = new ArrayList<HashMap<String, Integer>>();
			for(HashMap<String, Double> metric : scoresList)
			{
				ranksList.add(reachabilityRanks(metric, nodeScoreMap, degreeMap));
			}
			
			// Write the scores generated by each metric along with the degree and
			// original node score
			PrintWriter writer = new PrintWriter(new FileWriter(fullOutFile));
			
			writer.print("Gene\tGene id\tSrc\tNode score\tDegree\tScreen hits");
			if(twoScreens)
			{
				writer.print("\t2nd screen hits");
			}
			for(String name : metricNames)
			{
				writer.print("\t" + name);
			}
			for(String name : metricNames)
			{
				writer.print("\t" + name + "Rank");
			}
			writer.println();
			
			for(Vertex node : nodesToScore)
			{
				String nodeId = node.toString();
				StringBuffer buf = new StringBuffer();
				buf.append(Entrez.getSymbol(nodeId)).append("\t").append(nodeId).append("\t");
				
				if(sources.contains(nodeId))
				{
					buf.append("Y");
				}
				else
				{
					buf.append("N");
				}
				
				if(!nodeScoreMap.containsKey(nodeId) || !degreeMap.containsKey(nodeId))
				{
					throw new IllegalStateException("Cannot find node score " + 
							"and degree for " + nodeId);
				}
				buf.append("\t").append(nodeScoreMap.get(nodeId));
				buf.append("\t").append(degreeMap.get(nodeId)).append("\t");
				
				// Write the number of screen hits
				if(!screenHitsCounts.containsKey(nodeId))
				{
					throw new IllegalStateException("Cannot find number of screen hits for " + nodeId);
				}
				buf.append(screenHitsCounts.get(nodeId));
				if(twoScreens)
				{
					if(!screenHitsCounts2.containsKey(nodeId))
					{
						throw new IllegalStateException("Cannot find number of screen hits for " + nodeId);
					}
					buf.append("\t").append(screenHitsCounts2.get(nodeId));
				}
				
				// Finally, write the score for this node from each metric
				// and the rank from this metric
				for(HashMap<String, Double> scoreMetric : scoresList)
				{
					if(!scoreMetric.containsKey(nodeId))
					{
						throw new IllegalStateException("Cannot find all metric scores " + 
								"for " + nodeId);
					}
					buf.append("\t").append(scoreMetric.get(nodeId));
				}
				for(HashMap<String, Integer> rankMetric : ranksList)
				{
					if(!rankMetric.containsKey(nodeId))
					{
						throw new IllegalStateException("Cannot find all metric scores " + 
								"for " + nodeId);
					}
					buf.append("\t").append(rankMetric.get(nodeId));
				}
				writer.println(buf);
			}
			
			writer.close();
			
			// Write a summary of each metric that includes the number of screen
			// hits in the top predictions and the AUC
			writer = new PrintWriter(new FileWriter(summaryOutFile));
			
			writer.print("Name");
			for(String name : metricNames)
			{
				writer.print("\t" + name);
			}
			writer.println();
			
			// Print descriptions of each metric
			writer.print("Paths used");
			for(String name : metricNames)
			{
				writer.print("\t");
				if(name.contains("all"))
				{
					writer.print("all");
				}
				else if(name.contains("top"))
				{
					writer.print("top");
				}
				else
				{
					writer.print("N/A");
				}
			}
			writer.println();
			
			writer.print("Connectivity");
			for(String name : metricNames)
			{
				writer.print("\t");
				if(name.contains("St"))
				{
					writer.print("Source-target pairs");
				}
				else if(name.contains("Targ"))
				{
					writer.print("Targets");
				}
				else
				{
					writer.print("N/A");
				}
			}
			writer.println();
			
			writer.print("Scoring");
			for(String name : metricNames)
			{
				writer.print("\t");
				if(name.contains("Unw"))
				{
					writer.print("Unweighted");
				}
				else if(name.contains("Avg"))
				{
					writer.print("Weighted, avg");
				}
				else if(name.contains("Min"))
				{
					writer.print("Weighted, min");
				}
				else
				{
					writer.print("N/A");
				}
			}
			writer.println();
			
			// Calculate the AUC for each metric
			writer.print("AUC");
			for(HashMap<String, Integer> rankMetric : ranksList)
			{
				writer.print("\t" + areaUnderCurve(rankMetric, hitsInNodesToScore));
			}
			writer.println();
			
			// Calculate the number of screen hits appearing at various thresholds
			int[] rankThresholds = {10, 20, 50, 100, 150};
			for(int rt : rankThresholds)
			{
				writer.print("Screen hits in top " + rt);
				for(HashMap<String, Integer> rankMetric : ranksList)
				{
					// This version gives the number of screen hits
					writer.print("\t" + thresholdRankedPreds(rankMetric, hitsInNodesToScore, rt));
					writer.print(" [");
					// This version gives the number of screen hits and how many screens
					// those hits appear in
					writer.print(thresholdRankedPreds(rankMetric, hitsInNodesToScore, screenHitsCounts, rt));
					writer.print("]");
				}
				writer.println();
			}
			writer.close();
			
			// If an output files is provided for the gene pair scores,
			// calculate the scores and write them here
			if(usePairs)
			{
				// First the four unweighted metrics
				HashMap<String, Double> allStUnwPairs = reachabilityScore(satPaths,
						nodesToScore,
						false,
						false,
						false,
						true);
				HashMap<String, Double> topStUnwPairs = reachabilityScore(topPaths,
						nodesToScore,
						false,
						false,
						false,
						true);
				HashMap<String, Double> allTargUnwPairs = reachabilityScore(satPaths,
						nodesToScore,
						true,
						false,
						false,
						true);
				HashMap<String, Double> topTargUnwPairs = reachabilityScore(topPaths,
						nodesToScore,
						true,
						false,
						false,
						true);
				
				// Then the four weighted metrics that aggregate using the average
				HashMap<String, Double> allStAvgPairs = reachabilityScore(satPaths,
						nodesToScore,
						false,
						true,
						true,
						true);
				HashMap<String, Double> topStAvgPairs = reachabilityScore(topPaths,
						nodesToScore,
						false,
						true,
						true,
						true);
				HashMap<String, Double> allTargAvgPairs = reachabilityScore(satPaths,
						nodesToScore,
						true,
						true,
						true,
						true);
				HashMap<String, Double> topTargAvgPairs = reachabilityScore(topPaths,
						nodesToScore,
						true,
						true,
						true,
						true);
				
				// Then the four weighted metrics that aggregate using the min
				HashMap<String, Double> allStMinPairs = reachabilityScore(satPaths,
						nodesToScore,
						false,
						true,
						false,
						true);
				HashMap<String, Double> topStMinPairs = reachabilityScore(topPaths,
						nodesToScore,
						false,
						true,
						false,
						true);
				HashMap<String, Double> allTargMinPairs = reachabilityScore(satPaths,
						nodesToScore,
						true,
						true,
						false,
						true);
				HashMap<String, Double> topTargMinPairs = reachabilityScore(topPaths,
						nodesToScore,
						true,
						true,
						false,
						true);
				
				// Order the pairs so that the output is consistent across runs
				// with the same settings
				ArrayList<String> pairs = new ArrayList<String>(allStUnwPairs.keySet());
				Collections.sort(pairs);
				
				// The node score only and upper bound scores aren't meaningful
				// when calculating genetic interactions
				// Control for the node score and degree being used for tie-breaking
//				HashMap<String, Double> nodeScoreOnlyPairs = new HashMap<String, Double>();
//				for(String pair : pairs)
//				{
//					// 1 means that 100% of the targets are still connected after the
//					// gene is removed (therefore need to check node score to rank)
//					nodeScoreOnly.put(pair, 1.0);
//				}
				
				// Create a pairs score simlar to the upper bound metric for single genes
				// Give a score of 0 if both nodes in the pair are screen hits and a score
				// of 0.5 if only one node is a screen hit
//				HashMap<String, Double> upperBoundPairs = new HashMap<String, Double>();
//				for(String pair : pairs)
//				{
//					String[] nodeIds = pair.split("\t");
//					double score = 1;
//					if(screenHits.contains(nodeIds[0]))
//					{
//						score -= 0.5;
//					}
//					if(screenHits.contains(nodeIds[1]))
//					{
//						score -= 0.5;
//					}
//					upperBoundPairs.put(pair, score);
//				}
				
				
				// Compute the double deletion effects for each metric
				ArrayList<HashMap<String, Double>> pairsScoresList = new ArrayList<HashMap<String, Double>>();
				pairsScoresList.add(allStUnwPairs);
				pairsScoresList.add(topStUnwPairs);
				pairsScoresList.add(allTargUnwPairs);
				pairsScoresList.add(topTargUnwPairs);
				pairsScoresList.add(allStAvgPairs);
				pairsScoresList.add(topStAvgPairs);
				pairsScoresList.add(allTargAvgPairs);
				pairsScoresList.add(topTargAvgPairs);
				pairsScoresList.add(allStMinPairs);
				pairsScoresList.add(topStMinPairs);
				pairsScoresList.add(allTargMinPairs);
				pairsScoresList.add(topTargMinPairs);
//				pairsScoresList.add(nodeScoreOnlyPairs);
//				pairsScoresList.add(upperBoundPairs);
				
				// Remove the node score metric and upper bound metric
				scoresList.remove(scoresList.size()-1);
				scoresList.remove(scoresList.size()-1);
				ranksList.remove(ranksList.size()-1);
				ranksList.remove(ranksList.size()-1);
				metricNames.remove(metricNames.size()-1);
				metricNames.remove(metricNames.size()-1);
				if(scoresList.size() != pairsScoresList.size() ||
						ranksList.size() != pairsScoresList.size() ||
						metricNames.size() != pairsScoresList.size())
				{
					throw new IllegalStateException("Must use same metrics for single " +
							"and double deletions");
				}
				
				// Compute the predicted genetic interactions for each metric
				ArrayList<HashMap<String, Double>> genIntsList = new ArrayList<HashMap<String, Double>>();
				for(int s = 0; s < scoresList.size(); s++)
				{
					genIntsList.add(geneticInts(scoresList.get(s), pairsScoresList.get(s)));
				}
				
				// Rank the predicted genetic interactions for each metric
				ArrayList<HashMap<String, Integer>> pairsRanksList = new ArrayList<HashMap<String, Integer>>();
				for(HashMap<String, Double> metric : genIntsList)
				{
					pairsRanksList.add(reachabilityRanks(metric, nodeScoreMap, degreeMap, true));
				}
				
				writer = new PrintWriter(new FileWriter(pairsOutFile));

				writer.print("Gene A\tGene B\tGene id A\tGene id B\tSrc A\tSrc B\t" + 
						"Node score A\tNode score B\tDegree A\tDegree B\t" +
						"Screen hits A\tScreen hits B");
				if(twoScreens)
				{
					writer.print("\t2nd screen hits A\t2nd screen hits B");
				}
				for(String name : metricNames)
				{
					writer.print("\t" + name + "_int");
					writer.print("\t" + name + "_obs_AB");
					writer.print("\t" + name + "_exp_AB");
					writer.print("\t" + name + "_obs_A");
					writer.print("\t" + name + "_obs_B");
				}
				for(String name : metricNames)
				{
					writer.print("\t" + name + "Rank");
				}
				writer.println();
				
				for(String pair : pairs)
				{
					String[] nodeIds = pair.split("\t");
					String a = nodeIds[0];
					String b = nodeIds[1];

					StringBuffer buf = new StringBuffer();
					buf.append(Entrez.getSymbol(a)).append("\t").append(Entrez.getSymbol(b)).append("\t");
					buf.append(a).append("\t").append(b).append("\t");
					
					if(sources.contains(a))
					{
						buf.append("Y\t");
					}
					else
					{
						buf.append("N\t");
					}
					if(sources.contains(b))
					{
						buf.append("Y");
					}
					else
					{
						buf.append("N");
					}
					
					if(!nodeScoreMap.containsKey(a) || !degreeMap.containsKey(a))
					{
						throw new IllegalStateException("Cannot find node score " + 
								"and degree for " + a);
					}
					if(!nodeScoreMap.containsKey(b) || !degreeMap.containsKey(b))
					{
						throw new IllegalStateException("Cannot find node score " + 
								"and degree for " + b);
					}
					buf.append("\t").append(nodeScoreMap.get(a)).append("\t").append(nodeScoreMap.get(b));
					buf.append("\t").append(degreeMap.get(a)).append("\t").append(degreeMap.get(b));
					
					// Write the number of screen hits
					if(!screenHitsCounts.containsKey(a))
					{
						throw new IllegalStateException("Cannot find number of screen hits for " + a);
					}
					buf.append("\t").append(screenHitsCounts.get(a));
					if(!screenHitsCounts.containsKey(b))
					{
						throw new IllegalStateException("Cannot find number of screen hits for " + b);
					}
					buf.append("\t").append(screenHitsCounts.get(b));
					if(twoScreens)
					{
						if(!screenHitsCounts2.containsKey(a))
						{
							throw new IllegalStateException("Cannot find number of screen hits for " + a);
						}
						buf.append("\t").append(screenHitsCounts2.get(a));
						if(!screenHitsCounts2.containsKey(b))
						{
							throw new IllegalStateException("Cannot find number of screen hits for " + b);
						}
						buf.append("\t").append(screenHitsCounts2.get(b));
					}
					
					// Finally, write the scores for this pair from each metric
					// and the rank from this metric
					for(int s = 0; s < scoresList.size(); s++)
					{
						// The genetic interaction score
						if(!genIntsList.get(s).containsKey(pair))
						{
							throw new IllegalStateException("Cannot find all genetic " + 
									"interaction scores for " + pair);
						}
						buf.append("\t").append(genIntsList.get(s).get(pair));
						
						// The observed double deletion score
						if(!pairsScoresList.get(s).containsKey(pair))
						{
							throw new IllegalStateException("Cannot find all pairwise metric " + 
									"scores for " + pair);
						}
						buf.append("\t").append(pairsScoresList.get(s).get(pair));
						
						// The expected double deletion score and observed
						// single deletion scores				
						if(!scoresList.get(s).containsKey(a) || !scoresList.get(s).containsKey(b))
						{
							throw new IllegalStateException("Cannot find all metric scores " + 
									"for " + pair);
						}
						buf.append("\t").append(scoresList.get(s).get(a)*scoresList.get(s).get(b));
						buf.append("\t").append(scoresList.get(s).get(a));
						buf.append("\t").append(scoresList.get(s).get(b));
					}
					for(HashMap<String, Integer> rankMetric : pairsRanksList)
					{
						if(!rankMetric.containsKey(pair))
						{
							throw new IllegalStateException("Cannot find all metric ranks " + 
									"for " + pair);
						}
						buf.append("\t").append(rankMetric.get(pair));
					}
					writer.println(buf);
				}

				writer.close();
			}
		}
		catch(IOException e)
		{
			e.printStackTrace();
		}
	}
	
	
	/**
	 * Compare different methods of predicting genetic screens.
	 * Consider sources and nodes with high node scores (which may include targets).
	 * Copied from predictGeneticScreens and stripped to not take in known
	 * screen hits or calculate AUC.
	 * @param edgeFile
	 * @param srcFile
	 * @param targFile
	 * @param nodePriorsDir
	 * @param nodeScoresPrefix
	 * @param defaultPrior
	 * @param nodeScoresFile
	 * @param nodeScoreThresh
	 * @param numTopPaths
	 * @param satPathsFile all of the satisifed paths from the SDREM output
	 * @param fullOutFile the predicted deletion scores for all genes of interest
	 * @param pairsOutFile the predicted genetic interactions for all pairs of genes of
	 * interest.  Set as null or "" to ignore.
	 */
	public static void predictGeneticScreens(String edgeFile, String srcFile, String targFile,
			String nodePriorsFile, double defaultPrior,
			String nodeScoresDir,
			final String nodeScoresPrefix,
			double nodeScoreThresh,
			int numTopPaths,
			String satPathsFile,
			String fullOutFile,
			String pairsOutFile)
	{
		try
		{
			boolean usePairs = (pairsOutFile != null && !pairsOutFile.equals(""));
			
			// Load the satisfied paths
			Vertex.setNodePriors(nodePriorsFile, defaultPrior);
			Graph graph = new Graph();
			DataLoader.readEdgesFromEda(graph, edgeFile);
			DataLoader.readSources(graph, srcFile);
			DataLoader.readTargets(graph, targFile);
			ArrayList<Path> satPaths = graph.loadStoredPaths(StringPath.loadFile(satPathsFile, true));
			System.out.println("Loaded " + satPaths.size() + " satisfied paths to predict screen hits");

			// Store the top-ranked satisfied paths
			Collections.sort(satPaths, Path.getComparator("PathWeight"));
			Collections.reverse(satPaths); // Order largest to smallest
			// If numTopPaths is set to -1, calculate the number of top paths
			// dynamically by taking 5 * the number of targets
			if (numTopPaths < 0)
			{
				numTopPaths = 5 * graph.getTargets().size();
			}
			int thresh = Math.min(numTopPaths, satPaths.size());
			
			ArrayList<Path> topPaths = new ArrayList<Path>(thresh);
			for(int t = 0; t < thresh; t++)
			{
				topPaths.add(satPaths.get(t));
			}
			
			HashSet<String> sources = MapUtil.loadSet(srcFile, false);
			HashSet<String> sourcesAndInternal = new HashSet<String>(sources);

			// Find the node scores file
			File nsDir = new File(nodeScoresDir);
			File[] nsArray = nsDir.listFiles(new FilenameFilter()
			{
				public boolean accept(File dir, String name)
				{
					return (name.matches(nodeScoresPrefix + ".*"));
				}
			});
			if(nsArray.length != 1)
			{
				throw new IOException("Could not find unique node score file");
			}
			String nodeScoresFile = nsArray[0].getCanonicalPath();
			
			// Load the node score of all proteins, including sources and targets
			// Use the % of top paths for the node score
			HashMap<String, String> nodeScoreMapStr = MapUtil.loadMap(nodeScoresFile, 0, 5, true);
			HashMap<String, Double> nodeScoreMap = new HashMap<String, Double>(nodeScoreMapStr.size());
			
			// Determine which proteins have node scores >= the threshold
			for(String node : nodeScoreMapStr.keySet())
			{
				double score = Double.parseDouble(nodeScoreMapStr.get(node));
				nodeScoreMap.put(node, score);
				if(score >= nodeScoreThresh)
				{
					sourcesAndInternal.add(node);
				}
			}
			System.out.println("Scoring " + sourcesAndInternal.size() + " sources and high-scoring nodes");
			
			// Also load the degree of all proteins
			HashMap<String, String> degreeMapStr = MapUtil.loadMap(nodeScoresFile, 0, 4, true);
			HashMap<String, Integer> degreeMap = new HashMap<String, Integer>(degreeMapStr.size());
			for(String node : degreeMapStr.keySet())
			{
				degreeMap.put(node, Integer.parseInt(degreeMapStr.get(node)));
			}
			
			// Convert the protein names to Vertex objects
			ArrayList<Vertex> nodesToScore = new ArrayList<Vertex>();
			for(String nodeId : sourcesAndInternal)
			{
				Vertex node = Vertex.findVertex(nodeId);
				if(node == null)
				{
					throw new IllegalStateException("Cannot find vertex " + nodeId + " in the graph");
				}
				nodesToScore.add(node);
			}
			// Sort the nodes so that the output file is identical across runs
			// with the same settings
			Collections.sort(nodesToScore);
			
			// Calculate the predicted genetic screen scores using the different metrics
			// Name the variables using all/top-ranked satisfied paths, source-target pairs
			// or target connectivity only (single source mode), unweighted/weighted,
			// and average/min (only applicable for weighted)
			// First the four unweighted metrics
			HashMap<String, Double> allStUnw = reachabilityScore(satPaths,
					nodesToScore,
					false,
					false,
					false);
			HashMap<String, Double> topStUnw = reachabilityScore(topPaths,
					nodesToScore,
					false,
					false,
					false);
			HashMap<String, Double> allTargUnw = reachabilityScore(satPaths,
					nodesToScore,
					true,
					false,
					false);
			HashMap<String, Double> topTargUnw = reachabilityScore(topPaths,
					nodesToScore,
					true,
					false,
					false);
			
			// Then the four weighted metrics that aggregate using the average
			HashMap<String, Double> allStAvg = reachabilityScore(satPaths,
					nodesToScore,
					false,
					true,
					true);
			HashMap<String, Double> topStAvg = reachabilityScore(topPaths,
					nodesToScore,
					false,
					true,
					true);
			HashMap<String, Double> allTargAvg = reachabilityScore(satPaths,
					nodesToScore,
					true,
					true,
					true);
			HashMap<String, Double> topTargAvg = reachabilityScore(topPaths,
					nodesToScore,
					true,
					true,
					true);
			
			// Then the four weighted metrics that aggregate using the min
			HashMap<String, Double> allStMin = reachabilityScore(satPaths,
					nodesToScore,
					false,
					true,
					false);
			HashMap<String, Double> topStMin = reachabilityScore(topPaths,
					nodesToScore,
					false,
					true,
					false);
			HashMap<String, Double> allTargMin = reachabilityScore(satPaths,
					nodesToScore,
					true,
					true,
					false);
			HashMap<String, Double> topTargMin = reachabilityScore(topPaths,
					nodesToScore,
					true,
					true,
					false);
			
			// Control for the node score and degree being used for tie-breaking
			HashMap<String, Double> nodeScoreOnly = new HashMap<String, Double>();
			for(Vertex node : nodesToScore)
			{
				// 1 means that 100% of the targets are still connected after the
				// gene is removed (therefore need to check node score to rank)
				nodeScoreOnly.put(node.toString(), 1.0);
			}

			ArrayList<HashMap<String, Double>> scoresList = new ArrayList<HashMap<String, Double>>();
			scoresList.add(allStUnw);
			scoresList.add(topStUnw);
			scoresList.add(allTargUnw);
			scoresList.add(topTargUnw);
			scoresList.add(allStAvg);
			scoresList.add(topStAvg);
			scoresList.add(allTargAvg);
			scoresList.add(topTargAvg);
			scoresList.add(allStMin);
			scoresList.add(topStMin);
			scoresList.add(allTargMin);
			scoresList.add(topTargMin);
			scoresList.add(nodeScoreOnly);
			String[] metricNamesArr = {"allStUnw",
					"topStUnw",
					"allTargUnw",
					"topTargUnw",
					"allStAvg",
					"topStAvg",
					"allTargAvg",
					"topTargAvg",
					"allStMin",
					"topStMin",
					"allTargMin",
					"topTargMin",
					"nodeScore",};
			ArrayList<String> metricNames = new ArrayList<String>(Arrays.asList(metricNamesArr));
			// Rank the nodes for each metric
			ArrayList<HashMap<String, Integer>> ranksList = new ArrayList<HashMap<String, Integer>>();
			for(HashMap<String, Double> metric : scoresList)
			{
				ranksList.add(reachabilityRanks(metric, nodeScoreMap, degreeMap));
			}
			
			// Write the scores generated by each metric along with the degree and
			// original node score
			PrintWriter writer = new PrintWriter(new FileWriter(fullOutFile));
			
			writer.print("Gene\tSource\tNode score\tDegree");
			for(String name : metricNames)
			{
				writer.print("\t" + name);
			}
			for(String name : metricNames)
			{
				writer.print("\t" + name + "Rank");
			}
			writer.println();
			
			for(Vertex node : nodesToScore)
			{
				String nodeId = node.toString();
				StringBuffer buf = new StringBuffer();
				buf.append(nodeId).append("\t");
				
				if(sources.contains(nodeId))
				{
					buf.append("Y");
				}
				else
				{
					buf.append("N");
				}
				
				if(!nodeScoreMap.containsKey(nodeId) || !degreeMap.containsKey(nodeId))
				{
					throw new IllegalStateException("Cannot find node score " + 
							"and degree for " + nodeId);
				}
				buf.append("\t").append(nodeScoreMap.get(nodeId));
				buf.append("\t").append(degreeMap.get(nodeId));
				
				// Finally, write the score for this node from each metric
				// and the rank from this metric
				for(HashMap<String, Double> scoreMetric : scoresList)
				{
					if(!scoreMetric.containsKey(nodeId))
					{
						throw new IllegalStateException("Cannot find all metric scores " + 
								"for " + nodeId);
					}
					buf.append("\t").append(scoreMetric.get(nodeId));
				}
				for(HashMap<String, Integer> rankMetric : ranksList)
				{
					if(!rankMetric.containsKey(nodeId))
					{
						throw new IllegalStateException("Cannot find all metric scores " + 
								"for " + nodeId);
					}
					buf.append("\t").append(rankMetric.get(nodeId));
				}
				writer.println(buf);
			}
			
			writer.close();
			
			// If an output files is provided for the gene pair scores,
			// calculate the scores and write them here
			if(usePairs)
			{
				// First the four unweighted metrics
				HashMap<String, Double> allStUnwPairs = reachabilityScore(satPaths,
						nodesToScore,
						false,
						false,
						false,
						true);
				HashMap<String, Double> topStUnwPairs = reachabilityScore(topPaths,
						nodesToScore,
						false,
						false,
						false,
						true);
				HashMap<String, Double> allTargUnwPairs = reachabilityScore(satPaths,
						nodesToScore,
						true,
						false,
						false,
						true);
				HashMap<String, Double> topTargUnwPairs = reachabilityScore(topPaths,
						nodesToScore,
						true,
						false,
						false,
						true);
				
				// Then the four weighted metrics that aggregate using the average
				HashMap<String, Double> allStAvgPairs = reachabilityScore(satPaths,
						nodesToScore,
						false,
						true,
						true,
						true);
				HashMap<String, Double> topStAvgPairs = reachabilityScore(topPaths,
						nodesToScore,
						false,
						true,
						true,
						true);
				HashMap<String, Double> allTargAvgPairs = reachabilityScore(satPaths,
						nodesToScore,
						true,
						true,
						true,
						true);
				HashMap<String, Double> topTargAvgPairs = reachabilityScore(topPaths,
						nodesToScore,
						true,
						true,
						true,
						true);
				
				// Then the four weighted metrics that aggregate using the min
				HashMap<String, Double> allStMinPairs = reachabilityScore(satPaths,
						nodesToScore,
						false,
						true,
						false,
						true);
				HashMap<String, Double> topStMinPairs = reachabilityScore(topPaths,
						nodesToScore,
						false,
						true,
						false,
						true);
				HashMap<String, Double> allTargMinPairs = reachabilityScore(satPaths,
						nodesToScore,
						true,
						true,
						false,
						true);
				HashMap<String, Double> topTargMinPairs = reachabilityScore(topPaths,
						nodesToScore,
						true,
						true,
						false,
						true);
				
				// Order the pairs so that the output is consistent across runs
				// with the same settings
				ArrayList<String> pairs = new ArrayList<String>(allStUnwPairs.keySet());
				Collections.sort(pairs);	
				
				// Compute the double deletion effects for each metric
				ArrayList<HashMap<String, Double>> pairsScoresList = new ArrayList<HashMap<String, Double>>();
				pairsScoresList.add(allStUnwPairs);
				pairsScoresList.add(topStUnwPairs);
				pairsScoresList.add(allTargUnwPairs);
				pairsScoresList.add(topTargUnwPairs);
				pairsScoresList.add(allStAvgPairs);
				pairsScoresList.add(topStAvgPairs);
				pairsScoresList.add(allTargAvgPairs);
				pairsScoresList.add(topTargAvgPairs);
				pairsScoresList.add(allStMinPairs);
				pairsScoresList.add(topStMinPairs);
				pairsScoresList.add(allTargMinPairs);
				pairsScoresList.add(topTargMinPairs);
				
				// Remove the node score metric
				scoresList.remove(scoresList.size()-1);
				ranksList.remove(ranksList.size()-1);
				metricNames.remove(metricNames.size()-1);
				if(scoresList.size() != pairsScoresList.size() ||
						ranksList.size() != pairsScoresList.size() ||
						metricNames.size() != pairsScoresList.size())
				{
					throw new IllegalStateException("Must use same metrics for single " +
							"and double deletions");
				}
				
				// Compute the predicted genetic interactions for each metric
				ArrayList<HashMap<String, Double>> genIntsList = new ArrayList<HashMap<String, Double>>();
				for(int s = 0; s < scoresList.size(); s++)
				{
					genIntsList.add(geneticInts(scoresList.get(s), pairsScoresList.get(s)));
				}
				
				// Rank the predicted genetic interactions for each metric
				ArrayList<HashMap<String, Integer>> pairsRanksList = new ArrayList<HashMap<String, Integer>>();
				for(HashMap<String, Double> metric : genIntsList)
				{
					pairsRanksList.add(reachabilityRanks(metric, nodeScoreMap, degreeMap, true));
				}
				
				writer = new PrintWriter(new FileWriter(pairsOutFile));

				writer.print("Gene A\tGene B\tSource A\tSource B\t" + 
						"Node score A\tNode score B\tDegree A\tDegree B");
				for(String name : metricNames)
				{
					writer.print("\t" + name + "_int");
					writer.print("\t" + name + "_obs_AB");
					writer.print("\t" + name + "_exp_AB");
					writer.print("\t" + name + "_obs_A");
					writer.print("\t" + name + "_obs_B");
				}
				for(String name : metricNames)
				{
					writer.print("\t" + name + "Rank");
				}
				writer.println();
				
				for(String pair : pairs)
				{
					String[] nodeIds = pair.split("\t");
					String a = nodeIds[0];
					String b = nodeIds[1];

					StringBuffer buf = new StringBuffer();
					buf.append(a).append("\t").append(b).append("\t");
					
					if(sources.contains(a))
					{
						buf.append("Y\t");
					}
					else
					{
						buf.append("N\t");
					}
					if(sources.contains(b))
					{
						buf.append("Y");
					}
					else
					{
						buf.append("N");
					}
					
					if(!nodeScoreMap.containsKey(a) || !degreeMap.containsKey(a))
					{
						throw new IllegalStateException("Cannot find node score " + 
								"and degree for " + a);
					}
					if(!nodeScoreMap.containsKey(b) || !degreeMap.containsKey(b))
					{
						throw new IllegalStateException("Cannot find node score " + 
								"and degree for " + b);
					}
					buf.append("\t").append(nodeScoreMap.get(a)).append("\t").append(nodeScoreMap.get(b));
					buf.append("\t").append(degreeMap.get(a)).append("\t").append(degreeMap.get(b));
					
					// Finally, write the scores for this pair from each metric
					// and the rank from this metric
					for(int s = 0; s < scoresList.size(); s++)
					{
						// The genetic interaction score
						if(!genIntsList.get(s).containsKey(pair))
						{
							throw new IllegalStateException("Cannot find all genetic " + 
									"interaction scores for " + pair);
						}
						buf.append("\t").append(genIntsList.get(s).get(pair));
						
						// The observed double deletion score
						if(!pairsScoresList.get(s).containsKey(pair))
						{
							throw new IllegalStateException("Cannot find all pairwise metric " + 
									"scores for " + pair);
						}
						buf.append("\t").append(pairsScoresList.get(s).get(pair));
						
						// The expected double deletion score and observed
						// single deletion scores				
						if(!scoresList.get(s).containsKey(a) || !scoresList.get(s).containsKey(b))
						{
							throw new IllegalStateException("Cannot find all metric scores " + 
									"for " + pair);
						}
						buf.append("\t").append(scoresList.get(s).get(a)*scoresList.get(s).get(b));
						buf.append("\t").append(scoresList.get(s).get(a));
						buf.append("\t").append(scoresList.get(s).get(b));
					}
					for(HashMap<String, Integer> rankMetric : pairsRanksList)
					{
						if(!rankMetric.containsKey(pair))
						{
							throw new IllegalStateException("Cannot find all metric ranks " + 
									"for " + pair);
						}
						buf.append("\t").append(rankMetric.get(pair));
					}
					writer.println(buf);
				}

				writer.close();
			}
		}
		catch(IOException e)
		{
			e.printStackTrace();
		}
	}
	
	
	/**
	 * Calculate node scores based on the predicted impact to downstream targets
	 * if the node is removed.  Scores range from 1 (all targets still connected)
	 * to 0 (no "flow" to targets after deletion).  Does not calculate
	 * scores for pairs of nodes, only single node deletions.
	 * @param paths the set of paths to consider when assessing connectivity
	 * @param nodes the nodes of interest that should be deleted in silico and scored
	 * @param singleSource if true, only require that a path from a single soucre to
	 * a target must be present to activate the target.  Otherwise, treat source-target
	 * pairs independently (a target is activated distinctly by distinct sources).
	 * @param weighted if true, calculate weighted version of the scores instead of
	 * requiring complete disconnection
	 * @param takeAverage if true, average the reachability scores across targets or
	 * source-target pairs for each gene.  Otherwise take the min
	 * @return
	 */
	public static HashMap<String, Double> reachabilityScore(ArrayList<Path> paths,
			Collection<Vertex> nodes,
			boolean singleSource,
			boolean weighted,
			boolean takeAverage)
	{
		return reachabilityScore(paths, nodes, singleSource, weighted, takeAverage, false);
	}
	
	/**
	 * Calculate node scores based on the predicted impact to downstream targets
	 * if the node is removed.  Scores range from 1 (all targets still connected)
	 * to 0 (no "flow" to targets after deletion).
	 * @param paths the set of paths to consider when assessing connectivity
	 * @param nodes the nodes of interest that should be deleted in silico and scored
	 * @param singleSource if true, only require that a path from a single soucre to
	 * a target must be present to activate the target.  Otherwise, treat source-target
	 * pairs independently (a target is activated distinctly by distinct sources).
	 * @param weighted if true, calculate weighted version of the scores instead of
	 * requiring complete disconnection
	 * @param takeAverage if true, average the reachability scores across targets or
	 * source-target pairs for each gene.  Otherwise take the min
	 * @param usePairs if true, delete pairs of nodes at the same time and calculate
	 * the reachability score for double deletion
	 * @return
	 */
	public static HashMap<String, Double> reachabilityScore(ArrayList<Path> paths,
			Collection<Vertex> nodes,
			boolean singleSource,
			boolean weighted,
			boolean takeAverage,
			boolean usePairs)
	{
		// Build an index from targets or source-target pairs to paths connecting them
		HashMap<String, ArrayList<Path>> targToPaths = new HashMap<String, ArrayList<Path>>();
		for(Path curPath : paths)
		{
			String key;
			if(singleSource)
			{
				// Target only needs to be connected to a single source
				// to be "activated" so no need to consider s-t pairs
				key = curPath.getTarget().toString();
			}
			else
			{
				// Group paths by source-target pairs
				key = curPath.stKey();
			}
			
			ArrayList<Path> curList;
			if(targToPaths.containsKey(key))
			{
				curList = targToPaths.get(key);
			}
			else
			{
				curList = new ArrayList<Path>();
			}
			curList.add(curPath);
			targToPaths.put(key, curList);
		}
		
		// Scores for each node or node-node pair
		HashMap<String, Double> scoreMap = new HashMap<String, Double>();
		
		// Need to order the nodes in order to create the pairs
		ArrayList<Vertex> nodeList = new ArrayList<Vertex>(nodes);
		for(int n1 = 0; n1 < nodeList.size(); n1++)
		{
			// The node or pair of nodes that are being deleted from the network
			ArrayList<Vertex> delNodes = new ArrayList<Vertex>();
			Vertex vert1 = nodeList.get(n1);
			delNodes.add(vert1);

			// Loop through the second node to be paired, but if pairs are
			// not being used only calculate the score once and do not enter
			// the loop block again
			int n2 = n1 + 1;
			String key = vert1.toString(); // Default to single node scoreMap key
			while((usePairs && n2 < nodeList.size()) || (!usePairs && !scoreMap.containsKey(key)))
			{
				if(usePairs)
				{
					// Remove the second node in the pair from the previous iteration
					delNodes.clear();
					delNodes.add(vert1);
					
					Vertex vert2 = nodeList.get(n2++);
					delNodes.add(vert2);
					key = StrOps.sortAndCombine(vert1.toString(), vert2.toString());
					if(delNodes.size() > 2)
					{
						throw new IllegalStateException("Error, removing more than 2 nodes");
					}
				}

				// Weighted version must iterate through all paths whereas the
				// unweighted version can shortcut as soon as a path to the target
				// is found
				if(weighted)
				{
					ArrayList<Double> scores = new ArrayList<Double>(targToPaths.size());

					// See what fraction of the "flow" (the sum of path weights)
					// still remains when the node is removed.  Greater disruption
					// to a target or source-target pair leads to a lower (more severe) score.
					for(String targOrPair : targToPaths.keySet())
					{
						double disconnected = 0, all = 0;
						for(Path curPath : targToPaths.get(targOrPair))
						{
							all +=  curPath.getMetricValue("PathWeight");
							if(curPath.containsAnyVertex(delNodes))
							{
								disconnected += curPath.getMetricValue("PathWeight");
							}
						}
						scores.add(1 - disconnected / all);
					}

					// Aggregate the scores across all targets or s-t pairs
					if(takeAverage)
					{
						scoreMap.put(key, ArrayUtil.average(scores));
					}
					else
					{
						// Otherwise take the min score
						scoreMap.put(key, ArrayUtil.min(scores));
					}
				}
				else
				{
					// See if removing this node or nodes completely disconnects any targets 
					// or source-target pairs (depending on how paths were grouped above)
					int disconnectedPairs = 0;
					for(String targOrPair : targToPaths.keySet())
					{
						boolean foundConnection = false;
						Iterator<Path> curPaths = targToPaths.get(targOrPair).iterator();
						while(curPaths.hasNext() && !foundConnection)
						{
							// The pair will be still connected once we find a single
							// path that does not contain the deleted node(s)
							foundConnection = !curPaths.next().containsAnyVertex(delNodes);
						}

						if(!foundConnection)
						{
							disconnectedPairs++;
						}
					}

					// How many of the pairs were disconnected
					double ratio = 1 - (double)disconnectedPairs / targToPaths.size();
					scoreMap.put(key, ratio);
				}
			}  // End inner loop over vertices or calculation of single node score
		} // End outer loop over vertices
		
		return scoreMap;
	}
	
	/**
	 * Sort genes by their reachability scores using node score, degree, and node id
	 * for tie-breaking.  Uses ReachabilityComp to sort.  In the returned ranks,
	 * rank 1 corresponds to the best (most confident) prediction.  Only scores single
	 * genes.
	 * @param reachabilityScores a map from gene names to predicted phenotypic effects
	 * after a single deletion
	 * @param nodeScores from SDREM
	 * @param degrees
	 * @return a map from gene names to ranks
	 */
	public static HashMap<String, Integer> reachabilityRanks(HashMap<String, Double> reachabilityScores,
			HashMap<String, Double> nodeScores,
			HashMap<String, Integer> degrees)
	{
		return reachabilityRanks(reachabilityScores, nodeScores, degrees, false);
	}
	
	/**
	 * Sort genes or pairs of genes by their reachability scores using node score, degree, and node id
	 * for tie-breaking.  Uses ReachabilityComp or ReachabilityPairComp to sort.
	 * In the returned ranks, rank 1 corresponds to the best (most confident) prediction.
	 * Supports single gene deletions or pairs of genes
	 * @param reachabilityScores a map from gene names or gene-gene pairs
	 * to predicted phenotypic effects after a single/double deletion
	 * @param nodeScores from SDREM
	 * @param degrees
	 * @param usePairs if true, sort the gene-gene pairs using ReachabilityPairComp
	 * @return a map from gene names or gene-gene pairs to ranks
	 */
	public static HashMap<String, Integer> reachabilityRanks(HashMap<String, Double> reachabilityScores,
			HashMap<String, Double> nodeScores,
			HashMap<String, Integer> degrees,
			boolean usePairs)
	{
		ArrayList<String> sortedNodes = new ArrayList<String>(reachabilityScores.keySet());
		if(usePairs)
		{
			// In this case the reachability scores are actually predicted genetic interactions
			Collections.sort(sortedNodes, new ReachabilityPairComp(reachabilityScores, nodeScores, degrees));
		}
		else
		{
			Collections.sort(sortedNodes, new ReachabilityComp(reachabilityScores, nodeScores, degrees));
		}
		
		// Iterate through the ordered nodes and store the ranks
		HashMap<String, Integer> rankedNodes = new HashMap<String, Integer>();
		for(int r = 0; r < sortedNodes.size(); r++)
		{
			rankedNodes.put(sortedNodes.get(r), r + 1);
		}
		
		return rankedNodes;
	}
			
	/**
	 * Takes a list of ranked predictions and a set of the examples in the positive (true)
	 * class and calculates the area under the curve.  The ranks are specified such
	 * that 1 is the most confident prediction.
	 * Algorithm from http://www.csd.uwo.ca/faculty/ling/papers/ijcai03.pdf (Eq 1) and
	 * http://www.springerlink.com/content/nn141j42838n7u21/fulltext.pdf
	 * @param rankedPreds
	 * @param positiveExamples
	 * @return
	 */
	public static double areaUnderCurve(HashMap<String, Integer> rankedPreds,
			Set<String> positiveExamples)
	{
		// Assume that all ranked predictions that aren't in the positive class
		// are negative
		int nPreds = rankedPreds.size();
		double nPos = positiveExamples.size();
		int nNeg = nPreds - positiveExamples.size();
		
		int rankSum = 0;
		for(String posEx : positiveExamples)
		{
			if(!rankedPreds.containsKey(posEx))
			{
				throw new IllegalStateException(posEx + " is not ranked");
			}
			
			// The nodes are ranked with 1 being the most confident, but the equation
			// assumes the most confident is ranked last.  Transform ranks so that
			// 1 is worst and N is best
			rankSum += (1 + nPreds - rankedPreds.get(posEx));
		}
		
		return (rankSum - nPos*(nPos+1)/2)/(nPos*nNeg);
	}
	
	/**
	 * Takes a ranked list of predictions and returns the number of top-ranked
	 * predictions that are in the positive (true) class.
	 * @param rankedPreds ranked predictions, with 1 corresponding to the most
	 * confident prediction
	 * @param positiveExamples predictions that are present in this set are
	 * considered to be correct
	 * @param threshold the number of top-ranked predictions to consider
	 * @return
	 */
	public static int thresholdRankedPreds(HashMap<String, Integer> rankedPreds,
			Set<String> positiveExamples, int threshold)
	{
		int count = 0;
		for(String posEx : positiveExamples)
		{
			// Assumed predictions are ranked from 1 to N with 1 being the most confident
			if(rankedPreds.get(posEx) <= threshold)
			{
				count++;
			}
		}
		
		return count;
	}
	
	/**
	 * Takes a ranked list of predictions and returns the number of top-ranked
	 * predictions that are in the positive (true) class broken down to include
	 * the amount of evidence for each true prediction.  That is, if m true predictions
	 * are supported by exactly n pieces of evidence, the String returned will
	 * denote this by n(m).
	 * @param rankedPreds ranked predictions, with 1 corresponding to the most
	 * confident prediction
	 * @param positiveExamples predictions that are present in this set are
	 * considered to be correct
	 * @param exampleCounts the pieces of evidence that support each positive example.
	 * Must contain all elements in the set positiveExamples
	 * @param threshold the number of top-ranked predictions to consider
	 * @return a String summarizing the number of correct predictions and the
	 * amount of evidence for each correct prediction.
	 */
	public static String thresholdRankedPreds(HashMap<String, Integer> rankedPreds,
			Set<String> positiveExamples, HashMap<String, Integer> exampleCounts,
			int threshold)
	{
		HashMap<Integer, Integer> histogram = new HashMap<Integer, Integer>();
		for(String posEx : positiveExamples)
		{
			// Assumed predictions are ranked from 1 to N with 1 being the most confident
			if(rankedPreds.get(posEx) <= threshold)
			{
				if(!exampleCounts.containsKey(posEx))
				{
					throw new IllegalStateException("Cannot obtain the count for the example " + posEx);
				}
				int count = exampleCounts.get(posEx);
				
				int oldVal = 0;
				if(histogram.containsKey(count))
				{
					oldVal = histogram.get(count);
				}
				histogram.put(count, oldVal+1);
			}
		}
		
		ArrayList<Integer> sorted = new ArrayList<Integer>(histogram.keySet());
		Collections.sort(sorted);
		Collections.reverse(sorted);
		
		StringBuffer buf = new StringBuffer();
		for(int key : sorted)
		{
			buf.append(histogram.get(key)).append("(").append(key).append("):");
		}
		buf.deleteCharAt(buf.length()-1);
		
		return buf.toString();
	}
	
	/**
	 * Calculate the genetic interactions between pairs of genes given their
	 * observed single and double deletion phenotypes.  These phenotypes could
	 * be biological or derived computationally.  The genetic interaction is scored
	 * as:<br>
	 * eps = phenotype_observed_ab - phenotype_expected_ab<br>
	 * = phenotype_observed_ab - phenotype_observed_a * phenotype_observed_b<br>
	 * Therefore a negative genetic interaction means that the double knockout
	 * had a stronger than expected effect because smaller scores (e.g.
	 * colony size, growth fitness, or target connectivity) correspond to a more
	 * potent phenotype.
	 * @param singleScores
	 * @param doubleScores
	 * @return a map from gene-gene pairs to their genetic interaction score.
	 */
	public static HashMap<String, Double> geneticInts(HashMap<String, Double> singleScores,
			HashMap<String, Double> doubleScores)
	{
		HashMap<String, Double> geneticInts = new HashMap<String, Double>();
		
		for(String pair : doubleScores.keySet())
		{
			String[] nodes = pair.split("\t");
			String a = nodes[0];
			String b = nodes[1];
			
			if(!singleScores.containsKey(a) || !singleScores.containsKey(b))
			{
				throw new IllegalStateException("Cannot find single deletion scores for "
						+ a + " and " + b);
			}
			
			double expected = singleScores.get(a) * singleScores.get(b);
			geneticInts.put(pair, doubleScores.get(pair) - expected);
		}
		
		return geneticInts;
	}
	
	/**
	 * Find which TFs are differentially expressed in a time series experiment
	 * @param exprFile the DREM table that includes the expression values for
	 * genes that were used to build the DREM model (met the fold change and missing
	 * values criteria)
	 * @param targFile A headerless list of targets
	 * @param outFile
	 * @param noaPrefix The prefix of the Cytoscape noa (node attribute) files that
	 * specifiy which TFs are targets (exceed the activity score theshold) at each
	 * time point
	 * @param noaFiles The number of noa files, which are assumed to be 0-indexed
	 * @param timepoints The number of time points in the expression file (including the
	 * 0 time point)
	 * @param logThresh The log fold change expression threshold, which typically
	 * will be the same threshold used to build the DREM model
	 * @param orientedGraph An oriented network that has been loaded based on
	 * the stored path edges from the SDREM run
	 */
	public static void findDETFs(String exprFile, String targFile, String outFile,
			String noaPrefix, int noaFiles,
			int timepoints, double logThresh,
			EdgeOrientAlg orientedGraph)
	{
		try
		{
			// Offset for the gene name and spot id
			// Check the expression at the 0 time point (which will always be 0)
			// so that the labels are available for the TF activity analysis
			int offset = 2;
			
			// Parse the expression table from DREM to track which the time points
			// at which each gene exceed the log fold change threshold
			BufferedReader reader = new BufferedReader(new FileReader(exprFile));
			String[] timeLabels = new String[timepoints];
			String[] header = reader.readLine().split("\t");
			for(int t = 0; t < timepoints; t++)
			{
				timeLabels[t] = header[t+offset];
			}
			System.out.println("Checking time points " + timeLabels[0] +
					" through " + timeLabels[timepoints-1]);
			
			String line;
			HashMap<String, String> deTimes = new HashMap<String, String>();
			while((line = reader.readLine()) != null)
			{
				String[] parts = line.split("\t");
				String gene = parts[0];
				
				StringBuffer buf = new StringBuffer();
				for(int t = 0; t < timepoints; t++)
				{
					if(!parts[t+offset].equals(""))
					{
						double expr = Double.parseDouble(parts[t+offset]);
						if(Math.abs(expr) >= logThresh)
						{
							buf.append(timeLabels[t]).append(",");
						}
					}
				}
				buf.deleteCharAt(buf.length()-1);
				deTimes.put(gene, buf.toString());
			}
			reader.close();
			
			// Parse the noa files to build a map from TFs to a list of splits
			// where they are active
			System.out.println("Checking active splits " + timeLabels[0] +
					" through " + timeLabels[noaFiles-1]);
			ArrayList<HashSet<String>> activeTFs = new ArrayList<HashSet<String>>(noaFiles);
			HashSet<String> allActiveTFs = new HashSet<String>();
			for(int n = 0; n < noaFiles; n++)
			{
				HashSet<String> curActiveTFs = new HashSet<String>();
				
				reader = new BufferedReader(new FileReader(noaPrefix + n + ".noa"));
				reader.readLine(); // skip the header
				while((line = reader.readLine()) != null)
				{
					String parts[] = line.split(" ");
					if(parts[2].equalsIgnoreCase("Target"))
					{
						// Lookup the ORF name
						String orfName = MapUtil.getOrf(parts[0]);
						curActiveTFs.add(orfName);
						allActiveTFs.add(orfName);
					}
				}
				reader.close();
				
				activeTFs.add(curActiveTFs);
			}
			
			// Build the Strings that describe which splits the TF is active at
			HashMap<String, String> activeTimes = new HashMap<String, String>();
			for(String targ : allActiveTFs)
			{
				StringBuffer buf = new StringBuffer();
				for(int t = 0; t < noaFiles; t++)
				{
					if(activeTFs.get(t).contains(targ))
					{
						buf.append(timeLabels[t]).append(",");
					}
				}
				// We know each TF is active at at least one split
				buf.deleteCharAt(buf.length()-1);
				activeTimes.put(targ, buf.toString());
			}
			
			ArrayList<String> targList = MapUtil.loadList(targFile, false);
			if(!allActiveTFs.containsAll(targList))
			{
				throw new IllegalStateException("Target list and dynamic target " +
						"activities are inconsistent");
			}
			Collections.sort(targList);
			
			// Check all paths on the oriented graph to see if there is at least one path
			// to a given target that involves protein-DNA edges
			// Updated to calculate the fraction of paths that contain at least one protein-DNA edge
			ArrayList<Path> paths = orientedGraph.getSatisfiedPaths();
			// Store the number of paths to each target and the number of paths
			// with at least one protein-DNA edge to each target
			HashMap<String, Integer> pdiPaths = new HashMap<String, Integer>();
			HashMap<String, Integer> allPaths = new HashMap<String, Integer>();

			for(String targ : targList)
			{
				pdiPaths.put(targ, 0);
				allPaths.put(targ, 0);
			}

			for(Path path : paths)
			{
				String targ = path.getTarget().getName();

				if(!allActiveTFs.contains(targ))
				{
					throw new IllegalStateException(targ + " is not a recognized target name");
				}

				allPaths.put(targ, allPaths.get(targ) + 1);
				if(path.containsDirEdge())
				{
					pdiPaths.put(targ, pdiPaths.get(targ) + 1);
				}
			}
			
			PrintWriter writer = new PrintWriter(new FileWriter(outFile));
			writer.println("TF\tORF name\tExpr >= " + logThresh + "\tActive splits" +
					"\t# paths\t% paths with protein-DNA edge");
			
			for(String targ : targList)
			{
				writer.print(MapUtil.getStandard(targ) + "\t" + targ + "\t");
				if(deTimes.containsKey(targ))
				{
					writer.print(deTimes.get(targ) + "\t");
				}
				else
				{
					writer.print("None\t");
				}
				
				if(activeTimes.containsKey(targ))
				{
					writer.print(activeTimes.get(targ) + "\t");
				}
				else
				{
					writer.print("None\t");
				}
				
				int allCount = allPaths.get(targ);
				writer.print(allCount + "\t");
				
				if(allCount > 0)
				{
					writer.println(pdiPaths.get(targ)/(double)allCount);
				}
				else
				{
					writer.println("No paths");
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
	 * Loads the SDREM model from the node score file and determines how many negatives
	 * (non-gold standard nodes) exist among the possible internal nodes and targets.
	 * Uses that to calculate the false positive rate.
	 * @param actScoreFile the activity score file
	 * @param nodeScoreFile the node score file
	 * @param keggMembers a list of genes known to be involved in the HOG pathway
	 * based on the KEGG database.  Uses standard names.
	 * @param sciSigMembers a list of genes known to be involved in the HOG pathway
	 * based on the Science Signaling database.  Uses standard names.
	 * @param otherGSMembers members of the HOG pathway as depicted in Figure 1 of Krantz2009
	 * and taken from other HOG reviews.  Uses standard names.
	 * @param outFile
	 */
	public static void hogPredictionFPR(String actScoreFile,
			String nodeScoreFile,
			double nodeScoreThreshold,
			String keggMembers, String sciSigMembers,
			String otherGSMembers,
			String outFile)
	{
		try
		{
			// Load the gold standard
			HashSet<String> keggStandard = MapUtil.loadSet(keggMembers);
			// Translate from standard to ORF
			HashSet<String> goldStandard = new HashSet<String>();
			for(String gene : keggStandard)
			{
				goldStandard.add(MapUtil.getOrf(gene.toUpperCase().trim()));
			}

			HashSet<String> sciSigStandard = MapUtil.loadSet(sciSigMembers);
			// Translate from standard to ORF
			for(String gene : sciSigStandard)
			{
				goldStandard.add(MapUtil.getOrf(gene.toUpperCase().trim()));
			}

			HashSet<String> krantzStandard = MapUtil.loadSet(otherGSMembers);
			// Translate from standard to ORF
			for(String gene : krantzStandard)
			{
				goldStandard.add(MapUtil.getOrf(gene.toUpperCase().trim()));
			}
			
			// Load the SDREM model
			BufferedReader reader = new BufferedReader(new FileReader(nodeScoreFile));
			HashSet<String> potentialInternal = new HashSet<String>();
			HashSet<String> internal = new HashSet<String>();
			HashSet<String> targets = new HashSet<String>();
			String line = reader.readLine(); // Skip the header
			int nodeCnt = 0, srcCnt = 0;
			
			while((line = reader.readLine()) != null)
			{
				String[] parts = line.split("\t");
				String gene = parts[0];
				nodeCnt++;
				
				if(parts[2].equals("Y"))
				{
					srcCnt++;
				}
				
				// Allow for the possibility that a source is a target too
				if(parts[3].equals("Y"))
				{
					targets.add(gene);
				}
				// The non-source, non-target nodes
				else if(!parts[2].equals("Y"))
				{
					potentialInternal.add(gene);
					double nodeScore = Double.parseDouble(parts[5]);
					if(nodeScore >= nodeScoreThreshold)
					{
						internal.add(gene);
					}
				}
			}
			reader.close();
			
			// Next load the TFs that could have been selected as targets
			Set<String> potentialTargets = MapUtil.loadMap(actScoreFile, 0, 1, false).keySet();
			
			// Determine which potential internal nodes and targets are negatives w.r.t.
			// the gold standard and which are false positives
			HashSet<String> negativeInternal = MapUtil.subtract(potentialInternal, goldStandard);
			HashSet<String> negativeTargets = MapUtil.subtract(potentialTargets, goldStandard);
			HashSet<String> fpInternal = MapUtil.subtract(internal, goldStandard);
			HashSet<String> fpTargets = MapUtil.subtract(targets, goldStandard);
			
			PrintWriter writer = new PrintWriter(new FileWriter(outFile));
			// Write all the counts for sanity-checking
			writer.println("Sources: " + srcCnt);
			writer.println("Gold standard: " + goldStandard.size());
			writer.println("Total nodes: " + nodeCnt);
			writer.println("\nTargets: " + targets.size());
			writer.println("Targets false positives: " + fpTargets.size());
			writer.println("Targets negatives: " + negativeTargets.size());
			writer.println("Targets FPR: " + ((double) fpTargets.size())/negativeTargets.size());
			writer.println("\nInternal: " + internal.size());
			writer.println("Internal false positives: " + fpInternal.size());
			writer.println("Internal negatives: " + negativeInternal.size());
			writer.println("Internal FPR: " + ((double) fpInternal.size())/negativeInternal.size());
			writer.close();
		}
		catch(IOException e)
		{
			e.printStackTrace();
		}
	}
	
	/**
	 * A class that stores much of the state from a node score calculation
	 */
	public static class NodeScoreInfo
	{
		public HashMap<String, Integer> nodeMap;
		public double[][] finalScores;
		public HashMap<String, String> infoMap;

		public NodeScoreInfo(HashMap<String, Integer> nodeMap, double[][] finalScores, HashMap<String, String> infoMap)
		{
			this.nodeMap = nodeMap;
			this.finalScores = finalScores;
			this.infoMap = infoMap;
		}
	}
	
	/**
	 * Sort reachability scores (the predicted effect on target reachability from
	 * the network sources) of a single deletion.  Sorting using this comparator will
	 * put the most confident prediction first, there is no need to reverse the sorted
	 * list.  Node score, node degree, and node id is used to break ties in that order.
	 * Higher node scores or node degrees are preferred (placed earlier in the sorted list).
	 */
	public static class ReachabilityComp implements Comparator<String>
	{
		private HashMap<String, Double> reachabilityScore;
		private HashMap<String, Double> nodeScore;
		private HashMap<String, Integer> degree;
		
		public ReachabilityComp(HashMap<String, Double> reachabilityScoreMap,
				HashMap<String, Double> nodeScoreMap,
				HashMap<String, Integer> degreeMap)
		{
			reachabilityScore = reachabilityScoreMap;
			nodeScore = nodeScoreMap;
			degree = degreeMap;
		}
		
		// Break ties using node score, degree, node name (case-insensitive)
		// For the reachability score, a lower score is better but for the 
		// other tie-breaker metrics higher is better so switch the order of
		// the elements during the comparison so that the comparator will
		// correctly place better elements first
		public int compare(String v1, String v2)
		{
			if(!reachabilityScore.containsKey(v1))
			{
				throw new IllegalStateException("Missing reachability score for " + v1);
			}
			if(!reachabilityScore.containsKey(v2))
			{
				throw new IllegalStateException("Missing reachability score for " + v2);
			}
			
			double r1 = reachabilityScore.get(v1);
			double r2 = reachabilityScore.get(v2);
			if(r1 != r2)
			{
				return Double.compare(r1, r2);
			}
			else
			{
				// Try to break the tie using the node score
				if(!nodeScore.containsKey(v1))
				{
					throw new IllegalStateException("Missing node score for " + v1);
				}
				if(!nodeScore.containsKey(v2))
				{
					throw new IllegalStateException("Missing node score for " + v2);
				}
				
				double n1 = nodeScore.get(v1);
				double n2 = nodeScore.get(v2);
				if(n1 != n2)
				{
					return Double.compare(n2, n1);
				}
				else
				{
					// Try to break the tie using the node degree
					if(!degree.containsKey(v1))
					{
						throw new IllegalStateException("Missing degree for " + v1);
					}
					if(!degree.containsKey(v2))
					{
						throw new IllegalStateException("Missing degree for " + v2);
					}
					
					int d1 = degree.get(v1);
					int d2 = degree.get(v2);
					if(d1 != d2)
					{
				        if (d1 < d2)
				        {
				            return 1;
				        }
				        else // d1 > d2
				        {
				            return -1;
				        }
					}
					else
					{
						// Break the tie using the node names
						// This sorting String versions of numbers in the order
						// opposite of what would be expected (e.g. "3", "2", "1")
						// but is fine as long as it is consistent
						return v2.compareToIgnoreCase(v1);
					}
				}
			}
		}
	}
	
	/**
	 * Sort gene-gene pair reachability scores (the predicted effect on target reachability from
	 * the network sources) of a double deletions.  Sorting using this comparator will
	 * put the most confident prediction first, there is no need to reverse the sorted
	 * list.  Average node score, average node degree, and node ids are used to break ties
	 * in that order.  Higher node scores or node degrees are preferred (placed
	 * earlier in the sorted list).
	 */
	public static class ReachabilityPairComp implements Comparator<String>
	{
		private HashMap<String, Double> geneticScore;
		private HashMap<String, Double> nodeScore;
		private HashMap<String, Integer> degree;
		
		public ReachabilityPairComp(HashMap<String, Double> geneticScoreMap,
				HashMap<String, Double> nodeScoreMap,
				HashMap<String, Integer> degreeMap)
		{
			geneticScore = geneticScoreMap;
			nodeScore = nodeScoreMap;
			degree = degreeMap;
		}
		
		// Break ties using node score, degree, pair name (case-insensitive)
		// taking the average for the two nodes in the pair.
		// For the genetic interaction score prefer negative interactions (lower, i.e.
		// more negative, is better).  For the numeric tiebreakers
		// higher is better so switch the order of
		// the elements during the comparison so that the comparator will
		// correctly place better elements first
		public int compare(String pair1, String pair2)
		{
			if(!geneticScore.containsKey(pair1))
			{
				throw new IllegalStateException("Missing genetic interaction score for " + pair1);
			}
			if(!geneticScore.containsKey(pair2))
			{
				throw new IllegalStateException("Missing genetic interaction score for " + pair2);
			}
			
			double g1 = geneticScore.get(pair1);
			double g2 = geneticScore.get(pair2);
			if(g1 != g2)
			{
				// Sort stronger negative interactions first
				return Double.compare(g1, g2);
			}
			else
			{
				// Try to break the tie using the average node score
				String[] verts1 = pair1.split("\t");
				String[] verts2 = pair2.split("\t");
				
				if(!nodeScore.containsKey(verts1[0]))
				{
					throw new IllegalStateException("Missing node score for " + verts1[0]);
				}
				if(!nodeScore.containsKey(verts1[1]))
				{
					throw new IllegalStateException("Missing node score for " + verts1[1]);
				}
				if(!nodeScore.containsKey(verts2[0]))
				{
					throw new IllegalStateException("Missing node score for " + verts2[0]);
				}
				if(!nodeScore.containsKey(verts2[1]))
				{
					throw new IllegalStateException("Missing node score for " + verts2[1]);
				}
				
				double n1 = (nodeScore.get(verts1[0]) + nodeScore.get(verts1[1]))/2;
				double n2 = (nodeScore.get(verts2[0]) + nodeScore.get(verts2[1]))/2;
				if(n1 != n2)
				{
					return Double.compare(n2, n1);
				}
				else
				{
					// Try to break the tie using the node degree
					if(!degree.containsKey(verts1[0]))
					{
						throw new IllegalStateException("Missing degree for " + verts1[0]);
					}
					if(!degree.containsKey(verts1[1]))
					{
						throw new IllegalStateException("Missing degree for " + verts1[1]);
					}
					if(!degree.containsKey(verts2[0]))
					{
						throw new IllegalStateException("Missing degree for " + verts2[0]);
					}
					if(!degree.containsKey(verts2[1]))
					{
						throw new IllegalStateException("Missing degree for " + verts2[1]);
					}
					
					double d1 = (degree.get(verts1[0]) + degree.get(verts1[1]))/2;
					double d2 = (degree.get(verts2[0]) + degree.get(verts2[1]))/2;
					if(d1 != d2)
					{
				        if (d1 < d2)
				        {
				            return 1;
				        }
				        else // d1 > d2
				        {
				            return -1;
				        }
					}
					else
					{
						// Break the tie using the node names
						// This sorting String versions of numbers in the order
						// opposite of what would be expected (e.g. "3", "2", "1")
						// but is fine as long as it is consistent
						return pair2.compareToIgnoreCase(pair1);
					}
				}
			}
		}
	}
}
