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
import util.MapUtil;

/**
 * This class provides a variety of methods for generating
 * activity scores for targets.
 *
 */
public class TargetScores {

	public static final String SYN_FILE = "SGD_standardToOrf.txt";
	/** Map from standard names to ORF names */
	public static HashMap<String, String> synMap = null;
	public static final String REV_SYN_FILE = "SGD_orfToStandard.txt";
	/** Map from ORF names to standard names */
	public static HashMap<String, String> revSynMap = null;

	/**
	 * Uses the top ranked satisfied paths only by some criteria to obtain a
	 * distribution of metrics (path weight, edge use, etc.) for the real
	 * targets and an equal number of random targets.  Does not
	 * presently use this distribution to generate a score, but prints
	 * the distributions to a file and returns the targets' percentiles
	 * within the distribution
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
	public static HashMap<String, Double> rankedPathScore(int runs, int pathThreshold, int randTargsRatio,
			String rankingMetric, String edgeFile, String srcFile, String targFile,
			String runName, int pathLength, String outDir)
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
				HashMap<Vertex, Integer> randMap =
					DataLoader.readRealRandTargets(g, targFile, randTargsRatio);

				EdgeOrientAlg orient = new EdgeOrientAlg(g);

				orient.findPaths(pathLength);
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
	 * number of runs.
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
	 * @param outDir the output directory, the filename will be generated automatically.
	 * If the empty string "", then no output will be written
	 * @return a NodeScoreInfo object that summarizes the output
	 */
	public static NodeScoreInfo rankedPathNodeScore(int runs, int pathThreshold,
			String rankingMetric, String edgeFile, String srcFile, String targFile,
			String runName, int pathLength, String outDir)
	{
		if(revSynMap == null)
		{
			revSynMap = MapUtil.loadMap(REV_SYN_FILE, 0, 1, false);
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

				if(r > 0)
				{
					// Reset all state in between runs.
					// It is inconvenient to remove targets and paths and track
					// state correctly
					g = new Graph();

					DataLoader.readEdgesFromEda(g, edgeFile);
					DataLoader.readSources(g, srcFile);
					DataLoader.readTargets(g, targFile);
				}

				EdgeOrientAlg orient = new EdgeOrientAlg(g);

				orient.findPaths(pathLength);
				orient.findConflicts();				
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
	 * @param pathLength the max path length
	 * @param outDir the output directory, the filename will be generated automatically
	 */
	public static void weightRatioScore(int runs,
			String edgeFile, String srcFile, String targFile,
			String runName, int pathLength, String outDir)
	{
		weightRatioScore(runs, 1, edgeFile, srcFile, targFile, runName, pathLength, outDir);
	}

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
	 * @param runName will be used in the automatic filename creation
	 * @param pathLength the max path length
	 * @param outDir the output directory, the filename will be generated automatically
	 */
	public static void weightRatioScore(int runs, int randTargsRatio,
			String edgeFile, String srcFile, String targFile,
			String runName, int pathLength, String outDir)
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

				orient.findPaths(pathLength);
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
				synMap = MapUtil.loadMap(SYN_FILE, 1, 0, false);
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
						if(parts[t].equals("0"))
						{
							// Write 0 again if the TF doesn't bind this gene
							buf.append(parts[t]);
						}
						else
						{
							double val = Double.parseDouble(parts[t]);
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
						// Write the original value if this TF wasn't a target
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

}
