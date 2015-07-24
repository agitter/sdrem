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

import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.util.Calendar;
import java.util.Date;
import java.util.HashMap;
import java.util.Properties;

/**
 * Contains the main method for edge orientation.  Specifies the
 * input files.
 *
 */


public class SDREM {
	
	private static final String VERSION_MESSAGE = "SDREM version 1.1";

	public static void main(String[] args)
	{
		if(args.length < 1)
		{
			throw new RuntimeException("Not enough parameters\n" +
				"Usage: SDREM <properties_file>\n");
		}
		
		runSDREM(args[0]);
	}

	/**
	 * Run SDREM using the settings in the properties file
	 * @param propFile must end in .prop to prevent accidentally running with other files
	 */
	public static void runSDREM(String propFile)
	{
		if(!propFile.toLowerCase().endsWith((".prop"))
				&& !propFile.toLowerCase().endsWith((".props"))
				&& !propFile.toLowerCase().endsWith((".properties")))
		{
			throw new IllegalArgumentException("Expect properties file to end in .prop, .props, or " +
					".properties: " + propFile);
		}

		try
		{
			Properties defaults = new Properties();
			setDefaults(defaults);

			Properties props = new Properties(defaults);

			// Load the properties from the specified file
			FileInputStream propsIn = new FileInputStream(propFile);
			props.load(propsIn);
			propsIn.close();

			runSDREM(props.getProperty("edges.file"),
					props.getProperty("sources.file"),
					props.getProperty("model.dir"),
					props.getProperty("drem.settings.file"),
					props.getProperty("binding.priors.file"),
					props.getProperty("node.priors.file"),
					props.getProperty("stored.paths.dir"),
					Integer.parseInt(props.getProperty("iterations")),
					Integer.parseInt(props.getProperty("max.path.length")),
					Integer.parseInt(props.getProperty("path.enum.bound")),
					Integer.parseInt(props.getProperty("orientation.restarts")),
					Integer.parseInt(props.getProperty("drem.random.runs")),
					Integer.parseInt(props.getProperty("dist.tfs")),
					Integer.parseInt(props.getProperty("top.paths")),
					Double.parseDouble(props.getProperty("dist.thresh")),
					Double.parseDouble(props.getProperty("target.thresh")),
					Double.parseDouble(props.getProperty("node.thresh")),
					Double.parseDouble(props.getProperty("default.node.prior")),
					Double.parseDouble(props.getProperty("min.prior")),
					Float.parseFloat(props.getProperty("random.target.ratio")));
		}
		catch(IOException e)
		{
			e.printStackTrace();
			System.exit(1);
		}
	}

	// Taken from the runSDREM method in the sdrem-www project, but includes nodePriors file.
	// default node prior, random target ratio, number of top paths, and stored paths directory
	/**
	 * Run the iterative SDREM model using the parameters in the properties file
	 */
	public static void runSDREM(String edgeFile,
			String srcFile,
			String modelDir,
			String dremSettingsFile,
			String tfGenePriorsFile,
			String nodePriorsFile,
			String storedPathsDir,
			int iterations,
			int pathLength,
			int pathEnumBound,
			int orientationRestarts,
			int dremRandomizations,
			int distTfs,
			int topRankedPaths,
			double distThresh,
			double targetThresh,
			double nodeThresh,
			double defaultNodePrior,
			double minPrior,
			float randTargThreshold
			)
	{
		// Does not yet support resuming SDREM at an iteration > 1
		int firstItr = 1;
		int lastItr = firstItr + iterations - 1;

		if(!modelDir.endsWith("/"))
		{
			modelDir = modelDir + "/";
		}

		// Fixed settings (almost all of these are now set in the properties file)
		String metric = "PathWeight";

		try
		{
			// In each SDREM iteration, DREM constructs regulatory paths and the
			// active TFs on those paths are connected to the sources in the
			// network via network orientation
			for(int i = 0; i < iterations; i++)
			{
				int itr = i + firstItr;
				long start = System.currentTimeMillis();

				FileOutputStream fos = new FileOutputStream(modelDir + "itr" + itr + ".out");
				PrintStream ps = new PrintStream(fos);
				System.setOut(ps);

				System.out.println(VERSION_MESSAGE);

				if(i==0)
				{
					// Only need to set the node priors file once
					// but do it inside the loop so that the output message
					// is saved to a file
					Vertex.setNodePriors(nodePriorsFile, defaultNodePrior);
				}


				Calendar calendar = Calendar.getInstance();
				Date now = calendar.getTime();
				System.out.println(now);
				System.out.println("Starting iteration " + i + "\n");

				String dremOut = modelDir + itr;
				String targFile = dremOut + ".targets";
				String oldPrior = modelDir + tfGenePriorsFile + (itr - 1) + ".txt";
				String newPrior = modelDir + tfGenePriorsFile + itr + ".txt";
				String runName = "_itr" + itr;
				String dremDefaults = modelDir + dremSettingsFile;

				// Run DREM to generate the new activity scores
				DREMInterface.runDREM(dremDefaults, oldPrior, dremOut, dremRandomizations, distTfs, distThresh);

				calendar = Calendar.getInstance();
				now = calendar.getTime();
				System.out.println(now);
				System.out.println("DREM finished\n\n");


				HashMap<String, Double> targetScores = TargetScores.rankedPathScore(orientationRestarts, topRankedPaths,
						randTargThreshold, pathEnumBound, metric, edgeFile, srcFile, targFile, runName, pathLength,
						storedPathsDir, modelDir);

				TargetScores.NodeScoreInfo nodeResults = TargetScores.rankedPathNodeScore(orientationRestarts,
						topRankedPaths, metric, edgeFile, srcFile, targFile, runName, pathLength, pathEnumBound,
						storedPathsDir, modelDir);

				// Get the node scores, which are the % of top paths the node appears in
				// Only gets scores for non-source, non-target nodes
				HashMap<String, Double> nodeScores = TargetScores.findNodesOfInterest(nodeResults, 0);

				TargetScores.adjustPrior(targetScores, targetThresh, nodeScores, nodeThresh, oldPrior, newPrior, minPrior);

				// TODO the last orientation in rankedPathNodeScore could save this
				// Orient the network once more and save the orientations and
				// write a list of all edges on paths
				Graph g = new Graph();
				DataLoader.readEdgesFromEda(g, edgeFile);
				DataLoader.readSources(g, srcFile);
				DataLoader.readTargets(g, targFile);

				EdgeOrientAlg orient = new EdgeOrientAlg(g);
				// Find new paths only if a directory with previously stored paths
				// was not provided
				if(storedPathsDir == null || storedPathsDir.equals(""))
				{
					orient.findTopPaths(pathLength, pathEnumBound);
				}
				else
				{
					orient.loadStoredTopPaths(pathEnumBound, storedPathsDir);
				}
				orient.findConflicts();
				orient.randPlusSearchSln(orientationRestarts);
				orient.writeConflictOrientations(modelDir + "conflictOrientations_itr" + itr + ".txt");
				orient.writePathEdges(modelDir + "pathEdges_itr" + itr + ".txt");
				orient.writeSatisfiedStringPaths(modelDir + "satisfiedPaths_itr" + itr + ".txt.gz");


				calendar = Calendar.getInstance();
				now = calendar.getTime();
				System.out.println(now);
				long stop = System.currentTimeMillis();
				System.out.println("Total time (s): " + ((stop - start)/1000));
				ps.close();
			}

			// After all iterations have run, run meta-analysis on the results
			TargetScores.analyzeTargetChanges(lastItr,
					".targets", modelDir, modelDir + "targetsByIteration.txt",
					modelDir + "newTargets.txt",
					modelDir + "droppedTargets.txt");
		}
		catch(Exception e)
		{
			e.printStackTrace();
		}
	}

	/**
	 * Set the default SDREM properties
	 * @param defaults
	 */
	public static void setDefaults(Properties defaults)
	{
		defaults.setProperty("edges.file", "example/sampleEdges.txt");
		defaults.setProperty("sources.file", "example/sampleSources.txt");
		defaults.setProperty("model.dir", "example");
		defaults.setProperty("drem.settings.file", "dremSettings.txt");
		defaults.setProperty("binding.priors.file", "bindingPriors");
		defaults.setProperty("node.priors.file", ""); // Do not load specific node priors
		defaults.setProperty("stored.paths.dir", ""); // Do not load precomputed paths
		defaults.setProperty("iterations", "10");
		defaults.setProperty("max.path.length", "5");
		defaults.setProperty("path.enum.bound", "-1"); // -1 Signifies that all paths should be enumerated
		defaults.setProperty("orientation.restarts", "20");
		defaults.setProperty("drem.random.runs", "10");
		defaults.setProperty("dist.tfs", "50");
		defaults.setProperty("top.paths", "-1"); // -1 Signifies 5 * number of targets
		defaults.setProperty("dist.thresh", "0.5");
		defaults.setProperty("target.thresh", "0.8");
		defaults.setProperty("node.thresh", "0.01");
		defaults.setProperty("default.node.prior", "1.0");
		defaults.setProperty("min.prior", "0.01");
		defaults.setProperty("random.target.ratio", "1");
	}
}
