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
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.util.HashMap;
import java.util.Properties;

public class SDREM {
	public static void main(String[] args)
	{		
		try
		{			
			Properties defaults = new Properties();
			setDefaults(defaults);
			
			Properties props = new Properties(defaults);
			
			if(args.length < 1)
			{
				System.err.println("Not enough parameters\n" +
					"Usage: SDREM <properties_file>\n" + 
					"Using default properties");
			}
			else
			{
				// Load the properties from the specified file
				FileInputStream propsIn = new FileInputStream(args[0]);
				props.load(propsIn);
				propsIn.close();
			}

			runSDREM(props.getProperty("edges.file"),
					props.getProperty("sources.file"),
					props.getProperty("model.dir"),
					props.getProperty("drem.settings"),
					props.getProperty("binding.priors"),
					Integer.parseInt(props.getProperty("iterations")),
					Integer.parseInt(props.getProperty("max.path.length")),
					Integer.parseInt(props.getProperty("orientation.restarts")),
					Integer.parseInt(props.getProperty("drem.random.runs")),
					Integer.parseInt(props.getProperty("dist.tfs")),
					Double.parseDouble(props.getProperty("dist.thresh")),
					Double.parseDouble(props.getProperty("target.thresh")),
					Double.parseDouble(props.getProperty("node.thresh")),
					Double.parseDouble(props.getProperty("min.prior")));
		}
		catch(IOException e)
		{
			e.printStackTrace();
			System.exit(1);
		}
	}
	
	/**
	 * Set the default properties
	 * @param defaults
	 */
	public static void setDefaults(Properties defaults)
	{
		defaults.setProperty("edges.file", "example/sampleEdges.txt");
		defaults.setProperty("sources.file", "example/sampleSources.txt");
		defaults.setProperty("model.dir", "example");
		defaults.setProperty("drem.settings", "dremSettings.txt");
		defaults.setProperty("binding.priors", "priors");
		defaults.setProperty("iterations", "10");
		defaults.setProperty("max.path.length", "5");
		defaults.setProperty("orientation.restarts", "20");
		defaults.setProperty("drem.random.runs", "10");
		defaults.setProperty("dist.tfs", "50");
		defaults.setProperty("dist.thresh", "0.5");
		defaults.setProperty("target.thresh", "0.8");
		defaults.setProperty("node.thresh", "0.01");
		defaults.setProperty("min.prior", "0.01");
	}

	/**
	 * Run the iterative SDREM model using the parameters in the properties file
	 */
	public static void runSDREM(String edgeFile,
			String srcFile,
			String modelDir,
			String dremSettings,
			String tfGenePriors,
			int iterations,
			int pathLength,
			int orientationRestarts,
			int dremRandomizations,
			int distTfs,
			double distThresh,
			double targetThresh,
			double nodeThresh,
			double minPrior
			)
	{
		int firstItr = 1;
		int lastItr = firstItr + iterations - 1;
		
		if(!modelDir.endsWith("/"))
		{
			modelDir = modelDir + "/";
		}
		
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
				
				String dremOut = modelDir + itr;
				String targFile = dremOut + ".targets";
				String oldPrior = modelDir + tfGenePriors + (itr - 1) + ".txt";
				String newPrior = modelDir + tfGenePriors + itr + ".txt";
				String runName = "_itr" + itr;
				String dremDefaults = modelDir + dremSettings;

				// Fixed settings
				int threshold = -1;
				int ratio = 1;
				String metric = "PathWeight";
				String outDir = modelDir;
				
				// Run DREM to generate the new activity scores
				DREMInterface.runDREM(dremDefaults, oldPrior, dremOut, dremRandomizations, distTfs, distThresh);
				
				// Poll for the targets file to be created
				File targets = new File(dremOut + ".model.activities.flag");
				while(!targets.exists());
				
				System.out.println("DREM finished\n\n");

				
				HashMap<String, Double> targetScores = TargetScores.rankedPathScore(orientationRestarts, threshold,
						ratio, metric, edgeFile, srcFile, targFile, runName, pathLength, outDir);
				
				TargetScores.NodeScoreInfo nodeResults = TargetScores.rankedPathNodeScore(orientationRestarts,
						threshold, metric, edgeFile, srcFile, targFile, runName, pathLength, outDir);
				
				// Get the node scores, which are the % of top paths the node appears in
				// Only gets scores for non-source, non-target nodes
				HashMap<String, Double> nodeScores = TargetScores.findNodesOfInterest(nodeResults, 0);
				
				TargetScores.adjustPrior(targetScores, targetThresh, nodeScores, nodeThresh, oldPrior, newPrior, minPrior);

				// Orient the network once more and save the orientations and
				// write a list of all edges on paths
				Graph g = new Graph();
				DataLoader.readEdgesFromEda(g, edgeFile);
				DataLoader.readSources(g, srcFile);
				DataLoader.readTargets(g, targFile);
	
				EdgeOrientAlg orient = new EdgeOrientAlg(g);
				orient.findPaths(pathLength);
				orient.findConflicts();
				orient.randPlusSearchSln(orientationRestarts);
				orient.writeConflictOrientations(modelDir + "conflictOrientations_itr" + itr + ".txt");
				orient.writePathEdges(modelDir + "pathEdges_itr" + itr + ".txt");
				
				
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
}
