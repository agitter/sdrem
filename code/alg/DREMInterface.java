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
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Random;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.zip.GZIPOutputStream;

import util.Entrez;
import util.MapUtil;
import edu.cmu.cs.sb.drem.DREM_IO_Batch;

/**
 * This class serves as an interface between edge orientation and DREM.
 * Will likely want to include all of the source in the same project
 * or link one as a JAR of the other, but during development they
 * will remain separate.
 *
 */

public class DREMInterface {

	/**
	 * A list of seeds to be used by the random DREM runs
	 */
	public static ArrayBlockingQueue<Long> seeds;
	
	/**
	 * Run DREM programmatically.
	 * @param defaults A file that contains all input parameters, including the
	 * binding and expression data
	 * @param bindFile the TF-gene binding data with priors
	 * @param saveFile the prefix of the output filename
	 */
	public static void runDREM(String defaults, String bindFile, String saveFile)
	{
//		try
//		{
			// DREM_IO version is not thread-safe
//			DREM_IO.batchMode(defaults, bindFile, saveFile);
			new DREM_IO_Batch(defaults, bindFile, saveFile);
//		}
//		catch(IOException e)
//		{
//			e.printStackTrace();
//		}
	}
	
	
	/**
	 * Run DREM programmatically.  Instead of calculating TF activity scores directly,
	 * randomize the binding interactions and run DREM on the random data to
	 * obtain a distribution of the TF activity scores for top-ranked TFs.
	 * The target score is the percentile in this distribution times 5 (the usual
	 * path length).  The call to the DREM code will block until the DREM run
	 * completes.  The real and random DREM runs are executed in parallel using
	 * as many threads as there are procesors.
	 * @param defaults A file that contains all input parameters, including the
	 * binding and expression data
	 * @param bindFile the TF-gene binding data with priors
	 * @param saveFile the prefix of the output filename, which will also
	 * become the directory where the randomization data is stored
	 * @param randomizations the number of times to randomize the binding data
	 * and run DREM on the random data
	 * @param topTfs use at most this many TFs to form the activity score distribution
	 * @param percentileThreshold only TFs that fall at or above this percentile
	 * in the score distribution are considered to be targets
	 */
	public static void runDREM(String defaults, String bindFile, String saveFile,
			int randomizations, int topTfs, double percentileThreshold)
	{
		try
		{
			// Create a list of seeds for the binding data randomization
			// so that it can be reproduced.  Currently does not use a fixed
			// see to generate the other seeds so that the results are
			// pseudo-random.
			seeds = new ArrayBlockingQueue<Long>(randomizations);
			Random seedMaker = new Random();
//			Random seedMaker = new Random(55555);
			for(int r = 0; r < randomizations; r++)
			{
				seeds.add(seedMaker.nextLong());
			}
			
			// Create a directory for the randomization information
			File randDir = new File(saveFile);
			randDir.mkdir();
			
			// Setup the random DREM runs and DREM run on the real binding data
			Collection<Callable<Boolean>> dremRuns = new ArrayList<Callable<Boolean>>(randomizations + 1);
			for(int r = 0; r < randomizations; r++)
			{
				dremRuns.add(new RandomDREMRun(defaults, bindFile, saveFile, r));
			}
			dremRuns.add(new RealDREMRun(defaults, bindFile, saveFile));

			// Setup the ExcecutorService using one thread per processor on the machine
	        int numProcs = Runtime.getRuntime().availableProcessors();
	        System.out.println("Number of processors available to the Java Virtual Machine: " + numProcs);
	        ExecutorService executor = Executors.newFixedThreadPool(numProcs);
	        
	        // Submit the DREM runs and block until they complete
	        List<Future<Boolean>> results = executor.invokeAll(dremRuns);
	        executor.shutdown();
	        
	        // Verify that all of the runs completed successfully
	        for(Future<Boolean> result : results)
	        {
	        	if(!result.get())
	        	{
	        		throw new IllegalStateException("At least one DREM run did not run to completion");
	        	}
	        }
	        System.out.println("Real and random DREM runs completed");
	        
			
			// Use the TF activity scores from the randomized runs and the real
			// run to calculate a percentile for each TF in the activity
			// score distribution
			ArrayList<Double> dist = new ArrayList<Double>();
			for(int r = 0; r < randomizations; r++)
			{
				ArrayList<Double> curDist = new ArrayList<Double>();
				
				// Read the activity scores from this random run
				String randScores = saveFile + "/rand" + r + ".model.activities";
				BufferedReader reader = new BufferedReader(new FileReader(randScores));
				String line;
				while((line = reader.readLine()) != null)
				{
					String[] parts = line.split("\t");
					curDist.add(Double.parseDouble(parts[1]));
				}
				
				// Sort the activity scores from largest to smallest
				Collections.sort(curDist);
				Collections.reverse(curDist);
				
				// Take the top scores in the distribution
				int top = Math.min(topTfs, curDist.size());
				for(int t = 0; t < top; t++)
				{
					dist.add(curDist.get(t));
				}
			} // Load the random activity scores
			
			// Load the real activity scores and write those TFs
			// whose percentile is above the threshold
			String realScores = saveFile + ".model.activities";
			BufferedReader reader = new BufferedReader(new FileReader(realScores));
			
			String targScores = saveFile + ".targets";
			PrintWriter writer = new PrintWriter(new FileWriter(targScores));
			
			String targScoresStd = saveFile + ".targetsStd";
			PrintWriter writerStd = new PrintWriter(new FileWriter(targScoresStd));
			
			String line;
			while((line = reader.readLine()) != null)
			{
				String[] parts = line.split("\t");
				String curTf = parts[0];
				double curActivity = Double.parseDouble(parts[1]);
				
				// Calculate the percentile for this TF
				int greaterThan = 0;
				for(double randActivity : dist)
				{
					if(curActivity > randActivity)
					{
						greaterThan++;
					}
				}
				
				double percentile = ((double) greaterThan) / dist.size();
				if(percentile >= percentileThreshold)
				{
					// TODO have this use the actual path length?
					// Multiply the percentile by 5, the usual value for k,
					// to get the final target score
					writer.println(curTf + "\t" + (percentile * 5));
					writerStd.println(curTf + "\t" + percentile + "\t" +
							(percentile * 5) + "\t" + curActivity);
				}
			}
			reader.close();
			writer.close();
			writerStd.close();
		}
		catch(IOException e)
		{
			e.printStackTrace();
		}
		// TODO Reconsider how to handle these exceptions
		catch (InterruptedException e)
		{
			System.err.println("DREM run was interrupted");
			throw new RuntimeException(e);
		}
		catch (ExecutionException e)
		{
			System.err.println("DREM run exited with an exception");
			e.printStackTrace();
			throw new RuntimeException(e);
		}
	}
	
	
	/**
	 * Randomize the binding data in the specified priors files.  The values
	 * of the priors for a TF will not be changed, but the targets for each
	 * source will be iteratively randomly swapped using randomizePairs.
	 * Use a default 50 swaps per interaction.
	 * @param priorsIn
	 * @param priorsOut
	 * @see randomizePairs
	 */
	public static void randomizeBindingPriors(String priorsIn, String priorsOut)
	{
		try
		{
			String[][] interactions = readInteractions(priorsIn);
			randomizePairs(interactions, 50 * interactions.length);
			updateInteractions(priorsIn, priorsOut, interactions);
		}
		catch(IOException e)
		{
			e.printStackTrace();
		}
	}
	
	
	/**
	 * Randomize the binding data in the specified priors files.  The values
	 * of the priors for a TF will not be changed, but the targets for each
	 * source will be iteratively randomly swapped using randomizePairs.
	 * @param priorsIn
	 * @param priorsOut
	 * @param swaps
	 * @see randomizePairs
	 */
	public static void randomizeBindingPriors(String priorsIn, String priorsOut, long swaps)
	{
		try
		{
			String[][] interactions = readInteractions(priorsIn);
			randomizePairs(interactions, swaps);
			updateInteractions(priorsIn, priorsOut, interactions);
		}
		catch(IOException e)
		{
			e.printStackTrace();
		}
	}
	
	/**
	 * Swap the edges of the interactions for the specified number of iterations.  Will not
	 * make a swap if the edge already exists in the list.  The randomization
	 * is done in place.  Assumes directed interactions.
	 * @param interactions
	 * @param iterations the number of random swaps to make
	 */
	private static void randomizePairs(String[][] interactions, long iterations)
	{
		// Load all the edges into a HashSet for quick retrieval
		HashSet<String> allEdges = new HashSet<String>();
		for(int p = 0; p < interactions.length; p++)
		{
			// interactions[p][0] is the source
			// interactions[p][1] is the target
			String newPair = interactions[p][0] + "::" + interactions[p][1];
			allEdges.add(newPair);
		}

		// Get a seed from the thread-safe list
		Random rng = new Random(seeds.poll());
		
		long edgeConflicts = 0;
		long successfulSwaps = 0;
		do
		{
			int ind1 = rng.nextInt(interactions.length);
			int ind2 = rng.nextInt(interactions.length);

			// Use the array to get the two directed edges at these indices
			String[] edge1 = interactions[ind1];
			String[] edge2 = interactions[ind2];


			// Check to make sure we don't add an edge with an identical
			// source and target.
			// Also check that we didn't pick a pair of interactions with the same
			// source or same target
			if(!edge1[0].equals(edge2[1]) && !edge2[0].equals(edge1[1])
					&& !edge1[0].equals(edge2[0]) && !edge1[1].equals(edge2[1]))
			{
				// Now we know the swap will create a new directed edge that will
				// not be a self-edge
				// See if either of the new edges is already an existing edge
				String newEdge1 = edge1[0] + "::" + edge2[1];
				String newEdge2 = edge2[0] + "::" + edge1[1];

				// Now that the new edges are constructed we do the check
				if(allEdges.contains(newEdge1) || allEdges.contains(newEdge2))
				{
					// The edge already exists
					edgeConflicts++;
				}
				else
				{
					// Remove the old edges, add the new ones to the HashSet
					String oldPair1 = (new StringBuffer(edge1[0]).append("::").append(edge1[1])).toString();
					String oldPair2 = (new StringBuffer(edge2[0]).append("::").append(edge2[1])).toString();

					if(!allEdges.remove(oldPair1))
					{
						throw new IllegalStateException("Randomization error");
					}
					if(!allEdges.remove(oldPair2))
					{
						throw new IllegalStateException("Randomization error");
					}
					
					allEdges.add(newEdge1);
					allEdges.add(newEdge2);
					
					// Update the original array of interactions by swapping the targets
					String tmp = interactions[ind1][1];
					interactions[ind1][1] = interactions[ind2][1];
					interactions[ind2][1] = tmp;
					
					successfulSwaps++;
				}
			}
			// We created a pair between a gene and itself or chose the same sources or
			// targets
			else
			{
				edgeConflicts++;
			}
		}while(successfulSwaps < iterations);

		if(allEdges.size() != interactions.length)
		{
			throw new IllegalStateException("Randomization error");
		}

		// Need to repopulate the array with the random pairs
//		int curInd = 0;
//		for(String edge : allEdges)
//		{
//			interactions[curInd++] = edge.split("::");
//		}

		// If we have too many pairs an ArrayOutOfBoundsException will have already
		// been thrown, but we need to watch for not enough pairs.
//		if(curInd != interactions.length)
//		{
//			throw new RuntimeException("After randomizing only have " + curInd + " pairs instead of " + interactions.length);
//		}

		System.out.println("Made " + iterations + " swaps with " + edgeConflicts + " conflicts");
	}

	/**
	 * Read a binding file with priors and return a map from TFs to
	 * their prior (the max value with which they bind any gene)
	 * @param prior
	 * @return
	 * @throws IOException
	 */
	public static HashMap<String, Double> readPrior(String prior)
		throws IOException
	{
		BufferedReader reader = new BufferedReader(new FileReader(prior));
		
		// The first line gives the TF names
		String line = reader.readLine();
		// tfs[0] isn't a TF name
		String[] tfs = line.split("\t");
		
		double[] priors = new double[tfs.length];
		
		while((line = reader.readLine()) != null)
		{
			String[] values = line.split("\t");
			for(int i = 1; i < values.length; i++)
			{
				double curVal = Double.valueOf(values[i]);
				priors[i] = Math.max(curVal, priors[i]);
			}
		}
		
		HashMap<String, Double> priorMap = new HashMap<String, Double>();
		for(int i = 1; i < tfs.length; i++)
		{
			priorMap.put(tfs[i], priors[i]);
		}
		
		reader.close();
		
		return priorMap;
	}

	
	/**
	 * Read a binding file with priors and return an array of interactions
	 * where each entry in the array is another array that is [source][target]
	 * @param prior
	 * @return interactions String[][]
	 * @throws IOException
	 */
	public static String[][] readInteractions(String prior)
		throws IOException
	{
		BufferedReader reader = new BufferedReader(new FileReader(prior));
		
		// The first line gives the TF names
		String line = reader.readLine();
		// tfs[0] isn't a TF name
		String[] tfs = line.split("\t");
		
		ArrayList<String[]> intList = new ArrayList<String[]>();
		
		while((line = reader.readLine()) != null)
		{
			String[] values = line.split("\t");
			String gene = values[0];
			for(int i = 1; i < values.length; i++)
			{
				double curVal = Double.valueOf(values[i]);
				if(curVal > 0)
				{
					String[] curInt = {tfs[i], gene};
					intList.add(curInt);
				}
			}
		}
		
		reader.close();
		
		// Now store the interactions as a 2d array
		String[][] interactions = new String[intList.size()][2];
		for(int i = 0; i < intList.size(); i++)
		{
			interactions[i] = intList.get(i);
		}
		
		return interactions;
	}
	
	
	/**
	 * Update the binding interactions in origPriors with the randomized interactions
	 * @param origPriors
	 * @param updatedPriors
	 * @param interactions
	 * @throws IOException
	 */
	public static void updateInteractions(String origPriors, String updatedPriors,
			String[][] interactions) throws IOException
	{
		BufferedReader reader = new BufferedReader(new FileReader(origPriors));
		// Sore the TF names
		String[] tfs = reader.readLine().split("\t");
		// Store the gene names and a map for each gene
		ArrayList<String> genes = new ArrayList<String>();
		HashMap<String, HashSet<String>> genesBoundBy = new HashMap<String, HashSet<String>>();
		String line;
		while((line = reader.readLine()) != null)
		{
			String[] parts = line.split("\t");
			genes.add(parts[0]);
			HashSet<String> tmpSet = new HashSet<String>();
			genesBoundBy.put(parts[0], tmpSet);
		}
		reader.close();
		
		// Populate the genesBoundBy map with the set of TFs that bind each gene
		for(int i = 0; i < interactions.length; i++)
		{
			String curGene = interactions[i][1];
			if(!genesBoundBy.containsKey(curGene))
			{
				throw new IllegalStateException("Unrecognized gene" + curGene);
			}
			
			// Add the binding TF to the set
			String curTf = interactions[i][0];	
			genesBoundBy.get(curGene).add(curTf);
		}

		// Read the TF priors
		HashMap<String, Double> tfPriors = readPrior(origPriors);
		
		// Write the new binding file
		PrintWriter writer = new PrintWriter(new GZIPOutputStream(new FileOutputStream(updatedPriors)));
		writer.print(tfs[0]);
		for(int t = 1; t < tfs.length; t++)
		{
			writer.print("\t" + tfs[t]);
		}
		writer.println();
		
		// Write each gene and its new binding TFs and the correct prior
		for(String gene : genes)
		{
			writer.print(gene);
			
			HashSet<String> tfSet = genesBoundBy.get(gene);
			
			// Loop through the TFs and see which bind the gene
			for(int t = 1; t < tfs.length; t++)
			{
				String curTf = tfs[t];
				if(tfSet.contains(curTf))
				{
					// Assume all TFs have the correct priors because
					// the priors were read from the same binding file
					// used to create the TF list
					writer.print("\t" + tfPriors.get(curTf));
					
					// Remove the written interaction
					tfSet.remove(curTf);
				}
				else
				{
					writer.print("\t0");
				}
			} // Loop through all TFs
			writer.println();
			
			// Ensure all interactions were written
			if(tfSet.size() > 0)
			{
				throw new IllegalStateException("Unrecognized TF(s)");
			}
		} // Loop through all rows of genes
		
		writer.close();
	}
	
	// TODO typo, should be convert
	/**
	 * Generate a new model file that uses TF gene symbols instead of TF gene ids
	 * @param inModelFile
	 * @param outModelFile
	 */
	public static void covertModelIds(String inModelFile, String outModelFile)
	{
		try
		{
			BufferedReader reader = new BufferedReader(new FileReader(inModelFile));
			PrintWriter writer = new PrintWriter(new FileWriter(outModelFile));
			
			String line;
			while((line = reader.readLine()) != null)
			{
				String[] parts = line.split("\t");
				
				// No processing needed
				if(parts.length > 2 || parts[0].equalsIgnoreCase("INTERCEPT"))
				{
					writer.println(line);
				}
				else if(parts[0].equalsIgnoreCase("Num. Coefficients"))
				{
					// Write the next line as well
					writer.println(line);
					writer.println(reader.readLine());
				}
				else
				{
					// Convert the id
					writer.println(Entrez.getSymbol(parts[0]) + "\t" + parts[1]);
				}
			}
			
			reader.close();
			writer.close();
		}
		catch(IOException e)
		{
			e.printStackTrace();
		}
	}
	
	/**
	 * A class that runs DREM using real binding data.
	 *
	 */
	public static class RealDREMRun implements Callable
	{
		private String settings, bindFile, saveFile;
		
		public RealDREMRun(String settings, String bindFile, String saveFile)
		{
			this.settings = settings;
			this.bindFile = bindFile;
			this.saveFile = saveFile;
		}
		
		/**
		 * Run DREM using the settings provided in the constructor
		 * @return true if the DREM run completed without an Exception
		 */
		public Boolean call() throws Exception
		{
			System.out.println("Running DREM on real binding data");
			
			// Run DREM on the real data
			new DREM_IO_Batch(settings, bindFile, saveFile);
			
			// Return true if the run terminated without an exception
			return true;
		}
	}

	/**
	 * Randomize the TF-gene binding data and then run DREM
	 * on that randomized data.
	 *
	 */
	public static class RandomDREMRun implements Callable
	{
		private String settings, bindFile, saveFile;
		private int itr;
		
		public RandomDREMRun(String settings, String bindFile, String saveFile, int itr)
		{
			this.settings = settings;
			this.bindFile = bindFile;
			this.saveFile = saveFile;
			this.itr = itr;
		}
		
		/**
		 * Run DREM using the settings provided in the constructor
		 * @return true if the DREM run completed without an Exception
		 */
		public Boolean call() throws Exception
		{
			System.out.println("Running DREM on randomized binding data, iteration " + itr);
			
			String randBindFile = saveFile + "/rand" + itr + "binding.txt.gz";
			// Randomize the binding data
			randomizeBindingPriors(bindFile, randBindFile);
			// Run DREM on the randomized binding data
			new DREM_IO_Batch(settings, randBindFile, saveFile + "/rand" + itr);
			
			// Return true if the run terminated without an exception
			return true;
		}
	}
}
