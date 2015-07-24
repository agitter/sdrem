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

import util.Entrez;
import util.MapUtil;
import util.StrOps;


/**
 * Creates files that can be loaded into Cytoscape.
 */
public class CytoscapeInterface {

	public static void generateCytoscape(String srcFile, String targFile, String nodeFile,
			String pathEdges, int edgeMatches, String sifOut, String noaOut)
	{
		if(edgeMatches < 0 || edgeMatches > 2)
		{
			throw new IllegalArgumentException("Only 0, 1, or 2 proteins in a edge " + 
					"can be required to be a node of interest");
		}
		
		try
		{
			HashSet<String> sources = MapUtil.loadSet(srcFile, false);
			HashSet<String> targets = MapUtil.loadSet(targFile, false);
			HashSet<String> nodes = MapUtil.loadSet(nodeFile, false);

			HashSet<String> all = new HashSet<String>();
			all.addAll(sources);
			all.addAll(targets);
			all.addAll(nodes);

			PrintWriter noaWriter = new PrintWriter(new FileWriter(noaOut));
			noaWriter.println("Role (class=java.lang.String)");
			for(String protein : all)
			{
				noaWriter.print(MapUtil.getStandard(protein) + " = ");
				
				if(sources.contains(protein) && targets.contains(protein))
				{
					noaWriter.println("Both");
				}
				else if(sources.contains(protein))
				{
					noaWriter.println("Source");
				}
				else if(targets.contains(protein))
				{
					noaWriter.println("Target");
				}
				else
				{
					noaWriter.println("Other");
				}
			}
			noaWriter.close();
			
			
			PrintWriter sifWriter = new PrintWriter(new FileWriter(sifOut));
			BufferedReader reader = new BufferedReader(new FileReader(pathEdges));
			int count = 2;
			// Skip the heading
			String line = reader.readLine();
			while((line = reader.readLine()) != null)
			{
				String[] parts = line.split("\t");
				// parts[0] : source
				// parts[1] : interaction type
				// parts[2] : target
				// parts[3] : is oriented
				// parts[4] : weight
				String src = parts[0].trim().toUpperCase();
				String targ = parts[2].trim().toUpperCase();
				
				// See if enough interactions in the edge match the nodes of interest
				// to print the edge
				if((edgeMatches == 0) ||
					(edgeMatches == 1 && (all.contains(src) || all.contains(targ))) ||
					(edgeMatches == 2 && all.contains(src) && all.contains(targ)))
				{
					sifWriter.println(MapUtil.getStandard(src) + " " + parts[1] + " " +
							MapUtil.getStandard(targ));
					System.out.println(count);
				}
				count++;
			}
			reader.close();
			sifWriter.close();
		}
		catch(IOException e)
		{
			e.printStackTrace();
		}
	}
	
	
	/**
	 * This version of generateCytoscape takes a map that has the sets
	 * of source, target, and internal nodes already loaded.  It fixes the number
	 * of edges matches to be 2.
	 * @param modelNodes
	 * @param pathEdges
	 * @param sifOut The edges that are between two model members
	 * @param noaOut The roles of the nodes in the SDREM model
	 */
	public static void generateCytoscape(HashMap<String, Set<String>> modelNodes,
			String pathEdges, String sifOut, String noaOut)
	{
		// Fix edgeMatches
		int edgeMatches = 2;
		if(edgeMatches < 0 || edgeMatches > 2)
		{
			throw new IllegalArgumentException("Only 0, 1, or 2 proteins in a edge " + 
					"can be required to be a node of interest");
		}
		
		try
		{
			Set<String> sources = modelNodes.get("sources");
			Set<String> targets = modelNodes.get("targets");
			Set<String> internal = modelNodes.get("internal");

			HashSet<String> all = new HashSet<String>();
			all.addAll(sources);
			all.addAll(targets);
			all.addAll(internal);
			ArrayList<String> allList = new ArrayList<String>(all);
			Collections.sort(allList);

			PrintWriter noaWriter = new PrintWriter(new FileWriter(noaOut));
			noaWriter.println("Role (class=java.lang.String)");
			for(String protein : allList)
			{
				noaWriter.print(protein + " = ");
				
				if(sources.contains(protein) && targets.contains(protein))
				{
					noaWriter.println("Source/Target");
				}
				else if(sources.contains(protein))
				{
					noaWriter.println("Source");
				}
				else if(targets.contains(protein))
				{
					noaWriter.println("Target");
				}
				else
				{
					noaWriter.println("Internal");
				}
			}
			noaWriter.close();
						
			PrintWriter sifWriter = new PrintWriter(new FileWriter(sifOut));
			BufferedReader reader = new BufferedReader(new FileReader(pathEdges));
			// Skip the heading
			String line = reader.readLine();
			while((line = reader.readLine()) != null)
			{
				String[] parts = line.split("\t");
				// parts[0] : source
				// parts[1] : interaction type
				// parts[2] : target
				// parts[3] : is oriented
				// parts[4] : weight
				String src = parts[0].trim().toUpperCase();
				String targ = parts[2].trim().toUpperCase();
				
				// See if enough interactions in the edge match the nodes of interest
				// to print the edge
				if((edgeMatches == 0) ||
					(edgeMatches == 1 && (all.contains(src) || all.contains(targ))) ||
					(edgeMatches == 2 && all.contains(src) && all.contains(targ)))
				{
					// Don't try to find a synonym or id for the nodes
					sifWriter.println(src + " " + parts[1] + " " + targ);
				}
			}
			reader.close();
			sifWriter.close();
		}
		catch(IOException e)
		{
			e.printStackTrace();
		}
	}
	
	
	/**
	 * Similar to generateCytoscape except creates another noa file denoting screen hits.
	 * Also takes all input using Entrez gene ids and writes all output with gene symbols.
	 * Does not check that the sources and target are in any of the edges.
	 * @param srcFile a list of gene ids, no header.
	 * @param targFile the map from gene ids to target weights.  No header.  Gene ids
	 * are in column 0.
	 * @param nodeFile
	 * @param pathEdges has a header.
	 * @param screenFile the node priors file with gene ids in column 0 and no header.
	 * All genes in this file are screen hits.
	 * @param hvFile an optional list of host-virus interactions, no header.  Is not required
	 * to be consistent with the set of sources.  If provided, all edges are written to
	 * the sif file regardless of the edgeMatches value.
	 * @param edgeMatches
	 * @param sifOut
	 * @param typeNoaOut
	 * @param screenNoaOut
	 */
	public static void generateCytoscapeHuman(String srcFile, String targFile,
			String nodeFile, String pathEdges, String screenFile, String hvFile,
			int edgeMatches, String sifOut, String typeNoaOut, String screenNoaOut)
	{
		if(edgeMatches < 0 || edgeMatches > 2)
		{
			throw new IllegalArgumentException("Only 0, 1, or 2 proteins in a edge " + 
					"can be required to be a node of interest");
		}
		
		// Do we have host-virus PPI
		boolean haveHv = hvFile != null && !hvFile.equals("");
		
		try
		{
			HashSet<String> sources = MapUtil.loadSet(srcFile, false);
			Set<String> targets = MapUtil.loadMap(targFile, 0, 0, false).keySet();
			HashSet<String> internal = MapUtil.loadSet(nodeFile, false);

			HashSet<String> allHuman = new HashSet<String>();
			allHuman.addAll(sources);
			allHuman.addAll(targets);
			allHuman.addAll(internal);
			
			if(haveHv)
			{
				// Human proteins that directly interact with viral proteins.  These
				// may or may not be sources
				Set<String> h1 = MapUtil.loadMap(hvFile, 1, 1, false).keySet();
				allHuman.addAll(h1);
			}

			// Write the type of each node in the model
			PrintWriter noaWriter = new PrintWriter(new FileWriter(typeNoaOut));
			noaWriter.println("Role (class=java.lang.String)");
			for(String protein : allHuman)
			{
				noaWriter.print(Entrez.getSymbol(protein) + " = ");
				
				if(sources.contains(protein) && targets.contains(protein))
				{
					noaWriter.println("Both");
				}
				else if(sources.contains(protein))
				{
					noaWriter.println("Source");
				}
				else if(targets.contains(protein))
				{
					noaWriter.println("Target");
				}
				else
				{
					// Internal nodes and H1 proteins that are not members of the model
					noaWriter.println("Other");
				}
			}
			
			// Map from viral proteins to sets of human proteins they interact with
			HashMap<String, HashSet<String>> hvInts = null;
			
			if(haveHv)
			{
				hvInts = MapUtil.loadMultiMap(hvFile, 0, 1, false);
				for(String viralProt : hvInts.keySet())
				{
					noaWriter.println(viralProt + " = Viral");
				}
			}
			
			noaWriter.close();
			
			// Write which nodes in the model are screen hits
			Set<String> screenHits = MapUtil.loadMap(screenFile, 0, 0, false).keySet();
			noaWriter = new PrintWriter(new FileWriter(screenNoaOut));
			noaWriter.println("ScreenHit (class=java.lang.String)");
			for(String protein : allHuman)
			{
				noaWriter.print(Entrez.getSymbol(protein) + " = ");
				
				if(screenHits.contains(protein))
				{
					noaWriter.println("Y");
				}
				else
				{
					noaWriter.println("N");
				}
			}
			
			if(haveHv)
			{
				// Viral proteins are never screen hits
				for(String viralProt : hvInts.keySet())
				{
					noaWriter.println(viralProt + " = N");
				}
			}
			noaWriter.close();
			
			// Write the edges between the nodes in the model in accordance with
			// the number of matches desired
			PrintWriter sifWriter = new PrintWriter(new FileWriter(sifOut));
			BufferedReader reader = new BufferedReader(new FileReader(pathEdges));
			int lineNum = 2;
			// Skip the heading
			String line = reader.readLine();
			while((line = reader.readLine()) != null)
			{
				String[] parts = line.split("\t");
				// parts[0] : source
				// parts[1] : interaction type
				// parts[2] : target
				// parts[3] : is oriented
				// parts[4] : weight
				String src = parts[0].trim().toUpperCase();
				String targ = parts[2].trim().toUpperCase();
				
				// See if enough interactions in the edge match the nodes of interest
				// to print the edge
				if((edgeMatches == 0) ||
					(edgeMatches == 1 && (allHuman.contains(src) || allHuman.contains(targ))) ||
					(edgeMatches == 2 && allHuman.contains(src) && allHuman.contains(targ)))
				{
					sifWriter.println(Entrez.getSymbol(src) + " " + parts[1] + " " +
							Entrez.getSymbol(targ));
					System.out.println(lineNum);
				}
				lineNum++;
			}
			
			// If we have host-virus PPI, write them
			if(haveHv)
			{
				for(String viralProt : hvInts.keySet())
				{
					for(String h1 : hvInts.get(viralProt))
					{
						// Write pd for a directed edge even though it isn't a
						// protein-DNA edge
						sifWriter.println(viralProt + " pd " + Entrez.getSymbol(h1));
					}
				}
			}
			
			sifWriter.close();
		}
		catch(IOException e)
		{
			e.printStackTrace();
		}
	}
	
	public static void generateCytoscapeDynamic(String srcFile, String targFile, String nodeFile,
			String dynamicScoresFile, double scoreThresh, String noaOutPrefix)
	{
		try
		{
			// Load the dynamic scores to determine which TFs are active at each
			// time point
			HashMap<String, HashSet<String>> timeMap = new HashMap<String, HashSet<String>>();
			BufferedReader reader = new BufferedReader(new FileReader(dynamicScoresFile));
			String line;
			while((line = reader.readLine()) != null)
			{
				String[] parts = line.split("\t");
				String time = parts[2];
				HashSet<String> activeTargs;
				if(timeMap.containsKey(time))
				{
					activeTargs = timeMap.get(time);
				}
				else
				{
					activeTargs = new HashSet<String>();
					timeMap.put(time, activeTargs);
				}
				
				double score = Double.parseDouble(parts[1]);
				if(score >= scoreThresh)
				{
					activeTargs.add(MapUtil.getOrf(parts[0]));
				}
			}
			
			HashSet<String> sources = MapUtil.loadSet(srcFile, false);
			HashSet<String> targets = MapUtil.loadSet(targFile, false);
			HashSet<String> nodes = MapUtil.loadSet(nodeFile, false);

			HashSet<String> all = new HashSet<String>();
			all.addAll(sources);
			all.addAll(targets);
			all.addAll(nodes);
			for(String time : timeMap.keySet())
			{
				all.addAll(timeMap.get(time));
			}

			for(String time : timeMap.keySet())
			{
				HashSet<String> activeTargs = timeMap.get(time);
				
				PrintWriter noaWriter = new PrintWriter(new FileWriter(noaOutPrefix + time + ".noa"));
				noaWriter.println("Role (class=java.lang.String)");
				for(String protein : all)
				{
					noaWriter.print(MapUtil.getStandard(protein) + " = ");
					
					if(sources.contains(protein) && activeTargs.contains(protein))
					{
						noaWriter.println("Both");
					}
					else if(sources.contains(protein))
					{
						noaWriter.println("Source");
					}
					else if(activeTargs.contains(protein))
					{
						noaWriter.println("Target");
					}
					else
					{
						noaWriter.println("Other");
					}
				}
				noaWriter.close();
			}
		}
		catch(IOException e)
		{
			e.printStackTrace();
		}
	}
	
	
	// TODO document
	public static void multiConditionHuman(String[] modelNames, 
			String[] sourceFiles, String[] targetFiles, 
			String[] nodeScoresDirs, String nodeScoresPrefix,
			double[] nodeScoreThreshs, String[] pathEdgeFiles,
			String sifOut, String noaPrefix)
	{

		int conds = modelNames.length;
		if(conds != sourceFiles.length ||
				conds != targetFiles.length ||
				conds != nodeScoresDirs.length ||
				conds != nodeScoreThreshs.length ||
				conds != pathEdgeFiles.length)
		{
			throw new IllegalArgumentException("All input arrays must have the same length: " + conds);
		}
		
		try
		{
			// Load all of the predictions
			ArrayList<HashMap<String, Set<String>>> predictions = new ArrayList<HashMap<String, Set<String>>>();
			// Also track all proteins that are predicted in any method
			HashSet<String> allNodes = new HashSet<String>();
			
			for(int p = 0; p < conds; p++)
			{
				HashMap<String, Set<String>> newPreds = TargetScores.loadSDREM(sourceFiles[p], targetFiles[p],
						nodeScoresDirs[p], nodeScoresPrefix, nodeScoreThreshs[p]);
				predictions.add(newPreds);
				
				allNodes.addAll(newPreds.get("sources"));
				allNodes.addAll(newPreds.get("targets"));
				allNodes.addAll(newPreds.get("internal"));
			}
			
			// While looping through all conditions, store the set of edges to write
			HashSet<String> edgesToWrite = new HashSet<String>();
			
			// For each condition write the type of each node in the model
			for(int c = 0; c < conds; c++)
			{
				System.out.println("Condition " + (c+1) + ": " + modelNames[c]);
				
				// Track which nodes are in the model for this condition
				HashSet<String> inCurCond = new HashSet<String>();
				
				PrintWriter noaWriter = new PrintWriter(new FileWriter(noaPrefix + (c+1) + ".noa"));
				noaWriter.println("Condition" + (c+1) + " (class=java.lang.String)");
				for(String protein : allNodes)
				{
					noaWriter.print(Entrez.getSymbol(protein) + " = ");

					// Even if a node is a source at a target, treat it as a source only
					if(predictions.get(c).get("sources").contains(protein))
					{
						noaWriter.println("Source");
						inCurCond.add(protein);
					}
					else if(predictions.get(c).get("targets").contains(protein))
					{
						noaWriter.println("Target");
						inCurCond.add(protein);
					}
					else if(predictions.get(c).get("internal").contains(protein))
					{
						noaWriter.println("Internal");
						inCurCond.add(protein);
					}
					else
					{
						noaWriter.println("Absent");
					}
				}
				System.out.println(inCurCond.size() + " nodes in the SDREM model\n");
				noaWriter.close();
				
				// After writing the node annotations, store all of the edges
				// between the nodes that were included in this condition
				// Do not write them yet because one sif file will be written for a
				// all of the conditions
				BufferedReader reader = new BufferedReader(new FileReader(pathEdgeFiles[c]));

				// Skip the heading
				String line = reader.readLine();
				while((line = reader.readLine()) != null)
				{
					String[] parts = line.split("\t");
					// parts[0] : source
					// parts[1] : interaction type
					// parts[2] : target
					// parts[3] : is oriented
					// parts[4] : weight
					String src = parts[0].trim().toUpperCase();
					String targ = parts[2].trim().toUpperCase();
					
					// Only print edges for which both nodes are in the model
					if(inCurCond.contains(src) && inCurCond.contains(targ))
					{
						String srcSym = Entrez.getSymbol(src);
						if(srcSym == null)
						{
							System.err.println("Could not lookup gene id " + src);
							srcSym = "LOC" + src;
						}
						
						String targSym = Entrez.getSymbol(targ);
						if(srcSym == null)
						{
							System.err.println("Could not lookup gene id " + targ);
							targSym = "LOC" + targ;
						}
						
						edgesToWrite.add(srcSym + " " + parts[1] + " " + targSym);
					}
				}
				reader.close();
			}
			
		
			// Write the edges between the nodes in the models that were saved earlier
			// It is possible to have a directed edge in both directions
			PrintWriter sifWriter = new PrintWriter(new FileWriter(sifOut));
			for(String edge : edgesToWrite)
			{
				sifWriter.println(edge);
			}
			sifWriter.close();
		}
		catch(IOException e)
		{
			e.printStackTrace();
		}
	}
	
	/**
	 * Convert all protein names in a .sif file to title case
	 * @param inFile
	 * @param outFile
	 */
	public static void sifToTitle(String inFile, String outFile)
	{
		try
		{
			BufferedReader reader = new BufferedReader(new FileReader(inFile));
			PrintWriter writer = new PrintWriter(new FileWriter(outFile));
			
			String line;
			while((line = reader.readLine()) != null)
			{
				String[] parts = line.split(" ");
				if(parts.length != 3)
				{
					throw new IllegalStateException("Wrong file format");
				}
				
				writer.println(StrOps.titleCase(parts[0]) + " " +
						parts[1] + " " + 
						StrOps.titleCase(parts[2]));
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
	 * Convert all protein names in a .noa file to title case
	 * @param inFile
	 * @param outFile
	 */
	public static void noaToTitle(String inFile, String outFile)
	{
		try
		{
			BufferedReader reader = new BufferedReader(new FileReader(inFile));
			PrintWriter writer = new PrintWriter(new FileWriter(outFile));
			
			String line;
			// Rewrite the header
			writer.println(reader.readLine());
			while((line = reader.readLine()) != null)
			{
				String[] parts = line.split(" ");
				if(parts.length != 3)
				{
					throw new IllegalStateException("Wrong file format");
				}
				
				writer.println(StrOps.titleCase(parts[0]) + " " +
						parts[1] + " " + parts[2]);
			}
			
			reader.close();
			writer.close();
		}
		catch(IOException e)
		{
			e.printStackTrace();
		}
	}
}
