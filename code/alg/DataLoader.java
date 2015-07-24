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
import java.io.FileReader;
import java.io.FileWriter;
import java.io.FilenameFilter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Random;

import util.MapUtil;
import util.StrOps;

/**
 * Used to read in data from a variety of file formats and set
 * up the graph object.
 *
 */
public class DataLoader {

	/**
	 * Reads a list of directed and undirected edges from a file format
	 * similar to the EDA (Cytoscape edge attribute file) file format.
	 * Directed edges are denoted by (pd) and undirected by (pp).  The weights
	 * in the file with not be transformed in any way.  Lines take the form:<br>
	 * [vertex 1] \t [type] \t [vertex 2] \t [weight]
	 * @param g the Graph to which edges will be added
	 * @param edaFile
	 * @throws IOException
	 */
	public static void readEdgesFromEda(Graph g, String edaFile) throws IOException
	{
		BufferedReader reader = new BufferedReader(new FileReader(edaFile));
		String line;
		int count = 0;
		int initialVerts = Vertex.uniqueVertices();

		while((line = reader.readLine()) != null)
		{
			String[] parts = line.split("\t");

			// Create or retrieve the vertices
			Vertex v1 = Vertex.createVertex(parts[0]);
			Vertex v2 = Vertex.createVertex(parts[2]);

			// Parse the weight
			double weight = Double.parseDouble(parts[3]);

			// TODO may need a new case for post-translational modifications
			// Create the correct type of edge
			if(parts[1].equalsIgnoreCase("pd") || parts[1].equalsIgnoreCase("(pd)")
					|| parts[1].equalsIgnoreCase("ptm") || parts[1].equalsIgnoreCase("(ptm)"))
			{
				// pd is a directed protein-DNA edge from v1 to v2
				DirEdge e = new DirEdge(v1, v2, weight);
				g.addEdge(e);
			}
			else if(parts[1].equalsIgnoreCase("pp") || parts[1].equalsIgnoreCase("(pp)"))
			{
				// pp is an undirected protein-protein edge
				UndirEdge e = new UndirEdge(v1, v2, weight);
				g.addEdge(e);
			}
			else
			{
				throw new IOException(parts[1] + " is not a valid edge type");
			}
			count++;
		}

		System.out.println("Loaded " + count + " edges");
		System.out.println("Loaded " + (Vertex.uniqueVertices() - initialVerts)+ " new vertices " +
				"(" + Vertex.uniqueVertices() + " total)");
		reader.close();
	}

	/**
	 * Reads a list of directed edges from the EDA file (Cytoscape edge attribute file)
	 * output by the Physical Network Models (PNM) code.
	 * Directed edges are denoted by (pd) and undirected by (pp), but all edges should
	 * have been assigned a direction.  The same edges must have already been
	 * loaded in the graph, and the weights in the graph's edges are not modified.
	 * Any edge that is not in the PNM results is removed from the graph, and
	 * edges that are in the graph will be directed as specified.  Must be called
	 * before any edges in the graph are oriented.
	 * Lines take the form:<br>
	 * [vertex 1] [type] [vertex 2] = [direction]
	 * @param g the Graph from which edges will be removed
	 * @param edaFile the PNM output
	 * @throws IOException
	 */
	public static void readEdgesFromPNM(Graph g, String edaFile) throws IOException
	{
		// Track which edges that are already in the graph should be
		// removed because they are not in the PNM results
		HashSet<DirEdge> dirToRemove = new HashSet<DirEdge>(g.getDirEdges());
		HashSet<UndirEdge> undirToRemove = new HashSet<UndirEdge>(g.getUndirEdges());

		// Verify no undirected edges have been oriented
		for(UndirEdge e : undirToRemove)
		{
			if(e.isOriented())
			{
				throw new IllegalStateException("Undirected edges must not have been oriented");
			}
		}

		BufferedReader reader = new BufferedReader(new FileReader(edaFile));
		String line;
		int count = 0, pseudoCount = 0;

		// Skip the header
		reader.readLine();

		// For each directed edge in the PNM results, remove it from the
		// list of edges to remove (i.e. keep it in the graph)
		while((line = reader.readLine()) != null)
		{
			String[] parts = line.split(" ");

			// Skip over the false edges that are required for PNM
			if(parts[2].contains("TARG"))
			{
				pseudoCount++;
			}
			else
			{
				// Retrieve the vertices
				Vertex v1 = Vertex.createVertex(parts[0]);
				Vertex v2 = Vertex.createVertex(parts[2]);

				// Store the type and direction
				String type = parts[1];
				String direction = parts[4];

				// Keep the correct type of edge
				if(type.equalsIgnoreCase("pd") || type.equalsIgnoreCase("(pd)"))
				{
					if(direction != "+")
					{
						throw new IllegalStateException("Directed edge from " +
								v1 + " to " + v2 + " is oriented in the wrong direction");
					}

					// pd is a directed protein-DNA edge from v1 to v2
					ArrayList<Edge> searchEdges = g.getDirEdges(v1);

					// Find which of the edges out of v1 is the edge to v2
					DirEdge keepEdge = null;
					for(Edge e : searchEdges)
					{
						if(e.getTarget() == v2)
						{
							keepEdge = (DirEdge) e;
						}
					}

					if(keepEdge == null)
					{
						throw new IllegalStateException("Graph does not contain " +
								"a directed edge from " + v1 + " to " + " v2");
					}
					else if(keepEdge.getSource() != v1)
					{
						throw new IllegalStateException("Directed edge " + keepEdge +
								" should have source " + v1);
					}
					// If we found the edge, keep it in the graph
					else if(!dirToRemove.remove(keepEdge))
					{
						throw new IllegalStateException("Graph does not contain " +
								"the expected directed edge from " + v1 + " to " + " v2");
					}
				}
				else if(type.equalsIgnoreCase("pp") || type.equalsIgnoreCase("(pp)"))
				{
					// pp is an undirected protein-protein edge
					ArrayList<Edge> searchEdges = g.getUndirEdges(v1);

					// Find which of the edges out of v1 is the edge to v2
					UndirEdge keepEdge = null;
					for(Edge e : searchEdges)
					{
						if(e.containsVertex(v2))
						{
							keepEdge = (UndirEdge) e;
						}
					}

					if(keepEdge == null)
					{
						throw new IllegalStateException("Graph does not contain " +
								"an undirected edge between " + v1 + " and " + v2);
					}
					else if(!keepEdge.containsVertex(v1))
					{
						throw new IllegalStateException("Undirected edge " + keepEdge +
								" should contain vertex " + v1);
					}
					// If we found the edge, keep it in the graph...
					else if(!undirToRemove.remove(keepEdge))
					{
						throw new IllegalStateException("Graph does not contain " +
								"the expected undirected edge between " + v1 + " and " + v2);
					}

					// ...and orient the edge away from the new source
					if(direction.equals("+"))
					{
						keepEdge.setOrientation(keepEdge.findDirection(v1));
					}
					else if(direction.equals("-"))
					{
						keepEdge.setOrientation(keepEdge.findDirection(v2));
					}
					else
					{
						throw new IllegalStateException("Direction must be + or -");
					}
				}
				else
				{
					throw new IOException(type + " is not a valid edge type");
				}
				count++;
			}
		}

		System.out.println("Kept " + count + " directed edges from PNM");
		System.out.println("Found " + pseudoCount + " edges to fake targets in PNM results");
		reader.close();


		// Remove the edges that were not in the PNM results
		System.out.println("Removing " + (undirToRemove.size() + dirToRemove.size()) +
				" edges from the graph");
		for(DirEdge e : dirToRemove)
		{
			g.removeEdge(e);
		}
		for(UndirEdge e : undirToRemove)
		{
			g.removeEdge(e);
		}

		// Ensure that all undirected edges were oriented by the PNM results
		for(Edge e : g.getUndirEdges())
		{
			if(!e.isOriented())
			{
				throw new IllegalStateException(e + " must be oriented after loading " +
						" PNM results");
			}
		}
	}
	
	// TODO Will fail if both a PPI and protein-DNA edge between
	// the same vertices was present
	/**
	 * Reads a list of edges, which should all be directed, that were 
	 * written by EdgeOrientAlg.writePathEdges.  Only contains edges
	 * that were on a path, and does not have their weights. Lines take the form:<br>
	 * [vertex 1] \t [type] \t [vertex 2] \t [isOriented]
	 * @param g the Graph to which edges will be added
	 * @param pathEdgeFile
	 * @throws IOException
	 */
	public static void readPathEdges(Graph g, String pathEdgeFile) throws IOException
	{
		BufferedReader reader = new BufferedReader(new FileReader(pathEdgeFile));
		// Skip the header
		String line = reader.readLine();
		int count = 0;
		int initialVerts = Vertex.uniqueVertices();

		while((line = reader.readLine()) != null)
		{
			String[] parts = line.split("\t");

			// Create or retrieve the vertices
			Vertex v1 = Vertex.createVertex(parts[0]);
			Vertex v2 = Vertex.createVertex(parts[2]);

			// All edges should be oriented
			if(parts[3].equals("true"))
			{
				if(parts[1].equals("pd"))
				{
					// A directed edge from v1 to v2
					DirEdge e = new DirEdge(v1, v2);
					g.addEdge(e);
				}
				else if(parts[1].equals("pp"))
				{
					// A undirected edge that is then oriented from v1 to v2
					UndirEdge e = new UndirEdge(v1, v2);
					g.addEdge(e);
					e.setOrientation(Edge.FORWARD);
				}
				else
				{
					throw new IOException("Unrecongnized edge type: " + parts[1]);
				}
			}
			else
			{
				throw new IOException("All edges must be oriented");
			}
			count++;
		}

		System.out.println("Loaded " + count + " edges");
		System.out.println("Loaded " + (Vertex.uniqueVertices() - initialVerts)+ " new vertices " +
				"(" + Vertex.uniqueVertices() + " total)");
		reader.close();
	}

	/**
	 * Reads a file containing only target vertices and their target weights.
	 * The set of source vertices
	 * must be set separately.  Each line in the file has the form:<br>
	 * [target name] \t [target weight]
	 * <br>
	 * If a target weight is not specified for a particular target,
	 * it will be set to the default weight.  The target weight will overwrite
	 * the previous target weight for the vertex if one existed.
	 * @param g the Graph to which targets will be added
	 * @param targFile
	 * @throws IOException
	 */
	public static void readTargets(Graph g, String targFile) throws IOException
	{
		BufferedReader reader = new BufferedReader(new FileReader(targFile));
		String line;
		int count = 0;

		while((line = reader.readLine()) != null)
		{
			String[] parts = line.split("\t");

			Vertex v = Vertex.createVertex(parts[0]);

			// A target weight may or may not be present
			if(parts.length > 1)
			{
				double tWeight = Double.parseDouble(parts[1]);
				v.setTargetWeight(tWeight);
			}

			g.addTarget(v);
			count++;
		}

		System.out.println("Loaded " + count + " targets");
		reader.close();
	}

	/**
	 * Reads a file containing a list of source vertices.  Weights are
	 * assumed to have been specified elsewhere.  Each line has the form:<br>
	 * [source name]
	 * @param g the Graph to which sources will be added
	 * @param srcFile
	 * @throws IOException
	 */
	public static void readSources(Graph g, String srcFile) throws IOException
	{
		BufferedReader reader = new BufferedReader(new FileReader(srcFile));
		String line;
		int count = 0;

		while((line = reader.readLine()) != null)
		{
			Vertex v = Vertex.createVertex(line.trim());
			g.addSource(v);
			count++;
		}

		System.out.println("Loaded " + count + " sources");
		reader.close();
	}


	/**
	 * Randomly adds the specified number of sources and targets from the provided
	 * list of genes.  All sources and targets will be unique.
	 * @param g the Graph to which sources and targets will be added
	 * @param geneFile
	 * @param numSrc
	 * @param numTarg
	 * @throws IOException
	 */
	public static void readRandSrcsTargs(Graph g, String geneFile, int numSrc, int numTarg)
		throws IOException
	{
		BufferedReader reader = new BufferedReader(new FileReader(geneFile));
		String line;
		ArrayList<String> genes = new ArrayList<String>();

		// Read in and store all possible genes
		while((line = reader.readLine()) != null)
		{
			genes.add(line.trim());
		}

		if((numSrc + numTarg) > genes.size())
		{
			throw new IllegalArgumentException("There are not enough genes to " +
					"select " + numSrc + " sources and " + numTarg + " targets");
		}

		// Now randomly select the sources
		Random rand = new Random();
		for(int s = 0; s < numSrc; s++)
		{
			int index = rand.nextInt(genes.size());
			Vertex v = Vertex.createVertex(genes.get(index));
			g.addSource(v);

			// Remove the gene from the list so it will not be selected again
			genes.remove(index);
		}

		// Similarly, randomly select the targets
		for(int t = 0; t < numTarg; t++)
		{
			int index = rand.nextInt(genes.size());
			Vertex v = Vertex.createVertex(genes.get(index));
			g.addTarget(v);

			// Remove the gene from the list so it will not be selected again
			genes.remove(index);
		}

		reader.close();
	}





	/**
	 * Randomly adds the specifiec number of sources and targets from the provided
	 * list of genes.  All sources and targets will be unique.
	 * @param geneFile
	 * @param outDir the directory where the files will be written
	 * @param filename this sources file will be "sources_[filename]" and the targets
	 * file will be "targets_[filename]"
	 * @param numSrc
	 * @param numTarg
	 * @throws IOException
	 */
	public static void genRandSrcsTargs(String geneFile, String outDir, String filename, int numSrc, int numTarg)
		throws IOException
	{
		BufferedReader reader = new BufferedReader(new FileReader(geneFile));
		String line;
		ArrayList<String> genes = new ArrayList<String>();

		// Read in and store all possible genes
		while((line = reader.readLine()) != null)
		{
			genes.add(line.trim());
		}

		if((numSrc + numTarg) > genes.size())
		{
			throw new IllegalArgumentException("There are not enough genes to " +
					"select " + numSrc + " sources and " + numTarg + " targets");
		}

		// Now randomly select the sources and write them to the file
		PrintWriter writer = new PrintWriter(new FileWriter(outDir + "sources_" + filename));

		Random rand = new Random();
		for(int s = 0; s < numSrc; s++)
		{
			int index = rand.nextInt(genes.size());
			writer.println(genes.get(index));

			// Remove the gene from the list so it will not be selected again
			genes.remove(index);
		}
		writer.close();

		// Similarly, randomly select the targets and write them to the file
		writer = new PrintWriter(new FileWriter(outDir + "targets_" + filename));
		for(int t = 0; t < numTarg; t++)
		{
			int index = rand.nextInt(genes.size());
			writer.println(genes.get(index));

			// Remove the gene from the list so it will not be selected again
			genes.remove(index);
		}
		writer.close();

		reader.close();
	}

	/**
	 * Reads a file containing only target vertices and their target weights.
	 * The set of source vertices
	 * must be set separately.  Each line in the file has the form:<br>
	 * [target name] \t [target weight]
	 * <br>
	 * If a target weight is not specified for a particular target,
	 * it will be set to the default weight.  The target weight will overwrite
	 * the previous target weight for the vertex if one existed.
	 * If there are N real targets, N random targets will be added as well.  The
	 * random targets will be chosen from all genes in the edges of the graph and
	 * will not include the real targets of any sources that have already been added
	 * to the graph.  Thus, this should be called after calling readEdgesFromEda
	 * and readSources
	 * @param g the Graph to which targets will be added
	 * @param targFile
	 * @return a map from the random vertices to indices
	 * @throws IOException
	 */
	public static HashMap<Vertex, Integer> readRealRandTargets(Graph g, String targFile)
		throws IOException
	{
		return readRealRandTargets(g, targFile, 1);
	}

	/**
	 * Reads a file containing only target vertices and their target weights.
	 * The set of source vertices
	 * must be set separately.  Each line in the file has the form:<br>
	 * [target name] \t [target weight]
	 * <br>
	 * If a target weight is not specified for a particular target,
	 * it will be set to the default weight.  The target weight will overwrite
	 * the previous target weight for the vertex if one existed.
	 * If there are N real targets, (N * randTargsRatio) random targets will be added as well.  The
	 * random targets will be chosen from all genes in the edges of the graph and
	 * will not include the real targets or any sources that have already been added
	 * to the graph.  Thus, this should be called after calling readEdgesFromEda
	 * and readSources
	 * @param g the Graph to which targets will be added
	 * @param targFile
	 * @param randTargsRatio the ratio of random targets to real targets (i.e. 2 gives
	 * twice as many random targets as real targets).  Must be > 0.
	 * @return a map from the random vertices to indices
	 * @throws IOException
	 * @throws IllegalStateException if not enough random targets can be selected
	 */
	public static HashMap<Vertex, Integer> readRealRandTargets(Graph g, String targFile,
			float randTargsRatio) throws IOException
	{
		if(randTargsRatio <= 0)
		{
			throw new IllegalArgumentException("Must have random"
					+ " targets.  " + randTargsRatio + " is not valid.");
		}

		BufferedReader reader = new BufferedReader(new FileReader(targFile));
		String line;
		int count = 0;

		// Read the real targets
		while((line = reader.readLine()) != null)
		{
			String[] parts = line.split("\t");

			Vertex v = Vertex.createVertex(parts[0]);

			// A target weight may or may not be present
			if(parts.length > 1)
			{
				double tWeight = Double.parseDouble(parts[1]);
				v.setTargetWeight(tWeight);
			}

			g.addTarget(v);
			count++;
		}

		System.out.println("Loaded " + count + " real targets");
		reader.close();

		// Determine the set of genes to choose random targets from
		// TODO May want to only include TFs in this set
		HashSet<Vertex> allGenes = new HashSet<Vertex>();
		for(Edge e : g.getAllEdges())
		{
			allGenes.add(e.getSource());
			allGenes.add(e.getTarget());
		}
		allGenes.removeAll(g.getSources());
		allGenes.removeAll(g.getTargets());

		int randTargs = Math.round(count * randTargsRatio);

		if(allGenes.size() < randTargs)
		{
			throw new IllegalStateException("Cannot generate " + randTargs + " random targets");
		}

		// Convert to an array list because it's easier to select random
		// targets from
		ArrayList<Vertex> genes = new ArrayList<Vertex>(allGenes);
		Random rand = new Random();
		HashMap<Vertex, Integer> map = new HashMap<Vertex, Integer>();
		int mapInd = 0;

		// Randomly select the targets
		for(int t = 0; t < randTargs; t++)
		{
			int geneInd = rand.nextInt(genes.size());
			g.addTarget(genes.get(geneInd));

			// Update the map
			map.put(genes.get(geneInd), mapInd);
			mapInd++;

			// Remove the gene from the list so it will not be selected again
			genes.remove(geneInd);
		}

		System.out.println("Added " + randTargs + " random targets");

		return map;
	}
	
	/**
	 * Reads a file containing only target vertices and their target weights.
	 * The set of source vertices
	 * must be set separately.  Each line in the file has the form:<br>
	 * [target name] \t [target weight]
	 * <br>
	 * If a target weight is not specified for a particular target,
	 * it will be set to the default weight.  The target weight will overwrite
	 * the previous target weight for the vertex if one existed.
	 * If there are N real targets, (N * randTargsRatio) random targets will be added as well.
	 * The random targets will not include the real targets or any sources
	 * that have already been added to the graph.  Thus, this should be
	 * called after calling readEdgesFromEda and readSources
	 * @param g the Graph to which targets will be added
	 * @param targFile
	 * @param storedTargetsDir a directory containing files with the form [target].txt.gz.
	 * The targets named by these files are the only ones that can be used as random targets.
	 * @param randTargsRatio the ratio of random targets to real targets (i.e. 2 gives
	 * twice as many random targets as real targets).  Must be > 0.
	 * @return a map from the random vertices to indices
	 * @throws IOException
	 * @throws IllegalStateException if not enough random targets can be selected
	 * from the set of stored targets or if one of the stored targets in the in the network
	 */
	public static HashMap<Vertex, Integer> readRealRandStoredTargets(Graph g, String targFile,
			String storedTargetsDir, float randTargsRatio) throws IOException
	{
		if(randTargsRatio <= 0)
		{
			throw new IllegalArgumentException("Must have random"
					+ " targets.  " + randTargsRatio + " is not valid.");
		}

		BufferedReader reader = new BufferedReader(new FileReader(targFile));
		String line;
		int count = 0;

		// Read the real targets
		while((line = reader.readLine()) != null)
		{
			String[] parts = line.split("\t");

			Vertex v = Vertex.createVertex(parts[0]);

			// A target weight may or may not be present
			if(parts.length > 1)
			{
				double tWeight = Double.parseDouble(parts[1]);
				v.setTargetWeight(tWeight);
			}

			g.addTarget(v);
			count++;
		}

		System.out.println("Loaded " + count + " real targets");
		reader.close();

		// Determine the set of genes to choose random targets from
		File storedTargetsFile = new File(storedTargetsDir);
		File[] possibleTargsFiles = storedTargetsFile.listFiles(new FilenameFilter()
		{
			public boolean accept(File dir, String name)
			{
				return name.endsWith(".txt.gz");
			}
		});
		
		System.out.println("Found " + possibleTargsFiles.length + " stored targets");
		
		HashSet<Vertex> possibleTargs = new HashSet<Vertex>();
		for(File f : possibleTargsFiles)
		{
			String targName = StrOps.trimString(f.getName(), "", ".txt.gz");
			Vertex v = Vertex.findVertex(targName);
			
			if(v == null)
			{
				throw new IllegalStateException(targName + " is not in the interaction network");
			}
			possibleTargs.add(v);
		}
		possibleTargs.removeAll(g.getSources());
		possibleTargs.removeAll(g.getTargets());
		
		System.out.println(possibleTargs.size() + " possible random targets after removing souces and real targets");

		int randTargs = Math.round(count * randTargsRatio);

		if(possibleTargs.size() < randTargs)
		{
			throw new IllegalStateException("Cannot generate " + randTargs + " random targets");
		}

		// Convert to an array list because it's easier to select random
		// targets from
		ArrayList<Vertex> genes = new ArrayList<Vertex>(possibleTargs);
		Random rand = new Random();
		HashMap<Vertex, Integer> map = new HashMap<Vertex, Integer>();
		int mapInd = 0;

		// Randomly select the targets
		for(int t = 0; t < randTargs; t++)
		{
			int geneInd = rand.nextInt(genes.size());
			g.addTarget(genes.get(geneInd));

			// Update the map
			map.put(genes.get(geneInd), mapInd);
			mapInd++;

			// Remove the gene from the list so it will not be selected again
			genes.remove(geneInd);
		}

		System.out.println("Added " + randTargs + " random targets");

		return map;
	}


	/**
	 * Reads the file of targets and uses it to create a map of targets to indices
	 * @param targFile
	 * @return a map from vertices to indices
	 * @throws IOException
	 */
	public static HashMap<Vertex, Integer> createTargetMap(String targFile) throws IOException
	{
		BufferedReader reader = new BufferedReader(new FileReader(targFile));
		HashMap<Vertex, Integer> map = new HashMap<Vertex, Integer>();
		String line;
		int count = 0;

		while((line = reader.readLine()) != null)
		{
			String[] parts = line.split("\t");

			Vertex v = Vertex.createVertex(parts[0]);
			map.put(v, count);
			count++;
		}

		reader.close();

		return map;
	}
	
	/**
	 * Convert a sif file written by CytoscapeInterface.generateCytoscape
	 * to a path edges file (like those written by EdgeOrientAlg.writePathEdges)
	 * so that it can be loaded as a network.  Assumes yeast genes.
	 * @param sifFile
	 * @param outFile
	 */
	public static void convertSifToPathEdges(String sifFile, String outFile)
	{
		try
		{
			BufferedReader reader = new BufferedReader(new FileReader(sifFile));
			PrintWriter writer = new PrintWriter(new FileWriter(outFile));
			
			writer.println("Source\tType\tTarget\tOriented");
			
			String line;
			while((line = reader.readLine()) != null)
			{
				String[] parts = line.split(" ");
				writer.print(MapUtil.getOrf(parts[0]) + "\t");
				writer.print(parts[1] + "\t");
				writer.print(MapUtil.getOrf(parts[2]) + "\t");
				writer.println("true");
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
