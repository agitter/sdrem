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
import java.util.HashMap;
import java.util.HashSet;
import java.util.Random;

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

			// Create the correct type of edge
			if(parts[1].equalsIgnoreCase("pd") || parts[1].equalsIgnoreCase("(pd)"))
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
	 * will not include the real targets of any sources that have already been added
	 * to the graph.  Thus, this should be called after calling readEdgesFromEda
	 * and readSources
	 * @param g the Graph to which targets will be added
	 * @param targFile
	 * @param randTargsRatio the ratio of random targets to real targets (i.e. 2 gives
	 * twice as many random targets as real targets).  Must be > 0.
	 * @return a map from the random vertices to indices
	 * @throws IOException
	 */
	public static HashMap<Vertex, Integer> readRealRandTargets(Graph g, String targFile,
			int randTargsRatio)	throws IOException
	{
		if(randTargsRatio < 1)
		{
			throw new IllegalArgumentException("Must have at least as many random"
					+ " targets as real targets.  " + randTargsRatio + " is not valid.");
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
		HashSet<Vertex> allGenes = new HashSet<Vertex>();
		for(Edge e : g.getAllEdges())
		{
			allGenes.add(e.getSource());
			allGenes.add(e.getTarget());
		}
		allGenes.removeAll(g.getSources());
		allGenes.removeAll(g.getTargets());

		int randTargs = count * randTargsRatio;

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
}
