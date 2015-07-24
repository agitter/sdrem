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
import java.util.Set;
import java.util.Map;


/**
 * Loads actual sgdPathways.  Supports SGD and KEGG sgdPathways.
 */
public class PathwayLoader {

	public static String sgdDir = "../testData/sgd/12.0/data/";
	public static String keggDir = "../testData/kegg/";
	public static String sciSigDir = "../testData/scienceSignaling/";

	/** A counter appended to Vertex ids that keeps distinct pathways separate when
	 * searching for paths */
	private int pathCounter;

	// The following are SGD-specific
	/** A map from monomer and complex ids to lists of genes */
	private HashMap<String, HashSet<String>> monoCompMap;
	/** A map from enzymatic reactions to monomers in the monoCompMap */
	private HashMap<String, String> enzMap;
	/** A map from reactions to enzymatic reactions */
	private HashMap<String, HashSet<String>> reacMap;
	/** The set of all SGD pathways */
	public HashSet<SGDPathway> sgdPathways;

	// The following are KEGG-specific
	/** The set of all KEGG pathways */
	public HashSet<Graph> keggPathways;

	// The following are Science Signaling-specific
	/** The set of all Science Signaling pathways */
	public HashSet<Graph> sciPathways;

	/** Contains all KEGG and Science Signaling signaling pathways */
	public Graph signalPathways;

	/** Contains all SGD metabolic pathways (converted to gene pathways) */
	public Graph metabolicPathways;

	public static void main(String[] args)
	{
		PathwayLoader test = new PathwayLoader();

		try {
//			test.loadSgdPathways();
//			test.loadKeggPathways();
//			test.loadSciPathways();
			test.loadSignalPathways();
			test.loadMetabolicPathways();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	public PathwayLoader()
	{
		pathCounter = 0;
	}

	/**
	 * Loads a mapping of monomers to ORFs and complexes to ORFs.
	 * Uses the file proteins.dat
	 *
	 *@throws IOException
	 */
	public void loadMonomersComplexes() throws IOException
	{
		String data = "proteins.dat";

		// First load the list of non-ORF gene identifiers
		BufferedReader reader = new BufferedReader(new FileReader(sgdDir + "../../geneMap.txt"));
		HashMap<String, String> nonOrfMap = new HashMap<String, String>();
		String line;
		while((line = reader.readLine()) != null)
		{
			// Each line is "G3O-* \t <orf name>"
			String[] parts = line.trim().split("\t");
			nonOrfMap.put(parts[0], parts[1]);
		}
		reader.close();



		reader = new BufferedReader(new FileReader(sgdDir + data));
		monoCompMap = new HashMap<String,HashSet<String>>();

		// Skip the heading
		while((line = reader.readLine()).startsWith("#"));

		// Create a map from modified monomers to unmodified monomers
		// If the modified monomers' entries do not contain genes we can
		// use the modified monomers' genes
		HashMap<String, String> modifiedMap = new HashMap<String, String>();

		String id = "";
		String unmodified = "";
		boolean isMonomer = false;
		HashSet<String> geneList = new HashSet<String>();
		// During the first pass only create the monomer map.  It must be finished
		// before the complexes can be loaded
		do
		{
			line = line.trim();
			// Store the id
			if(line.startsWith("UNIQUE-ID - "))
			{
				if(!id.equals(""))
					throw new IOException("Two ids found for a single entry: " + id +
							" and " + line.substring(12).trim());

				// Remove "UNIQUE-ID - "
				id = line.substring(12).trim();
			}
			// See if this is a monomer or something else
			else if(line.contains("TYPES - Polypeptides") || 
					line.contains("TYPES - Apo-GcvH"))
			{
				isMonomer = true;
			}
			// Store the gene.  There should only be one per entry
			else if(line.startsWith("GENE - "))
			{
				if(geneList.size() > 0)
					throw new IOException("Two genes found for a single entry: " + id);

				// Remove "GENE - "
				String gene = line.substring(7).trim();

				// Handle non-ORF names
				if(!gene.startsWith("Q") && !gene.startsWith("Y"))
				{
					String orfName = nonOrfMap.get(gene);
					if(orfName == null || orfName.equals(""))
					{
						throw new IOException(id + " mapping is " + gene);
					}
					else
					{
						gene = orfName;
					}
				}

				geneList.add(gene);
			}
			// This monomer is a modified form of another monomer so
			// we use that unmodified monomer's gene mapping
			else if(line.startsWith("UNMODIFIED-FORM - "))
			{
				// We only need the unmodified form if there aren't genes listed
				// for this entry
				if(!unmodified.equals("") && geneList.size() == 0 && isMonomer)
					throw new IOException("Two unmodified forms found for a single entry "
							+ id);

				// Remove "UNMODIFIED-FORM - "
				unmodified = line.substring(18).trim();
			}
			// We reached the end of the entry so store the mapping if it's
			// a monomer and complete
			else if(line.startsWith("//"))
			{
				if(isMonomer)
				{
					if(id.equals(""))
						throw new IOException("Cannot add entry, no id");

					if(geneList.size() == 0)
					{
						// See if this monomer is the modified form
						// of another
						if(!unmodified.equals(""))
						{
							// The other monomer may not have been added yet
							// so just track the mapping now and add
							// the genes later
							modifiedMap.put(id, unmodified);
						}
						// See if we can use the id as the gene name
						else if(id.startsWith("Y") && id.endsWith("-MONOMER"))
						{
							geneList.add(id.substring(0, id.indexOf("-MONOMER")));
							monoCompMap.put(id, geneList);
						}
						else
						{
							throw new IOException("No genes provided for monomer " + id);
						}
					}
					// The gene list is not empty
					else
					{
						monoCompMap.put(id, geneList);
					}
				}

				// Reset all entry variables
				id = "";
				unmodified = "";
				isMonomer = false;
				geneList = new HashSet<String>();
			}
		}
		while((line = reader.readLine()) != null);
		reader.close();

		// Before moving on to the complexes, try adding
		// the modified monomers
		for(Map.Entry<String, String> pair : modifiedMap.entrySet())
		{
			id = pair.getKey();
			unmodified = pair.getValue();

			// Get the gene list for the unmodified monomer
			geneList = monoCompMap.get(unmodified);

			if(geneList == null || geneList.size() == 0)
			{
				throw new IOException("Could not find gene mapping for unmodifed " +
						"monomer " + unmodified + " of monomer " + id);
			}

			monoCompMap.put(id, geneList);
		}






		// Now pass through the file again to load the complexes
		reader = new BufferedReader(new FileReader(sgdDir + data));

		// Skip the heading
		while((line = reader.readLine()).startsWith("#"));

		// Reset the entry variables
		id = "";
		boolean isComplex = false;
		geneList = new HashSet<String>();
		// Create the complex map using the previously created monomer map
		do
		{
			line = line.trim();
			// Store the id
			if(line.startsWith("UNIQUE-ID - "))
			{
				if(!id.equals(""))
				{
					throw new IOException("Two ids found for a single entry: " +
							id + " " + line.substring(12));
				}

				// Remove "UNIQUE-ID - "
				id = line.substring(12).trim();
			}
			// See if this is a complex or something else
			// The "Co-E-Clth" complexes are not added here because they
			// do not map to monomers
			else if(line.contains("TYPES - Protein-Complexes"))
			{
				isComplex = true;
			}
			// Find the monomer name and retrieve the associated gene.
			// There could be one or more monomers per entry.  Some complexes
			// contain mulitple copies of a single monomer but we only add
			// each one once
			else if(line.startsWith("COMPONENTS - "))
			{
				// Remove "COMPONENTS - "
				String monomer = line.substring(13).trim();

				Set<String> genes = monoCompMap.get(monomer);

				if(genes == null || genes.size() == 0)
				{
					throw new IOException("No genes for monomer " + monomer +
							" in complex " + id);
				}

				geneList.addAll(genes);
			}
			// We reached the end of the entry so store the mapping if it's
			// a complex and complete
			else if(line.startsWith("//"))
			{
				if(isComplex)
				{
					if(id.equals("") || geneList.size() == 0)
						throw new IOException("Cannot add entry " + id);

					monoCompMap.put(id, geneList);
				}

				// Reset all entry variables
				id = "";
				isComplex = false;
				geneList = new HashSet<String>();
			}
		}
		while((line = reader.readLine()) != null);

		reader.close();
	}




	/**
	 * Loads a mapping enzyme reactions to monomers.  Requires
	 * that the monomer mapping has already been loaded
	 * Uses the file enzrxns.dat
	 *
	 *@throws IOException
	 */
	public void loadEnzymes() throws IOException
	{
		String data = "enzrxns.dat";
		if(monoCompMap == null)
		{
			loadMonomersComplexes();
		}


		BufferedReader reader = new BufferedReader(new FileReader(sgdDir + data));
		enzMap = new HashMap<String,String>();

		// Skip the heading
		String line;
		while((line = reader.readLine()).startsWith("#"));

		String id = "", monomer = "";
		boolean isEnzyme = false;
		do
		{
			line = line.trim();
			// Store the id
			if(line.startsWith("UNIQUE-ID - "))
			{
				if(!id.equals(""))
					throw new IOException("Two ids found for a single entry: " + id +
							" and " + line.substring(12).trim());

				// Remove "UNIQUE-ID - "
				id = line.substring(12).trim();
			}
			// See if this is an enzyme reaction or something else
			else if(line.contains("TYPES - Enzymatic-Reactions"))
			{
				isEnzyme = true;
			}
			// Store the monomer.  There should only be one per entry
			// and it must match a monomer in the monomer map
			else if(line.startsWith("ENZYME - "))
			{
				if(!monomer.equals(""))
					throw new IOException("Two monomers found for a single entry: " + id);

				// Remove "ENZYME - "
				monomer = line.substring(9).trim();

				// Ensure the monomer is in the map
				if(!monoCompMap.containsKey(monomer))
				{
					throw new IOException(monomer + " is not a valid enzyme for " + id);
				}
			}
			// We reached the end of the entry so store the mapping if it's complete
			else if(line.startsWith("//"))
			{
				if(isEnzyme)
				{
					if(id.equals(""))
						throw new IOException("Cannot add entry, no id");

					// TODO remove once SGD addresses my concern
					if(id.equals("ENZRXN3O-946"))
					{
						System.out.println("Skipping id " + id + " because of " + 
								"bad formatting in the file " + data);
					}
					// TODO remove once SGD addresses my concern
					else if(id.equals("ENZRXN3O-1458"))
					{
						enzMap.put("ENZRXN3O-1458-PLACEHOLDER-1", "YER073W-MONOMER");
						enzMap.put("ENZRXN3O-1458-PLACEHOLDER-2", "YOR374W-MONOMER");
					}
					else if(monomer.equals(""))
					{
						//throw new IOException("No enzymes provided for reaction " + id);
						System.err.println("No enzymes provided for reaction " + id);
					}
					else
					{
						enzMap.put(id, monomer);
					}
				}

				// Reset all entry variables
				id = "";
				monomer = "";
				isEnzyme = false;
			}
		}
		while((line = reader.readLine()) != null);
		reader.close();
	}



	/**
	 * Loads a mapping of reactions to enzymatic reactions.  Requires
	 * that the enzymatic reaction mapping has already been loaded
	 * Uses the file reactions.dat
	 *
	 *@throws IOException
	 */
	public void loadReactions() throws IOException
	{
		String data = "reactions.dat";
		if(enzMap == null || monoCompMap == null)
		{
			loadEnzymes();
		}


		BufferedReader reader = new BufferedReader(new FileReader(sgdDir + data));
		reacMap = new HashMap<String,HashSet<String>>();

		// Skip the heading
		String line;
		while((line = reader.readLine()).startsWith("#"));

		String id = "";
		HashSet<String> enzReacs = new HashSet<String>();
		do
		{
			// Unlike the other entities, we do not check the type
			// for reactions.  Rather we use the presence of
			// the ENZYMATIC-REACTION attribute to see if the
			// entry should be added to the map

			line = line.trim();
			// Store the id
			if(line.startsWith("UNIQUE-ID - "))
			{
				if(!id.equals(""))
					throw new IOException("Two ids found for a single entry: " + id +
							" and " + line.substring(12).trim());

				// Remove "UNIQUE-ID - "
				id = line.substring(12).trim();
			}
			// Store the enzymatic reaction.  It must appears in the previously
			// constructed map of enzymatic reactions
			else if(line.startsWith("ENZYMATIC-REACTION - "))
			{
				// Remove "ENZYMATIC-REACTION - "
				String enzReac = line.substring(21).trim();

				// TODO remove once SGD curators respond
				// These weren't added to enzMap
				if(!enzReac.equals("ENZRXN3O-946") && !enzReac.equals("ENZRXN3O-1458"))
				{

					// Ensure the enzymatic reaction is in the map
					if(!enzMap.containsKey(enzReac))
					{
						throw new IOException(enzReac + " is not a " +
								"valid enzymatic reaction for " + id);
					}

					enzReacs.add(enzReac);
				}
			}
			// We reached the end of the entry so store the mapping if it's complete
			else if(line.startsWith("//"))
			{
				// TODO remove once SGD curators respond
				if(id.equals("ALDHDEHYDROG-RXN"))
				{
					enzReacs.add("ENZRXN3O-1458-PLACEHOLDER-1");
					enzReacs.add("ENZRXN3O-1458-PLACEHOLDER-2");
				}

				if(enzReacs.size() > 0)
				{
					if(id.equals(""))
						throw new IOException("Cannot add entry, no id");

					reacMap.put(id, enzReacs);
				}
//				else
//				System.err.println(id);

				// Reset all entry variables
				id = "";
				enzReacs = new HashSet<String>();
			}
		}
		while((line = reader.readLine()) != null);
		reader.close();
	}


	/**
	 * Loads the sgdPathways.  Requires
	 * that the reaction, enzymatic reaction, and monomer mappings have
	 * already been loaded.
	 * Uses the file pathways.dat
	 *
	 *@throws IOException
	 */
	public void loadSgdPathways() throws IOException
	{
		String data = "pathways.dat";
		if(reacMap == null || enzMap == null || monoCompMap == null)
		{
			loadReactions();
		}


		BufferedReader reader = new BufferedReader(new FileReader(sgdDir + data));
		sgdPathways = new HashSet<SGDPathway>();

		// Skip the heading
		String line;
		while((line = reader.readLine()).startsWith("#"));

		String id = "", name = "";
		// A bad path is any path for which one of the reactions
		// is not present in the reaction map.  At this point
		// we choose to completely discard the pathway if this is the case
		boolean badPath = false;
		SGDPathway curPath = new SGDPathway();
		do
		{
			line = line.trim();
			// Store the id
			if(line.startsWith("UNIQUE-ID - "))
			{
				if(!id.equals(""))
					throw new IOException("Two ids found for a single entry: " + id +
							" and " + line.substring(12).trim());

				// Remove "UNIQUE-ID - "
				id = line.substring(12).trim();
			}
			else if(line.startsWith("COMMON-NAME - "))
			{
				if(!name.equals(""))
					throw new IOException("Two common names found for a single " +
							"entry: " + name +
							" and " + line.substring(14).trim());

				// Remove "COMMON-NAME - "
				name = line.substring(14).trim();
			}
			// Add the reaction to the pathway
			else if(line.startsWith("REACTION-LAYOUT - "))
			{
				// Remove "REACTION-LAYOUT - (" and the trailing ")"
				line = line.substring(19, line.length()-1).trim();

				// The line will be of the form:
				// <reaction> (:LEFT-PRIMARIES <pri> <pri>) (:DIRECTION <dir>) (:RIGHT-PRIMARIES <pri>)

				// Now use the (: as a delimiter to split the the components
				String[] parts = line.split("\\(:");

				// parts[0] is the reaction
				String reaction = parts[0].trim();

				// Reactions may be surrounded by | |
				// Strip these if present
				if(reaction.startsWith("|") && reaction.endsWith("|"))
				{					
					reaction = reaction.substring(1, reaction.length() - 1);
				}

				// Ensure the reaction is in the map
				if(reacMap.containsKey(reaction))
				{
					// Use the reaction to get the set of genes on the pathway
					// edge from the other mappings
					HashSet<String> enzReacs = reacMap.get(reaction);
					HashSet<String> genes = new HashSet<String>();

					for(String enz : enzReacs)
					{
						// Get the monomer/complex and then use that
						// to get the set of genes
						genes.addAll(monoCompMap.get(enzMap.get(enz)));
					}

					// Get the direction from parts[2] so that the source and
					// target can be determined
					// If the direction is nil skip this pathway
					// Remove the trailing )
					parts[2] = parts[2].trim();
					parts[2] = parts[2].substring(0, parts[2].length() - 1);
					String[] dir = parts[2].trim().split(" ");
//					if(dir.length < 2)
//					{
//					dir = null;
//					}
					if(dir[1].equals(":L2R") || dir[1].equals(":R2L"))
					{
						// parts[1] is of the form LEFT-PRIMARIES <pri> <pri>)
						// Remove the trailing )
						parts[1] = parts[1].trim();
						parts[1] = parts[1].substring(0, parts[1].length() - 1);
						String[] leftPri = parts[1].trim().split(" ");
						// Make sure leftPri[0] is LEFT-PRIMARIES
						if(!leftPri[0].equals("LEFT-PRIMARIES"))
						{
							throw new IOException("Parse error " + id);
						}

						// Add the actual primaries to a set to eliminate
						// duplicates
						HashSet<String> leftPriSet = new HashSet<String>();
						for(int l = 1; l < leftPri.length; l++)
						{
							leftPriSet.add(leftPri[l].trim());
						}


						// parts[3] is of the form RIGHT-PRIMARIES <pri> <pri>)
						// Remove the trailing )
						parts[3] = parts[3].trim();
						parts[3] = parts[3].substring(0, parts[3].length() - 1);
						String[] rightPri = parts[3].trim().split(" ");
						// Make sure rightPri[0] is RIGHT-PRIMARIES
						if(!rightPri[0].equals("RIGHT-PRIMARIES"))
						{
							throw new IOException("Parse error " + id);
						}

						// Add the actual primaries to a set to eliminate
						// duplicates
						HashSet<String> rightPriSet = new HashSet<String>();
						for(int r = 1;r < rightPri.length; r++)
						{
							rightPriSet.add(rightPri[r].trim());
						}


						// Add this reaction to the pathway
						for(String lPri : leftPriSet)
						{
							for(String rPri : rightPriSet)
							{
								if(dir[1].equals(":L2R"))
								{
									curPath.addEdge(lPri, rPri, genes);
								}
								else if(dir[1].equals(":R2L"))
								{
									curPath.addEdge(rPri, lPri, genes);
								}
								else
								{
									// Should be no way to reach this
									throw new IOException("Parse error " + id);
								}
							}
						}
					}
					else
					{
						// The direction wasn't :L2R or :R2L
						badPath = true;
					}

				}
				// If any reaction is not present we will not add this pathway
				else
				{
					badPath = true;
				}
			}
			// We reached the end of the entry so store the pathway if it's complete
			else if(line.startsWith("//"))
			{
				// Don't add the pathway if no edges were found
				// or if it is a bad path
				if(!badPath && curPath.numEdges() > 0)
				{
					if(id.equals(""))
						throw new IOException("Cannot add entry, no id");

					if(name.equals(""))
						throw new IOException("Cannot add entry, no name");

					curPath.setId(id);
					curPath.setName(name);
					sgdPathways.add(curPath);
				}

				// Reset all entry variables
				id = "";
				name = "";
				badPath = false;
				curPath = new SGDPathway();
			}
		}
		while((line = reader.readLine()) != null);
		reader.close();
	}


	/**
	 * Loads metabolic pathways from the SGD database, then merges them
	 * into one unified pathway Graph.  Vertices in the different
	 * original pathways are still distinct nodes in the merged Graph.
	 * @return Graph a reference to the signaling pathways
	 */
	public Graph loadMetabolicPathways() throws IOException
	{
		loadSgdPathways();

		metabolicPathways = new Graph();
		for(SGDPathway sgdPath : sgdPathways)
		{
			Graph genePath = sgdPath.toGenePathway(pathCounter);
			pathCounter++;
			metabolicPathways.mergeGraph(genePath);
		}		
		return metabolicPathways;
	}

	/**
	 * Loads the MAPK and Phosphatidylinositol signaling pathways
	 * for yeast from the KEGG database
	 * @throws IOException
	 */
	public void loadKeggPathways() throws IOException
	{
		String[] files = {"sce04011.xml", "sce04070.xml"};
		//String[] files = {"test.xml"};
		//String[] files = {"sce04070.xml"};


		// Split the relation subtypes into directed and undirected keywords
		String[] dirRels = {"activation", "inhibition", "expression", 
				"repression", "indirect effect", "phosphorylation", 
				"dephosphorylation", "glycosylation", "ubiquitination",
		"methylation"};
		String[] undirRels = {"compound", "hidden compound",
				"binding/association", "dissociation"};
		// TODO place these into one of the categories once they are encountered
		String[] unknownRels = {"state change", "missing interaction"};

		HashSet<String> dirRelSet = new HashSet<String>();
		for(String rel : dirRels)
		{
			dirRelSet.add(rel);
		}

		HashSet<String> undirRelSet = new HashSet<String>();
		for(String rel : undirRels)
		{
			undirRelSet.add(rel);
		}

		keggPathways = new HashSet<Graph>();
		for(int f = 0; f < files.length; f++)
		{
			System.out.println("Loading " + files[f]);

			BufferedReader reader = new BufferedReader(new FileReader(keggDir + files[f]));
			// Store all edges in the pathway.  Since a graph and pathway are both
			// represented as a collection of edges, use a Graph object.
			// All vertices will be added as sources because subpathways can
			// start at any point
			Graph pathway = new Graph();
			// Map of entry ids to a list of vertices
			HashMap<String, ArrayList<Vertex>> entryMap = new HashMap<String, ArrayList<Vertex>>();
			// Set of all vertices in the pathway
			HashSet<Vertex> allVertices = new HashSet<Vertex>();

			int entries = 0, relations = 0;

			String line;
			while((line = reader.readLine()) != null)
			{
				line = line.trim();

				// Only entry and relation lines are needed in the XML file.
				// It is assumed that all entries come before any relations.
				if(line.startsWith("<entry"))
				{
					// We are only interested in entries that are yeast genes
					// (not sgdPathways or KEGG ortholog identifiers)
					if(line.contains("type=\"gene\""))
					{
						// Get the id
						int start = line.indexOf("id=\"") + 4;
						int end = line.indexOf("\"", start);
						String id = line.substring(start, end);

						// Force all entry ids to be unique
						if(entryMap.containsKey(id))
						{
							throw new IllegalStateException("Entry id " + id +
							" already exists");
						}

						// Get the gene name(s)
						start = line.indexOf("name=\"") + 6;
						end = line.indexOf("\"", start);
						String nameStr = line.substring(start, end);
						String[] names = nameStr.split(" ");

						// Create the list of vertices for this id
						ArrayList<Vertex> vertices = new ArrayList<Vertex>();
						for(int n = 0; n < names.length; n++)
						{
							if(names[n].startsWith("sce:"))
							{
								// TODO need a good way to handle multiple ids per
								// gene
//								String standName = standardName(names[n].substring(4));
//								Vertex v = Vertex.createVertex(standName + "_id" + id + "pa" + pathCounter);
								Vertex v = Vertex.createVertex(names[n].substring(4)
										+ "_id" + id + "pa" + pathCounter);
								vertices.add(v);
								allVertices.add(v);
							}
							else
							{
								throw new IllegalStateException(names[n] + " is not a valid gene name");
							}
						}

						entryMap.put(id, vertices);
						entries++;
					}
					// We also need to handle group of genes that operate together
					// as a complex
					else if(line.contains("type=\"group\""))
					{
						// Get the id
						int start = line.indexOf("id=\"") + 4;
						int end = line.indexOf("\"", start);
						String id = line.substring(start, end);

						// Add all genes that belong to the group
						ArrayList<Vertex> vertices = new ArrayList<Vertex>();
						// Find the ids of components in the group
						boolean groupDone = false;
						// Check !groupDone first so that it short circuits
						// appropriately
						while(!groupDone && (line = reader.readLine()) != null)
						{
							line = line.trim();
							if(line.equals("</entry>"))
							{
								groupDone = true;
							}
							else if(line.startsWith("<component"))
							{
								// Get the component id
								start = line.indexOf("id=\"") + 4;
								end = line.indexOf("\"", start);
								String comp = line.substring(start, end);

								// Add all genes for this id
								vertices.addAll(entryMap.get(comp));
							}
						} // end search for component ids
						entryMap.put(id, vertices);

						entries++;
					} // end group parsing
				} // end entry parsing
				else if(line.startsWith("<relation"))
				{
					// First get the entry ids
					int start = line.indexOf("entry1=\"") + 8;
					int end = line.indexOf("\"", start);
					String id1 = line.substring(start, end);

					start = line.indexOf("entry2=\"") + 8;
					end = line.indexOf("\"", start);
					String id2 = line.substring(start, end);

					// Ensure the entry ids correspond to any genes
					if(!entryMap.containsKey(id1))
					{
						throw new IllegalStateException("Could not find id " + id1);
					}
					if(!entryMap.containsKey(id2))
					{
						throw new IllegalStateException("Could not find id " + id2);
					}


					// Then determine what the relationship is
					// http://www.genome.jp/kegg/xml/docs/ gives possible
					// relation subtypes
					boolean relDone = false, isDir = false, isUndir = false;
					// Must check !relDone first so that the conditional
					// short circuits.  Otherwise a line will be skipped
					while(!relDone && (line = reader.readLine()) != null)
					{
						line = line.trim();

						if(line.equals("</relation>"))
						{
							relDone = true;
						}
						else if(line.startsWith("<subtype"))
						{
							start = line.indexOf("name=\"") + 6;
							end = line.indexOf("\"", start);
							String name = line.substring(start, end);

							isDir |= dirRelSet.contains(name);
							isUndir |= undirRelSet.contains(name);
						}
					}

					// Make sure we didn't find a relation that is directed and
					// undirected
					if(isDir && isUndir)
					{
						throw new IllegalStateException("Relation " + id1 +
								":" + id2 + " is directed and undirected");
					}

					// Iterate through the lists of vertices corresponding to each id
					// adding pairwise edges
					for(Vertex v1 : entryMap.get(id1))
					{
						for(Vertex v2 : entryMap.get(id2))
						{
							// Don't allow self-interactions
							if(!v1.getName().equals(v2.getName()))
							{
								// undirected edge is default

								if(isDir)
								{
									pathway.addEdge(new DirEdge(v1, v2));
								}
								else
								{
									pathway.addEdge(new UndirEdge(v1, v2));
								}
							}
						}
					}

					relations++;
				} // end relation parsing
			} // end reading file

			// Add all vertices as sources
			pathway.addSources(allVertices);
			keggPathways.add(pathway);
			pathCounter++;

			reader.close();

			System.out.println("Entries: " + entries + "  Relations: " + relations);


			// TODO remove testing
			// Find paths of depth 5 regardless of whether
			// or not they terminate at a target
//			ArrayList<Path> paths = pathway.findPaths(3, true);
//			for(Path p : paths)
//			{
//			System.out.println(p);
//			}

//			PrintWriter writer = new PrintWriter(new FileWriter(keggDir + files[f] + ".out"));
//			pathway.writeEdges(writer);
//			writer.close();

		} // end iterating through all files

	}


	/**
	 * Load the Filamentous Growth, HOG, and Pheromone Signaling pathways
	 * from the Science Signaling Database of Cell Signaling.
	 * The data files were hand-curated from the website images.
	 *
	 * @throws IOException
	 */
	public void loadSciPathways() throws IOException
	{
		String[] files = {"FilamentousGrowth.txt", "HOG.txt", "Pheromone.txt"};
//		String[] files = {"HOG.txt"};


		sciPathways = new HashSet<Graph>();
		for(int f = 0; f < files.length; f++)
		{
			System.out.println("Loading " + files[f]);

			BufferedReader reader = new BufferedReader(new FileReader(sciSigDir + files[f]));
			// Store all edges in the pathway.  Since a graph and pathway are both
			// represented as a collection of edges, use a Graph object.
			// All vertices will be added as sources because subpathways can
			// start at any point
			Graph pathway = new Graph();
			// Set of all vertices in the pathway
			HashSet<Vertex> allVertices = new HashSet<Vertex>();

			String line;
			while((line = reader.readLine()) != null)
			{
				String[] parts = line.trim().split("\t");
				// parts[1] contains the interaction type
				// -> and -| are directed
				// -- is undirected

				Vertex v1 = Vertex.createVertex(parts[0].trim() + "_pa" + pathCounter);
				Vertex v2 = Vertex.createVertex(parts[2].trim() + "_pa" + pathCounter);
//				Vertex v1 = Vertex.createVertex(standardName(parts[0].trim()));
//				Vertex v2 = Vertex.createVertex(standardName(parts[2].trim()));


				if(parts[1].trim().equals("--"))
				{
					pathway.addEdge(new UndirEdge(v1,v2));
				}
				else if(parts[1].trim().equals("->") ||
						parts[1].trim().equals("-|"))
				{
					pathway.addEdge(new DirEdge(v1,v2));
				}
				else
				{
					throw new IllegalStateException(parts[1] + " is not a valid edge type");
				}

				allVertices.add(v1);
				allVertices.add(v2);
			} // end reading file
			reader.close();

			pathway.addSources(allVertices);
			sciPathways.add(pathway);
			pathCounter++;

			// TODO remove testing
			// Find paths of depth 5 regardless of whether
			// or not they terminate at a target
//			ArrayList<Path> paths = pathway.findPaths(3, true);
//			for(Path p : paths)
//			{
//			System.out.println(p);
//			}

		} // end loop through all files
	}


	/**
	 * Loads signaling pathways from the KEGG and Science
	 * Signaling databases, then merges them into one unified
	 * pathway Graph.  Vertices in the different original pathways
	 * are still distinct nodes in the merged Graph.
	 * @return Graph a reference to the signaling pathways
	 */
	public Graph loadSignalPathways() throws IOException
	{
		loadKeggPathways();
		loadSciPathways();

		signalPathways = new Graph();
		for(Graph keggPath : keggPathways)
		{
			signalPathways.mergeGraph(keggPath);
		}
		for(Graph sciPath : sciPathways)
		{
			signalPathways.mergeGraph(sciPath);
		}

		// TODO remove testing
//		ArrayList<Path> paths = signalPathways.findPaths(3, true);
//		for(Path p : paths)
//		{
//		System.out.println(p);
//		}

//		System.out.println("There are " + signalPathways.getSources().size() +
//		" vertices in the signaling network");

		return signalPathways;
	}

	// TODO move to a util file
	// ORF to standard
	private static HashMap<String, String> synonyms;
	public String standardName(String orfName) throws IOException
	{
		if(synonyms == null)
		{
			synonyms = loadSynonyms();
		}

		if(synonyms.containsKey(orfName))
		{
			return synonyms.get(orfName);
		}
		else
		{
			return orfName;
		}
	}

	// TODO move to a util file
	public HashMap<String, String> loadSynonyms() throws IOException
	{
		HashMap<String, String> synonyms = new HashMap<String, String>();
		BufferedReader synReader = new BufferedReader(new FileReader(keggDir + "SGD_orfToStandard.txt"));

		String line;
		while((line = synReader.readLine()) != null)
		{
			String[] parts = line.split("\t");
			// parts[] is ORF name, standard name
			// The best match comes first in the list, and some standard names
			// will appear later as aliases
			if(!synonyms.containsKey(parts[0]))
			{
				synonyms.put(parts[0], parts[1]);
			}
		}
		synReader.close();

		return synonyms;
	}

	// TODO merge with the version that loads all pathways
	/**
	 * Load the Filamentous Growth, HOG, and/or Pheromone Signaling pathways
	 * from the Science Signaling Database of Cell Signaling.
	 * The data files were hand-curated from the website images.
	 *
	 * @throws IOException
	 */
	public Graph loadSciPathway(String file) throws IOException
	{
//		String[] files = {"FilamentousGrowth.txt", "HOG.txt", "Pheromone.txt"};


		sciPathways = new HashSet<Graph>();

		System.out.println("Loading " + file);

		BufferedReader reader = new BufferedReader(new FileReader(sciSigDir + file));
		// Store all edges in the pathway.  Since a graph and pathway are both
		// represented as a collection of edges, use a Graph object.
		// All vertices will be added as sources because subpathways can
		// start at any point
		Graph pathway = new Graph();
		// Set of all vertices in the pathway
		HashSet<Vertex> allVertices = new HashSet<Vertex>();

		String line;
		while((line = reader.readLine()) != null)
		{
			String[] parts = line.trim().split("\t");
			// parts[1] contains the interaction type
			// -> and -| are directed
			// -- is undirected

			Vertex v1 = Vertex.createVertex(parts[0].trim() + "_pa" + pathCounter);
			Vertex v2 = Vertex.createVertex(parts[2].trim() + "_pa" + pathCounter);
//			Vertex v1 = Vertex.createVertex(standardName(parts[0].trim()));
//			Vertex v2 = Vertex.createVertex(standardName(parts[2].trim()));


			if(parts[1].trim().equals("--"))
			{
				pathway.addEdge(new UndirEdge(v1,v2));
			}
			else if(parts[1].trim().equals("->") ||
					parts[1].trim().equals("-|"))
			{
				pathway.addEdge(new DirEdge(v1,v2));
			}
			else
			{
				throw new IllegalStateException(parts[1] + " is not a valid edge type");
			}

			allVertices.add(v1);
			allVertices.add(v2);
		} // end reading file
		reader.close();

		pathway.addSources(allVertices);
//		sciPathways.add(pathway);
		pathCounter++;

		return pathway;
	}
	
	
	// TODO Merge with the version that loads all pathways
	/**
	 * Loads the MAPK and Phosphatidylinositol signaling pathways
	 * for yeast from the KEGG database
	 * @throws IOException
	 */
	public Graph loadKeggPathway(String file) throws IOException
	{
//		String[] files = {"sce04011.xml", "sce04070.xml"};


		// Split the relation subtypes into directed and undirected keywords
		String[] dirRels = {"activation", "inhibition", "expression", 
				"repression", "indirect effect", "phosphorylation", 
				"dephosphorylation", "glycosylation", "ubiquitination",
		"methylation"};
		String[] undirRels = {"compound", "hidden compound",
				"binding/association", "dissociation"};
		// TODO place these into one of the categories once they are encountered
		String[] unknownRels = {"state change", "missing interaction"};

		HashSet<String> dirRelSet = new HashSet<String>();
		for(String rel : dirRels)
		{
			dirRelSet.add(rel);
		}

		HashSet<String> undirRelSet = new HashSet<String>();
		for(String rel : undirRels)
		{
			undirRelSet.add(rel);
		}

		keggPathways = new HashSet<Graph>();

			System.out.println("Loading " + file);

			BufferedReader reader = new BufferedReader(new FileReader(keggDir + file));
			// Store all edges in the pathway.  Since a graph and pathway are both
			// represented as a collection of edges, use a Graph object.
			// All vertices will be added as sources because subpathways can
			// start at any point
			Graph pathway = new Graph();
			// Map of entry ids to a list of vertices
			HashMap<String, ArrayList<Vertex>> entryMap = new HashMap<String, ArrayList<Vertex>>();
			// Set of all vertices in the pathway
			HashSet<Vertex> allVertices = new HashSet<Vertex>();

			int entries = 0, relations = 0;

			String line;
			while((line = reader.readLine()) != null)
			{
				line = line.trim();

				// Only entry and relation lines are needed in the XML file.
				// It is assumed that all entries come before any relations.
				if(line.startsWith("<entry"))
				{
					// We are only interested in entries that are yeast genes
					// (not sgdPathways or KEGG ortholog identifiers)
					if(line.contains("type=\"gene\""))
					{
						// Get the id
						int start = line.indexOf("id=\"") + 4;
						int end = line.indexOf("\"", start);
						String id = line.substring(start, end);

						// Force all entry ids to be unique
						if(entryMap.containsKey(id))
						{
							throw new IllegalStateException("Entry id " + id +
							" already exists");
						}

						// Get the gene name(s)
						start = line.indexOf("name=\"") + 6;
						end = line.indexOf("\"", start);
						String nameStr = line.substring(start, end);
						String[] names = nameStr.split(" ");

						// Create the list of vertices for this id
						ArrayList<Vertex> vertices = new ArrayList<Vertex>();
						for(int n = 0; n < names.length; n++)
						{
							if(names[n].startsWith("sce:"))
							{
								// TODO need a good way to handle multiple ids per
								// gene
//								String standName = standardName(names[n].substring(4));
//								Vertex v = Vertex.createVertex(standName + "_id" + id + "pa" + pathCounter);
								
								// General way
//								Vertex v = Vertex.createVertex(names[n].substring(4)
//										+ "_id" + id + "pa" + pathCounter);
								Vertex v = Vertex.createVertex(names[n].substring(4)
										+ "_pa" + pathCounter);
								
								vertices.add(v);
								allVertices.add(v);
							}
							else
							{
								throw new IllegalStateException(names[n] + " is not a valid gene name");
							}
						}

						entryMap.put(id, vertices);
						entries++;
					}
					// We also need to handle group of genes that operate together
					// as a complex
					else if(line.contains("type=\"group\""))
					{
						// Get the id
						int start = line.indexOf("id=\"") + 4;
						int end = line.indexOf("\"", start);
						String id = line.substring(start, end);

						// Add all genes that belong to the group
						ArrayList<Vertex> vertices = new ArrayList<Vertex>();
						// Find the ids of components in the group
						boolean groupDone = false;
						// Check !groupDone first so that it short circuits
						// appropriately
						while(!groupDone && (line = reader.readLine()) != null)
						{
							line = line.trim();
							if(line.equals("</entry>"))
							{
								groupDone = true;
							}
							else if(line.startsWith("<component"))
							{
								// Get the component id
								start = line.indexOf("id=\"") + 4;
								end = line.indexOf("\"", start);
								String comp = line.substring(start, end);

								// Add all genes for this id
								vertices.addAll(entryMap.get(comp));
							}
						} // end search for component ids
						entryMap.put(id, vertices);

						entries++;
					} // end group parsing
				} // end entry parsing
				else if(line.startsWith("<relation"))
				{
					// First get the entry ids
					int start = line.indexOf("entry1=\"") + 8;
					int end = line.indexOf("\"", start);
					String id1 = line.substring(start, end);

					start = line.indexOf("entry2=\"") + 8;
					end = line.indexOf("\"", start);
					String id2 = line.substring(start, end);

					// Ensure the entry ids correspond to any genes
					if(!entryMap.containsKey(id1))
					{
						throw new IllegalStateException("Could not find id " + id1);
					}
					if(!entryMap.containsKey(id2))
					{
						throw new IllegalStateException("Could not find id " + id2);
					}


					// Then determine what the relationship is
					// http://www.genome.jp/kegg/xml/docs/ gives possible
					// relation subtypes
					boolean relDone = false, isDir = false, isUndir = false;
					// Must check !relDone first so that the conditional
					// short circuits.  Otherwise a line will be skipped
					while(!relDone && (line = reader.readLine()) != null)
					{
						line = line.trim();

						if(line.equals("</relation>"))
						{
							relDone = true;
						}
						else if(line.startsWith("<subtype"))
						{
							start = line.indexOf("name=\"") + 6;
							end = line.indexOf("\"", start);
							String name = line.substring(start, end);

							isDir |= dirRelSet.contains(name);
							isUndir |= undirRelSet.contains(name);
						}
					}

					// Make sure we didn't find a relation that is directed and
					// undirected
					if(isDir && isUndir)
					{
						throw new IllegalStateException("Relation " + id1 +
								":" + id2 + " is directed and undirected");
					}

					// Iterate through the lists of vertices corresponding to each id
					// adding pairwise edges
					for(Vertex v1 : entryMap.get(id1))
					{
						for(Vertex v2 : entryMap.get(id2))
						{
							// Don't allow self-interactions
							if(!v1.getName().equals(v2.getName()))
							{
								// undirected edge is default

								if(isDir)
								{
									boolean exists = false;
									for(DirEdge e : pathway.getDirEdges())
									{
										if(e.getSource().getName().equals(v1.getName()) &&
												e.getTarget().getName().equals(v2.getName()))
										{
											exists = true;
										}
									}
									if(!exists)
									{
										pathway.addEdge(new DirEdge(v1, v2));
									}
								}
								else
								{
									boolean exists = false;
									for(UndirEdge e : pathway.getUndirEdges())
									{
										if((e.getSource().getName().equals(v1.getName()) &&
												e.getTarget().getName().equals(v2.getName())) ||
											(e.getSource().getName().equals(v2.getName()) &&
												e.getTarget().getName().equals(v1.getName())))
										{
											exists = true;
										}
									}
									if(!exists)
									{
										pathway.addEdge(new UndirEdge(v1, v2));
									}
								}
							}
						}
					}

					relations++;
				} // end relation parsing
			} // end reading file

			// Add all vertices as sources
			pathway.addSources(allVertices);
//			keggPathways.add(pathway);
			pathCounter++;

			reader.close();

			System.out.println("Entries: " + entries + "  Relations: " + relations);
			
			return pathway;

	}
}
