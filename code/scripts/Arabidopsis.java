package scripts;

import java.io.*;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;

import util.MapUtil;
import util.StatUtil;

public class Arabidopsis
{

	public static final String ARAB_DIR = "../../arabidopsis/";
	public static final String SYN_FILE = ARAB_DIR + "arabidopsis_genelist.txt";
	/** Map from symbols to ORF names */
	public static HashMap<String, String> synMap = null;
	/** Map from ORF names to symbols */
	public static HashMap<String, String> revSynMap = null;
	
	public static void main(String[] args)
	{
//		prepWangData(ARAB_DIR + "Wang2010/Wang2010Expression_ProbeIds.txt",
//				ARAB_DIR + "Wang2010/GPL198_probeToOrf.txt",
//				ARAB_DIR + "Wang2010/Wang2010Expression_Orfs.txt");
		
//		prepBindingGrid(ARAB_DIR + "proteinDna/arabidopsis_tf_gene_agris.txt",
//				ARAB_DIR + "proteinDna/arabidopsis_tf_gene_agris_priors.txt");
	
//		sdrem092211();
//		sdrem052112();
		
		// Use the real PPI network instead of the network with fake edges for TFs
//		checkSources(ARAB_DIR + "hpaEffectorsGrouped.txt",
//				ARAB_DIR + "ppi/allPPI.txt");
	}

	public static void prepWangData(String exprFile, String platFile, String outFile)
	{
		try
		{
			HashMap<String, String> probeMap = MapUtil.loadMap(platFile, 0, 1, true);
			BufferedReader reader = new BufferedReader(new FileReader(exprFile));
			PrintWriter writer = new PrintWriter(new FileWriter(outFile));
			
			String line = reader.readLine();
			String[] header = line.split("\t");
			writer.print("Gene");
			for(int h = 1; h < header.length; h++)
			{
				writer.print("\t" + header[h]);
			}
			writer.println();
			
			int skipped = 0;
			while((line = reader.readLine()) != null)
			{
				String[] parts = line.split("\t");
				String probe = parts[0].toUpperCase();

				if(probeMap.containsKey(probe) && !probeMap.get(probe).equals(""))
				{
					writer.print(probeMap.get(probe));
				
					for(int p = 1; p < parts.length; p++)
					{
						writer.print("\t" + parts[p]);
					}
					writer.println();
				}
				else
				{
					skipped++;
				}
			}
			
			writer.close();
			reader.close();
			
			System.out.println("Skipping " + skipped + " bad rows");
		}
		catch(IOException e)
		{
			e.printStackTrace();
		}
	}

	public static void prepBindingGrid(String sparseFile, String gridFile)
	{
		try
		{
			BufferedReader reader = new BufferedReader(new FileReader(sparseFile));
			HashMap<String, Integer> colMap = new HashMap<String, Integer>();
			HashMap<String, Integer> rowMap = new HashMap<String, Integer>();
			
			int colCounter = 0, rowCounter = 0;
			
			// Skip header
			String line = reader.readLine();
			while((line = reader.readLine()) != null)
			{
				String[] parts = line.split("\t");
				for(String tf : parts[0].split("\\|"))
				{
					if(!colMap.containsKey(tf))
					{
						colMap.put(tf, colCounter);
						colCounter++;
					}
				}
				
				for(String gene : parts[1].split("\\|"))
				{
					if(!rowMap.containsKey(gene))
					{
						rowMap.put(gene, rowCounter);
						rowCounter++;
					}
				}
			}
			
			reader.close();
			
			boolean[][] binding = new boolean[colCounter][rowCounter];
			reader = new BufferedReader(new FileReader(sparseFile));
			
			// Skip header
			line = reader.readLine();
			while((line = reader.readLine()) != null)
			{
				String[] parts = line.split("\t");
				for(String tf : parts[0].split("\\|"))
				{
					for(String gene : parts[1].split("\\|"))
					{
						binding[colMap.get(tf)][rowMap.get(gene)] = true;
					}
				}
			}
			reader.close();
			
			String[] tfs = new String[colCounter];
			String[] genes = new String[rowCounter];
			
			for(String tf : colMap.keySet())
			{
				tfs[colMap.get(tf)] = tf;
			}
			
			for(String gene : rowMap.keySet())
			{
				genes[rowMap.get(gene)] = gene;
			}
			
			PrintWriter writer = new PrintWriter(new FileWriter(gridFile));
			
			writer.print("Gene");
			for(int t = 0; t < tfs.length; t++)
			{
				writer.print("\t" + tfs[t]);
			}
			writer.println();
			
			for(int g = 0; g < rowCounter; g++)
			{
				writer.print(genes[g]);
				for(int t = 0; t < colCounter; t++)
				{
					if(binding[t][g])
					{
						writer.print("\t0.5");
					}
					else
					{
						writer.print("\t0");
					}
				}
				writer.println();
			}
			
			writer.close();
		}
		catch(IOException e)
		{
			e.printStackTrace();
		}
	}
	/**
	 * Calculate the overlap between all model predictions (sources, targets, internal
	 * nodes) and the screen
	 * @param srcFile sources are pathogen proteins
	 * @param targFile
	 * @param nodeScoresFile
	 * @param nodeScoreThresh cutoff used to select significant internal nodes
	 * @param tfListFile a list of all proteins that are TFs (i.e. have predicted binding data)
	 * @param outFile
	 */
	public static void evaluatePredictions(String srcFile, String targFile,
			String nodeScoresFile, double nodeScoreThresh, String tfListFile,
			String outFile)
	{
		try
		{
			HashSet<String> srcList = MapUtil.loadSet(srcFile, false);
			// Load the targets but ignore their target weights
			Set<String> targList = MapUtil.loadMap(targFile, 0, 1, false).keySet();

			// Load the node score of all proteins, including sources and targets
			// Use the % of top paths for the node score
			HashMap<String, String> nodeScoreMap = MapUtil.loadMap(nodeScoresFile, 0, 5, true);

			// Determine which proteins have node scores >= the threshold
			HashSet<String> internalNodes = new HashSet<String>();
			for(String node : nodeScoreMap.keySet())
			{
				double score = Double.parseDouble(nodeScoreMap.get(node));
				if(score >= nodeScoreThresh)
				{
					internalNodes.add(node);
				}
			}

			// Sources, targets, and proteins with high node scores are all
			// considered to be important in the response
			HashSet<String> inModel = new HashSet<String>();
			inModel.addAll(srcList);
			inModel.addAll(targList);
			inModel.addAll(internalNodes);

			// Load the set of all known TFs
			HashSet<String> allTfs = MapUtil.loadSet(tfListFile, false);


			PrintWriter writer = new PrintWriter(new FileWriter(outFile));
			writer.println("Gene symbol\tGene id\tRole\tHave protein-DNA binding data");


			// Loop through all nodes in the model and write which sets they belong to
			for(String protein : inModel)
			{
				StringBuffer buf = new StringBuffer();
				String sym =getSymbol(protein);

				buf.append(sym).append("\t").append(protein).append("\t");

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
				if(allTfs.contains(protein))
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

			writer.close();
		}
		catch(IOException e)
		{
			e.printStackTrace();
		}
	}
	
	public static String getOrf(String symbolOrig)
	{
		if(synMap == null)
		{
			loadMap();
		}

		String symbol = symbolOrig.toUpperCase().replace(' ', '_');
		if(synMap.containsKey(symbol))
		{
			return synMap.get(symbol);
		}
		
		return symbol;
	}
	
	public static String getSymbol(String orf)
	{
		if(revSynMap == null)
		{
			loadMap();
		}
		
		String orfNoIso = stripIsoform(orf);
		if(revSynMap.containsKey(orfNoIso))
		{
			return revSynMap.get(orfNoIso);
		}
		
		return orfNoIso;
	}
	
	public static void loadMap()
	{
		try
		{
			synMap = new HashMap<String, String>();
			revSynMap = new HashMap<String, String>();
			
			BufferedReader reader = new BufferedReader(new FileReader(SYN_FILE));
			// Skip the header
			String line = reader.readLine();
			
			while((line = reader.readLine()) != null)
			{
				String[] parts = line.toUpperCase().split("\t");
				String orf = stripIsoform(parts[0]);
				String symbol = parts[1].replace(' ', '_');
				
				// Use the first occurrence of an id for the mapping if there is
				// ambiguity
				if(!synMap.containsKey(symbol))
				{
					synMap.put(symbol, orf);
				}
				
				if(!revSynMap.containsKey(orf))
				{
					revSynMap.put(orf, symbol);
				}
			}
			
			System.out.println("Loaded " + synMap.size() + " symbol -> ORF mappings");
			System.out.println("Loaded " + revSynMap.size() + " ORF -> symbol mappings");
			
			reader.close();
		}
		catch(IOException e)
		{
			e.printStackTrace();
		}
	}
	
	public static String stripIsoform(String fullId)
	{
		String id = fullId.toUpperCase();
		
		int ind = id.lastIndexOf('.');
		if(ind < 0)
		{
			return id;
		}
		
		return id.substring(0, ind);
	}
	
	public static void sdrem092211()
	{
		String modelDir = ARAB_DIR + "sdrem092211/";
		
		evaluatePredictions(ARAB_DIR + "hpaEffectors.txt",
				modelDir + "10.targets",
				modelDir + "nodeScores_itr10_PathWeight_20_1000.txt",
				0.01,
				ARAB_DIR + "proteinDna/tfList.txt",
				modelDir + "itr10Eval.txt");

		gridColsToSymbols(modelDir + "tfActivityPriors_round9.txt",
				modelDir + "tfActivityPriors_round9_symbols.txt");
		convertModelIds(modelDir + "10.model",
				modelDir + "10_symbols.model");

		String cytoscapeDir = modelDir + "cytoscape/";
		generateCytoscape(ARAB_DIR + "hpaEffectors.txt",
				modelDir + "10.targets",
				modelDir + "internalNodes10.txt",
				modelDir + "pathEdges_itr10.txt",
				2,
				cytoscapeDir + "10match" + 2 + ".sif",
				cytoscapeDir + "10role.noa");
	}
	
	// This SDREM run is a continuation of sdrem051712, which is why iteration
	// 3 is being used as the final iteration (7 iterations had already been run)
	public static void sdrem052112()
	{
		String modelDir = ARAB_DIR + "sdrem052112/";
		
		evaluatePredictions(ARAB_DIR + "hpaEffectorsGrouped.txt",
				modelDir + "3.targets",
				modelDir + "nodeScores_itr3_PathWeight_20_1000.txt",
				0.01,
				ARAB_DIR + "proteinDna/tfList.txt",
				modelDir + "itr3Eval.txt");

		gridColsToSymbols(modelDir + "tfActivityPriors_round2.txt",
				modelDir + "tfActivityPriors_round2_symbols.txt");
		convertModelIds(modelDir + "3.model",
				modelDir + "3_symbols.model");

		String cytoscapeDir = modelDir + "cytoscape/";
		generateCytoscape(ARAB_DIR + "hpaEffectorsGrouped.txt",
				modelDir + "3.targets",
				modelDir + "internalNodes3.txt",
				modelDir + "pathEdges_itr3.txt",
				2,
				cytoscapeDir + "3match" + 2 + ".sif",
				cytoscapeDir + "3role.noa");
		
		// Calculate a p-value for the overlap between the predictions and
		// the functional validation from Mukhtar2011
		// overlap: 6
		// predictions: 83
		// validated: 15
		// Arabidopsis proteins: 5052
		double pval = StatUtil.overlapSignificance(6, 15, 83, 5052);
		System.out.println("Overlap with Mukhtar2011 validation (p-val): " + pval);
	}
	
	/**
	 * Rewrite a human TF binding grid that uses gene ORF names for TF columns to use
	 * gene symbols for TF column names.
	 * @param inFile
	 * @param outFile
	 */
	public static void gridColsToSymbols(String inFile, String outFile)
	{
		try
		{
			BufferedReader reader = new BufferedReader(new FileReader(inFile));
			String[] header = reader.readLine().split("\t");
			
			PrintWriter writer = new PrintWriter(new FileWriter(outFile));
			// Write the gene name header column
			writer.print(header[0]);
			
			// Write the TF names
			for(int h = 1; h < header.length; h++)
			{
				writer.print("\t" + getSymbol(header[h]));
			}
			writer.println();
			
			// Print the rest of the file
			String line;
			while((line = reader.readLine()) != null)
			{
				writer.println(line);
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
	 * Generate a new model file that uses TF gene symbols instead of TF ORFs
	 * @param inModelFile
	 * @param outModelFile
	 */
	public static void convertModelIds(String inModelFile, String outModelFile)
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
					writer.println(getSymbol(parts[0]) + "\t" + parts[1]);
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
	 * Creates files that can be loaded into Cytoscape.
	 */
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
				HashSet<String> sourcesOrig = MapUtil.loadSet(srcFile, false);
				HashSet<String> sources = new HashSet<String>();
				for(String s : sourcesOrig)
				{
					sources.add(s.toUpperCase());
				}
				
				Set<String> targets = MapUtil.loadMap(targFile, 0, 1, false).keySet();
				HashSet<String> nodes = MapUtil.loadSet(nodeFile, false);

				HashSet<String> all = new HashSet<String>();
				all.addAll(sources);
				all.addAll(targets);
				all.addAll(nodes);

				PrintWriter noaWriter = new PrintWriter(new FileWriter(noaOut));
				noaWriter.println("Role (class=java.lang.String)");
				for(String protein : all)
				{
					noaWriter.print(getSymbol(protein) + " = ");

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
						sifWriter.println(getSymbol(src) + " " + parts[1] + " " +
								getSymbol(targ));
						System.out.println(count);
					}
					count++;
				}

				sifWriter.close();
			}
			catch(IOException e)
			{
				e.printStackTrace();
			}
		}
		
	/**
	 * Check how many sources are present in the interaction network.
	 * Ignores edge directionality.
	 * @param sourceFile no header
	 * @param edgeFile no header
	 */
	public static void checkSources(String sourceFile, String edgeFile)
	{
		// Want the sources to be capitalized
		Set<String> sources = MapUtil.loadMap(sourceFile, 0, 0, false).keySet();
		HashSet<String> networkProts = MapUtil.union(MapUtil.loadMap(edgeFile, 0, 1, false).keySet(),
				MapUtil.loadMap(edgeFile, 2, 1, false).keySet());
		
		HashSet<String> missing = MapUtil.subtract(sources, networkProts);
		System.out.println(missing.size() + " of " + sources.size() + 
				" are missing from the interaction network:");
		for(String mis : missing)
		{
			System.out.println(mis);
		}
	}
}
