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

package util;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;

public class Entrez {

	// Human only for now
	private static final String ENTREZ_GENE_INFO = "../../human/entrez/Homo_sapiens.gene_info";
	private static final String ENTREZ_GENE_HISTORY = "../../human/entrez/Homo_sapiens.gene_history";
	/** Map from discontinued Entrez Gene id to current id */
	private static HashMap<String, String> geneHistory = null;
	/** Map from Entrez Gene id to official gene symbol */
	private static HashMap<String, String> gene2Symbol = null;
	/** Map from official gene symobls and gene name synonyms to Entrez Gene ids.
	 * If a string matches an official symbol and a synonym for a different symbol,
	 * the id corresponding to the official symbol is used.  If multiple synonyms
	 * match the same string, the lower gene id is returned (due to the way the
	 * gene_info file is sorted) */
	private static HashMap<String, String> name2Gene = null;
	
	public static void main(String[] args)
	{
//		System.out.println(isDiscontinuedId("llkj"));
//		System.out.println(isDiscontinuedId("168"));
//		System.out.println(isDiscontinuedId("139"));
//		System.out.println(getCurrentId("llkj"));
//		System.out.println(getCurrentId("168"));
//		System.out.println(getCurrentId("139"));
//		System.out.println(getId("ArhGAP11B"));
//		System.out.println(getId("aaa1"));
//		System.out.println(getId("XYXYXYXYXYYXXY"));
//		System.out.println(getSymbol("404744"));
//		System.out.println(getSymbol("5193"));
//		System.out.println(getSymbol("ABCDEF"));
//		System.out.println(getSymbol(getId("TP53")));
//		System.out.println(getId(getSymbol("7157")));
//		System.out.println(getId("7157"));
//		System.out.println(getId("LOC7157"));
//		System.out.println(getId("TP53-ISOFORM"));
//		System.out.println(getId("TP53-FORM"));
//		System.out.println(getId("TP53-isoform100"));
//		System.out.println(getId("TP53-formN"));
//		System.out.println(getId("TP5-3"));
//		System.out.println(getId("TP-53"));
//		System.out.println(getId("TP53-"));
//		System.out.println(getId("-TP53"));
//		System.out.println(getId("TP--53"));
//		System.out.println(getId("-T-P--53-"));
//		System.out.println(getId("TP53--"));
//		System.out.println(getId("--TP53"));
//		System.out.println(getId("168"));
	}
	
	/**
	 * Load the Entrez Gene maps
	 */
	public static void loadMaps()
	{
		try
		{
			// Allow '-' as a valid value
			geneHistory = MapUtil.loadMap(ENTREZ_GENE_HISTORY, 2, 1, true, " ");
			gene2Symbol = MapUtil.loadMap(ENTREZ_GENE_INFO, 1, 2, false);

			// First load a map from official symbols to gene ids
			// Note that not all gene symbols are unique
			name2Gene = MapUtil.loadMap(ENTREZ_GENE_INFO, 2, 1, false);

			// Then add synonyms and cross references if they do not conflict
			// with an existing name or symbol
			BufferedReader reader = new BufferedReader(new FileReader(ENTREZ_GENE_INFO));
			String line;
			int collisions = 0;
			while((line = reader.readLine()) != null)
			{
				String[] parts = line.split("\t");
				String id = parts[1].trim();
				
				String allSyn = parts[4].trim().toUpperCase();

				if(!allSyn.equals("-") && !allSyn.equals(""))
				{
					String[] synonyms = allSyn.split("\\|");
					for(String syn : synonyms)
					{
						if(name2Gene.containsKey(syn))
						{
							collisions++;
						}
						else
						{
							name2Gene.put(syn, id);
						}
					}
				}
				
				// Add cross references to other databases
				String allXref = parts[5].trim().toUpperCase();

				if(!allXref.equals("-") && !allXref.equals(""))
				{
					String[] xrefs = allXref.split("\\|");
					for(String xref : xrefs)
					{
						if(name2Gene.containsKey(xref))
						{
							collisions++;
						}
						else
						{
							name2Gene.put(xref, id);
						}
					}
				}
			}

			reader.close();
//			System.out.println(collisions + " synonyms were not loaded");
		}
		catch(IOException e)
		{
			throw new IllegalStateException("Maps cannot be loaded", e);
		}
	}
	
	public static void clearMaps()
	{
		geneHistory = null;
		gene2Symbol = null;
		name2Gene = null;
	}
	
	public static void mapsNeeded()
	{
		if(geneHistory == null || gene2Symbol == null || name2Gene == null)
		{
			loadMaps();
		}
	}
	
	/**
	 * Return the Entrez Gene id.  If there are multiple matches for the
	 * name and it is an official symbol, an arbitrary matching id is
	 * selected.  Also checks if the name is actually a valid id and if
	 * the name is of the form LOCXXXXXXX, where XXXXXXX is a valid id.
	 * If a match is still not found, tries removing isoform suffixes
	 * '-ISOFORM*' and '-FORM*' where * is the wildcard.
	 * As a last resort, try removing '-' as long as it isn't surrounded
	 * by numbers.  Supports looking up cross references to other databases if
	 * the database name is prepended, e.g. 'Ensembl:ENSG00000114771'.
	 * @param name An official symbol or synonym
	 * @return the id or null if there is no match
	 */
	public static String getId(String name)
	{
		mapsNeeded();
		
		name = name.trim().toUpperCase();
		if(name2Gene.containsKey(name))
		{
			String id = name2Gene.get(name);
			String[] matches = id.split("\\|");
			id = matches[0];
			return getCurrentId(id);
		}
		// Check if the given string is already a valid id
		else if(isGeneId(name))
		{
			return getCurrentId(name);
		}
		// Check if a valid id is contained in the name
		else if(name.startsWith("LOC"))
		{
			String nameSuffix = name.substring(3);
			if(isGeneId(nameSuffix))
			{
				return getCurrentId(nameSuffix);
			}
		}
		// Try removing suffixes that denote an isoform
		else if(name.contains("-ISOFORM"))
		{
			int isoInd = name.lastIndexOf("-ISOFORM");
			return getId(name.substring(0, isoInd));
		}
		// Try removing suffixes that denote an isoform
		else if(name.contains("-FORM"))
		{
			int isoInd = name.lastIndexOf("-FORM");
			return getId(name.substring(0, isoInd));
		}
		// Try removing -
		else if(name.contains("-"))
		{
			String modName = name;
			int dashInd = name.indexOf('-');
			while(dashInd >= 0)
			{
				
				if(dashInd == 0)
				{
					modName = modName.substring(dashInd+1);
					
					// Adjust for the character that was removed
					dashInd = modName.indexOf('-', dashInd - 1);
				}
				else if(dashInd == modName.length()-1)
				{
					modName = modName.substring(0, dashInd);
					dashInd = -1;
				}
				// Don't want to join two adjacent numbers
				// that were previously separated
				else if(!(modName.substring(dashInd-1, dashInd).matches("\\d+") &&
						modName.substring(dashInd+1, dashInd+2).matches("\\d+")))
				{
					modName = modName.substring(0, dashInd) + modName.substring(dashInd+1);

					// Adjust for the character that was removed
					dashInd = modName.indexOf('-', dashInd - 1);
				}
				else
				{
					dashInd = modName.indexOf('-', dashInd + 1);
				}
			}
			
			// Don't want to recursively call getId in case there were
			// other dashes that weren't removed
			if(name2Gene.containsKey(modName))
			{
				String id = name2Gene.get(modName);
				String[] matches = id.split("\\|");
				id = matches[0];
				return getCurrentId(id);
			}
		}
		
		// An id could not be found
		return null;
	}
	
	/**
	 * Determine if this Entrez Gene id has been discontinued
	 * @param id
	 * @return
	 */
	public static boolean isDiscontinuedId(String id)
	{
		mapsNeeded();
		return geneHistory.containsKey(id.trim());
	}
	
	/**
	 * Determine if this string is a Entrez Gene id (current or discontinued)
	 * @param id
	 * @return
	 */
	public static boolean isGeneId(String id)
	{
		mapsNeeded();
		return gene2Symbol.containsKey(id) || geneHistory.containsKey(id);
	}
	
	/**
	 * Obtain the current Entrez Gene id if this id has been
	 * discontinued.  If the id has been discontinued but not replaced
	 * by a current id, return the discontinued id.
	 * @param id
	 * @return
	 */
	public static String getCurrentId(String id)
	{
		mapsNeeded();
		id = id.trim();
		
		if(!isDiscontinuedId(id))
		{
			return id;
		}
		else
		{
			String newId = geneHistory.get(id);
			if(newId.equals("-"))
			{
				return id;
			}
			else
			{
				return newId;
			}
		}
	}
	
	/**
	 * Use the Entrez Gene id to lookup the official gene symbol.
	 * Translates the id to the current id first if it is discontinued.
	 * If the id is discontinued but does not map to a current id,
	 * does not try to lookup the old gene symbol of the discontinued id.
	 * @param id
	 * @return the symbol or null if the id could not be found
	 */
	public static String getSymbol(String id)
	{
		mapsNeeded();
		
		id = getCurrentId(id);
		
		if(gene2Symbol.containsKey(id))
		{
			return gene2Symbol.get(id);
		}
		else
		{
			return null;
		}
	}
}
