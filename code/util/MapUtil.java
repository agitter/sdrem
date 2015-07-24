package util;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Set;

public class MapUtil
{
	public static final String SYN_FILE = "../testData/binding/SGD_standardToOrf.txt";
	/** Map from yeast standard names to ORF names */
	public static HashMap<String, String> synMap = null;
	public static final String REV_SYN_FILE = "../testData/binding/SGD_orfToStandard.txt";
	/** Map from yeast ORF names to standard names */
	public static HashMap<String, String> revSynMap = null;
	
	/**
	 * Return the ORF name for the standard name of this yeast gene.
	 * @param standard
	 * @return the ORF name if it was found in the mapping, or the
	 * original standard name if it was not found
	 */
	public static String getOrf(String standard)
	{
		standard = standard.trim().toUpperCase();
		
		if(synMap == null)
		{
			// Don't want the code to crash if the hard-coded file isn't found
			File mapFile = new File(SYN_FILE);
			if(mapFile.exists())
			{
				synMap = MapUtil.loadMap(SYN_FILE, 1, 0, false);
			}
			else
			{
				synMap = new HashMap<String, String>();
			}				
		}
		
		if(synMap.containsKey(standard))
		{
			String orf = synMap.get(standard);
			if(orf.contains("|"))
			{
				String[] orfs = orf.split("\\|");
				return orfs[0];
			}
			else
			{
				return orf;
			}
		}
		else
		{
			return standard;
		}
	}
	
	/**
	 * Return the standard name for the ORF name of this yeast gene.
	 * @param orf
	 * @return the standard name if it was found in the mapping, or the
	 * original ORF name if it was not found
	 */
	public static String getStandard(String orf)
	{
		orf = orf.trim().toUpperCase();
		
		if(revSynMap == null)
		{
			// Don't want the code to crash if the hard-coded file isn't found
			File mapFile = new File(REV_SYN_FILE);
			if(mapFile.exists())
			{
				revSynMap = MapUtil.loadMap(REV_SYN_FILE, 0, 1, false);
			}
			else
			{
				revSynMap = new HashMap<String, String>();
			}
		}
		
		if(revSynMap.containsKey(orf))
		{
			String stand = revSynMap.get(orf);
			if(stand.contains("|"))
			{
				String[] stands = stand.split("\\|");
				return stands[0];
			}
			else
			{
				return stand;
			}
		}
		else
		{
			return orf;
		}
	}
	
	
	/**
	 * Create and return a mapping of entries in the keyInd column of the file
	 * to entries in the valInd column.  If a key appears more than once, the multiple
	 * values will be collected in a HashSet.  Puts all strings in upper case.
	 * Supports value columns that contain multiple values separated by '|'.
	 * Ignores values of '-' and ''.
	 * @param file tab-separated file
	 * @param keyInd
	 * @param valInd
	 * @param header true if there is a header row
	 * @return
	 */
	public static HashMap<String, HashSet<String>> loadMultiMap(String file, int keyInd, int valInd, boolean header)
	{
		try
		{
			BufferedReader reader = new BufferedReader(new FileReader(file));
			HashMap<String, HashSet<String>> map = new HashMap<String, HashSet<String>>();
	
			if(header)
			{
				reader.readLine();
			}
	
			String line;
			while((line = reader.readLine()) != null)
			{
				String[] parts = line.split("\t");
				if(parts.length > Math.max(keyInd, valInd))
				{
					String key = parts[keyInd].trim().toUpperCase();
					String val = parts[valInd].trim().toUpperCase();
	
					if(!val.equals("-") && !val.equals(""))
					{
						HashSet<String> existing;
						if(map.containsKey(key))
						{
							existing = map.get(key);
						}
						else
						{
							existing = new HashSet<String>();
						}
	
						// There may be multiple values
						String[] valParts = val.split("\\|");
						for(String valPart : valParts)
						{
							existing.add(valPart);
						}
	
						map.put(key, existing);
					}
				}
				else
				{
					System.err.println("Not enough columns: " + line);
				}
			}			
	
			reader.close();
	
			return map;
		}
		catch(IOException e)
		{
			e.printStackTrace();
			return null;
		}
	}

	/**
	 * Create and return a mapping of entries in the keyInd column of the file
	 * to entries in the valInd column.  If a key appears more than once, the multiple
	 * values will be collected in a HashSet.  Puts all strings in upper case.
	 * Assumes no header. Supports value columns that contain
	 * multiple values separated by '|'.  Ignores values of '-' and ''.
	 * @param file tab-separated file
	 * @param keyInd
	 * @param valInd
	 * @return
	 */
	public static HashMap<String, HashSet<String>> loadMultiMap(String file, int keyInd, int valInd)
	{
		return loadMultiMap(file, keyInd, valInd, false);
	}

	/**
	 * Create and return a mapping of entries in the keyInd column of the file
	 * to entries in the valInd column.  If a key appears more than once, the multiple
	 * values will be listed together separated by '|'.  Puts all strings in upper case.
	 * Ignores values of '-' and ''.
	 * @param file tab-separated file
	 * @param keyInd
	 * @param valInd
	 * @param header true if there is a header row
	 * @return
	 */
	public static HashMap<String, String> loadMap(String file, int keyInd, int valInd, boolean header)
	{
		return loadMap(file, keyInd, valInd, header, "-");
	}
		

	/**
	 * Create and return a mapping of entries in the keyInd column of the file
	 * to entries in the valInd column.  If a key appears more than once, the multiple
	 * values will be listed together separated by '|'.  Puts all strings in upper case.
	 * Ignores values of ignoreVal and ''.
	 * @param file tab-separated file
	 * @param keyInd
	 * @param valInd
	 * @param header true if there is a header row
	 * @param ignoreVal do not load keys if their value is this string
	 * @return
	 */
	public static HashMap<String, String> loadMap(String file, int keyInd, int valInd, boolean header, String ignoreVal)
	{
		try
		{
			BufferedReader reader = new BufferedReader(new FileReader(file));
			HashMap<String, String> map = new HashMap<String, String>();

			if(header)
			{
				reader.readLine();
			}
	
			String line;
			while((line = reader.readLine()) != null)
			{
				String[] parts = line.split("\t");
				if(parts.length > Math.max(keyInd, valInd))
				{
					String key = parts[keyInd].trim().toUpperCase();
					String val = parts[valInd].trim().toUpperCase();
	
					if(!val.equals(ignoreVal) && !val.equals(""))
					{
						if(map.containsKey(key))
						{
							String existing = map.get(key);
							map.put(key, existing + "|" + val);
						}
						else
						{
							map.put(key, val);
						}
					}
				}
				else
				{
					System.err.println("Not enough columns: " + line);
				}
			}			
	
			reader.close();
	
			return map;
		}
		catch(IOException e)
		{
			e.printStackTrace();
			return null;
		}
	}

	/**
	 * Create and return a mapping of entries in the keyInd column of the file
	 * to entries in the valInd column.  If a key appears more than once, the multiple
	 * values will be listed together separated by '|'.  Puts all strings in upper case.
	 * Assumes no header.  Ignores values of '-' and ''.
	 * @param file tab-separated file
	 * @param keyInd
	 * @param valInd
	 * @return
	 */
	public static HashMap<String, String> loadMap(String file, int keyInd, int valInd)
	{
		return loadMap(file, keyInd, valInd, false);
	}
	
	
	/**
	 * Create and return a set of all Strings in the file.  Assumes no header.
	 * @param file
	 * @return
	 */
	public static HashSet<String> loadSet(String file)
	{
		return loadSet(file, false);
	}
	
	/**
	 * Create and return a set of all Strings in the file.
	 * @param file
	 * @param header true if there is a header line
	 * @return
	 */
	public static HashSet<String> loadSet(String file, boolean header)
	{
		try
		{
			BufferedReader reader = new BufferedReader(new FileReader(file));
			HashSet<String> set = new HashSet<String>();
			
			String line;
			if(header)
			{
				line = reader.readLine();
			}
			
			while((line = reader.readLine()) != null)
			{
				set.add(line.trim());
			}
						
			return set;
		}
		catch(IOException e)
		{
			e.printStackTrace();
			return null;
		}
	}
	
	/**
	 * Load multiple sets by creating a map from Strings to sets.  The elements
	 * in the set are given in the 1st column, the name of the set to which the element
	 * belongs is in the 2nd column of the input file.  Assumes no header line.
	 * @param file
	 * @return
	 */
	public static HashMap<String, HashSet<String>> loadSets(String file)
	{
		return loadSets(file, false);
	}
	
	/**
	 * Load multiple sets by creating a map from Strings to sets.  The elements
	 * in the set are given in the 1st column, the name of the set to which the element
	 * belongs is in the 2nd column of the input file
	 * @param file
	 * @param header
	 * @return
	 */
	public static HashMap<String, HashSet<String>> loadSets(String file, boolean header)
	{
		try
		{
			BufferedReader reader = new BufferedReader(new FileReader(file));
			HashMap<String, HashSet<String>> map = new HashMap<String, HashSet<String>>();
			
			String line;
			if(header)
			{
				line = reader.readLine();
			}
			
			while((line = reader.readLine()) != null)
			{
				String[] parts = line.trim().split("\t");
				String key = parts[1];
				String val = parts[0];
				
				HashSet<String> set;
				if(map.containsKey(key))
				{
					set = map.get(key);
				}
				else
				{
					set = new HashSet<String>();
				}
				set.add(val);
				map.put(key, set);
			}
			
			return map;
		}
		catch(IOException e)
		{
			e.printStackTrace();
			return null;
		}
	}
	
	
	/**
	 * Create and return a list of all Strings in the file.
	 * Assumes no header line.
	 * @param file
	 * @return
	 */
	public static ArrayList<String> loadList(String file)
	{
		return loadList(file, false);
	}
	
	/**
	 * Create and return a list of all Strings in the file.
	 * @param file
	 * @param header true if there is a header line
	 * @return
	 */
	public static ArrayList<String> loadList(String file, boolean header)
	{
		try
		{
			BufferedReader reader = new BufferedReader(new FileReader(file));
			ArrayList<String> list = new ArrayList<String>();
			
			String line;
			if(header)
			{
				line = reader.readLine();
			}
			
			while((line = reader.readLine()) != null)
			{
				list.add(line.trim());
			}
						
			return list;
		}
		catch(IOException e)
		{
			e.printStackTrace();
			return null;
		}
	}
	
	// Does this belong here?
	/**
	 * Return a new set containing the intersection of the two parameter sets
	 */
	public static <T> HashSet<T> intersection(Set<T> set1, Set<T> set2)
	{
		HashSet<T> intersect = new HashSet<T>(set1);
		intersect.retainAll(set2);
		return intersect;
	}
	
	// Does this belong here?
	/**
	 * Return a new set containing the intersection of the two parameter sets
	 */
	public static <T> HashSet<T> union(Set<T> set1, Set<T> set2)
	{
		HashSet<T> uni = new HashSet<T>(set1);
		uni.addAll(set2);
		return uni;
	}
	
	// Does this belong here?
	/**
	 * Return a new set containing the members of set 1 that are not present
	 * in set2
	 */
	public static <T> HashSet<T> subtract(Set<T> set1, Set<T> set2)
	{
		HashSet<T> sub = new HashSet<T>(set1);
		sub.removeAll(set2);
		return sub;
	}
	
	/**
	 * Return true if the two collections contain the same elements
	 * @param col1
	 * @param col2
	 * @return
	 */
	public static <T> boolean sameElements(Collection<T> col1, Collection<T> col2)
	{
		return col1.containsAll(col2) && col2.containsAll(col1);
	}
}
