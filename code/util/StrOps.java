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
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.Collection;
import java.util.HashSet;

public class StrOps
{

	/**
	 * Create the directory and any necessary parent directories if they do not exist
	 * @param directory
	 * @return a reference to the directory
	 * @throws IOException
	 */
	public static File createDirectory(String directory) throws IOException
	{
		File newDirectory = new File(directory);
		if (!newDirectory.exists())
		{
			// Throw an exception if the directory cannot be created
			if (!newDirectory.mkdirs())
			{
				throw new IOException("Could not create the directory " + directory);
			}
		}
		return newDirectory;
	}

	/**
	 * Returns the original String after removing the prefix (and everything before
	 * it) and the suffix (and everything after it).  If the prefix and/or suffix
	 * are "", they are ignored.  If the prefix and/or suffix are not "" and cannot
	 * both be found in the String or the prefix comes after the suffix, the original
	 * String is returned.  The first match of the prefix and last match of the suffix
	 * are taken
	 * @param original
	 * @param prefix
	 * @param suffix
	 * @return
	 */
	public static String trimString(String original, String prefix, String suffix)
	{
		if(prefix.equals("") && suffix.equals(""))
		{
			return original;
		}
		else if(prefix.equals(""))
		{
			int back = original.lastIndexOf(suffix);
			if(back <= 0)
			{
				return original;
			}
			else
			{
				return original.substring(0, back);
			}
		}
		else if(suffix.equals(""))
		{
			int front = original.indexOf(prefix);
			if(front < 0)
			{
				return original;
			}
			else
			{
				return original.substring(front + prefix.length());
			}
		}
		else
		{
			int front = original.indexOf(prefix);
			int back = original.lastIndexOf(suffix);
			if(front < 0 || back <= 0 || (back < front + prefix.length()))
			{
				return original;
			}
			else
			{
				return original.substring(front + prefix.length(), back);
			}
		}
	}
	
	/**
	 * Return [alphabetically first string]\t[alphabetically second string].
	 * Case insensitive
	 * @param str1
	 * @param str2
	 */
	public static String sortAndCombine(String str1, String str2)
	{
		if(str1.compareToIgnoreCase(str2) <= 0)
		{
			return str1 + "\t" + str2;
		}
		else
		{
			return str2 + "\t" + str1;
		}
	}
	
	/**
	 * Capitalize only the first letter
	 * @param str
	 * @return
	 */
	public static String titleCase(String str)
	{
		String first = str.toUpperCase().substring(0, 1);
		String remaining = str.toLowerCase().substring(1);
		return first + remaining;
	}

	/**
	 * Merge two files by writing a new file where each line is
	 * [file1Line][spacer][file2Line].  If one file is longer
	 * than the other writes the lines from the longer file
	 * to the merged file without the spacer.
	 * @param file1
	 * @param file2
	 * @param outFile
	 * @param spacer
	 * @throws IOException
	 */
	public static void mergeFilesHoriz(String file1, String file2, String outFile,
			String spacer) throws IOException
	{
		BufferedReader reader1 = new BufferedReader(new FileReader(file1));
		BufferedReader reader2 = new BufferedReader(new FileReader(file2));
		PrintWriter writer = new PrintWriter(new FileWriter(outFile));
		
		boolean done1 = false, done2 = false;
		String line1 = null, line2 = null;
		while(!done1 && !done2)
		{
			line1 = reader1.readLine();
			line2 = reader2.readLine();
			
			if(line1 == null)
			{
				done1 = true;
			}
			if(line2 == null)
			{
				done2 = true;
			}
			
			// Write both lines
			if(!done1 && !done2)
			{
				writer.println(line1 + spacer + line2);
			}
			// Only write from file1
			else if(!done1)
			{
				writer.println(line1);
			}
			// Only write from file2
			else if(!done2)
			{
				writer.println(line2);
			}
		}
		
		// One of the files has been written complete, see if the other still
		// has lines
		if(!done1)
		{
			while((line1 = reader1.readLine()) != null)
			{
				writer.println(line1);
			}
		}
		if(!done2)
		{
			while((line2 = reader2.readLine()) != null)
			{
				writer.println(line2);
			}
		}
		
		writer.close();
	}

	/**
	 * Return the first occurrence of key in parts[].  Case-sensitive.
	 * @param parts
	 * @param key
	 * @return
	 */
	public static int findString(String[] parts, String key)
	{
		for(int i = 0; i < parts.length; i++)
		{
			if(parts[i].equals(key))
			{
				return i;
			}
		}
		return -1;
	}
	
	/**
	 * Return the first occurrence of key in the first line of a
	 * tab-delimited file.  Case-sensitive.
	 * @param filename
	 * @param key
	 * @return
	 * @throws IOException
	 */
	public static int findStringHeader(String filename, String key)
		throws IOException
	{
		BufferedReader reader = new BufferedReader(new FileReader(filename));
		int val = findString(reader.readLine().split("\t"), key);
		reader.close();
		return val;
	}
	
	/**
	 * Write a 2-column file with the standard gene name and ORF name
	 * for all genes in the orfNameFile, which has no header
	 * @param orfNameFile no header
	 * @param outFile
	 */
	public static void writeStandard(String orfNameFile, String outFile)
	{
		try
		{
			HashSet<String> orfSet = MapUtil.loadSet(orfNameFile, false);
			
			PrintWriter writer = new PrintWriter(new FileWriter(outFile));
			for(String orf : orfSet)
			{
				writer.println(MapUtil.getStandard(orf) + "\t" + orf);
			}
			writer.close();
		}
		catch(IOException e)
		{
			e.printStackTrace();
		}
	}
	
	/**
	 * Return the number of lines in a file
	 * @param file
	 * @return
	 * @throws IOException
	 */
	public static int countLines(String file) throws IOException
	{
		BufferedReader reader = new BufferedReader(new FileReader(file));
		int lines = 0;
		while(reader.readLine() != null)
		{
			lines++;
		}
		return lines;
	}
	
	/**
	 * Write each String in the collection to a file, one String per line.
	 * @param strings
	 * @param filename
	 * @throws IOException
	 */
	public static void writeCollection(Collection<String> strings, String filename)
		throws IOException
	{
		PrintWriter writer = new PrintWriter(new FileWriter(filename));
		for(String s : strings)
		{
			writer.println(s);
		}
		writer.close();
	}
}
