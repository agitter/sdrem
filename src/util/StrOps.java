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

import java.io.File;
import java.io.IOException;

public class StrOps {

	/**
	 * Create the directory if it does not exist
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
	 * be found in the String or the prefix comes after the suffix, the original
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

}
