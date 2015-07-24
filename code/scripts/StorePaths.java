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

package scripts;

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.util.Properties;

import alg.EdgeOrientAlg;

/**
 * Store and/or filter stored paths so that paths do not need to be
 * enumerated multiple times during an SDREM run.
 *
 */
public class StorePaths
{

	/**
	 * @param args
	 */
	public static void main(String[] args)
	{
		if(args.length != 1)
		{
			throw new IllegalArgumentException("Usage: StorePaths <properties_file>");
		}
		
		Properties defaults = new Properties();
		setDefaults(defaults);
		
		Properties props = new Properties(defaults);
		
		// Load the properties from the specified file
		try
		{
			FileInputStream propsIn = new FileInputStream(args[0]);
			props.load(propsIn);
			propsIn.close();
			
			// For each property that must be specified (there is no default),
			// verify it is provided before proceeding
			String mode = props.getProperty("mode").toLowerCase();
			if(mode.contains("store"))
			{
				String sourcesFile = checkAndGet(props, "sources.file");
				String targetsFile = checkAndGet(props, "targets.file");
				String edgesFile = checkAndGet(props, "edges.file");
				String priorsFile = props.getProperty("node.priors.file");
				int pathLength = Integer.parseInt(props.getProperty("max.path.length"));
				double defPrior = Double.parseDouble(props.getProperty("default.node.prior"));
				String stored = checkAndGet(props, "stored.paths.dir");
				
				// Split the stored paths directory into the parent path and the output directory
				// name
				File storedFile = new File(stored);
				String parent = storedFile.getParent();
				String storedName = storedFile.getName();
				
				EdgeOrientAlg.storePaths(sourcesFile,
					targetsFile,
					edgesFile,
					priorsFile,
					pathLength,
					defPrior,
					storedName,
					parent);
			}
			
			// Allow for both storing and filtering
			if(mode.contains("filter"))
			{
				int paths = Integer.parseInt(props.getProperty("path.enum.bound"));
				String stored = checkAndGet(props, "stored.paths.dir");
				String filtered = checkAndGet(props, "filtered.paths.dir");
				
				EdgeOrientAlg.filterStoredPathsMt(paths, stored, filtered);
			}
			
			if(!mode.contains("store") && !mode.contains("filter"))
			{
				throw new IllegalArgumentException("mode must be store, filter, or store/filter");
			}
		}
		catch(IOException e)
		{
			e.printStackTrace();
		}
	}
	
	/**
	 * Set the default properties
	 * @param defaults
	 */
	public static void setDefaults(Properties defaults)
	{
		defaults.setProperty("mode", "store/filter");
		
		// Options for storing paths
		defaults.setProperty("sources.file", "");
		defaults.setProperty("targets.file", "");
		defaults.setProperty("edges.file", "");
		defaults.setProperty("node.priors.file", "");
		defaults.setProperty("max.path.length", "5");
		defaults.setProperty("default.node.prior", "1.0");
		defaults.setProperty("stored.paths.dir", "");
		
		// Options for filtering paths
		defaults.setProperty("path.enum.bound", "100000");
		defaults.setProperty("filtered.paths.dir", "");
	}

	/**
	 * Obtain the property value or throw an exception
	 * if the property value is the empty string
	 * @param props
	 * @param name
	 * @return
	 * @throws IllegalArgumentException
	 */
	public static String checkAndGet(Properties props, String name) throws IllegalArgumentException
	{
		String str = props.getProperty(name).trim();
		
		if(str.equals(""))
		{
			throw new IllegalArgumentException("Must specify " + name);
		}
		
		return str;
	}

}
