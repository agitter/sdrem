package edu.cmu.cs.sb.drem;

import edu.cmu.cs.sb.core.*;
import edu.umd.cs.piccolo.nodes.PPath;
import javax.swing.*;
import javax.swing.event.*;
import java.awt.*;
import java.awt.event.*;
import java.util.*;
import java.text.NumberFormat;
import java.io.*;

import javax.imageio.*;
import java.awt.image.*;
import edu.umd.cs.piccolo.nodes.PImage;
import util.MapUtil;

/**
 * Class to encapsulate window used to specify a file to save a DREM model
 */ 
public class DREMGui_SaveModel extends JPanel
{
	final static String SYN_FILE = "SGD_standardToOrf.txt";
	static HashMap synMap = null;
	
	final static Color bgColor = Color.white;
	final static Color fgColor = Color.black;
	final DREM_Timeiohmm theTimeiohmm;
	final DREM_Timeiohmm.Treenode treecopy;
	final JFileChooser theChooser;

	/**
	 * Class constructor
	 */
	public DREMGui_SaveModel(final DREM_Timeiohmm theTimeiohmm,final DREM_Timeiohmm.Treenode treecopy,
			final JFrame theFrame,final DREMGui theDREMGui)
	{
		this.theTimeiohmm = theTimeiohmm;
		this.treecopy = treecopy;

		setLayout(new BoxLayout(this, BoxLayout.PAGE_AXIS));
		setForeground(fgColor);
		theChooser = new JFileChooser();
		add(theChooser);
		theChooser.setDialogType(JFileChooser.SAVE_DIALOG);
		theChooser.addActionListener(new ActionListener() 
		{
			public void actionPerformed(ActionEvent e) 
			{
				// set label's icon to the current image
				String state = (String)e.getActionCommand();

				if (state.equals(JFileChooser.CANCEL_SELECTION))
				{
					theFrame.setVisible(false);
					theFrame.dispose();
				}
				else if (state.equals(JFileChooser.APPROVE_SELECTION))
				{
					File f = theChooser.getSelectedFile();

					try
					{
						PrintWriter pw = new PrintWriter(new FileOutputStream(f));
						pw.print(theTimeiohmm.saveString(treecopy));
						pw.println("COLORS");
						pw.print(theDREMGui.saveColors());
						pw.close();
						
						saveActivityScoresDynamic(f);
						saveActivityScores(f);
					}
					catch (final IOException fex)
					{
						javax.swing.SwingUtilities.invokeLater(new Runnable() 
						{
							public void run() 
							{
								JOptionPane.showMessageDialog(null, fex.getMessage(), 
										"Exception thrown", JOptionPane.ERROR_MESSAGE);
							}
						});
						fex.printStackTrace(System.out);
					}
					theDREMGui.bsavedchange = true;
					theFrame.setVisible(false);
					theFrame.dispose();
				}
			}
		});			      
	}
	
	
	/**
	 * Uses the filename to create three files:
	 * [filename].model - the model file
	 * [filename].model.activities - the TF activities
	 * [filename].model.activitiesStd - the TF activities with their standard names
	 * @param filename
	 */
	public void batchSave(String filename)
	{
		theChooser.setSelectedFile(new File(filename + ".model"));
		theChooser.approveSelection();
	}
	
	/**
	 * A temporary means for retrieving the max TF activity scores
	 * @param modelFile
	 */
	private void saveActivityScores(File modelFile) throws IOException
	{
		if(synMap == null)
		{
			synMap = MapUtil.loadMap(SYN_FILE, 1, 0, false);
		}
		
		String filename = modelFile.getAbsolutePath() + ".activities";
		PrintWriter writer = new PrintWriter(new FileWriter(filename));
		
		String filenameStd = modelFile.getAbsolutePath() + ".activitiesStd";
		PrintWriter writerStd = new PrintWriter(new FileWriter(filenameStd));
		
		for(int t = 0; t < theTimeiohmm.tfNames.length; t++)
		{
			double score = theTimeiohmm.dMaxTFActivityScore[t];
			// Write the true activity score of all TFs
//			score = score / (1 + score);
//			if(score >= 0.9)
//			if(score > 10)
//			{
				String tf = theTimeiohmm.tfNames[t].toUpperCase();
				String stdName = tf;
				if(synMap.containsKey(tf))
				{
					String newName = (String) synMap.get(tf);
					if(newName.contains("|"))
					{
						newName = newName.substring(0, newName.indexOf('|'));
					}
					tf = newName;
				}
				
				writer.println(tf + "\t" + score);
				writerStd.println(stdName + "\t" + score);
//			}
		}
		writer.close();
		writerStd.close();
		
		String filenameFlag = modelFile.getAbsolutePath() + ".activities.flag";
		PrintWriter writerFlag = new PrintWriter(new FileWriter(filenameFlag));
		writerFlag.println("DONE");
		writerFlag.close();
	}
	
	/**
	 * A temporary means for retrieving the max TF activity scores at
	 * each time point.  Must be called before saveActivityScores if at all
	 * because it does write its own flag.
	 * @param modelFile
	 */
	private void saveActivityScoresDynamic(File modelFile) throws IOException
	{
		if(synMap == null)
		{
			synMap = MapUtil.loadMap(SYN_FILE, 1, 0, false);
		}
		
		String filename = modelFile.getAbsolutePath() + ".activitiesDynamic";
		PrintWriter writer = new PrintWriter(new FileWriter(filename));
		
		writer.print(theTimeiohmm.saveActivityScores(treecopy));
		
		writer.close();
	}
}