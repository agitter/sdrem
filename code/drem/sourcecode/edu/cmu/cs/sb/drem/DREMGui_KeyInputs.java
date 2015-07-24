package edu.cmu.cs.sb.drem;

import edu.cmu.cs.sb.core.*;
import edu.umd.cs.piccolo.nodes.PPath;
import javax.swing.*;
import javax.swing.event.*;
import java.awt.*;
import java.awt.event.*;
import java.util.*;
import java.text.NumberFormat;
import javax.swing.*;

/**
 * Class encapsulates the window used to specify the criteria that determines
 * transcription factor labels appear on the main interface
 */
public class DREMGui_KeyInputs extends JPanel implements ActionListener, ChangeListener, ItemListener
{

	final static Color bgColor = Color.white;
	final static Color fgColor = Color.black;
	JSlider theSlider;
	JSlider theSliderPercent;
	JButton hideButton;
	JButton colorButton;
	Hashtable theDictionary;
	Hashtable theDictionaryPercent;
	DREMGui theDREMGui;
	JLabel pvalLabel;
	JLabel percentLabel;

	ButtonGroup group = new ButtonGroup();
	JRadioButton splitButton;
	JRadioButton edgeSplitButton;
	JRadioButton edgeFullButton;
	JRadioButton activityButton;
	JRadioButton firstActivityButton;
	JFrame theFrame;

	NumberFormat nf2;


	/**
	 * Class constructor - builds interface window
	 */
	public DREMGui_KeyInputs(JFrame theFrame, DREMGui theDREMGui)
	{
		this.theFrame = theFrame;
		nf2 = NumberFormat.getInstance(Locale.ENGLISH);
		nf2.setMinimumFractionDigits(2);
		nf2.setMaximumFractionDigits(2);  

		setLayout(new BoxLayout(this, BoxLayout.PAGE_AXIS));
		setBackground(bgColor);
		setForeground(fgColor);
		int ninitval = (int) Math.round(-100*Math.log(theDREMGui.dkeyinputpvalue)/Math.log(10))/10;	

		//pval = 10^-X
		//log p-val/log 10 = -X
		//X = -log p-val/log 10
		theDREMGui.dkeyinputpvalue = Math.pow(10,-ninitval/10.0);

		JLabel theTopLabel = new JLabel("  Only display key TFs with a score less than 10^-X where X is:");
		JPanel topPanel = new JPanel();
		topPanel.setAlignmentX(Component.LEFT_ALIGNMENT);
		topPanel.add(theTopLabel);
		topPanel.setBackground(new Color((float) 0.0,(float) 1.0,(float) 0.0,(float) 0.4));

		add(topPanel);
		theSlider = new JSlider(0,120,ninitval);
		theDictionary = new Hashtable();
		for (int nindex = 0; nindex <= 12; nindex++)
		{
			theDictionary.put(new Integer(nindex*10),new JLabel(""+nindex));
		}
		theSlider.setLabelTable(theDictionary);
		this.theDREMGui = theDREMGui;
		theSlider.setMajorTickSpacing(10);
		theSlider.setMinorTickSpacing(5);
		theSlider.setPaintTicks(true);
		theSlider.setPaintLabels(true);
		theSlider.addChangeListener(this);
		theSlider.setPaintTicks(true);
		theSlider.setAlignmentX(Component.LEFT_ALIGNMENT);
		add(theSlider);

		JPanel labelPanel = new JPanel();
		pvalLabel = new JLabel("X = "+ninitval/10.0+
				";  score threshold is "+ doubleToSz(Math.pow(10,-ninitval/10.0)));
		labelPanel.setAlignmentX(Component.LEFT_ALIGNMENT);
		labelPanel.add(pvalLabel);
		labelPanel.setBackground(Color.white);
		add(labelPanel);

		JLabel binLabel;
		binLabel = new JLabel(" Compute and display key TF significance based on: ");

		JPanel binPanel = new JPanel();
		binPanel.setBackground(new Color((float) 0.0,(float) 1.0,(float) 0.0,(float) 0.4));
		binPanel.add(binLabel);
		binPanel.setAlignmentX(Component.LEFT_ALIGNMENT);
		add(binPanel);

		binLabel.setAlignmentX(Component.LEFT_ALIGNMENT);
		theDREMGui.blowbelow = true;
		theDREMGui.blowabove =true; 

		splitButton = new JRadioButton("Split Significance");
		edgeSplitButton = new JRadioButton("Path Significance Conditional on Split");
		edgeFullButton = new JRadioButton("Path Significance Overall");
		activityButton = new JRadioButton("Activity Score");
		firstActivityButton = new JRadioButton("Activity Score (first appearance per path)");


		group = new ButtonGroup();

		if (theDREMGui.nKeyInputType == 0)
		{
			splitButton.setSelected(true);
		}
		else if (theDREMGui.nKeyInputType == 1)
		{
			edgeSplitButton.setSelected(true);
		} 
		else if (theDREMGui.nKeyInputType == 2)
		{
			edgeFullButton.setSelected(true);
		}
		else if(theDREMGui.nKeyInputType == 778)
		{
			firstActivityButton.setSelected(true);
		}
		// Default to activity score
		else // if (theDREMGUI.nKeyInputType == 777)
		{
			activityButton.setSelected(true);
		}
		group.add(edgeSplitButton);
		group.add(edgeFullButton);
		group.add(splitButton);
		group.add(activityButton);
		group.add(firstActivityButton);

		splitButton.addItemListener(this);
		edgeSplitButton.addItemListener(this);
		edgeFullButton.addItemListener(this);
		activityButton.addItemListener(this);
		firstActivityButton.addItemListener(this);

		add(edgeSplitButton);
		add(edgeFullButton);
		add(splitButton);
		add(activityButton);
		add(firstActivityButton);

		splitButton.setAlignmentX(Component.LEFT_ALIGNMENT);
		splitButton.setBackground(Color.white);
		edgeSplitButton.setAlignmentX(Component.LEFT_ALIGNMENT);
		edgeSplitButton.setBackground(Color.white);
		edgeFullButton.setAlignmentX(Component.LEFT_ALIGNMENT);
		edgeFullButton.setBackground(Color.white);
		activityButton.setAlignmentX(Component.LEFT_ALIGNMENT);
		activityButton.setBackground(Color.white);
		firstActivityButton.setAlignmentX(Component.LEFT_ALIGNMENT);
		firstActivityButton.setBackground(Color.white);

		JLabel theTopLabel2 = new JLabel("For Path Significance Conditional on Split - Minimum Split %: ");
		JPanel topPanel2 = new JPanel();
		topPanel2.setAlignmentX(Component.LEFT_ALIGNMENT);
		topPanel2.add(theTopLabel2);
		topPanel2.setBackground(new Color((float) 0.0,(float) 1.0,(float) 0.0,(float) 0.4));
		add(topPanel2);


		if (DREM_Timeiohmm.BDEBUG)
		{
			System.out.println("!!!!"+theDREMGui.dsplitpercent);
		}
		int ninitvalpercent = (int) Math.round(theDREMGui.dsplitpercent*100);

		//pval = 10^-X
		//log p-val/log 10 = -X
		//X = -log p-val/log 10
		theDREMGui.dsplitpercent = ninitvalpercent/100.0;

		theSliderPercent = new JSlider(0,100,ninitvalpercent);
		theDictionaryPercent = new Hashtable();
		for (int nindex = 0; nindex <= 10; nindex++)
		{
			theDictionaryPercent.put(new Integer(nindex*10),new JLabel(""+nindex*10));
		}
		theSliderPercent.setLabelTable(theDictionaryPercent);
		theSliderPercent.setMajorTickSpacing(10);
		theSliderPercent.setMinorTickSpacing(5);
		theSliderPercent.setPaintTicks(true);
		theSliderPercent.setPaintLabels(true);
		theSliderPercent.addChangeListener(this);
		theSliderPercent.setPaintTicks(true);
		theSliderPercent.setAlignmentX(Component.LEFT_ALIGNMENT);
		add(theSliderPercent);


		JPanel labelPanel2 = new JPanel();
		percentLabel = new JLabel(" Minimum Split % is "+nf2.format(100*theDREMGui.dsplitpercent));
		labelPanel2.setAlignmentX(Component.LEFT_ALIGNMENT);
		labelPanel2.add(percentLabel);
		labelPanel2.setBackground(Color.white);
		labelPanel2.setMaximumSize(new Dimension(Integer.MAX_VALUE,20));
		add(labelPanel2);

		hideButton = new JButton("Hide Key TF Labels");
		hideButton.setActionCommand("hide");
		hideButton.addActionListener(this);

		colorButton = new JButton("Change Labels Color");
		colorButton.setActionCommand("color");
		colorButton.setMinimumSize(new Dimension(800,20));
		colorButton.addActionListener(this);
		colorButton.setForeground(theDREMGui.keyInputLabelColor);

		JPanel buttonPanel = new JPanel();

		buttonPanel.setAlignmentX(Component.LEFT_ALIGNMENT);
		buttonPanel.setBackground(Color.white);
		buttonPanel.add(hideButton); 
		buttonPanel.add(colorButton);

		JButton helpButton = new JButton(Util.createImageIcon("Help16.gif"));
		helpButton.addActionListener(this);
		helpButton.setActionCommand("help");
		buttonPanel.add(helpButton);

		buttonPanel.setMaximumSize(new Dimension(Integer.MAX_VALUE,20));
		add(buttonPanel);
	}

	/**
	 * Responds to changes in the significance threshold at which transcription factor 
	 * labels should appears 
	 */
	public void stateChanged(ChangeEvent e) 
	{
		JSlider source = (JSlider)e.getSource();

		if (!source.getValueIsAdjusting()) 
		{
			if (source == theSlider)
			{
				theDREMGui.dkeyinputpvalue = Math.pow(10,-source.getValue()/10.0);
				pvalLabel.setText("X = "+source.getValue()/10.0+
						"; score threshold is "+doubleToSz(theDREMGui.dkeyinputpvalue));
				int nsize = theDREMGui.hidesigList.size();
				for (int nindex = 0; nindex < nsize; nindex++)
				{
					DREMGui.SigInfoRec theSigInfoRec = 
						(DREMGui.SigInfoRec) theDREMGui.hidesigList.get(nindex);	      
					if (theSigInfoRec.ntype != 0 && theSigInfoRec.ntype != 777 && theSigInfoRec.ntype != 778)
					{
						double dprewidth =theSigInfoRec.theSigText.getWidth();
						theDREMGui.setSigText(theSigInfoRec.theSigTF, theSigInfoRec.theSigText, theSigInfoRec.ntype, theSigInfoRec.node);
						theSigInfoRec.border.translate(dprewidth-theSigInfoRec.theSigText.getWidth(),0);
					}
					else
					{	      
						theDREMGui.setSigText(theSigInfoRec.theSigTF, theSigInfoRec.theSigText, theSigInfoRec.ntype, theSigInfoRec.node);
					}
				}
			}
			else if (source == theSliderPercent)
			{
				theDREMGui.dsplitpercent = source.getValue()/100.0;
				percentLabel.setText(" Minimum Split % is "
						+nf2.format(100*theDREMGui.dsplitpercent));
				int nsize = theDREMGui.hidesigList.size();
				for (int nindex = 0; nindex < nsize; nindex++)
				{
					DREMGui.SigInfoRec theSigInfoRec = 
						(DREMGui.SigInfoRec) theDREMGui.hidesigList.get(nindex);	      
					if (theSigInfoRec.ntype == 1)  //split path type
					{
						double dprewidth =theSigInfoRec.theSigText.getWidth();
						theDREMGui.setSigText(theSigInfoRec.theSigTF, theSigInfoRec.theSigText, theSigInfoRec.ntype, theSigInfoRec.node);
						theSigInfoRec.border.translate(dprewidth-theSigInfoRec.theSigText.getWidth(),0);
					}
				}
			}
		}
	}

	/**
	 * Converts a double value to a string as formatted as displayed on this interface
	 */
	public static String doubleToSz(double dval)
	{
		String szexp;
		double dtempval = dval;
		int nexp = 0;

		NumberFormat nf3 = NumberFormat.getInstance(Locale.ENGLISH);
		nf3.setMinimumFractionDigits(3);
		nf3.setMaximumFractionDigits(3);     

		NumberFormat nf4 = NumberFormat.getInstance(Locale.ENGLISH);
		nf4.setMinimumFractionDigits(2);
		nf4.setMaximumFractionDigits(2); 

		if (dval <= 0)
		{
			szexp = "0.000";
		}
		else
		{
			while ((dtempval<0.9995)&&(dtempval>0))
			{
				nexp--;
				dtempval = dtempval*10;
			}
			dtempval = dval*Math.pow(10,-nexp);
			if (nexp < -3)
				szexp = nf4.format(dtempval)+" * 10^"+nexp;
			else
				szexp = nf3.format(dval);
		}

		return szexp;
	}

	/**
	 * Responds to changes with the radio buttons
	 */    
	public void itemStateChanged(ItemEvent e) 
	{
		// TODO these should be global constants
		if (splitButton.isSelected())
		{
			theDREMGui.nKeyInputType = 0;
		}
		else if (edgeSplitButton.isSelected())
		{
			theDREMGui.nKeyInputType = 1;
		} 
		else if (edgeFullButton.isSelected())
		{
			theDREMGui.nKeyInputType = 2;
		}
		else if(activityButton.isSelected())
		{
			theDREMGui.nKeyInputType = 777;
		}
		else if(firstActivityButton.isSelected())
		{
			theDREMGui.nKeyInputType = 778;
		}

		int nsize = theDREMGui.hidesigList.size();
		for (int nindex = 0; nindex < nsize; nindex++)
		{
			DREMGui.SigInfoRec theSigInfoRec = 
				(DREMGui.SigInfoRec) theDREMGui.hidesigList.get(nindex);

			if ( theSigInfoRec.ntype != 0 && theSigInfoRec.ntype != 777 && theSigInfoRec.ntype != 778)
			{
				double dprewidth =theSigInfoRec.theSigText.getWidth();
				theDREMGui.setSigText(theSigInfoRec.theSigTF, theSigInfoRec.theSigText, theSigInfoRec.ntype, theSigInfoRec.node);
				theSigInfoRec.border.translate(dprewidth-theSigInfoRec.theSigText.getWidth(),0);
			}
			else
			{      
				theDREMGui.setSigText(theSigInfoRec.theSigTF, theSigInfoRec.theSigText,  theSigInfoRec.ntype, theSigInfoRec.node);
			}

		}
	}

	/**
	 * Responds to buttons being pressed on the interface window
	 */
	public void actionPerformed(ActionEvent e) 
	{
		String szcommand = e.getActionCommand();


		if (szcommand.equals("help"))
		{
			String szMessage =
				"This window controls the Key Input labels appearing on the DREM output interface.   "+
				"Consult section 4.4 of the user manual for more details on this window. " +
				"TF activity scores do not use pvalues so the TFs shown will be those with " +
				"activity scores greater than or equal to the inverse of the threshold.";
			Util.renderDialog(theFrame,szMessage,-350,-100);
		}
		else if (szcommand.equals("color"))
		{
			Color newColor = JColorChooser.showDialog(
					this,
					"Choose Color",
					theDREMGui.keyInputLabelColor);

			if (newColor != null)
			{
				theDREMGui.keyInputLabelColor=newColor;
				colorButton.setForeground(newColor);
				int nsize = theDREMGui.hidesigList.size();

				for (int nindex = 0; nindex < nsize; nindex++)
				{
					DREMGui.SigInfoRec theSigInfoRec = 
						(DREMGui.SigInfoRec) theDREMGui.hidesigList.get(nindex);	      

					theSigInfoRec.theSigText.setTextPaint(newColor);
				}
			}
		}
		else if (szcommand.equals("hide"))
		{
			int nsize = theDREMGui.hidesigList.size();

			if (theDREMGui.bshowkeyinputs)
			{
				hideButton.setText("Show Key TF Labels");
				theDREMGui.bshowkeyinputs = false;

				for (int nindex = 0; nindex < nsize; nindex++)
				{
					DREMGui.SigInfoRec theSigInfoRec = 
						(DREMGui.SigInfoRec) theDREMGui.hidesigList.get(nindex);	      

					PPath rect = (PPath) theSigInfoRec.theSigText.getParent();
					rect.setVisible(false);
					rect.setPickable(false);
					theSigInfoRec.theSigText.setVisible(false);
					theSigInfoRec.theSigText.setPickable(false);
				}
			}
			else
			{
				hideButton.setText("Hide Key TF Labels");
				theDREMGui.bshowkeyinputs = true;
				for (int nindex = 0; nindex < nsize; nindex++)
				{
					DREMGui.SigInfoRec theSigInfoRec = 
						(DREMGui.SigInfoRec) theDREMGui.hidesigList.get(nindex);
					PPath rect = (PPath) theSigInfoRec.theSigText.getParent();
					rect.setVisible(true);
					rect.setPickable(true);
					theSigInfoRec.theSigText.setVisible(true);
					theSigInfoRec.theSigText.setPickable(true);
				}
			}
		}
	}
}