package edu.cmu.cs.sb.drem;
import edu.cmu.cs.sb.core.*;

import javax.swing.ToolTipManager;
import java.awt.geom.*;
import edu.umd.cs.piccolo.PCanvas;
import edu.umd.cs.piccolo.PCamera;
import edu.umd.cs.piccolo.PNode;
import edu.umd.cs.piccolo.PLayer;
import edu.umd.cs.piccolo.event.PZoomEventHandler;
import edu.umd.cs.piccolo.event.PBasicInputEventHandler;
import edu.umd.cs.piccolo.event.PInputEvent;
import edu.umd.cs.piccolo.nodes.PPath;
import edu.umd.cs.piccolox.PFrame;
import edu.umd.cs.piccolo.nodes.PText;
import edu.umd.cs.piccolo.nodes.PImage;
import edu.umd.cs.piccolox.nodes.PLine;
import edu.umd.cs.piccolo.util.*;
import java.awt.geom.*;
import java.util.*;
import javax.swing.*;
import javax.swing.event.*;
import java.awt.*;
import java.awt.event.*;
import java.text.*;
import java.io.*;
import javax.imageio.*;
import java.awt.image.*;

/**
 * Class for the main interface window of a DREM regulatory map
 */
public class DREMGui extends PFrame implements ComponentListener
{

	boolean battachlabels = true;
	static int SCREENWIDTH = 800;
	static int SCREENHEIGHT = 600;
	static int LEFTBUFFER = 25;
	static int RIGHTBUFFER =75;
	static double TOPBUFFER = 10;
	static double BOTTOMBUFFER = 80;
	double REALHEIGHT =0; //not constant
	static int TICKLENGTH = 15;
	static double INITKEYLEVEL = .001;
	static int BUTTONRECWIDTH = 110;
	double INITLEFT;
	boolean binit = true;
	static Object datasetlock = new Object();
	static Color buttonColor = new Color(255,246,143);
	DREM_Timeiohmm theTimeiohmm;
	double[][] data;
	int[][] pmavalues;
	String[] genenames;
	PCanvas canvas;
	PLine[] plArray;
	PLayer theLayer;
	ArrayList hideList = new ArrayList();//array list of nodes and edges whose visibility can be toggled
	ArrayList hidesigList;
	ArrayList hidegolabelsList;
	ArrayList hidetfsetlabelsList;
	ArrayList hidepredictlabelsList;
	ArrayList hidegenesetlabelsList;
	TreeSet circleSet = new TreeSet(new CircleRecCompare());
	Hashtable htColors;
	Hashtable htHidden = new Hashtable();
	HashSet htRequired = new HashSet();
	Hashtable htNodes = new Hashtable();
	Hashtable htPathToLineSet = new Hashtable();
	Hashtable htTextVisible;
	JFrame filterStaticFrame;
	JFrame filterGOFrame;
	JFrame defineFrame;
	JFrame predictFrame;    
	JFrame keyinputFrame;
	JFrame yscaleFrame;
	JFrame saveDREMFrame;
	JFrame saveModelFrame;
	SelectedNodeRec theSelectedRec = new SelectedNodeRec();
	DREMGui_GOFilter theGOFilter;
	DREMGui_InterfaceOptions theYScalegui;
	DREM_Timeiohmm.Treenode treecopy;

	PNode genetableButton;
	PNode setButton;
	PNode gotableButton;
	PNode staticButton;
	PNode gofilterButton;
	PNode scaleButton;
	PNode predictButton;
	PNode hideButton;
	PNode timeButton;
	PNode siginButton;
	PNode saveDREMButton;
	PNode saveModelButton;
	PText hideText;
	PImage helpButton;
	int nheight;
	int nwidth;
	int ncolortime=0;

	int ncircleID= 0;   
	double dheightunits;
	double dminheightunits;
	double dwidthunits;
	double dnodek = 1; 
	double[] dwidthunitsInterval;
	double[] dwidthunitsCum;
	double dmax;
	double dmin;
	double dkeyinputpvalue = INITKEYLEVEL;
	double dsplitpercent=0;
	Color SPLITCOLOR = Color.green;

	boolean bapplygenesetlabels=false;;
	boolean bapplytfsetlabels=false;
	boolean bapplygolabels=false;

	boolean bholdedge = false;
	boolean bfiltergo=false;
	boolean bfiltergeneset=false;
	boolean bfilterinput=false;
	boolean brealXaxis = false;
	PText filterText;

	Color keyInputLabelColor = Color.black;
	Color goLabelColor = Color.black;
	Color predictLabelColor = Color.black;
	Color genesetLabelColor = Color.orange;
	Color tfLabelColor = Color.red;

	boolean binvalidreal = false;//whether sampling rate input is valid


	boolean[] bGOVisible;
	boolean[] bSetVisible;
	boolean[] bPathVisible;
	boolean[] bTFVisible;
	boolean bglobalVisible;

	boolean blowestbelow = true;
	boolean blowestabove = true;
	boolean bmiddlebelow = true;
	boolean bmiddleabove = true;
	boolean bhighestbelow = true;
	boolean bhighestabove = true;
	boolean blowbelow = true;
	boolean blowabove = true;

	boolean bglobalnode = true;
	boolean bshowpredict = false;
	boolean bshowkeyinputs = true;
	boolean bshowgolabels = true;
	boolean bshowtfsetlabels = true;
	boolean bshowgenesetlabels = true;
	int nKeyInputType;
	NumberFormat nf3;

	int numcolor = 0;
	int nparentcolorindex = 0;
	Hashtable htColorIDtoLinesList = new Hashtable();
	double dscaley;
	double dscalex = 0.5;

	Hashtable htColorIDtoColor = new Hashtable();
	Hashtable htLineIDtoColorID = new Hashtable();
	Hashtable htColorIDtoCircleList = new Hashtable();

	int[] storedbestpath;
	boolean bsavedchange;

	double dinitmin;
	double dinitmax;
	Color[] lineColorsA;

	Color edgeColorsTriples[][] = {
			{new Color((float) 32/255,(float) 178/255,(float) 170/255,(float) 1),
				new Color((float) 255/255,(float) 0/255,(float) 0/255,(float) 1),//red
				new Color((float) 255/255,(float) 0/255,(float) 255/255,(float) 1)},//pink

				{new Color((float) 255/255,(float) 128/255,(float) 0/255,(float) 1),//orange
					new Color((float) 128/255,(float) 128/255,(float) 128/255,(float) 1),//gray 
					new Color((float) 205/255,(float) 51/255,(float) 51/255,(float) 1)},//brown3

					{new Color((float) 0/255,(float) 153/255,(float) 153/255,(float) 1),
						new Color((float) 139/255,(float) 115/255,(float) 85/255,(float) 1),//burlywood4
						new Color((float) 204/255,(float) 0/255,(float) 204/255,(float) 1)},//beige

						{new Color((float) 46/255,(float) 139/255,(float) 87/255,(float) 1),//sea green
							new Color((float) 128/255,(float)0/255,(float) 255/255,(float) 1),
							new Color((float) 128/255,(float) 0/255,(float) 128/255,(float) 1)},//purple},//light green

							{new Color((float) 205/255,(float) 91/255,(float) 69/255,(float) 1),//coral
								new Color((float) 14/255,(float) 59/255,(float) 59/255,(float) 1),//rosy brown3
								new Color((float) 255/255,(float) 128/255,(float) 255/255,(float) 1)},

								{new Color((float) 85/255,(float) 107/255,(float) 47/255,(float) 1),//dark olive green
									new Color((float) 128/255,(float) 128/255,(float) 255/255,(float) 1),//blue
									new Color((float) 255/255,(float) 0/255,(float) 128/255,(float) 1)},

									{new Color((float) 162/255,(float) 162/255,(float) 104/255,(float) 1),//gray
										new Color((float) 0/255,(float) 102/255,(float) 102/255,(float) 1),
										new Color((float) 255/255,(float) 102/255,(float) 102/255,(float) 1)},//pink

										{new Color((float) 153/255,(float) 153/255,(float) 255/255,(float) 1),
											new Color((float) 153/255,(float) 102/255,(float)0/255,(float) 1),
											new Color((float) 0/255,(float) 102/255,(float) 102/255,(float) 1)},

											{new Color((float) 255/255,(float) 0/255,(float) 102/255,(float) 1),
												new Color((float) 51/255,(float) 51/255,(float) 51/255,(float) 1),
												new Color((float) 102/255,(float) 255/255,(float) 102/255,(float) 1)},

												{new Color((float) 0/255,(float) 102/255,(float) 153/255,(float) 1),
													new Color((float) 155/255,(float) 0/255,(float) 0/255,(float) 1),
													new Color((float) 102/255,(float) 0/255,(float) 102/255,(float) 1)},

													{new Color((float) 102/255,(float) 0/255,(float) 102/255,(float) 1),
														new Color((float) 146/255,(float) 92/255,(float) 62/255,(float) 1),
														new Color((float) 255/255,(float) 204/255,(float) 204/255,(float) 1)},//pink

														{new Color((float) 209/255,(float) 133/255,(float) 0/255,(float) 1),//gray
															new Color((float) 73/255,(float) 39/255,(float) 84/255,(float) 1),
															new Color((float) 241/255,(float) 104/255,(float) 104/255,(float) 1)},//pink

															{new Color((float) 170/255,(float) 15/255,(float) 163/255,(float) 1),
																new Color((float) 205/255,(float) 102/255,(float) 29/255,(float) 1),//chocoloate
																new Color((float) 205/255,(float) 92/255,(float) 92/255,(float) 1)}};//,indian red


	/**
	 * Class constructor
	 */
	public DREMGui(DREM_Timeiohmm theTimeiohmm,DREM_Timeiohmm.Treenode treecopy,
			boolean brealXaxisDEF,double dYaxisDEF, double dXaxisDEF,
			int nKeyInputTypeDEF,double dKeyInputXDEF, double dpercentDEF, String szFinal,
			double dnodekDEF)
	{
		super("DREM - Main Interface "+szFinal, false, null);
		this.treecopy = treecopy;
		storedbestpath = new int[theTimeiohmm.storedbestpath.length];
		bsavedchange = theTimeiohmm.bsavedchange;
		for (int nindex = 0; nindex < storedbestpath.length; nindex++)
		{
			storedbestpath[nindex] = theTimeiohmm.storedbestpath[nindex];
		}

		synchronized(datasetlock)
		{
			if (DREM_Timeiohmm.BDEBUG)
			{
				System.out.println("made it!!!!!!!!!");
			}
			dscaley = dYaxisDEF;
			dscalex = dXaxisDEF;
			dnodek = dnodekDEF;
			brealXaxis = brealXaxisDEF;
			nKeyInputType = nKeyInputTypeDEF;
			dsplitpercent = dpercentDEF/100;
			//rounding to the nearest tenth
			dkeyinputpvalue = Math.pow(10,-(Math.round(10*dKeyInputXDEF))/10.0);  //2.51 --> 25.1 --> 2.5

			this.theTimeiohmm = theTimeiohmm;
			data =  theTimeiohmm.theDataSet.data;
			pmavalues =theTimeiohmm.theDataSet.pmavalues;
			genenames = theTimeiohmm.theDataSet.genenames;
			bPathVisible = new boolean[data.length];
			bTFVisible = new boolean[data.length];
			bGOVisible = new boolean[data.length];
			bSetVisible = new boolean[data.length];

			for (int nrow = 0; nrow < data.length; nrow++)
			{
				bTFVisible[nrow] = true;
				bPathVisible[nrow] = true;
				bGOVisible[nrow] = true;
				bSetVisible[nrow] = true;
			}
			plArray = new PLine[data.length];
			bglobalVisible = true;

			htColors = new Hashtable();

			addComponentListener(this);
			nf3 = NumberFormat.getInstance(Locale.ENGLISH);
			nf3.setMinimumFractionDigits(3);
			nf3.setMaximumFractionDigits(3);

			dmax =  Math.abs(data[0][0]);
			dmin =  Math.abs(data[0][0]);
			double[] dmaxRow = new double[data.length];
			double[] dminRow = new double[data.length];

			if (DREM_Timeiohmm.BDEBUG)
			{
				System.out.println(dmax+"\t"+dmin);
			}

			for (int nrow = 0; nrow < data.length; nrow++)
			{         
				dmaxRow[nrow] = 0;
				dminRow[nrow] = 0;
				for (int ncol = 0; ncol < data[0].length; ncol++)
				{
					if (pmavalues[nrow][ncol] != 0)
					{
						if (data[nrow][ncol] < dminRow[nrow])
						{
							dminRow[nrow] = data[nrow][ncol];
						}

						if (data[nrow][ncol] > dmaxRow[nrow])
						{
							dmaxRow[nrow] = data[nrow][ncol];
						}
					}
				} 
			}

			Arrays.sort(dminRow);
			Arrays.sort(dmaxRow);

			dmin = dminRow[0];
			dmax = dmaxRow[dmaxRow.length-1];
			dinitmin = dmin;
			dinitmax = dmax;
			if (DREM_Timeiohmm.BDEBUG)
			{
				System.out.println(dinitmin+" "+dinitmax);
			}

			//potential synchronization issue with getting height
			nheight = getHeight();
			REALHEIGHT =  (nheight - TOPBUFFER - BOTTOMBUFFER);
			dminheightunits = REALHEIGHT/(dmax-dmin);
			treecopy.dsigma = (dmax-dmin)*0.02;
			binit = true;

			if (DREM_Timeiohmm.BDEBUG)
			{
				System.out.println("calling notifyall");
			}

			datasetlock.notifyAll();
		}
	}

	/**
	 * Sets the colors of the individual genes
	 */
	public void setGeneColors()
	{
		Color currColor;

		if (ncolortime == 0)
		{
			for (int nrow = 0; nrow < plArray.length; nrow++)
			{
				plArray[nrow].setStrokePaint(lineColorsA[nrow]);
			}
		}
		else
		{
			int nval = (int) Math.pow(theTimeiohmm.nmaxchild,data[0].length-ncolortime);

			for (int nrow = 0; nrow < plArray.length; nrow++)
			{
				int nid = nval*(storedbestpath[nrow]/nval);
				Integer colorIDobj = ((Integer) htLineIDtoColorID.get(ncolortime+";"+nid));
				if ((colorIDobj==null)||(htColorIDtoColor==null))
				{
					if (DREM_Timeiohmm.BDEBUG)
					{
						System.out.println(ncolortime+";"+nid+"\t"+colorIDobj+"\t"+htColorIDtoColor);
					}
				}

				currColor = (Color) htColorIDtoColor.get(colorIDobj);
				plArray[nrow].setStrokePaint(currColor);
			}
		}
	}

	/**
	 * Tool tip text updater.
	 */
	private static class ToolTipTextUpdater
	extends PBasicInputEventHandler
	{
		/** Canvas. */
		private PCanvas canvas;

		/** Node. */
		private PNode node;
		String szText;


		/**
		 * Create a new tool tip text updater for the specified
		 * canvas and node.
		 *
		 * @param canvas canvas
		 * @param node node
		 */
		public ToolTipTextUpdater(final PCanvas canvas,
				final PNode node, final String szText)
		{
			this.canvas = canvas;
			this.node = node;
			this.szText = szText;
		}


		/** @see PBasicInputEventHandler */
		public void mouseEntered(final PInputEvent e)
		{
			if (node.getVisible())
			{
				canvas.setToolTipText(szText);
			}
			else
			{
				canvas.setToolTipText(null);
			}

		}

		/** @see PBasicInputEventHandler */
		public void mouseExited(final PInputEvent e)
		{
			canvas.setToolTipText(null);
		}
	}

	/**
	 * Closes the frame windows
	 */
	public void closeWindows()
	{
		if (filterStaticFrame != null)
		{
			filterStaticFrame.setVisible(false);
			filterStaticFrame.dispose();
			filterStaticFrame = null;
		} 
		if (filterGOFrame != null)
		{
			filterGOFrame.setVisible(false);
			filterGOFrame.dispose();
			filterGOFrame = null;
		} 
		if (defineFrame != null)
		{
			defineFrame.setVisible(false);
			defineFrame.dispose();
			defineFrame = null;
		} 
		if (predictFrame != null)
		{
			predictFrame.setVisible(false);
			predictFrame.dispose();
			predictFrame = null;
		} 
		if (keyinputFrame != null)
		{
			keyinputFrame.setVisible(false);
			keyinputFrame.dispose();
			keyinputFrame = null;
		} 
		if (yscaleFrame != null)
		{
			yscaleFrame.setVisible(false);
			yscaleFrame.dispose();
			yscaleFrame = null;
		} 
	}



	/**
	 * Saves the edge color selections to an ouput file
	 */
	public String saveColors()
	{
		StringBuffer szbuf =  new StringBuffer();
		float[] f = new float[4];
		for (int nindex = 0; nindex < numcolor; nindex++)
		{
			Color currColor = (Color) htColorIDtoColor.get(new Integer(nindex));
			currColor.getRGBComponents(f);
			szbuf.append(f[0]+"\t"+f[1]+"\t"+f[2]+"\t"+f[3]+"\n");
		}

		return szbuf.toString();
	}


	/**
	 * Undo gene display based on the GO selection
	 */
	public void unselectGO()
	{
		for (int ngene = 0; ngene < bGOVisible.length; ngene++)
		{
			bGOVisible[ngene] = true;
			boolean bvisible = bglobalVisible&&bSetVisible[ngene]&&
			bPathVisible[ngene]&&bTFVisible[ngene];
			plArray[ngene].setVisible(bvisible); 
			plArray[ngene].setPickable(bvisible); 
		}
		bfiltergo = false;
		setFilterText();
	}

	/**
	 * Changes the gene visibility based on the selected GO category szSelectedGO
	 */
	public void selectGO(String szSelectedGO)
	{
		GoAnnotations tga = theTimeiohmm.theDataSet.tga;

		tga.szSelectedGO = szSelectedGO;
		String[] genenames = theTimeiohmm.theDataSet.genenames;
		for (int ngene = 0; ngene < genenames.length; ngene++)
		{
			HashSet hsGO = tga.labelsForID(genenames[ngene]);
			if (hsGO.contains(szSelectedGO))
			{
				bGOVisible[ngene] = true;
				boolean bvisible = bglobalVisible&& bSetVisible[ngene]&&
				bPathVisible[ngene]&&bTFVisible[ngene];
				plArray[ngene].setVisible(bvisible); 
				plArray[ngene].setPickable(bvisible); 
			}
			else
			{
				bGOVisible[ngene] = false;
				plArray[ngene].setVisible(false); 
				plArray[ngene].setPickable(false); 
			}
		}
		bfiltergo = true;
		setFilterText();
	}

	/**
	 * Generates the string in filterText which displays information as to how the genes were selecte
	 */
	public void setFilterText()
	{
		String sz = "Genes selected based on ";

		int ncount = 0;
		if (bfiltergo) ncount++;
		if (bfiltergeneset) ncount++;
		if (bfilterinput) ncount++;

		boolean bpath = theSelectedRec.theCircleID!=null;
		if (bpath)
		{
			sz += "a path constraint";

			if (ncount == 1)
				sz += " and ";
			else if (ncount >= 2)
				sz += ", ";
		}

		if (bfiltergo)
		{
			GoAnnotations tga = theTimeiohmm.theDataSet.tga;
			String szGO = "the GO category "+((GoAnnotations.Rec) tga.htGO.get(tga.szSelectedGO)).sztermName;
			if (bfiltergeneset)
			{
				if (bfilterinput)
				{
					sz+= "a gene set, TF input, and "+szGO;
				}
				else
				{
					if (bpath)
					{
						sz+= "a gene set, and "+szGO;
					}
					else
					{
						sz+= "a gene set and "+szGO;
					}
				}
			}
			else
			{
				if (bfilterinput)
				{
					if (bpath)
					{
						sz += "TF input, and "+ szGO;  
					} 
					else
					{
						sz += "TF input and "+ szGO;  
					}

				}
				else
				{
					sz += szGO;            
				}
			}
		}
		else
		{
			if (bfiltergeneset)
			{
				if (bfilterinput)
				{
					if (bpath)
					{
						sz+= "a gene set, and TF input";
					}
					else
					{
						sz+= "a gene set and TF input";
					}
				}
				else
				{
					sz+= "a gene set";
				}
			}
			else
			{
				if (bfilterinput)
				{
					sz+= "TF input";
				}
				else
				{
					if (!bpath)
					{
						sz = "";
					}
				}
			}
		}

		filterText.setText(sz);
	}

	/**
	 * Controls the display of the tex showing the filter information
	 */
	public void renderFilterText()
	{
		filterText = new PText();
		setFilterText();
		filterText.translate(LEFTBUFFER+15, 2);
		filterText.setPickable(false);
		filterText.setFont(new Font("times",Font.PLAIN,14));
		if (!bglobalnode)
		{
			filterText.setVisible(false);
			filterText.setPickable(false);
		}
		canvas.getCamera().addChild(filterText);
	}


	/**
	 * Empty method
	 */
	public void componentHidden(ComponentEvent e) {} 

	/**
	 * Empty method
	 */
	public void componentMoved(ComponentEvent e) {}  

	/**
	 * Calls drawmain
	 */
	public void componentShown(ComponentEvent e) 
	{
		drawmain();
	}

	/**
	 * Calls drawmain
	 */
	public void componentResized(ComponentEvent e)
	{
		drawmain();    
	}

	/**
	 * Sets the screen size
	 */         
	public void beforeInitialize()
	{
		if (DREM_Timeiohmm.BDEBUG)
		{
			System.out.println("enter beforeinitialize");
		}
		setSize(SCREENWIDTH,SCREENHEIGHT);
		if (DREM_Timeiohmm.BDEBUG)
		{
			System.out.println("leave beforeinitialize");
		}
	}

	/**
	 * calls drawmain
	 */
	public void initialize()
	{
		drawmain();
	}

	/**
	 * Responsible for laying out the main interface window
	 */
	public void drawmain()
	{
		synchronized(datasetlock)
		{
			while (!binit)
			{
				try
				{
					datasetlock.wait();
				} 
				catch (InterruptedException e) 
				{
				}
			}

			hideList = new ArrayList();   
			ncircleID = 0;
			numcolor = 0;
			nparentcolorindex = 0;
			htColorIDtoLinesList = new Hashtable();
			circleSet = new TreeSet(new CircleRecCompare());

			htColors = new Hashtable();
			htNodes = new Hashtable();
			htPathToLineSet = new Hashtable();
			if (hidesigList != null)
			{
				htTextVisible = new Hashtable();
				int nsize = hidesigList.size();
				for (int nindex = 0; nindex < nsize; nindex++)
				{
					SigInfoRec theSigInfoRec =(SigInfoRec) hidesigList.get(nindex);
					htTextVisible.put(theSigInfoRec.ndepth+";"+theSigInfoRec.npathscore+
					";"+theSigInfoRec.ntype, Boolean.valueOf(theSigInfoRec.theSigText.getVisible()));
				}
			}
			else 
			{
				htTextVisible = null;
			}
			hidesigList = new ArrayList();
			hidegolabelsList = new ArrayList();
			hidetfsetlabelsList = new ArrayList();
			hidegenesetlabelsList = new ArrayList();
			hidepredictlabelsList = new ArrayList();

			canvas = getCanvas();
			ToolTipManager.sharedInstance().registerComponent(canvas);
			theLayer = canvas.getLayer();
			PCamera theCamera = canvas.getCamera();
			theLayer.removeAllChildren();
			theCamera.removeAllChildren();

			dmin= dinitmin/dscaley;
			dmax= dinitmax/dscaley;
			double dmaxmindiff = (dmax-dmin);

			nheight = getHeight();
			REALHEIGHT =  (nheight - TOPBUFFER - BOTTOMBUFFER);
			dheightunits =REALHEIGHT/dmaxmindiff;
			nwidth = getWidth();      

			int numintervals = data[0].length - 1;
			dwidthunits = (nwidth-LEFTBUFFER-RIGHTBUFFER)*dscalex;

			dwidthunitsInterval = new double[numintervals];
			dwidthunitsCum = new double[numintervals+1];
			dwidthunitsCum[0] = 0;

			boolean bincreasing = true;

			try
			{
				double dprev = Util.removeUnits(theTimeiohmm.theDataSet.dsamplemins[0]);
				for (int ni = 0; ni < dwidthunitsInterval.length; ni++)
				{
					double dnext = Util.removeUnits(theTimeiohmm.theDataSet.dsamplemins[ni+1]);
					dwidthunitsInterval[ni] = dnext-dprev;
					if (dwidthunitsInterval[ni] < 0)
					{
						bincreasing = false;
						brealXaxis = false;
						binvalidreal = true;
					}
					dwidthunitsCum[ni+1] = dwidthunitsCum[ni] + dwidthunitsInterval[ni];
					dprev = dnext;
				}
				double dtotal = dwidthunitsCum[dwidthunitsCum.length-1]; 
				for (int ni = 1; ni < dwidthunitsCum.length; ni++)
				{
					dwidthunitsInterval[ni-1] = dwidthunits * dwidthunitsInterval[ni-1]/dtotal;
					dwidthunitsCum[ni]= dwidthunits * dwidthunitsCum[ni]/dtotal;
				}
			}
			catch (IllegalArgumentException iae)
			{
				binvalidreal = true;
				brealXaxis = false;
				for (int ni = 1; ni < dwidthunitsCum.length; ni++)
				{
					dwidthunitsInterval[ni-1] = dwidthunits/numintervals;
					dwidthunitsCum[ni]= ni*dwidthunits/numintervals;
				}
			}


			if ((!brealXaxis)||(!bincreasing))
			{
				for (int ni = 1; ni < dwidthunitsCum.length; ni++)
				{
					dwidthunitsInterval[ni-1] = dwidthunits/numintervals;
					dwidthunitsCum[ni]= ni*dwidthunits/numintervals;
				}
			}

			PLine plxaxis = new PLine();
			double doriginy = REALHEIGHT + dheightunits*dmin;
			plxaxis.addPoint(0, LEFTBUFFER,doriginy);
			plxaxis.addPoint(1, dscalex*nwidth-RIGHTBUFFER/2, doriginy);
			plxaxis.setStroke(new BasicStroke(4));
			plxaxis.setStrokePaint(Color.gray);
			theLayer.addChild(plxaxis);

			for (int ncol = 0; ncol < theTimeiohmm.theDataSet.dsamplemins.length; ncol++)
			{    
				PLine pltick = new PLine();
				double dxp =  dwidthunitsCum[ncol]+LEFTBUFFER;
				pltick.addPoint(0, dxp, doriginy);
				pltick.addPoint(1, dxp, doriginy+TICKLENGTH);
				pltick.setStroke(new BasicStroke(4));
				pltick.setStrokePaint(Color.gray);
				theLayer.addChild(pltick);

				PText tickLabel = new PText(""+theTimeiohmm.theDataSet.dsamplemins[ncol]);
				tickLabel.translate(dxp+4,doriginy+4);
				theLayer.addChild(tickLabel);
			}

			int npower=(int) Math.floor(Math.log(Math.max(Math.abs(dmin),dmax))/Math.log(10));
			NumberFormat nftick = NumberFormat.getInstance(Locale.ENGLISH);
			if (npower < 0)
			{
				nftick.setMaximumFractionDigits(-npower);
				nftick.setMinimumFractionDigits(-npower);
			}

			double drealincrementy = Math.pow(10,npower);
			double dincrementy =drealincrementy *dheightunits;
			double drealy = -drealincrementy;
			double dticky = doriginy+dincrementy;

			while (drealy < dinitmax)
			{
				dticky -= dincrementy;
				drealy += drealincrementy;
				PLine pltick = new PLine();
				pltick.addPoint(0, LEFTBUFFER-TICKLENGTH, dticky);
				pltick.addPoint(1, LEFTBUFFER, dticky);
				pltick.setStroke(new BasicStroke(4));
				pltick.setStrokePaint(Color.gray);
				theLayer.addChild(pltick);
				String sztickval = nftick.format(drealy);
				PText tickLabel = new PText(sztickval);
				int nextra;

				tickLabel.translate(LEFTBUFFER-TICKLENGTH-4*sztickval.length()+6,dticky-17);
				theLayer.addChild(tickLabel);
			}

			PLine plyaxis = new PLine();
			plyaxis.addPoint(0, LEFTBUFFER,dticky);

			dticky = doriginy;
			drealy =0;

			while (drealy > dinitmin)
			{
				dticky += dincrementy;
				drealy -= drealincrementy;
				PLine pltick = new PLine();
				pltick.addPoint(0, LEFTBUFFER-TICKLENGTH, dticky);
				pltick.addPoint(1, LEFTBUFFER, dticky);
				pltick.setStroke(new BasicStroke(4));
				pltick.setStrokePaint(Color.gray);
				theLayer.addChild(pltick);
				String sztickval = nftick.format(drealy);
				PText tickLabel = new PText(sztickval);
				tickLabel.translate(LEFTBUFFER-TICKLENGTH-3*sztickval.length()+3,dticky-17);
				theLayer.addChild(tickLabel);
			}

			plyaxis.addPoint(1, LEFTBUFFER, dticky);
			plyaxis.setStroke(new BasicStroke(4));
			plyaxis.setStrokePaint(Color.gray);
			theLayer.addChild(plyaxis);

			Random lineColors = new Random(74312);
			if (lineColorsA == null)
			{
				lineColorsA = new Color[data.length];
			}

			for (int nrow = 0; nrow < data.length; nrow++)
			{ 
				Color lineColor;

				if (lineColorsA[nrow] == null)
				{

					lineColor = new Color(lineColors.nextInt(226),
							lineColors.nextInt(226),
							lineColors.nextInt(226));
					lineColorsA[nrow] = lineColor;

					String szr = genenames[nrow];	 
				}

				PLine pl = new PLine();

				pl.addInputEventListener(new ToolTipTextUpdater(canvas, pl,genenames[nrow]));

				int nadded = 0;
				for (int ncol = 0; ncol < data[nrow].length; ncol++)
				{    
					if (pmavalues[nrow][ncol] != 0)
					{
						double dxp =  dwidthunitsCum[ncol]+LEFTBUFFER;
						double dyp = REALHEIGHT - dheightunits*(data[nrow][ncol]-dmin);
						pl.addPoint(nadded,dxp,dyp);
						nadded++;
					}           
				}

				ArrayList geneList = (ArrayList) htPathToLineSet.get(new Integer(storedbestpath[nrow]));

				if (geneList == null)
				{
					geneList = new ArrayList();
				}

				geneList.add(new Integer(nrow));
				htPathToLineSet.put( new Integer(storedbestpath[nrow]) ,geneList);
				plArray[nrow] = pl;
				boolean bvisible = bglobalVisible&&bPathVisible[nrow]&&bTFVisible[nrow]&&
				bGOVisible[nrow]&&bSetVisible[nrow];	       

				pl.setVisible(bvisible);
				pl.setPickable(bvisible);
				theLayer.addChild(pl);
			}

			drawNodes(treecopy,0,0,0,0,0,null,null,-1,nparentcolorindex,treecopy.nminparentlevel);
			setGeneColors();
			Iterator circleItr = circleSet.iterator();
			while (circleItr.hasNext())
			{
				CircleRec theCircleRec = (CircleRec) circleItr.next();
				theLayer.addChild(theCircleRec.circle);

				if (!bglobalnode)
				{
					theCircleRec.circle.setVisible(false);
					theCircleRec.circle.setPickable(false);
				}
			}  

			int npredictsize =  hidepredictlabelsList.size();
			for (int nindex = 0; nindex < npredictsize; nindex++)
			{
				theLayer.addChild((PText) hidepredictlabelsList.get(nindex));
			}

			int ntfsetsize =  hidetfsetlabelsList.size();
			for (int nindex = 0; nindex < ntfsetsize; nindex++)
			{
				theLayer.addChild((PText) hidetfsetlabelsList.get(nindex));
			}

			int ngenesetsize =  hidegenesetlabelsList.size();
			for (int nindex = 0; nindex < ngenesetsize; nindex++)
			{
				theLayer.addChild((PText) hidegenesetlabelsList.get(nindex));
			}

			int ngosize =  hidegolabelsList.size();
			for (int nindex = 0; nindex < ngosize; nindex++)
			{
				theLayer.addChild((PText) hidegolabelsList.get(nindex));
			}

			genetableButton = PPath.createRectangle((float) 0.0,(float) 0.0,(float) BUTTONRECWIDTH,(float) 18.0);
			PText genetableText = new PText("Gene Table");
			genetableText.setFont(new Font("times",Font.PLAIN,12));
			genetableText.translate(26,2);
			genetableButton.setPaint(buttonColor);
			genetableButton.addChild(genetableText);   

			genetableButton.addInputEventListener(new PBasicInputEventHandler() 
			{
				public void mousePressed(PInputEvent event) 
				{
					javax.swing.SwingUtilities.invokeLater(new Runnable() 
					{
						public void run() 
						{
							JFrame frame = new JFrame("Table of Selected Genes");
							frame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
							frame.setLocation(20,50);
							boolean[] bdisplay = new boolean[plArray.length];
							for (int ngene = 0; ngene < plArray.length; ngene++)
							{
								bdisplay[ngene] =  bPathVisible[ngene]&&bTFVisible[ngene]&&
								bGOVisible[ngene]&&bSetVisible[ngene];
							}

							DREMGui_GeneTable newContentPane = 
								new DREMGui_GeneTable(frame,theTimeiohmm.theDataSet,theTimeiohmm,
										theTimeiohmm.bindingSignGene,
										theTimeiohmm.bindingGeneindex,
										theTimeiohmm.bindingSignTF,
										theTimeiohmm.bindingTFindex,
										theTimeiohmm.tfNames,bdisplay);
							//need to bindingpvalGeneIndex reference
							newContentPane.setOpaque(true); //content panes must be opaque
							frame.setContentPane(newContentPane);
							//Display the window.
							frame.pack();
							frame.setVisible(true);
						}
					});
				}
			});

			final DREMGui ftheDREMGui = this;

			setButton = PPath.createRectangle((float) 0.0,(float) 0.0,(float) BUTTONRECWIDTH,(float) 18.0);

			PText setText = new PText("Select by Gene Set");
			setText.setFont(new Font("times",Font.PLAIN,12));
			setText.translate(4,2);
			setButton.setPaint(buttonColor);
			setButton.addChild(setText);   

			setButton.addInputEventListener(new PBasicInputEventHandler() 
			{
				public void mousePressed(PInputEvent event) 
				{
					javax.swing.SwingUtilities.invokeLater(new Runnable() 
					{
						public void run() 
						{
							if (defineFrame == null)
							{
								defineFrame = new JFrame("Select Genes Based on Defined Gene Set");
								defineFrame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
								defineFrame.setLocation(400,200);
								DREMGui_DefineGeneSet newContentPane = 
									new DREMGui_DefineGeneSet(defineFrame,theTimeiohmm.theDataSet.tga,
											ftheDREMGui,ftheDREMGui.treecopy);
								newContentPane.setOpaque(true); 
								//content panes must be opaque
								defineFrame.setContentPane(newContentPane);

								//Display the window.
								defineFrame.pack();
							}
							else
							{
								defineFrame.setExtendedState(Frame.NORMAL);
							}
							defineFrame.setVisible(true);
						}
					});
				}
			});

			gotableButton = PPath.createRectangle((float) 0.0,(float) 0.0,(float) BUTTONRECWIDTH,(float) 18.0);
			PText gotableText = new PText("GO Table");
			gotableText.setFont(new Font("times",Font.PLAIN,12));
			gotableText.translate(30,2);
			gotableButton.setPaint(buttonColor);
			gotableButton.addChild(gotableText);   

			gotableButton.addInputEventListener(new PBasicInputEventHandler() 
			{
				public void mousePressed(PInputEvent event) 
				{
					javax.swing.SwingUtilities.invokeLater(new Runnable() 
					{
						public void run() 
						{
							JFrame frame = new JFrame("GO Enrichment for Selected Genes");
							frame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
							frame.setLocation(20,50);
							frame.setSize(826,509);
							frame.setVisible(true);

							boolean[] bdisplay = new boolean[plArray.length];
							for (int ngene = 0; ngene < plArray.length; ngene++)
							{
								bdisplay[ngene] =  bPathVisible[ngene]&&bTFVisible[ngene]&&
								bGOVisible[ngene]&&bSetVisible[ngene];
							}
							DREMGui_GOTable newContentPane = new DREMGui_GOTable(ftheDREMGui, frame,
									theTimeiohmm.theDataSet,theTimeiohmm.bindingSignGene,
									theTimeiohmm.bindingGeneindex,theTimeiohmm.bindingSignTF,
									theTimeiohmm.bindingTFindex,theTimeiohmm.tfNames,bdisplay);
							//need to add geneIndex here

							newContentPane.setOpaque(true); //content panes must be opaque
							frame.setContentPane(newContentPane);
							//Display the window.
							frame.pack();
						}
					});
				}
			});


			staticButton = PPath.createRectangle((float) 0.0,(float) 0.0,(float) BUTTONRECWIDTH,(float) 18.0);
			PText staticText = new PText("Select by TFs");
			staticText.setFont(new Font("times",Font.PLAIN,12));
			staticText.translate(14,2);
			staticButton.setPaint(buttonColor);
			staticButton.addChild(staticText);   

			staticButton.addInputEventListener(new PBasicInputEventHandler() 
			{
				public void mousePressed(PInputEvent event) 
				{
					javax.swing.SwingUtilities.invokeLater(new Runnable() 
					{
						public void run() 
						{
							if (filterStaticFrame == null)
							{
								filterStaticFrame = new JFrame("Select Genes based on TF Input Constraints");
								filterStaticFrame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
								filterStaticFrame.setLocation(500,200);
								DREMGui_FilterStatic newContentPane = new DREMGui_FilterStatic(filterStaticFrame,ftheDREMGui,
										theTimeiohmm.theDataSet.tga,ftheDREMGui.treecopy);
								newContentPane.setOpaque(true); 
								//content panes must be opaque
								filterStaticFrame.setContentPane(newContentPane);

								//Display the window.
								filterStaticFrame.pack();
							}
							else
							{
								filterStaticFrame.setExtendedState(Frame.NORMAL);
							}
							filterStaticFrame.setVisible(true);
						}
					});
				}
			});

			gofilterButton = PPath.createRectangle((float) 0.0,(float) 0.0,(float) BUTTONRECWIDTH,(float) 18.0);

			gofilterButton.addInputEventListener(new PBasicInputEventHandler() 
			{
				public void mousePressed(PInputEvent event) 
				{
					javax.swing.SwingUtilities.invokeLater(new Runnable() 
					{
						public void run() 
						{
							if (filterGOFrame == null)
							{
								filterGOFrame = new JFrame("Select Genes based on GO Category");
								filterGOFrame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
								filterGOFrame.setLocation(500,250);

								theGOFilter = new DREMGui_GOFilter(filterGOFrame,ftheDREMGui,
										ftheDREMGui.theTimeiohmm.theDataSet.tga.theRecIDdrem,ftheDREMGui.treecopy);
								theGOFilter.setOpaque(true);
								filterGOFrame.setContentPane(theGOFilter);

								//Display the window.
								filterGOFrame.pack();
							}
							else
							{
								filterGOFrame.setExtendedState(Frame.NORMAL);
							}
							filterGOFrame.setVisible(true);
						}

					});
				}
			});


			PText gofilterText = new PText("Select by GO");
			gofilterText.setFont(new Font("times",Font.PLAIN,12));
			gofilterText.translate(22,2);
			gofilterButton.setPaint(buttonColor);
			gofilterButton.addChild(gofilterText);   

			scaleButton = PPath.createRectangle((float) 0.0,(float) 0.0,(float) BUTTONRECWIDTH,(float) 18.0);

			PText scaleText = new PText("Interface Options");
			scaleText.setFont(new Font("times",Font.PLAIN,12));
			scaleText.translate(10,2);
			scaleButton.setPaint(buttonColor);
			scaleButton.addChild(scaleText);   

			scaleButton.addInputEventListener(new PBasicInputEventHandler() 
			{
				public void mousePressed(PInputEvent event) 
				{

					javax.swing.SwingUtilities.invokeLater(new Runnable() 
					{
						public void run() 
						{
							if (yscaleFrame == null)
							{
								yscaleFrame = new JFrame("Interface Options");
								yscaleFrame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
								yscaleFrame.setLocation(650,150);
								theYScalegui = new DREMGui_InterfaceOptions(yscaleFrame,ftheDREMGui);
								theYScalegui.setOpaque(true); 
								//content panes must be opaque
								yscaleFrame.setContentPane(theYScalegui);
								//Display the window.
								yscaleFrame.pack();
							}
							else
							{
								yscaleFrame.setExtendedState(Frame.NORMAL);
							}
							yscaleFrame.setVisible(true);
						}
					});
				}
			});

			predictButton = PPath.createRectangle((float) 0.0,(float) 0.0,(float) BUTTONRECWIDTH,(float) 18.0);

			PText predictText = new PText("Predict");
			predictText.setFont(new Font("times",Font.PLAIN,12));
			predictText.translate(35,2);
			predictButton.setPaint(buttonColor);
			predictButton.addChild(predictText);   

			predictButton.addInputEventListener(new PBasicInputEventHandler() 
			{
				public void mousePressed(PInputEvent event) 
				{
					javax.swing.SwingUtilities.invokeLater(new Runnable() 
					{
						public void run() 
						{
							if (predictFrame == null)
							{
								predictFrame = new JFrame("Predict Time Series based on TF Input");
								predictFrame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
								predictFrame.setLocation(400,300);
								DREMGui_Predict newContentPane = new DREMGui_Predict(predictFrame,
										ftheDREMGui,ftheDREMGui.treecopy);
								newContentPane.setOpaque(true); 
								//content panes must be opaque
								predictFrame.setContentPane(newContentPane);

								//Display the window.
								predictFrame.pack();
							}
							else
							{
								predictFrame.setExtendedState(Frame.NORMAL);
							}
							predictFrame.setVisible(true);
						}
					});
				}
			});

			hideButton = PPath.createRectangle((float) 0.0,(float) 0.0,(float) BUTTONRECWIDTH,(float) 18.0);

			if (bglobalnode)
			{
				hideText = new PText("Hide Nodes");
			}
			else
			{
				hideText = new PText("Show Nodes");
			}
			hideText.setFont(new Font("times",Font.PLAIN,12));
			hideText.translate(23,2);
			hideButton.setPaint(buttonColor);
			hideButton.addChild(hideText);   

			hideButton.addInputEventListener(new PBasicInputEventHandler() 
			{
				public void mousePressed(PInputEvent event) 
				{

					int nsize = hideList.size();

					if (bglobalnode)
					{
						hideText.setText("Show Nodes");
						for (int nindex = 0; nindex < nsize; nindex++)
						{
							PNode theNode = ((PNode) hideList.get(nindex));
							theNode.setVisible(false);
							theNode.setPickable(false);
						}

						if (battachlabels)
						{
							hidelabels();
						}

						bglobalnode = false;
					}
					else
					{
						hideText.setText("Hide Nodes");
						for (int nindex = 0; nindex < nsize; nindex++)
						{
							PNode theNode = ((PNode) hideList.get(nindex));
							theNode.setVisible(true);
							theNode.setPickable(true);
						}

						if (battachlabels)
						{
							showlabels();
						}
						bglobalnode = true;
					}

				}
			});


			saveDREMButton = PPath.createRectangle((float) 0.0,(float) 0.0,(float) BUTTONRECWIDTH,(float) 18.0);


			PText saveDREMText;
			saveDREMText = new PText(" Save Image ");
			saveDREMText.setFont(new Font("times",Font.PLAIN,12));
			saveDREMText.translate(18,2);
			saveDREMButton.setPaint(buttonColor);
			saveDREMButton.addChild(saveDREMText);   

			saveDREMButton.addInputEventListener(new PBasicInputEventHandler() 
			{
				public void mousePressed(PInputEvent event) 
				{
					javax.swing.SwingUtilities.invokeLater(new Runnable() 
					{
						public void run() 
						{
							if (saveDREMFrame == null)
							{
								saveDREMFrame = new JFrame("Save as Image");
								saveDREMFrame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
								saveDREMFrame.setLocation(400,300);
								DREMGui_SaveDREM newContentPane = new DREMGui_SaveDREM(ftheDREMGui,saveDREMFrame);
								newContentPane.setOpaque(true); 
								//content panes must be opaque
								saveDREMFrame.setContentPane(newContentPane);
								//Display the window.
								saveDREMFrame.pack();
							}
							else
							{
								saveDREMFrame.setExtendedState(Frame.NORMAL);
							}
							saveDREMFrame.setVisible(true);
						}
					});
				}
			});

			saveModelButton = PPath.createRectangle((float) 0.0,(float) 0.0,(float) BUTTONRECWIDTH,(float) 18.0);

			PText saveModelText;
			saveModelText = new PText(" Save Model ");
			saveModelText.setFont(new Font("times",Font.PLAIN,12));
			saveModelText.translate(18,2);
			saveModelButton.setPaint(buttonColor);
			saveModelButton.addChild(saveModelText);  

			saveModelButton.addInputEventListener(new PBasicInputEventHandler() 
			{
				public void mousePressed(PInputEvent event) 
				{
					javax.swing.SwingUtilities.invokeLater(new Runnable() 
					{
						public void run() 
						{
							if (saveModelFrame == null)
							{
								saveModelFrame = new JFrame("Save Model to File");
								saveModelFrame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
								saveModelFrame.setLocation(400,300);
								DREMGui_SaveModel newContentPane = 
									new DREMGui_SaveModel(ftheDREMGui.theTimeiohmm,treecopy,
											saveModelFrame,ftheDREMGui);
								newContentPane.setOpaque(true); 
								//content panes must be opaque
								saveModelFrame.setContentPane(newContentPane);
								//Display the window.
								saveModelFrame.pack();
							}
							else
							{
								saveModelFrame.setExtendedState(Frame.NORMAL);
							}
							saveModelFrame.setVisible(true);
						}
					});
				}
			});

			timeButton = PPath.createRectangle((float) 0.0,(float) 0.0,(float) BUTTONRECWIDTH,(float) 18.0);

			PText timeText;
			if (bglobalVisible)
			{
				timeText = new PText("Hide Time Series");
			}
			else
			{
				timeText = new PText("Show Time Series");
			}
			timeText.setFont(new Font("times",Font.PLAIN,12));
			timeText.translate(7,2);
			timeButton.setPaint(buttonColor);
			timeButton.addChild(timeText);   

			final PText ftimeText = timeText;
			final PNode[] fplArray = plArray;
			timeButton.addInputEventListener(new PBasicInputEventHandler() 
			{
				public void mousePressed(PInputEvent event) 
				{
					if (bglobalVisible)
					{
						ftimeText.setText("Show Time Series");
						for (int nindex = 0; nindex < fplArray.length; nindex++)
						{
							fplArray[nindex].setVisible(false);
							fplArray[nindex].setPickable(false);
						}
						bglobalVisible = false;
					}
					else
					{
						ftimeText.setText("Hide Time Series");
						for (int nindex = 0; nindex < fplArray.length; nindex++)
						{
							boolean bvisible = bPathVisible[nindex]&&bTFVisible[nindex]&&
							bGOVisible[nindex]&&bSetVisible[nindex];
							fplArray[nindex].setVisible(bvisible);
							fplArray[nindex].setPickable(bvisible);
						}
						bglobalVisible = true;
					}
				}
			});



			siginButton = PPath.createRectangle((float) 0.0,(float) 0.0,(float) BUTTONRECWIDTH,(float) 18.0);

			PText siginText = new PText("Key TFs Labels");
			siginText.setFont(new Font("times",Font.PLAIN,12));
			siginText.translate(10,2);
			siginButton.setPaint(buttonColor);
			siginButton.addChild(siginText);  

			nheight = getHeight();
			int ninset = 10;
			int nspacing = ninset + BUTTONRECWIDTH;

			nwidth = getWidth();
			INITLEFT = nwidth /2.0- 3*nspacing;
			helpButton = new PImage(Util.getImageURL("Help24.gif"));
			final DREMGui thisFrame = this;

			helpButton.addInputEventListener(new PBasicInputEventHandler() 
			{
				public void mousePressed(PInputEvent event) 
				{
					javax.swing.SwingUtilities.invokeLater(new Runnable() 
					{
						public void run() 
						{
							JDialog helpDialog = new JDialog(thisFrame, "Help", false);
							Container theHelpDialogPane = helpDialog.getContentPane();

							helpDialog.setBackground(Color.white);
							theHelpDialogPane.setBackground(Color.white);
							String szMessage = 
								"This is the main interface window of DREM. "+
								"The window shows a plot of the expression profiles of genes that were "+
								"not filtered along with the map learned by DREM, and associated TF labels.  "+
								"Green nodes are split nodes, and nodes size is proportional to the standard deviation of their assocoated "+
								"gaussian emission distribution.  "+
								"Left clicking on an edge of the map shows only genes assigned to that path through the model.  Right clicking "+
								"on an edge on TF label box brings up info about the TFs regulating genes assigned to the path.\n\n"+
								"Along the bottom there are 12 buttons.  These buttons functions as follows: \n"+
								"Predict - displays for a given TF input vector the probability of transitioning to each state \n"+
								"Interface Options - displays menu to adjust interface options \n"+
								"GO Table - displays a GO enrichment analysis for the currently selected genes\n"+
								"Gene Table - displays a gene table for the currently selected genes\n"+
								"Key TFs Labels - adjust the criteria and threshold for display TF labels on the map\n"+
								"Select by TFs - select genes by which TFs they regulated by \n"+
								"Select by GO - selects genes based on the GO category to which they are annotated\n"+
								"Select by Gene Set - selects genes based on a defined set\n"+
								"Hide Nodes/Show Nodes - hide/shows the nodes \n"+
								"Hide Time Series/Show Time Series - hide/shows the time series expression patterns on the interface\n"+
								"Save Model - exports the model to a file can then be loaded later under the 'Saved Model File' option\n"+
								"Save Image - saves the interface window to a graphics file\n";

							JTextArea textArea = new JTextArea(szMessage,10,60);
							textArea.setLineWrap(true);
							textArea.setWrapStyleWord(true);
							textArea.setBackground(Color.white);
							textArea.setEditable(false);

							JScrollPane jsp2 = new JScrollPane(textArea);

							theHelpDialogPane.add(jsp2);
							theHelpDialogPane.setSize(820,600);
							theHelpDialogPane.validate();

							helpDialog.setLocation(thisFrame.getX()+50,thisFrame.getY()+25);

							helpDialog.setSize(820,600);
							helpDialog.setDefaultCloseOperation(DISPOSE_ON_CLOSE);
							helpDialog.setVisible(true);
						}
					});
				}
			});   

			theCamera.addChild(saveDREMButton);
			theCamera.addChild(saveModelButton);
			theCamera.addChild(helpButton);
			theCamera.addChild(staticButton);
			theCamera.addChild(setButton);
			theCamera.addChild(gofilterButton);
			theCamera.addChild(scaleButton);
			theCamera.addChild(predictButton);
			theCamera.addChild(hideButton);
			theCamera.addChild(timeButton); 
			theCamera.addChild(siginButton);
			theCamera.addChild(genetableButton);
			theCamera.addChild(gotableButton);

			scaleButton.translate(INITLEFT,nheight-55);
			predictButton.translate(INITLEFT,nheight-77);

			genetableButton.translate(INITLEFT+nspacing,nheight-55);
			gotableButton.translate(INITLEFT+nspacing,nheight-77);

			siginButton.translate(INITLEFT+2*nspacing,nheight-77);
			staticButton.translate(INITLEFT+2*nspacing,nheight-55);

			setButton.translate(INITLEFT+3*nspacing,nheight-55);
			gofilterButton.translate(INITLEFT+3*nspacing,nheight-77);

			helpButton.translate(INITLEFT+6*nspacing,nheight-59);
			saveDREMButton.translate(INITLEFT+5*nspacing,nheight-55);
			saveModelButton.translate(INITLEFT+5*nspacing,nheight-77);

			hideButton.translate(INITLEFT+4*nspacing,nheight-77);
			timeButton.translate(INITLEFT+4.0*nspacing,nheight-55);

			siginButton.addInputEventListener(new PBasicInputEventHandler() 
			{
				public void mousePressed(PInputEvent event) 
				{
					javax.swing.SwingUtilities.invokeLater(new Runnable() 
					{
						public void run() 
						{
							if (keyinputFrame == null)
							{
								keyinputFrame = new JFrame("Key Transcription Factors (TFs) Labels");
								keyinputFrame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
								keyinputFrame.setLocation(600,350);
								DREMGui_KeyInputs newContentPane = new DREMGui_KeyInputs(keyinputFrame,ftheDREMGui);
								newContentPane.setOpaque(true); 
								//content panes must be opaque
								keyinputFrame.setContentPane(newContentPane);

								//Display the window.
								keyinputFrame.pack();
							}
							else
							{
								keyinputFrame.setExtendedState(Frame.NORMAL);
							}
							keyinputFrame.setVisible(true);
						}
					});
				}
			});

			int nsize = hidesigList.size();
			for (int nindex = 0; nindex < nsize; nindex++)
			{
				SigInfoRec theSigInfoRec =  (SigInfoRec) hidesigList.get(nindex);

				theLayer.addChild(theSigInfoRec.border);

				double dprewidth =theSigInfoRec.theSigText.getWidth();
				setSigText(theSigInfoRec.theSigTF, theSigInfoRec.theSigText, theSigInfoRec.ntype, theSigInfoRec.node);

				theSigInfoRec.border.translate(dprewidth-theSigInfoRec.theSigText.getWidth(),0);
			}
			renderFilterText();

		}  //synchronized
	}


	/**
	 * Hides all annotation labels
	 */
	public void hidelabels()
	{
		int ngolabelsize =  hidegolabelsList.size();
		for (int nindex = 0; nindex < ngolabelsize; nindex++)
		{
			PText theText = (PText) hidegolabelsList.get(nindex);
			theText.setVisible(false);
			theText.setPickable(false);
		}

		int ngenesetlabelsize =  hidegenesetlabelsList.size();
		for (int nindex = 0; nindex < ngenesetlabelsize; nindex++)
		{
			PText theText = (PText) hidegenesetlabelsList.get(nindex);
			theText.setVisible(false);
			theText.setPickable(false);
		}

		int ntflabelsize =  hidetfsetlabelsList.size();
		for (int nindex = 0; nindex < ntflabelsize; nindex++)
		{
			PText theText = (PText) hidetfsetlabelsList.get(nindex);
			theText.setVisible(false);
			theText.setPickable(false);
		}

		int npredictsize =  hidepredictlabelsList.size();
		for (int nindex = 0; nindex < npredictsize; nindex++)
		{
			PText theText = (PText) hidepredictlabelsList.get(nindex);
			theText.setVisible(false);
			theText.setPickable(false);
		}

		int nsizesig = hidesigList.size();
		for (int nindex = 0; nindex < nsizesig; nindex++)
		{
			DREMGui.SigInfoRec theSigInfoRec = 
				(DREMGui.SigInfoRec) hidesigList.get(nindex);	      
			PPath rect = (PPath) theSigInfoRec.theSigText.getParent();
			rect.setVisible(false);
			rect.setPickable(false);
			theSigInfoRec.theSigText.setVisible(false);
			theSigInfoRec.theSigText.setPickable(false);
		}

		filterText.setVisible(false);
		filterText.setPickable(false);
	}


	/**
	 * Shows all annotation labels
	 */
	public void showlabels()
	{
		if ((bshowgolabels)&&(bapplygolabels))
		{
			int ngolabelsize =  hidegolabelsList.size();
			for (int nindex = 0; nindex < ngolabelsize; nindex++)
			{
				PText theText = (PText) hidegolabelsList.get(nindex);
				theText.setVisible(true);
			}
		}


		if ((bshowgenesetlabels)&&(bapplygenesetlabels))
		{
			int ngenesetlabelsize =  hidegenesetlabelsList.size();
			for (int nindex = 0; nindex < ngenesetlabelsize; nindex++)
			{
				PText theText = (PText) hidegenesetlabelsList.get(nindex);
				theText.setVisible(true);
			}
		}

		if ((bshowtfsetlabels)&&(bapplytfsetlabels))
		{
			int ntflabelsize =  hidetfsetlabelsList.size();
			for (int nindex = 0; nindex < ntflabelsize; nindex++)
			{
				PText theText = (PText) hidetfsetlabelsList.get(nindex);
				theText.setVisible(true);
			}
		}

		if (bshowpredict)
		{
			int npredictsize =  hidepredictlabelsList.size();
			for (int nindex = 0; nindex < npredictsize; nindex++)
			{
				PText theText = (PText) hidepredictlabelsList.get(nindex);
				theText.setVisible(true);
			}
		}

		if (bshowkeyinputs)
		{
			int nsizesig = hidesigList.size();
			for (int nindex = 0; nindex < nsizesig; nindex++)
			{
				DREMGui.SigInfoRec theSigInfoRec = 
					(DREMGui.SigInfoRec) hidesigList.get(nindex);	      
				PPath rect = (PPath) theSigInfoRec.theSigText.getParent();
				rect.setVisible(true);
				rect.setPickable(true);
				theSigInfoRec.theSigText.setVisible(true);
				theSigInfoRec.theSigText.setPickable(true);
			}
		}

		filterText.setVisible(true);
		filterText.setPickable(false);
	}

	/**
	 * Information about circle on the interface
	 */
	static class CircleID
	{
		int ndepth;
		int nminparentlevel;
		int nscore;
		int nid;
		int nprevminparentlevel;

		CircleID(int ndepth, int nminparentlevel, int nprevminparentlevel, int nscore,int nid)
		{
			this.nminparentlevel = nminparentlevel;
			this.nprevminparentlevel = nprevminparentlevel;
			this.ndepth = ndepth;
			this.nscore = nscore;
			this.nid = nid;
		}

	}

	/**
	 * Record for the node currently selected
	 */
	static class SelectedNodeRec
	{
		PNode selectedNode;
		boolean bcircle;
		CircleID theCircleID=null;
	}

	/**
	 * Sets the text for the significant transcription factors.
	 */
	public void setSigText(TreeSet tsSigTF, PText theSigText, int ntype)
	{
		StringBuffer szNamesbuf = new StringBuffer();
		// If the current type of label to be displayed (nKeyInputType) is not equal to
		// the the type of label this text corresponds to, we need to erase theSigText
		if (ntype == nKeyInputType)
		{
			Iterator itr = tsSigTF.iterator();
			HashSet htAdded = new HashSet();  

			if (itr.hasNext())
			{
				DREM_Timeiohmm.SigTFRecv2 theTFRec =  (DREM_Timeiohmm.SigTFRecv2) itr.next();

				if ((theTFRec.dpval <= dkeyinputpvalue)&&((nKeyInputType != 1)
						||(theTFRec.dpercent >= dsplitpercent)))
				{
					szNamesbuf.append(" "+theTFRec.szname+" \n");
					htAdded.add(theTFRec.szname);
				}

				while (itr.hasNext())
				{
					theTFRec =  (DREM_Timeiohmm.SigTFRecv2) itr.next();
					if ((theTFRec.dpval <= dkeyinputpvalue)&&
							((nKeyInputType != 1)||(theTFRec.dpercent >= dsplitpercent))&&
							(!htAdded.contains(theTFRec.szname)))
					{
						szNamesbuf.append(" "+theTFRec.szname+" \n");
						htAdded.add(theTFRec.szname);
					}
				}
			}
		}
		String szNames = szNamesbuf.toString();
		theSigText.setText(szNames);
		PPath rect = (PPath) theSigText.getParent();

		rect.reset();
		if (!szNames.equals(""))
		{ 
			rect.setPathToRectangle((float) 0.0,(float) 0.0,
					(float) (theSigText.getWidth()),(float) theSigText.getHeight());
		}   

	}
	
	// TODO should be able to merge the two versions of setSigText
	/**
	 * Sets the text for the significant transcription factors.  This version
	 * is specific to when activity scores are used to determine which
	 * TFs are significant and will only label a node with a TF name
	 * if it is the first time on the path that the TF appears, but is
	 * called by all types.  For the other types, the simpler version of
	 * setSigtext will be called here automatically.
	 * 
	 * @param ptr The tree node for which TF labels are being drawn.
	 */
	public void setSigText(TreeSet tsSigTF, PText theSigText, int ntype,
			DREM_Timeiohmm.Treenode ptr)
	{
		// Use the simple version if we are not only printing
		// significant TFs that first appear at this node on the path
		if(ntype != 778)
		{
			setSigText(tsSigTF, theSigText, ntype);
		}
		// Only applicable to activity scores when showing TFs only
		// the first time they appear per path
		else
		{
			StringBuffer szNamesbuf = new StringBuffer();
			// If the current type of label to be displayed (nKeyInputType) is not equal to
			// the the type of label this text corresponds to, we need to erase theSigText
			if (ntype == nKeyInputType)
			{
				// First traverse the ancestors of this node to determine
				// which TFs have already appeared at previous nodes
				// on the path
				HashMap ancestorScores = ptr.getAncestorActivityScores();

				Iterator itr = tsSigTF.iterator();
				HashSet htAdded = new HashSet();  

				DREM_Timeiohmm.SigTFRecv2 theTFRec;
				while (itr.hasNext())
				{
					theTFRec = (DREM_Timeiohmm.SigTFRecv2) itr.next();

					// If the activity value is above the threshold...
					// (activity values have been transformed to resemble p-values here)
					if (theTFRec.dpval <= dkeyinputpvalue)
					{
						// ... then check if it has already appeared in a path before using it
						String tfName = theTFRec.szname;
						// Check if the path out of the split has been added as " [#]"
						// to the TF name
						if(tfName.endsWith("]"))
						{
							tfName = tfName.substring(0,tfName.lastIndexOf(" ["));
						}

						// Only add the TF to this node's set of labels if
						// it wasn't already added to a set of labels at
						// an ancestor node on the path
						if(!ancestorScores.containsKey(tfName) ||
								((Double) ancestorScores.get(tfName)).doubleValue() > dkeyinputpvalue)
						{
							// Add the name with the path index
							szNamesbuf.append(" "+theTFRec.szname+" \n");
							htAdded.add(theTFRec.szname);
						}
					}
				}

			}
			String szNames = szNamesbuf.toString();
			theSigText.setText(szNames);
			PPath rect = (PPath) theSigText.getParent();

			rect.reset();
			if (!szNames.equals(""))
			{ 
				rect.setPathToRectangle((float) 0.0,(float) 0.0,
						(float) (theSigText.getWidth()),(float) theSigText.getHeight());
			}   
		}
	}

	/**
	 * Record for information about significant transcription factors
	 */
	public static class SigInfoRec
	{
		PText theSigText;
		TreeSet theSigTF;
		int ntype;
		PNode border;
		int npathscore;
		int ndepth;
		DREM_Timeiohmm.Treenode node;

		SigInfoRec(PNode border, PText theSigText,TreeSet tsSigTF,int ntype,
				int npathscore,int ndepth,DREM_Timeiohmm.Treenode node)
		{
			this.theSigText = theSigText;
			this.theSigTF = tsSigTF;
			this.ntype = ntype;
			this.border = border;
			this.npathscore = npathscore;
			this.ndepth = ndepth;
			this.node = node;
		}

	}

	/**
	 * Draws nodes on the interface screen
	 */
	public PBasicInputEventHandler drawNodes(DREM_Timeiohmm.Treenode ptr,final int ndepth,double parentx,
			double parenty,final int ncurrscore, int nchild,
			Color currColor,Color prevColor,int nprevcolorID, 
			int ncurrparentcolorindex, int nprevminparentlevel)
	{

		PBasicInputEventHandler pbiehLine = null;

		if (ptr != null)
		{
			double dnodex;
			double dnodey;
			double ddiameter;

			dnodex =  dwidthunitsCum[ndepth]+LEFTBUFFER;
			dnodey = REALHEIGHT - dheightunits*(ptr.dmean-dmin);
			ddiameter = Math.sqrt(dnodek)*dminheightunits*Math.sqrt(ptr.dsigma)/2;
			PPath line = null;
			int nlocalcolorID=numcolor - 1;;

			boolean bnewcolor = true;
			if ((ptr.parent != null)&&(ptr.parent.parentptrA != null)
					&&(ptr.parent.parentptrA.length != 1))
			{
				currColor = prevColor;
				nlocalcolorID = nprevcolorID;
				bnewcolor = false;
			}

			if (currColor == null)
			{
				currColor = (Color) htColorIDtoColor.get(new Integer(numcolor));
				if (currColor == null)
				{
					if (numcolor < theTimeiohmm.savedColors.size())
					{
						currColor = (Color) theTimeiohmm.savedColors.get(numcolor);
					}
					else
					{
						if (nchild == 0)
						{
							if (nchild < edgeColorsTriples[0].length)
							{
								currColor = 
									edgeColorsTriples[nparentcolorindex%edgeColorsTriples.length][nchild];
							}
							else
							{
								currColor = 
									edgeColorsTriples[(nparentcolorindex+nchild-edgeColorsTriples[0].length+1)
									                  %edgeColorsTriples[0].length][edgeColorsTriples[0].length-1];
							}

							if (DREM_Timeiohmm.BDEBUG)
							{
								System.out.println(nparentcolorindex+"\t"+nchild+"\t"+currColor);
							}

							if (ptr.parent != null)
							{
								ptr.parent.nparentcolorindex = nparentcolorindex;
							}
							nparentcolorindex++;
							ncurrparentcolorindex = nparentcolorindex;
						}
						else
						{
							if (nchild < edgeColorsTriples[0].length)
							{
								currColor = 
									edgeColorsTriples[ptr.parent.nparentcolorindex%edgeColorsTriples.length][nchild];
							}
							else
							{
								currColor = 
									edgeColorsTriples[(ptr.parent.nparentcolorindex+nchild-edgeColorsTriples[0].length+1)
									                  %edgeColorsTriples[0].length][edgeColorsTriples[0].length-1];
							}

							if (DREM_Timeiohmm.BDEBUG)
							{
								System.out.println(ptr.parent.nparentcolorindex+"\t"+nchild+"\t"+currColor);
							}
						}
					}
				}
				nlocalcolorID = numcolor;
				numcolor++;
			}	 

			if (bnewcolor)
			{
				htColorIDtoColor.put(new Integer(nlocalcolorID),currColor);
			}

			htLineIDtoColorID.put(ndepth+";"+ncurrscore, new Integer(nlocalcolorID));

			PBasicInputEventHandler borderTextHide = new PBasicInputEventHandler()
			{
				public void mousePressed(PInputEvent event) 
				{
					if (event.getButton() == MouseEvent.BUTTON1)
					{  
						PNode thenode = event.getPickedNode();
						if (thenode.getVisible())
						{
							thenode.setVisible(false);
							PNode par = thenode.getParent();
							if(par instanceof PPath)
							{
								// Set the color to be semi-transparent
								((PPath) par).setStrokePaint(new Color(0,0,0,75));
							}
						}
						else
						{
							thenode.setVisible(true);
							PNode par = thenode.getParent();
							if(par instanceof PPath)
							{
								// Set the color to be opaque black
								((PPath) par).setStrokePaint(new Color(0,0,0));
							}
						}
					}
				}
			};

			if (ndepth == 0)
			{
				final int fnchild = nchild;
				final DREM_Timeiohmm.Treenode fptr = ptr;
				final DREMGui theDREMGui = this;
				pbiehLine = new PBasicInputEventHandler() 
				{
					public void mousePressed(PInputEvent event) 
					{
						if (event.getButton() == MouseEvent.BUTTON3)
						{  
							javax.swing.SwingUtilities.invokeLater(new Runnable() 
							{
								public void run() 
								{
									JFrame frame = new JFrame("Path Table "+nKeyInputType);
									frame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
									frame.setLocation(200,200);

									DREMGui_EdgeTable newContentPane = new DREMGui_EdgeTable(theDREMGui,
											frame, theTimeiohmm,
											fptr,ncurrscore,ndepth,
											fnchild,nKeyInputType,true);
									newContentPane.setOpaque(true); //content panes must be opaque
									frame.setContentPane(newContentPane);

									//Display the window.
									frame.pack();
									frame.setVisible(true);
								}  
							});
						}
					}
				};
			}
			else if (ndepth > 0)
			{
				line = PPath.createLine((float) parentx,(float) parenty,(float) dnodex,(float) dnodey);
				PPath line2 = PPath.createLine((float) parentx,(float) parenty,(float) dnodex,(float) dnodey);

				line2.setStroke(new BasicStroke(8));
				line2.setStrokePaint(Color.black);
				canvas.getLayer().addChild(line2);

				line.setStroke(new BasicStroke(5));
				line.setStrokePaint(currColor);
				ArrayList linesList = (ArrayList) htColorIDtoLinesList.get(new Integer(nlocalcolorID));
				if (linesList == null)
				{
					linesList = new ArrayList();
					linesList.add(line);
					htColorIDtoLinesList.put(new Integer(nlocalcolorID), linesList);
				}
				else
				{
					linesList.add(line);
				}

				if (!bglobalnode)
				{
					line.setVisible(false);
					line.setPickable(false);
					line2.setVisible(false);
					line2.setPickable(false);
				}
				canvas.getLayer().addChild(line);

				hideList.add(line);
				hideList.add(line2);

				final int fnchild = nchild;
				final DREM_Timeiohmm.Treenode fptr = ptr;
				final DREMGui theDREMGui = this;
				pbiehLine = new PBasicInputEventHandler() 
				{
					public void mousePressed(PInputEvent event) 
					{
						if (event.getButton() == MouseEvent.BUTTON3)
						{  
							javax.swing.SwingUtilities.invokeLater(new Runnable() 
							{
								public void run() 
								{
									JFrame frame = new JFrame("Path Table "+nKeyInputType);
									frame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
									frame.setLocation(200,200);

									DREMGui_EdgeTable newContentPane = new DREMGui_EdgeTable(theDREMGui,
											frame, theTimeiohmm,
											fptr.parent,ncurrscore,ndepth,
											fnchild,nKeyInputType,false);
									newContentPane.setOpaque(true); //content panes must be opaque
									frame.setContentPane(newContentPane);

									//Display the window.
									frame.pack();
									frame.setVisible(true);
								}  
							});
						}
					}
				};
				line.addInputEventListener(pbiehLine);
			}

			PNode circle = PPath.createEllipse((float)(dnodex-ddiameter/2.0),(float)(dnodey-ddiameter/2.0),
					(float) ddiameter,(float) ddiameter);
			String szTip = "Mean "+ nf3.format(ptr.dmean) +"; Std. Dev. "+nf3.format(ptr.dsigma);
			ToolTipTextUpdater tipEvent = new ToolTipTextUpdater(canvas, circle,szTip);
			circle.addInputEventListener(tipEvent);

			if (ptr.numchildren <= 1)
			{
				ArrayList circleList = (ArrayList) htColorIDtoCircleList.get(new Integer(nlocalcolorID));
				if (circleList == null)
				{
					circleList = new ArrayList();
					circleList.add(circle);
					htColorIDtoCircleList.put(new Integer(nlocalcolorID), circleList);
				}
				else
				{
					circleList.add(circle);
				}

				if ((ptr.parentptrA == null)||(ptr.parentptrA.length == 1))
				{
					circle.setPaint(currColor);
				}
				else
				{
					circle.setPaint(prevColor);
				}
			}
			else
			{
				circle.setPaint(SPLITCOLOR);
			}

			PBasicInputEventHandler[] pbiehLineUse = new PBasicInputEventHandler[ptr.numchildren];

			// Recursivly draw the children
			for (int nnextchild = 0; nnextchild < ptr.numchildren; nnextchild++)
			{
				int npathscore = ncurrscore+
				nnextchild*(int) Math.pow(theTimeiohmm.nmaxchild,data[0].length-ndepth-1);
				if (ptr.numchildren == 1)
				{
					pbiehLineUse[nnextchild] =  
						drawNodes(ptr.nextptr[nnextchild],ndepth+1,dnodex,dnodey,
								npathscore,nnextchild,currColor,prevColor,nprevcolorID,
								ncurrparentcolorindex,ptr.nminparentlevel );
				}
				else
				{
					pbiehLineUse[nnextchild] = 
						drawNodes(ptr.nextptr[nnextchild],ndepth+1,dnodex,dnodey,
								npathscore,nnextchild,null,currColor,nlocalcolorID,ncurrparentcolorindex,ptr.nminparentlevel );
				}
			}

			final PNode fcircle = circle;
			final DREM_Timeiohmm.Treenode fptr = ptr;

			if ((ptr.parentptrA == null)||(ptr.parentptrA.length == 1))
			{
				htNodes.put(circle,new CircleID(ndepth, ptr.nminparentlevel,
						ptr.nminparentlevel, ncurrscore,ncircleID));
			}
			else
			{
				if (DREM_Timeiohmm.BDEBUG)
				{
					System.out.println(ptr.nminparentlevel+" "+nprevminparentlevel+"!!!!");
				}

				htNodes.put(circle,new CircleID(ndepth, ptr.nminparentlevel,
						nprevminparentlevel, ncurrscore,ncircleID));
			}

			NumberFormat nf1 = NumberFormat.getInstance(Locale.ENGLISH);

			if (ptr.dpredictweight < 1)
			{
				nf1.setMaximumIntegerDigits(0);
				nf1.setMinimumFractionDigits(2);
				nf1.setMaximumFractionDigits(2);
			}
			else
			{
				nf1.setMinimumFractionDigits(1);
				nf1.setMaximumFractionDigits(1);
			}

			ptr.thepredictText = new PText(nf1.format(ptr.dpredictweight));
			ptr.thepredictText.setVisible(bshowpredict&&(bglobalnode||!battachlabels));
			ptr.thepredictText.translate(dnodex-ddiameter/2.0,dnodey-ddiameter/3.0);
			ptr.thepredictText.setFont(new Font("times",Font.BOLD,10));
			ptr.thepredictText.setTextPaint(predictLabelColor);
			ptr.thepredictText.setPickable(false);

			ptr.goText = new PText(ptr.szgolabel);
			ptr.goText.setVisible(bshowgolabels&&bapplygolabels&&(bglobalnode||!battachlabels));
			ptr.goText.setFont(new Font("times",Font.BOLD,14));
			ptr.goText.setTextPaint(Color.black);
			ptr.goText.setPickable(false);
			ptr.goText.translate(dnodex+ddiameter/2.0,dnodey-ddiameter/2.0);
			hidegolabelsList.add(ptr.goText);

			ptr.genesetText = new PText(ptr.szgenesetlabel);
			ptr.genesetText.setVisible(bshowgenesetlabels&&bapplygenesetlabels
					&&(bglobalnode||!battachlabels));
			ptr.genesetText.setFont(new Font("times",Font.BOLD,14));
			ptr.genesetText.setTextPaint(genesetLabelColor);
			ptr.genesetText.setPickable(false);
			ptr.genesetText.translate(dnodex+ddiameter/2.0,dnodey-ddiameter/2.0);
			hidegenesetlabelsList.add(ptr.genesetText);

			ptr.tfsetText = new PText(ptr.sztfsetlabel);
			ptr.tfsetText.setVisible(bshowtfsetlabels&&bapplytfsetlabels
					&&(bglobalnode||!battachlabels));
			ptr.tfsetText.setFont(new Font("times",Font.BOLD,14));
			ptr.tfsetText.setTextPaint(tfLabelColor);
			ptr.tfsetText.setPickable(false);
			ptr.tfsetText.translate(dnodex+ddiameter/2.0,dnodey-ddiameter/2.0);
			hidetfsetlabelsList.add(ptr.tfsetText);
			hidepredictlabelsList.add(ptr.thepredictText);

			for (int nchildindex = 0; nchildindex < ptr.numchildren; nchildindex++)
			{
				PText thesigEdgeFullText = new PText();
				PPath border = new PPath();

				if ((!bglobalnode)&&(battachlabels)||(!bshowkeyinputs))
				{
					border.setVisible(false);
					border.setPickable(false);
				}
				border.addChild(thesigEdgeFullText);

				double dnexty = REALHEIGHT - dheightunits*(ptr.nextptr[nchildindex].dmean-dmin);
				double dradius = Math.sqrt(dnodek)*dminheightunits*Math.sqrt(ptr.nextptr[nchildindex].dsigma)/4;

				if (ptr.tsSigTFFull != null)
				{
					setSigText(ptr.tsSigTFFull[nchildindex], thesigEdgeFullText,3,ptr);
					border.translate(dnodex+dwidthunitsInterval[ndepth]-dradius,dnexty-dradius); 
					thesigEdgeFullText.setFont(new Font("Arial",Font.BOLD,14));
					thesigEdgeFullText.setTextPaint(keyInputLabelColor);
					int nextscore = ncurrscore+
					nchildindex*(int) Math.pow(theTimeiohmm.nmaxchild,data[0].length-ndepth-1);
					hidesigList.add(new SigInfoRec(border,thesigEdgeFullText,
							ptr.tsSigTFFull[nchildindex],2,nextscore,ndepth,ptr));
					thesigEdgeFullText.addInputEventListener(borderTextHide);
					thesigEdgeFullText.addInputEventListener(pbiehLineUse[nchildindex]);

					if (htTextVisible != null)
					{
						boolean bshow =  ((Boolean)htTextVisible.get(ndepth+";"+nextscore+";"+2)).booleanValue();
						thesigEdgeFullText.setVisible(bshow);
						thesigEdgeFullText.setPickable(bshowkeyinputs);
					}
				}
			}

			PNode borderCircle;
			if (ptr.numchildren >= 2)
			{
				PText thesigText = new PText();
				borderCircle = new PPath();
				if ((!bglobalnode)&&(battachlabels)||(!bshowkeyinputs))
				{
					borderCircle.setVisible(false);
					borderCircle.setPickable(false);
				}

				borderCircle.addChild(thesigText);
				borderCircle.translate(dnodex+ddiameter/2,dnodey-ddiameter/2);
				setSigText(ptr.tsSigTF, thesigText,0,ptr);

				if (htTextVisible != null)
				{
					boolean bshow = ((Boolean)htTextVisible.get(ndepth+";"+ncurrscore+";"+0)).booleanValue();
					thesigText.setVisible(bshow);
					thesigText.setPickable(bshowkeyinputs);
				}

				thesigText.setFont(new Font("Arial",Font.BOLD,14));
				thesigText.setTextPaint(keyInputLabelColor);
				hidesigList.add(new SigInfoRec(borderCircle,thesigText,ptr.tsSigTF,0,ncurrscore,ndepth,ptr));

				// Repeat for the TF activity scores
				PText thesigTextAct = new PText();
				PNode borderCircleAct = new PPath();
				if ((!bglobalnode)&&(battachlabels)||(!bshowkeyinputs))
				{
					borderCircleAct.setVisible(false);
					borderCircleAct.setPickable(false);
				}

				borderCircleAct.addChild(thesigTextAct);
				borderCircleAct.translate(dnodex+ddiameter/2,dnodey-ddiameter/2);
				setSigText(ptr.tsSigTFActivity, thesigTextAct, 777,ptr);

				if (htTextVisible != null)
				{
					boolean bshow = ((Boolean)htTextVisible.get(ndepth+";"+ncurrscore+";"+777)).booleanValue();
					thesigTextAct.setVisible(bshow);
					thesigTextAct.setPickable(bshowkeyinputs);
				}

				thesigTextAct.setFont(new Font("Arial",Font.BOLD,14));
				thesigTextAct.setTextPaint(keyInputLabelColor);
				hidesigList.add(new SigInfoRec(borderCircleAct,thesigTextAct,ptr.tsSigTFActivity,777,ncurrscore,ndepth,ptr));
				
				// TODO can this be shared with the other text and circle for the activity scores?
				// Repeat for the TF activity scores that are shown
				// only the first time the TF is active on a path
				PText thesigTextFirstAct = new PText();
				PNode borderCircleFirstAct = new PPath();
				if ((!bglobalnode)&&(battachlabels)||(!bshowkeyinputs))
				{
					borderCircleFirstAct.setVisible(false);
					borderCircleFirstAct.setPickable(false);
				}

				borderCircleFirstAct.addChild(thesigTextFirstAct);
				borderCircleFirstAct.translate(dnodex+ddiameter/2,dnodey-ddiameter/2);
				setSigText(ptr.tsSigTFActivity, thesigTextFirstAct, 778, ptr);

				if (htTextVisible != null)
				{
					boolean bshow = ((Boolean)htTextVisible.get(ndepth+";"+ncurrscore+";"+778)).booleanValue();
					thesigTextFirstAct.setVisible(bshow);
					thesigTextFirstAct.setPickable(bshowkeyinputs);
				}

				thesigTextFirstAct.setFont(new Font("Arial",Font.BOLD,14));
				thesigTextFirstAct.setTextPaint(keyInputLabelColor);
				hidesigList.add(new SigInfoRec(borderCircleFirstAct,thesigTextFirstAct,ptr.tsSigTFActivity,778,ncurrscore,ndepth,ptr));
				
				
				for (int nchildindex = 0; nchildindex < ptr.numchildren; nchildindex++)
				{
					PText thesigEdgeSplitText = new PText();
					PNode border = new PPath();

					if ((!bglobalnode)&&(battachlabels)||(!bshowkeyinputs))
					{
						border.setVisible(false);
						border.setPickable(false);
					}

					border.addChild(thesigEdgeSplitText);
					setSigText(ptr.tsSigTFEdgeSplit[nchildindex], thesigEdgeSplitText, 3, ptr);

					double dnexty = REALHEIGHT - dheightunits*(ptr.nextptr[nchildindex].dmean-dmin);			
					double dradius = Math.sqrt(dnodek)*(dminheightunits)*Math.sqrt(ptr.nextptr[nchildindex].dsigma)/4;

					border.translate(dnodex+dwidthunitsInterval[ndepth]-dradius,dnexty-dradius);			  

					thesigEdgeSplitText.setFont(new Font("Arial",Font.BOLD,14));

					thesigEdgeSplitText.setTextPaint(keyInputLabelColor);
					thesigEdgeSplitText.addInputEventListener(pbiehLineUse[nchildindex]);
					thesigEdgeSplitText.addInputEventListener(borderTextHide);
					int nextscore = ncurrscore+nchildindex*(int) Math.pow(theTimeiohmm.nmaxchild,data[0].length-ndepth-1);
					if (htTextVisible != null)
					{
						boolean bshow = ((Boolean)htTextVisible.get(ndepth+";"+nextscore+";"+1)).booleanValue();
						thesigEdgeSplitText.setVisible(bshow);
						thesigEdgeSplitText.setPickable(bshowkeyinputs);
					} 

					hidesigList.add(new SigInfoRec(border, thesigEdgeSplitText,
							ptr.tsSigTFEdgeSplit[nchildindex],1,nextscore,ndepth,ptr));
				}

				final DREMGui ftheDREMGui = this;	   
				PBasicInputEventHandler pbiehCircle1 = new PBasicInputEventHandler()
				{
					public void mousePressed(PInputEvent event) 
					{
						if (event.getButton() == MouseEvent.BUTTON3)
						{
							javax.swing.SwingUtilities.invokeLater(new Runnable() 
							{
								public void run() 
								{
									JFrame frame = new JFrame("Split Table");
									frame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
									frame.setLocation(200,200);
									CircleID theCircleID = (CircleID) htNodes.get(fcircle);

									if (fptr.numchildren == 2)
									{
										DREMGui_SplitTable newContentPane = 
											new DREMGui_SplitTable(ftheDREMGui,frame, theTimeiohmm,fptr,theCircleID.nscore);
										newContentPane.setOpaque(true); //content panes must be opaque
										frame.setContentPane(newContentPane);
									}
									else
									{
										JTabbedPane tabbedPane = new JTabbedPane();

										for (int ntable = 0; ntable < fptr.numchildren; ntable++)
										{
											DREMGui_SplitTable newContentPane = 
												new DREMGui_SplitTable(ftheDREMGui,frame, theTimeiohmm,
														fptr,theCircleID.nscore,ntable,fptr.orderA[ntable]);
											String szLabel;
											String szToolTip;
											if (fptr.numchildren == 3)
											{
												if (ntable == 0)
												{
													szLabel = "Low vs. Others ";
													szToolTip = "Low";
												}
												else if (ntable == 1)
												{
													szLabel = "Middle vs. Others ";
													szToolTip = "Middle";
												}
												else
												{
													szLabel = "High vs. Others ";
													szToolTip = "High";
												}
											}
											else
											{
												szLabel = "Child "+ntable+" vs. Others ";
												szToolTip = "Child "+ ntable;
											}
											tabbedPane.addTab(szLabel,null,newContentPane,szToolTip);
										}
										tabbedPane.setOpaque(true);
										frame.setContentPane(tabbedPane);
									}

									//Display the window.
									frame.pack();
									frame.setVisible(true);
								}
							});
						}
					}
				};

				circle.addInputEventListener(pbiehCircle1);
				borderCircle.addInputEventListener(pbiehCircle1);
				thesigText.addInputEventListener(borderTextHide);
				borderCircleAct.addInputEventListener(pbiehCircle1);
				thesigTextAct.addInputEventListener(borderTextHide);
				borderCircleFirstAct.addInputEventListener(pbiehCircle1);
				thesigTextFirstAct.addInputEventListener(borderTextHide);
			}
			else
			{
				circle.addInputEventListener(pbiehLine);
			}

			if ((theSelectedRec.theCircleID != null) && (theSelectedRec.theCircleID.nid == ncircleID))
			{
				if (theSelectedRec.bcircle)
				{
					htColors.put(circle,(Color) circle.getPaint());
					circle.setPaint(Color.yellow);
					theSelectedRec.selectedNode = circle;
				}
				else
				{
					htColors.put(line,(Color) line.getStrokePaint());
					line.setStrokePaint(Color.yellow);
					theSelectedRec.selectedNode = line;
				}
			}   

			class myPBasicInputEventHandler extends PBasicInputEventHandler
			{
				boolean bcircle;

				myPBasicInputEventHandler(boolean bcircle)
				{
					this.bcircle = bcircle;
				}

				public void mousePressed(PInputEvent event) 
				{
					if (event.getButton() == MouseEvent.BUTTON1)
					{ 
						boolean badd = true;
						PNode pickednode =  event.getPickedNode();
						if (theSelectedRec.selectedNode != null)
						{
							CircleID theCircleID = theSelectedRec.theCircleID;
							htRequired.remove(theCircleID);
							//removes the current circle ID from the set of required circle IDs
							Color oldcolor = (Color) htColors.get(theSelectedRec.selectedNode);
							if (theSelectedRec.bcircle)
							{
								theSelectedRec.selectedNode.setPaint(oldcolor);
							}
							else
							{
								((PPath)theSelectedRec.selectedNode).setStrokePaint(oldcolor);
							}

							//sets paint back to original color 
							Enumeration enumPaths = htHidden.keys();

							while (enumPaths.hasMoreElements())
							{ 
								Integer intPath = (Integer) enumPaths.nextElement();
								int npath = intPath.intValue();

								int nval;
								if (theSelectedRec.bcircle)
								{
									nval = (int) Math.pow(theTimeiohmm.nmaxchild,data[0].length-
											theCircleID.nminparentlevel);
								}
								else
								{
									nval = (int) Math.pow(theTimeiohmm.nmaxchild,data[0].length-theCircleID.nprevminparentlevel);
								}

								if (DREM_Timeiohmm.BDEBUG)
								{
									System.out.println("B"+theCircleID.nscore+" "+(nval*(npath/nval)+" "+data[0].length+" "+
											theCircleID.nminparentlevel+" "+theCircleID.nprevminparentlevel+" "+
											nval+" "+npath+" "+theCircleID.ndepth));
								}

								if ((theCircleID.nscore/nval) != (npath/nval))
								{

									int nblocking = ((Integer) htHidden.get(intPath)).intValue();
									nblocking--;
									if (DREM_Timeiohmm.BDEBUG)
									{
										System.out.println("hereB "+nblocking);
									}

									htHidden.put(new Integer(npath), new Integer(nblocking));
									if (nblocking == 0)
									{
										ArrayList lines = (ArrayList) htPathToLineSet.get(intPath);
										int nsize = lines.size();
										for (int nindex = 0; nindex < nsize; nindex++)
										{
											int ngeneindex = ((Integer) lines.get(nindex)).intValue();
											bPathVisible[ngeneindex] = true;
											boolean bvisible = bglobalVisible&& bTFVisible[ngeneindex]
											                                               &&bGOVisible[ngeneindex]&&bSetVisible[ngeneindex];
											plArray[ngeneindex].setVisible(bvisible);
											plArray[ngeneindex].setPickable(bvisible);
										}	    
									}
								}
							}

							if (pickednode.equals(theSelectedRec.selectedNode))
							{
								//current node will no longer be yellow
								theSelectedRec.selectedNode = null;
								theSelectedRec.theCircleID = null;
								badd = false;
								setFilterText();
							}
						}

						if (badd)
						{
							CircleID theCircleID = (CircleID) htNodes.get(fcircle);
							htRequired.add(theCircleID);
							if (bcircle)
							{
								htColors.put(pickednode,(Color) pickednode.getPaint());
								pickednode.setPaint(Color.yellow);
							}
							else
							{
								htColors.put(pickednode,((PPath) pickednode).getStrokePaint());
								((PPath) pickednode).setStrokePaint(Color.yellow);
							}

							theSelectedRec.selectedNode = pickednode;
							theSelectedRec.bcircle = bcircle;
							theSelectedRec.theCircleID = theCircleID;
							setFilterText();

							Enumeration els = htPathToLineSet.keys();
							while (els.hasMoreElements())
							{
								Integer intPath = (Integer) els.nextElement();

								int nintPath = intPath.intValue();

								int nval;

								if (bcircle)
								{
									nval = (int) Math.pow(theTimeiohmm.nmaxchild,data[0].length-
											theCircleID.nminparentlevel);
								}
								else
								{
									nval = (int) Math.pow(theTimeiohmm.nmaxchild,data[0].length-theCircleID.nprevminparentlevel);
								}

								if (DREM_Timeiohmm.BDEBUG)
								{
									System.out.println("A"+theCircleID.nscore+" "
											+(nval*(nintPath/nval))+" "+data[0].length+" "+
											theCircleID.nminparentlevel+" "+theCircleID.nprevminparentlevel+" "+nval+" "+nintPath);
								}

								if ((theCircleID.nscore/nval) != (nintPath/nval))
								{	
									Integer intBlocking = (Integer) htHidden.get(intPath);
									if (DREM_Timeiohmm.BDEBUG)
									{
										System.out.println("hereA "+intBlocking);
									}

									if ((intBlocking == null)||(intBlocking.intValue()==0))
									{
										ArrayList lines = (ArrayList) htPathToLineSet.get(intPath);
										int nsize = lines.size();
										for (int nindex = 0; nindex < nsize; nindex++)
										{
											int ngeneindex = ((Integer) lines.get(nindex)).intValue();
											bPathVisible[ngeneindex] = false;
											plArray[ngeneindex].setVisible(false);
											plArray[ngeneindex].setPickable(false);
										}
										htHidden.put(intPath, new Integer(1));
									}
									else
									{
										int nBlocking = intBlocking.intValue();
										htHidden.put(intPath, new Integer(nBlocking+1));
									}
								}
							}

							if (!bholdedge)
							{
								if (bcircle)
								{
									ncolortime = Math.min(theCircleID.ndepth+1,data[0].length-1);
								}
								else
								{
									ncolortime = theCircleID.ndepth;
								}
								if (theYScalegui != null)
								{
									theYScalegui.theColorSlider.setValue(ncolortime);
								}

								setGeneColors();
							}
						}                          
					}
				}
			};

			PBasicInputEventHandler pbiehc = new myPBasicInputEventHandler(true);
			circle.addInputEventListener(pbiehc);

			if (line != null)
			{
				PBasicInputEventHandler pbiehl = new myPBasicInputEventHandler(false);
				line.addInputEventListener(pbiehl);
			}
			hideList.add(circle);

			circleSet.add(new CircleRec(circle, ddiameter,ncircleID));
			ncircleID++;
		}      

		return pbiehLine;
	}


	/**
	 * Record for the circle nodes on the interface
	 */
	static class CircleRec
	{

		double ddiameter;
		int nid;
		PNode circle;

		CircleRec(PNode circle, double ddiameter, int nid)
		{
			this.circle = circle;
			this.ddiameter = ddiameter;
			this.nid = nid;
		}
	}


	/**
	 * Comparator for circlerec
	 */
	static class CircleRecCompare implements Comparator
	{
		/**
		 * Nodes with greater diameter get lower priorty, then nodes with lower ID
		 */
		public int compare(Object c1, Object c2)
		{
			CircleRec cr1 = (CircleRec) c1;
			CircleRec cr2 = (CircleRec) c2;

			if (cr1.ddiameter > cr2.ddiameter)
				return -1;
			else if (cr1.ddiameter < cr2.ddiameter)
				return 1;
			else if (cr1.nid < cr2.nid)
				return -1;
			else if (cr1.nid > cr2.nid)
				return 1;
			else 
				return 0;

		}
	}
	
	// TODO more permanent solution
	/**
	 * Save the model and TF activities when running in non-interactive
	 * batch mode
	 * @param filename use this name to create the output files
	 */
	public void batchSave(String filename)
	{
		if (saveModelFrame == null)
		{
			saveModelFrame = new JFrame("Save Model to File");
			saveModelFrame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
			saveModelFrame.setLocation(400,300);
			// Use "this" here, which is not exactly what is done in drawmain
			DREMGui_SaveModel newContentPane = 
				new DREMGui_SaveModel(this.theTimeiohmm,treecopy,
						saveModelFrame,this);
//			newContentPane.setOpaque(true); 
			//content panes must be opaque
			saveModelFrame.setContentPane(newContentPane);
			//Display the window.
//			saveModelFrame.pack();
		}
		
		Container c = saveModelFrame.getContentPane();
		// The content pane shoul be a DREMGUI_SaveModel
		((DREMGui_SaveModel) c).batchSave(filename);
	}

}
