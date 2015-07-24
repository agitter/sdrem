package edu.cmu.cs.sb.drem;

import edu.cmu.cs.sb.core.*;
import javax.swing.*;
import java.util.*;
import java.text.*;
import java.io.*;
import java.awt.*;
import java.awt.event.*;
import java.util.zip.*;
import edu.umd.cs.piccolo.nodes.PText;

/**
 * This class implements the core methods for learning the DREM maps
 */
public class DREM_Timeiohmm 
{
	int nglobaltime;
	double nodepenalty =10;
	boolean bmodeldiff = true;
	static final boolean BDEBUG = false;
	static final boolean BDEBUGMODEL = false;
	boolean ballowmergeval =false;
	boolean bsavedchange = false;
	boolean BEQUALSTD = false;
	boolean BVITERBI =false;
	boolean BREGDREM =false;
	boolean bhasmerge = false;
	double MINPROB = .00000000000000000001;
	double DEFAULTSIGMA = .5;
	double MAXFUTUREITR = 30;
	double MERGEMIN = -.0015;
	double DELAYMIN = -.0015;

	double MERGEMINDIFF = 0;
	double DELAYMINDIFF = 0;
	int nrandomseed;
	DREMGui progressDREMGui;
	int nglobaliteration;
	int ninitsearchval;

	double RESPLITMIN = -.0015;
	double RESPLITMINDIFF = 0;

	boolean brealXaxisDEF;
	double dYaxisDEF;
	double dXaxisDEF;
	double dnodekDEF;
	int nKeyInputTypeDEF;
	double dKeyInputXDEF;
	double dpercentDEF;
	JButton currentButton;

	JLabel statusLabel;
	JLabel statusLabel15;
	JLabel statuscountLabel;

	static final String SCORESDIR = "TOPAALoc";
	static final String SZDELIM = "|;,";
	DREM_NaiveBayes filteredClassifier;

	double dprevouterbestlog;
	/**last trainlikelihood*/
	double dtrainlike;
	boolean bfilterbinding = false;
	double EPSILON =.01; 
	double EPSILONDIFF = .01;
	double BEPSILON =0.001;
	int nmaxchild;
	int ntotalcombined;//num filtered and non-filtered

	double[][] CONSTANTA;

	int MINPATH =5;
	static final double RIDGE = 1;
	static final double MINSIGMA = .0001;

	String[] testgenenames;
	/** The subset of theDataSet.data that belongs to the training set??? 
	 * The expression data: row are genes, columns are time points in experiments*/
	double[][] traindata;
	/** The subset of bindingSignGene that belongs to the training set??? 
	 * Each row is a gene and gives the interaction sign (-1,0,1) for TFs
	 * regulating the gene */
	int[][] trainSign;
	/** The subset of bindingValGene that belongs to the training set. 
	 * Each row is a gene and gives the interaction value for TFs
	 * regulating the gene */
	double[][] trainVal;

	String[] traingenenames;
	/** The subset of theDataSet.data that belongs to the training set??? 
	 * The expression data: row are genes, columns are time points in experiments*/
	double[][] testdata;
	/** The subset of bindingSignGene that belongs to the test set??? 
	 * Each row is a gene and gives the interaction sign (-1,0,1) for TFs
	 * regulating the gene */
	int[][] testSign;
	/** The subset of bindingValGene that belongs to the test set. 
	 * Each row is a gene and gives the interaction value for TFs
	 * regulating the gene */
	double[][] testVal;
	
	/** The subset of theDataSet.pmavalues that belongs to the training set
	 * 0 if data value is missing, non-zero if present??? */
	int[][] trainpma;
	/** The subset of theDataSet.pmavalues that belongs to the test set
	 * 0 if data value is missing, non-zero if present??? */
	int[][] testpma;

	/** The subset of bindingGeneIndex that belongs to the training set??? 
	 * Each row is a gene and gives the index of TFs regulating the gene */
	int[][] trainIndex;
	/** The subset of bindingGeneIndex that belongs to the test set??? 
	 * Each row is a gene and gives the index of TFs regulating the gene */
	int[][] testIndex;

	int[][] trainSignTF;
	double[][] trainValTF;
	/** The subset of bindingTFIndex that belongs to the training set??? 
	 * Each row is a TF and gives the index of genes the TF regulates */
	int[][] trainTFIndex;
	int[][] testSignTF;
	double[][] testValTF;
	/** The subset of bindingTFIndex that belongs to the test set??? 
	 * Each row is a TF and gives the index of genes the TF regulates */
	int[][] testTFIndex;


	/** Each row is a TF and gives the interaction value for the genes it regulates*/
	double[][] bindingValTF; 
	/** Each row is a gene and gives the interaction value for TFs regulating the gene*/
	double[][] bindingValGene;
	/**each row is a TF and gives the interaction sign (-1,0,1) for the genes it regulates*/
	int[][] bindingSignTF; 
	/**each row is a gene and gives the interaction sign (-1,0,1) for TFs regulating the gene*/
	int[][] bindingSignGene;
	/**each row is a TF and gives the index of genes the TF regulates */
	int[][] bindingTFindex;
	/**each row is a gene and gives the index of TFs regulating the gene */
	int[][] bindingGeneindex;

	String[] tfNames;
	/**A sorted array of all unique signs for TF-gene binding interactions???*/
	int[] dBindingSigns;
	HashMap htBinding;
	/**A set of all unique TF-gene interaction signs???*/
	HashSet hsUniqueInput;

	Treenode treeptr;
	NumberFormat nf2;
	NumberFormat nf3;
	double dbestlog;
	Treenode bestTree = null;
	DREM_DataSet theDataSet;


	double[][][] theInstancesTF;
	int[][][] theInstancesTFIndex;
	/** Duplicated copies of trainVal, where the first index + 2
	 * gives the number of copies of each gene.
	 * Within each duplicated set,
	 * each row is a gene and gives the interaction value for TFs regulating the gene???*/
	double[][][] theInstances;
	/** Duplicated copies of trainValIndex, where the first index + 2
	 * gives the number of copies of each gene.
	 * Within each duplicated set,
	 * each row is a gene and gives the index of TFs regulating the gene??? */
	int[][][] theInstancesIndex;
	int[][] ylabels;
	JLabel statusLabel2;
	JLabel statusLabel3;
	int numtotalPath = 1;
	double dbesttrainlike;

	/**Number of boolean items for the static vector.  Equal to the number of TFs*/
	int numbits;
	double TRAINRATIO = .75;

	Random theRandom;

	int ntrain;
	int ntest;
	int numcols;
	int nholdout;

	double HOLDOUTRATIO =0;
	double[][] holdoutdata;
	int[][] holdoutpma;
	int[][] holdoutSign;
	int[][] holdoutIndex;
	double[][] holdoutVal;

	boolean[] bholdout;
	int[] storedbestpath;
	String szinitfileval;
	String szexpname;
	String sznodepenaltyval;
	String szconvergenceval;
	ArrayList savedColors = new ArrayList();
	double dgloballike;
	private int ntotalID;

	Hashtable htBackPts = new Hashtable();

	/**The probability that binding is functional, i.e. that the genes
	 * a TF binds are regulated by the TF.  Used to calculated
	 * TF activity scores.  Can be set in the DREM defaults file. */
	double dProbBindingFunctional = 0.8;
	/** All TF activiy priors are adjusted by ACTIVITY_EPSILON
	 * to avoid divide by 0 errors when the TF activity prior is 1. */
	static final double ACTIVITY_EPSILON = 1E-8;
	/** The prior belief that a TF is active in the condition
	 * from which the expression data was obtained.  The belief
	 * is derived from the binding values of the TF's binding
	 * interactions.  Priors range from 0 to 1-ACTIVITY_EPSILON (inclusive). */
	double[] dTFActivityPrior;
	/**An array of activity scores for all TFs based on how the
	 * genes the TF regulates transition to the next states compared
	 * to how all genes into the split transistion.  The max is the
	 * max score over this node and its descendent nodes.*/
	double[] dMaxTFActivityScore;
	

	/**
	 * Class constructor - provides the execution control on the method
	 */
	public DREM_Timeiohmm(DREM_DataSet theDataSet,String szbinding, 
			String sznumchildval, String szepsilonval,String szprunepathval,
			String szdelaypathval, String szmergepathval, 
			String szepsilonvaldiff,String szprunepathvaldiff,
			String szdelaypathvaldiff, String szmergepathvaldiff,
			String szseedval,
			boolean bstaticcheckval,boolean ballowmergeval, int ninitsearchval,
			String szinitfileval,  JLabel statusLabel, JLabel statusLabel15,
			JLabel statusLabel2,JLabel statusLabel3, JLabel statuscountLabel,
			JButton endSearchButton,boolean bstaticsearchval,
			final boolean  brealXaxisDEF,final double dYaxisDEF,final double dXaxisDEF,
			final double dnodekDEF,
			final int nKeyInputTypeDEF,final double dKeyInputXDEF,final double dpercentDEF,
			final JButton currentButton, String szexpname,String sznodepenaltyval,
			boolean bpenalizedmodel, String szconvergenceval, double dProbBindingFunctional) throws Exception
			{
		this.szconvergenceval = szconvergenceval;
		this.sznodepenaltyval = sznodepenaltyval;
		this.statusLabel = statusLabel;
		this.statusLabel15 = statusLabel15;
		this.statusLabel2 = statusLabel2;
		this.statusLabel3 = statusLabel3;
		this.statuscountLabel = statuscountLabel;
		this.szexpname = szexpname;
		this.ballowmergeval = ballowmergeval;
		this.brealXaxisDEF =brealXaxisDEF;
		this.dYaxisDEF = dYaxisDEF;
		this.dXaxisDEF = dXaxisDEF;
		this.dnodekDEF = dnodekDEF;
		this.nKeyInputTypeDEF = nKeyInputTypeDEF;
		this.dKeyInputXDEF = dKeyInputXDEF;
		this.dpercentDEF = dpercentDEF;
		this.currentButton = currentButton;
		this.ninitsearchval = ninitsearchval;
		this.dProbBindingFunctional = dProbBindingFunctional;


		BREGDREM = (!bstaticsearchval);
		if (BDEBUG)
		{
			System.out.println(szepsilonval+"&&&&"+szseedval);
		}
		this.szinitfileval = szinitfileval;
		this.bfilterbinding = bstaticcheckval; 

		nf2 = NumberFormat.getInstance(Locale.ENGLISH);
		nf2.setMinimumFractionDigits(2);
		nf2.setMaximumFractionDigits(2);

		nf3 = NumberFormat.getInstance(Locale.ENGLISH);
		nf3.setMinimumFractionDigits(3);
		nf3.setMaximumFractionDigits(3);

		this.statusLabel2 = statusLabel2;
		this.statusLabel3 = statusLabel3;

		EPSILON = Double.parseDouble(szepsilonval)/100;
		DELAYMIN = Double.parseDouble(szdelaypathval)/100;
		RESPLITMIN  = Double.parseDouble(szprunepathval)/100;
		MERGEMIN  = Double.parseDouble(szmergepathval)/100;
		BEPSILON = Double.parseDouble(szconvergenceval)/100;

		EPSILONDIFF = Double.parseDouble(szepsilonvaldiff);
		DELAYMINDIFF = Double.parseDouble(szdelaypathvaldiff);
		RESPLITMINDIFF  = Double.parseDouble(szprunepathvaldiff);
		MERGEMINDIFF  = Double.parseDouble(szmergepathvaldiff);

		nrandomseed = Integer.parseInt(szseedval);
		theRandom = new Random(nrandomseed);
		nodepenalty = Double.parseDouble(sznodepenaltyval);
		this.bmodeldiff = bpenalizedmodel;

		nmaxchild = Integer.parseInt(sznumchildval);

		this.theDataSet = theDataSet;       

		CONSTANTA = new double[nmaxchild+1][];
		for (int nindex = 1; nindex <= nmaxchild; nindex++)
		{
			CONSTANTA[nindex] = new double[nindex];
			for (int njindex = 0; njindex < nindex; njindex++)
				CONSTANTA[nindex][njindex] = 1.0/nindex;
		}

		if (szbinding.equals(""))
		{
			BREGDREM = true;
		}
		else
		{
			readTFGeneData(szbinding);
			buildFilteredClassifier();
		}

		//compute mean and std at each time point of genes assigned to the profile

		splitdata();  
		treeptr = new Treenode();

		boolean binitfile = (!szinitfileval.equals(""));
		if ((binitfile)&&(ninitsearchval==0))
		{	   
			readInitTree(treeptr);
			computeOrders(treeptr);
			combineTrainAndTest();
			viterbi(treeptr);
			if (bindingSignGene != null)
			{
				computeStats(treeptr,treeptr);
				computeminparentlevel(treeptr);
			}
			bsavedchange = true;
		}
		else
		{     
			bsavedchange = false;

			if ((binitfile)&&(ninitsearchval==1))
			{
				readInitTree(treeptr);
			}
			else 
			{
				double[] dsigma = new double[traindata[0].length];
				double[] dmeans = new double[traindata[0].length];
				computeDataStats(traindata,trainpma, dsigma,dmeans);
				buildEmptyTree(0,treeptr,dsigma,dmeans);
			}

			searchstage1();
			int numintervals = traindata[0].length-1;
			int[] path = new int[numintervals];

			int nremoved = 0;


			if (endSearchButton != null)
			{
				endSearchButton.setEnabled(false);
			}

			if (statusLabel3 != null)
			{
				statusLabel3.setText(" ");
			}

			if (statusLabel2 != null)
			{
				statusLabel2.setText(" ");    
			}

			int NUMRESPLITS;
			if (bmodeldiff)
			{
				NUMRESPLITS = 0;
			}
			else
			{
				NUMRESPLITS = 1;
			}

			for (int ni = 1; ni <= NUMRESPLITS; ni++)
			{
				if (statusLabel3 != null)
				{
					statusLabel3.setText(" Deleting Paths: "+nremoved+" removed so far ");
				}

				splitdata();
				if (BVITERBI)
				{
					dbestlog = trainhmmV(treeptr,true);
				}
				else
				{
					dbestlog = trainhmm(treeptr,true);
				}

				dbesttrainlike = dtrainlike;

				if (!bmodeldiff)
				{
					if (BVITERBI)
					{
						dbestlog = testhmmV(treeptr);
					}
					else
					{
						dbestlog = testhmm(treeptr);
					}
				}

				if (BDEBUG)
				{
					System.out.println("dbestlog = "+dbestlog);
				}

				boolean bdeleteagain = false;
				do
				{
					if (BDEBUG)
					{
						System.out.println("trying to delete after resplit");
					}

					traverse(bestTree, 0,false);
					dprevouterbestlog = dbestlog; 
					dbestlog = Double.NEGATIVE_INFINITY;
					bdeleteagain = traverseanddelete(path,treeptr,treeptr,true,
							RESPLITMIN,RESPLITMINDIFF);
					if (dbestlog == Double.NEGATIVE_INFINITY)
						dbestlog = dprevouterbestlog;

					if (BDEBUG)
					{
						System.out.println("****** delete  best log likelihood "+dbestlog);
					}

					traverse(bestTree, 0,true);
					if (bestTree != treeptr)
					{
						nremoved++;
						numtotalPath--;
						if (statuscountLabel != null)
						{
							statuscountLabel.setText(" Number of paths in model so far: "+numtotalPath);
							statusLabel3.setText(" Deleting Paths: "+nremoved+" removed so far");
							statusLabel2.setText(" "); 
						}
						treeptr = bestTree;
					}

					if (BDEBUG)
					{
						System.out.println("CCC\tdelete\t"+dbesttrainlike+"\t"+dbestlog+"\t"+dprevouterbestlog);
					}
				}
				while (bdeleteagain);
			}

			boolean bagaindelay;
			boolean bagaindelayouter;

			int numdelay = 0;

			do
			{
				bagaindelayouter = false;
				for (int ndesiredlevel = 1; ndesiredlevel < numintervals; ndesiredlevel++)
				{
					do
					{
						if (!bmodeldiff)
						{
							if (BVITERBI)
							{
								dbestlog = testhmmV(treeptr);
							}
							else
							{
								dbestlog = testhmm(treeptr);
							}
						}

						if (BDEBUG)
						{
							System.out.println("trying to delay "+ndesiredlevel+" "+dbestlog);
						}

						dprevouterbestlog = dbestlog; 

						dbestlog = Double.NEGATIVE_INFINITY;
						bagaindelay = traverseanddelay(path,treeptr,ndesiredlevel,treeptr,false);

						if (dbestlog == Double.NEGATIVE_INFINITY)
						{
							dbestlog = dprevouterbestlog;
						}
						else
						{ 
							numdelay++;
							if (statusLabel3 != null)
							{
								statusLabel3.setText(" Number of delayed is "+numdelay+" so far");
							}
							treeptr = bestTree;  

							double dimprovedelay = (dbestlog-dprevouterbestlog)/Math.abs(dbestlog);
							String szimprovedelay = nf3.format(100*dimprovedelay) + "%";
							String szimprovedelayDIFF = nf3.format(dbestlog-dprevouterbestlog);

							if (statusLabel2 != null)
							{
								statusLabel2.setText(" Delay improvement: "+szimprovedelayDIFF+" ("+szimprovedelay+")");
							}
						}

						if (BDEBUG)
						{
							System.out.println("****** delay  best log likelihood "+dbestlog);
						}
						traverse(bestTree, 0,false);

						treeptr = bestTree;

						if (BDEBUG)
						{
							System.out.println("YYY\tdelay\t"+dbesttrainlike+"\t"+dbestlog+"\t"+dprevouterbestlog);
						}

						if (bagaindelay)
						{
							bagaindelayouter = true;
						}
					}
					while (bagaindelay);
				}
			}
			while (bagaindelayouter);	 

			if (ballowmergeval)
			{
				traverseandmerge(path);
			}
			treeptr = bestTree;

			traverse(bestTree, 0,false);
			if (BDEBUG)
			{
				System.out.println("Test set "+dbestlog);
				System.out.println("calling combine tranandtest");
			}

			combineTrainAndTest();

			if (BVITERBI)
			{
				trainhmmV(treeptr,true);
			}
			else
			{
				double dtemplog = trainhmm(treeptr,true);
				nglobaltime++;
				int numnodes = countNodes(treeptr,nglobaltime);
				if (BDEBUGMODEL)
				{ 
					System.out.println(nodepenalty+"\t final train is\t"+dtemplog+"\t"+
							numnodes+"\t"+(dtemplog+numnodes*nodepenalty)+"\t"+ dgloballike);
				}
			}
			traverse(bestTree, 0,false);

			boolean bagain;

			do 
			{
				bagain = false;

				computeOrders(treeptr); 
				viterbi(treeptr);

				int[] bestpath = new int[path.length];
				for (int nindex = 0; nindex < bestpath.length; nindex++)
				{
					bestpath[nindex] = -1;
				}

				MinPathRec theMinPathRec = traverseanddeleteMinPath(path,bestpath,treeptr);
				if (BDEBUG)
				{
					System.out.println("min is "+theMinPathRec.nval+"\t"+"level is "+theMinPathRec.nlevel);
					for (int nindex = 0; nindex < bestpath.length; nindex++)
					{
						System.out.println(theMinPathRec.bestpath[nindex]);
					}
				}

				if (theMinPathRec.nval < MINPATH)
				{
					deleteMinPath(theMinPathRec.bestpath,theMinPathRec.nlevel,treeptr); 
					bagain = true;
					traverse(treeptr, 0,true);

					if (BVITERBI)
					{
						trainhmmV(treeptr,true);
					}
					else
					{
						trainhmm(treeptr,true);
					}

					if (BDEBUG)
					{
						System.out.println("after retrain");
					}

					traverse(treeptr, 0,true);

					nremoved++;
					numtotalPath--;
					if (statuscountLabel != null)
					{
						statuscountLabel.setText(" Number of paths in model so far: "+numtotalPath);
						statusLabel3.setText(" Deleting Paths: "+nremoved+" removed so far");
						statusLabel2.setText(" ");
					}
				}
			} while (bagain);

			if (bindingSignGene !=null)
			{
				computeStats(treeptr,treeptr);
				computeminparentlevel(treeptr);
			}
		}
			}


	/**
	 * Sets up the stored data so that the training data includes all data except any external held out validation
	 * data, meaning the data which might otherwise used for model selection is just for parameter learning
	 */
	public void combineTrainAndTest()
	{
		nholdout = (int) (theDataSet.data.length * HOLDOUTRATIO);
		ntrain  = (theDataSet.data.length-nholdout);
		numcols = theDataSet.data[0].length;
		traindata = new double[ntrain][numcols];
		trainpma = new int[ntrain][numcols];
		trainIndex = new int[ntrain][];
		trainTFIndex = new int[numbits][];
		trainVal = new double[ntrain][];
		trainValTF = new double[numbits][];
		traingenenames = new String[ntrain];

		int ntrainindex = 0;

		for (int nindex = 0; nindex < theDataSet.data.length; nindex++)
		{
			if (!bholdout[nindex])
			{
				for (int ncol = 0; ncol < numcols; ncol++)
				{
					traindata[ntrainindex][ncol] = theDataSet.data[nindex][ncol];
					trainpma[ntrainindex][ncol]  = theDataSet.pmavalues[nindex][ncol]; 
				}

				trainVal[ntrainindex] = bindingValGene[nindex];
				trainIndex[ntrainindex] = bindingGeneindex[nindex];
				traingenenames[ntrainindex] = theDataSet.genenames[nindex];
				ntrainindex++;  
			}
		}
		makeTFindex(trainVal,trainIndex,trainValTF,trainTFIndex);
		trainSign = sign(trainVal);
		trainSignTF = sign(trainValTF);

		loadInstances();
	}

	/**
	 * Splits the non-held out validation data into a training set for learning the model parameters
	 * and a test set for model selection
	 */
	public void splitdata()
	{

		//splits data into training and testsets

		nholdout = (int) (theDataSet.data.length * HOLDOUTRATIO);
		if (bmodeldiff)
		{
			ntrain  = (theDataSet.data.length-nholdout);
			ntest   = 0;
		}
		else
		{
			ntrain  = (int)((theDataSet.data.length-nholdout)*TRAINRATIO);
			ntest   = theDataSet.data.length - ntrain - nholdout;
		}

		numcols = theDataSet.data[0].length;


		//build arrays with time-series training data
		//and binding p-value
		testgenenames = new String[ntest];

		testdata  = new double[ntest][numcols];
		testpma = new int[ntest][numcols];
		traindata = new double[ntrain][numcols];
		trainpma = new int[ntrain][numcols];
		testVal  = new double[ntest][];
		trainVal = new double[ntrain][];
		trainValTF = new double[numbits][];
		testValTF = new double[numbits][];
		trainIndex = new int[ntrain][];
		testIndex = new int[ntest][];
		trainTFIndex = new int[numbits][];
		testTFIndex = new int[numbits][];


		double[] foldrandom = new double[theDataSet.data.length];
		double[] foldrandomcopy = new double[theDataSet.data.length];

		// If a set of genes to be held out does not exist yet, create one
		if (bholdout == null)
		{
			holdoutdata = new double[nholdout][numcols];
			bholdout = new boolean[theDataSet.data.length];
			holdoutVal = new double[nholdout][numbits];
			holdoutpma  = new int[nholdout][numcols];
			holdoutIndex = new int[nholdout][];

			// Assign random values and store a copy of them
			for (int nindex = 0; nindex < foldrandom.length; nindex++)
			{
				foldrandom[nindex] =  theRandom.nextDouble();
				foldrandomcopy[nindex] = foldrandom[nindex];
			}
			Arrays.sort(foldrandomcopy);
			double dcutoff = foldrandomcopy[holdoutVal.length];

			int nholdoutindex = 0;

			for (int nindex = 0; nindex < foldrandom.length; nindex++)
			{
				// The cutoff was chosen to ensure that the desired
				// number of genes are held out
				if (foldrandom[nindex] < dcutoff)
				{
					int nnonzero = 0;
					holdoutVal[nholdoutindex] = bindingValGene[nindex];
					holdoutIndex[nholdoutindex] = bindingGeneindex[nindex];

					for (int ncol = 0; ncol < numcols; ncol++)
					{
						holdoutdata[nholdoutindex][ncol] =  theDataSet.data[nindex][ncol]; 
						holdoutpma[nholdoutindex][ncol] = theDataSet.pmavalues[nindex][ncol];
					}

					bholdout[nindex] = true;
					nholdoutindex++;
				}
				else
				{
					bholdout[nindex] = false;
				}
			}
			
			holdoutSign = sign(holdoutVal);
		}

		/*---------------------------------------------*/
		//random drawing a set of ntrain elements from data.length elements
		int[] includetrain = new int[ntrain];
		boolean[] btrain = new boolean[theDataSet.data.length];
		if (BDEBUG)
		{
			System.out.println(ntrain+" "+theDataSet.data.length);
		}

		int nincludeindex = 0;
		int nbtrainindex= 0;
		int nonholdout = 0;
		while (nincludeindex < ntrain)
		{
			if (!bholdout[nbtrainindex])
			{
				includetrain[nincludeindex] = nbtrainindex;
				btrain[nbtrainindex] = true;
				nincludeindex++;
				nonholdout++;
			}
			nbtrainindex++;
		}
		if (BDEBUG)
		{
			System.out.println("done");
		}

		//drawing nrandom elements from a set of numtotalannotated elements
		//where each element is equally likely
		for (; nbtrainindex < btrain.length; nbtrainindex++)
		{
			if (!bholdout[nbtrainindex])
			{
				if (theRandom.nextDouble() < ((double) ntrain/(double) (nonholdout+1)))
				{
					int nreplaceindex =(int) Math.floor(ntrain*theRandom.nextDouble());
					btrain[nbtrainindex] = true;
					btrain[includetrain[nreplaceindex]] = false;
					includetrain[nreplaceindex] = nbtrainindex;
				}
				else
				{
					btrain[nbtrainindex] = false;
				}
				nonholdout++;
			}
		}
		/*------------------------------------------------------*/

		int ntrainindex = 0;
		int ntestindex = 0;  
		for (int nindex = 0; nindex < btrain.length; nindex++)
		{
			if (!bholdout[nindex])
			{
				if (btrain[nindex])
				{
					for (int ncol = 0; ncol < numcols; ncol++)
					{
						traindata[ntrainindex][ncol] = theDataSet.data[nindex][ncol]; 
						trainpma[ntrainindex][ncol] = theDataSet.pmavalues[nindex][ncol];
					}

					trainVal[ntrainindex] = bindingValGene[nindex];
					trainIndex[ntrainindex] = bindingGeneindex[nindex];
					ntrainindex++;

				} 
				else
				{
					testgenenames[ntestindex] = theDataSet.genenames[nindex];

					for (int ncol = 0; ncol < numcols; ncol++)
					{
						testdata[ntestindex][ncol]=theDataSet.data[nindex][ncol]; 
						testpma[ntestindex][ncol]  = theDataSet.pmavalues[nindex][ncol];
					}

					testVal[ntestindex] = bindingValGene[nindex];
					testIndex[ntestindex] = bindingGeneindex[nindex];
					ntestindex++;

				}
			}
		}

		makeTFindex(trainVal,trainIndex,trainValTF,trainTFIndex);
		makeTFindex(testVal,testIndex,testValTF,testTFIndex);

		testSign  = sign(testVal);
		trainSign = sign(trainVal);
		trainSignTF = sign(trainValTF);
		testSignTF = sign(testValTF);
		
		loadInstances();
	}

	
	/**
	 * Generates pvalTF and pvalTFindex based on pval and pvalindex so that the rows correspond to each TF
	 * and each entry in a row pvalTFindex corresponds to an index of a gene the TF regulates and the corresponding 
	 * entry in pvalTF corresponds to the interaction value.
	 */
	public void makeTFindex(double[][] pval, int[][] pvalindex, double[][] pvalTF, int[][] pvalTFindex)
	{
		// Use counts to determine how large to make the array at each row
		// of pvalTF and pvalTFindex.  counts will contain the number
		// of genes the TF interacts with
		int[] counts = new int[pvalTF.length];
		for (int i = 0; i < pvalindex.length; i++)
		{
			for (int j = 0; j < pvalindex[i].length; j++)
			{
				counts[pvalindex[i][j]]++;
			}
		}

		for (int i= 0; i < pvalTFindex.length; i++)
		{
			pvalTF[i] = new double[counts[i]];
			pvalTFindex[i] = new int[counts[i]];
		}

		//reusing counts as number found so far
		for (int i = 0; i < counts.length; i++)
		{
			counts[i] = 0;
		}

		for (int ngene = 0; ngene < pvalindex.length; ngene++)
		{
			for (int nregindex = 0; nregindex < pvalindex[ngene].length; nregindex++)
			{
				int nTF = pvalindex[ngene][nregindex];
				int nTFregindex = counts[pvalindex[ngene][nregindex]];
				pvalTF[nTF][nTFregindex] = pval[ngene][nregindex] ;
				pvalTFindex[nTF][nTFregindex] = ngene;
				counts[nTF]++;
			}
		}
	}

	/**
	 * This class is responsible for manipulating the TF-gene data so that
	 * it will be in format readily usable by a weighted classifier.
	 * For a split with n children, the classifier needs n copies of each instance
	 * because each instance will be assumed to have every label, but with different
	 * weights for different labels.
	 */
	void loadInstances()
	{
		//first dimension determines the number of instances we will have
		//each second dimension corresponds to an instance
		//the third dimension is the values in the instances

		theInstances = new double[nmaxchild-1][][];
		ylabels = new int[nmaxchild-1][];
		theInstancesIndex = new int[nmaxchild-1][][];
		theInstancesTF = new double[nmaxchild-1][][];
		theInstancesTFIndex = new int[nmaxchild-1][][];

		// numchild is the number of child nodes there will be coming out
		// of the split (i.e. classes).
		for (int numchild = 2; numchild <= nmaxchild; numchild++)
		{
			theInstances[numchild-2] = new double[numchild*trainSign.length][];
			ylabels[numchild-2] = new int[numchild*trainSign.length];
			theInstancesIndex[numchild-2]  = new int[numchild*trainSign.length][];
			theInstancesTF[numchild-2] = new double[numchild*numbits][];
			theInstancesTFIndex[numchild-2] = new int[numchild*numbits][];

			int nspot = 0;
			for (int ngene = 0; ngene < trainVal.length; ngene++)
			{
				// Make numchild copies of the binding data and assign the class labels
				for (int nchild = 0; nchild < numchild; nchild++)
				{
					theInstances[numchild-2][nspot] = trainVal[ngene];
					ylabels[numchild-2][nspot] = nchild;
					theInstancesIndex[numchild-2][nspot] = trainIndex[ngene];
					nspot++;
				}
			}

			makeTFindex(theInstances[numchild-2],theInstancesIndex[numchild-2],
					theInstancesTF[numchild-2],theInstancesTFIndex[numchild-2]);

		}
	}


	/**
	 * Record of index and regulatory value 
	 */
	private static class TFRegRec
	{
		int ntfindex;
		double dregval;

		TFRegRec(int ntfindex, double dregval)
		{
			this.ntfindex = ntfindex;
			this.dregval = dregval;
		}
	}

	///////////////////////////////////////////////////////
	/**
	 * Reads in the tf-gene interaction data currently in the format TF, gene, and if specified interaction value
	 * where each column is delimitted by tabs
	 */
	private void parseThreeColFormat(BufferedReader br) throws IOException
	{
		HashMap htGeneToTFArray = new HashMap();
		HashMap htTFtoInteger = new HashMap();
		int ntf = 0;
		String szLine;
		while ((szLine = br.readLine())!=null)
		{
			StringTokenizer st = new StringTokenizer(szLine,"\t");
			String szTF = st.nextToken();
			String szGene = st.nextToken();
			double dinput;
			if (st.hasMoreTokens())
			{
				String szToken = st.nextToken();
				try
				{
					dinput = Double.parseDouble(szToken);
					if(Double.isInfinite(dinput) || Double.isNaN(dinput))
					{
						throw new IllegalArgumentException(szToken +" is not a valid score for a TF-gene interaction");
					}
				}
				catch(NumberFormatException nfex)
				{
					throw new IllegalArgumentException(szToken +" is not a valid score for a TF-gene interaction");
				}
			}
			else
			{
				// If not binding value is given, the default value is 1.0
				dinput = 1;
			}

			Integer objInt = (Integer) htTFtoInteger.get(szTF);
			int ncurrTF;
			if (objInt == null)
			{
				ncurrTF = ntf;
				htTFtoInteger.put(szTF, new Integer(ntf));
				ntf++;
			}
			else
			{
				ncurrTF = ((Integer) objInt).intValue();
			}
			ArrayList al = (ArrayList) htGeneToTFArray.get(szGene);
			if (al == null)
			{
				al = new ArrayList();
			}
			al.add(new TFRegRec(ncurrTF,dinput));
			htGeneToTFArray.put(szGene,al);
		}
		br.close();

		numbits = htTFtoInteger.size();
		tfNames = new String[numbits];
		Iterator itrTFSet = htTFtoInteger.keySet().iterator();

		while (itrTFSet.hasNext())
		{
			String szTF = (String) itrTFSet.next();
			tfNames[((Integer) htTFtoInteger.get(szTF)).intValue()] = szTF;
		}


		Iterator itrGeneSet = htGeneToTFArray.keySet().iterator();
		while (itrGeneSet.hasNext())
		{
			String szid = (String) itrGeneSet.next();
			ArrayList al = (ArrayList) htGeneToTFArray.get(szid);
			int numnonzero = al.size();


			TFRegRec[] thePairsRecs = new TFRegRec[numnonzero];
			for (int nbit = 0; nbit < numnonzero; nbit++)
			{
				thePairsRecs[nbit] = (TFRegRec) al.get(nbit);
			}
			Arrays.sort(thePairsRecs, new TFRegRecCompare());
			int[] nonzeroindex = new int[numnonzero];
			double[] nonzerovals = new double[numnonzero];
			for (int i = 0; i < numnonzero; i++)
			{
				TFRegRec theTFRec = thePairsRecs[i];
				nonzeroindex[i] = theTFRec.ntfindex;
				nonzerovals[i] = theTFRec.dregval;
			}
			loadBinding(szid, nonzeroindex, nonzerovals, numnonzero);
		}
	}

	//////////////////////////////////////////////////////////////////////////////////
	/**
	 * Comparator object sorts strictly based on ntfindex
	 */
	private class TFRegRecCompare implements Comparator
	{
		public int compare(Object o1, Object o2)
		{
			TFRegRec tfr1 = (TFRegRec) o1;
			TFRegRec tfr2 = (TFRegRec) o2;

			if (tfr1.ntfindex < tfr2.ntfindex)
			{
				return -1;
			}
			else if (tfr1.ntfindex > tfr2.ntfindex)
			{
				return 1;
			}
			else
			{
				return 0;
			}
		}
	}


	/////////////////////////////////////////////////////////////////////////////////
	/**
	 * Reads in the TF-gene input data assuming the header is already read and 
	 * the data is in grid format
	 */
	private void parseGridFormat(BufferedReader br) throws IOException
	{
		int[] nonzeroindex = new int[numbits];
		double[] nonzerovals = new double[numbits];
		String szLine;
		StringTokenizer st;

		while ((szLine = br.readLine()) != null)
		{
			int numnonzero = 0;

			if (!szLine.trim().equals(""))
			{
				st = new StringTokenizer(szLine,"\t");
				int numTokens = st.countTokens();
				if (numbits != (numTokens-1))
				{
					throw new IllegalArgumentException(
							"Found a line with "+(numTokens-1)+" entries, expecting "+numbits);
				}
				String szid = st.nextToken();

				for (int nindex = 0; nindex < numbits; nindex++)
				{
					String sznum = st.nextToken();
					try
					{   
						double tempval = Double.parseDouble(sznum);
						if(Double.isInfinite(tempval) || Double.isNaN(tempval))
						{
							throw new IllegalArgumentException(sznum+" is not a valid value for a TF-gene interaction!");
						}
						else if (tempval != 0)
						{
							nonzeroindex[numnonzero] = nindex;
							nonzerovals[numnonzero] = tempval;
							numnonzero++;
						}
					}
					catch (NumberFormatException ex)
					{
						if ((numbits == 2) && (nindex == 0))
						{
							throw new IllegalArgumentException("If TF-gene data is in column format, then the first two "+
							"columns must have the headers 'TF' and 'Gene'");
						}
						else
						{
							throw new IllegalArgumentException(sznum+" is not a valid value for a TF-gene interaction!");
						}
					}
				}

				//storing each vector of p-values in hashtable with gene identifier
				loadBinding(szid, nonzeroindex, nonzerovals,numnonzero);
			}
		}
		br.close();
	}

	/////////////////////////////////////////// 
	/**
	 * Responsible for loading the TF-gene interaction input data from the file into the variable fields
	 */
	public void readTFGeneData(String szbinding) throws IOException
	{
		htBinding = new HashMap();
		hsUniqueInput = new HashSet();
		BufferedReader br = null;

		try
		{
			br = new BufferedReader(new InputStreamReader(
					new GZIPInputStream(new FileInputStream(szbinding))));
		}
		catch (IOException ex)
		{
			br = new BufferedReader(new FileReader(szbinding));
		}

		String szLine = br.readLine();
		StringTokenizer st = new StringTokenizer(szLine,"\t");
		String szh1="";

		if (szLine == null)
		{
			throw new IllegalArgumentException("Empty TF-gene interaction input file found!");
		}
		else if (szLine.startsWith("\t"))
		{
			numbits = st.countTokens();
		}
		else
		{
			numbits = st.countTokens()-1;
			szh1 = st.nextToken();
		}


		tfNames = new String[numbits];

		for (int ntfindex = 0; ntfindex < numbits; ntfindex++)
		{
			tfNames[ntfindex] = st.nextToken();
		}


		boolean bthreecol = (((numbits ==2))&&(szh1.equalsIgnoreCase("TF"))
				&&(tfNames[0].equalsIgnoreCase("GENE")));

		if (bthreecol)
		{
			parseThreeColFormat(br);
		}
		else
		{
			parseGridFormat(br);
		}
		

		// 0 is always a possible binding interaction value even if it doesn't appear
		// in the data ???
		hsUniqueInput.add(new Integer(0));
		if (bfilterbinding)
		{
			//filter those genes without binding data
			int nbinding = 0;  //count of the number of genes with binding info

			//first determine how many genes we have binding data for
			Object[] hitA  = new Object[theDataSet.data.length];
			for (int nrow = 0; nrow < theDataSet.data.length; nrow++)
			{
				Object obj = getBindingObject(theDataSet.genenames[nrow],theDataSet.probenames[nrow]); 

				if (obj != null)
				{
					hitA[nrow] = obj;
					nbinding++;
				}
				else
				{
					hitA[nrow] = null;
				}
			}

			if (BDEBUG)
			{
				System.out.println("nbinding = "+nbinding);
			}

			bindingSignGene= new int[nbinding][];
			bindingGeneindex = new int[nbinding][];

			int nbindingindex = 0;
			// stores whether or not should keep gene
			boolean[] bbindingdata = new boolean[theDataSet.data.length];
			

			// going through all data values storing into 
			// binding interaction values and signs
			// flagging for filtering those records without binding
			for (int nrow = 0; nrow < theDataSet.data.length; nrow++)
			{				
				if (hitA[nrow] != null)
				{
					BindingGeneRec theBindingGeneRec = (BindingGeneRec) hitA[nrow]; 

					bindingValGene[nbindingindex] = theBindingGeneRec.bindingValGene;
					bindingSignGene[nbindingindex] = sign(theBindingGeneRec.bindingValGene);
					bindingGeneindex[nbindingindex] = theBindingGeneRec.bindingGeneindex; 
					for (int nindex =0; nindex < bindingSignGene[nbindingindex].length; nindex++)
					{
						hsUniqueInput.add(new Integer(bindingSignGene[nbindingindex][nindex]));
					}
					nbindingindex++;
					bbindingdata[nrow] = true;
				}
				else
				{
					bbindingdata[nrow] = false;
				}
			}

			theDataSet = (DREM_DataSet) theDataSet.filtergenesgeneral(bbindingdata, nbinding, true);
		}
		else
		{
			//not filtering genes with missing binding values instead setting to 0
			//transfering binding values in a hashmap to an array
			//0 values if do not have a binding value for that gene
			bindingValGene = new double[theDataSet.data.length][];
			bindingSignGene = new int[theDataSet.data.length][];
			bindingGeneindex = new int[theDataSet.data.length][];

			for (int nrow = 0; nrow < theDataSet.data.length; nrow++)
			{
				Object obj = getBindingObject(theDataSet.genenames[nrow], theDataSet.probenames[nrow]);

				if (obj != null)
				{
					BindingGeneRec theBindingGeneRec = (BindingGeneRec) obj;
					bindingValGene[nrow] = theBindingGeneRec.bindingValGene;
					// Store the sign of each interaction value instead of the value itself
					bindingSignGene[nrow] = sign(theBindingGeneRec.bindingValGene);
					bindingGeneindex[nrow] = theBindingGeneRec.bindingGeneindex; 
					for (int nindex =0; nindex < bindingSignGene[nrow].length; nindex++)
					{
						hsUniqueInput.add(new Integer(bindingSignGene[nrow][nindex]));
					}
				}
				else
				{
					bindingValGene[nrow] = new double[0];
					bindingSignGene[nrow] = new int[0];
					bindingGeneindex[nrow] = new int[0];
				}
			}
		}

		System.out.println("Number of selected genes is "+bindingSignGene.length); 
		bindingTFindex = new int[numbits][];
		bindingValTF = new double[numbits][];
		makeTFindex(bindingValGene, bindingGeneindex, bindingValTF, bindingTFindex);
		bindingSignTF = sign(bindingValTF);
		
		Iterator itrKeys = hsUniqueInput.iterator();
		int nel = 0;

		dBindingSigns = new int[hsUniqueInput.size()];
		while (itrKeys.hasNext())
		{
			Integer key = (Integer) itrKeys.next();
			dBindingSigns[nel] = key.intValue();
			nel++;
		}
		Arrays.sort(dBindingSigns);
		
		// The prior belief for a TF's activity is the maximum
		// binding value it takes in either direction (activation or repression,
		// i.e. + or -)
		// Also initialize the max TF activity score
		dMaxTFActivityScore = new double[numbits];
		dTFActivityPrior = new double[numbits];
		double dmaxPrior = 0;
		for(int ntf = 0; ntf < bindingValTF.length; ntf++)
		{
			dMaxTFActivityScore[ntf] = 0;
			dTFActivityPrior[ntf] = 0;
			for(int ngene = 0; ngene < bindingValTF[ntf].length; ngene++)
			{
				dTFActivityPrior[ntf] =
					Math.max(dTFActivityPrior[ntf], Math.abs(bindingValTF[ntf][ngene]));
			}
			
			dmaxPrior = Math.max(dmaxPrior, dTFActivityPrior[ntf]);
		}
		
		// Activity priors must be in the range [0,1].  Because the absolute
		// value of the binding values was used, no prior will be < 0.
		// However, if any value is > 1, the priors must be normalized.
		if(dmaxPrior > 1)
		{
			for(int ntf = 0; ntf < bindingValTF.length; ntf++)
			{
				dTFActivityPrior[ntf] /= dmaxPrior;
			}
		}
		
		// To avoid dividing by 0 when a TF activity prior is 1, adjust
		// all priors by ACTIVITY_EPSILON.  Priors are allowed to be 0.
		for(int ntf = 0; ntf < bindingValTF.length; ntf++)
		{
			if(dTFActivityPrior[ntf] > ACTIVITY_EPSILON)
			{
				dTFActivityPrior[ntf] -= ACTIVITY_EPSILON;
			}
			else
			{
				dTFActivityPrior[ntf] = 0;
			}
		}
	}
	
	/**
	 * Finds the sign of the input value, where the sign is -1, 0, or 1.
	 * Returns 0 if the input is neither > 0 or < 0.
	 * @param val The input value
	 * @return The sign of val
	 */
	static int sign(double val)
	{
		if(val < 0)
		{
			return -1;
		}
		else if(val > 0)
		{
			return 1;
		}
		else
		{
			return 0;
		}
	}
	
	/**
	 * Finds the sign of the input value for every element in the array,
	 * where the sign is -1, 0, or 1.
	 * Returns 0 if the input is neither > 0 or < 0.
	 * @param valArray The input values
	 * @return A new array where the sign function has been applied to every value in
	 * the input array
	 */
	static int[] sign(double[] valArray)
	{
		int[] signArray = new int[valArray.length];
		for(int nindex = 0; nindex < valArray.length; nindex++)
		{
			signArray[nindex] = sign(valArray[nindex]);
		}
		
		return signArray;
	}
	
	/**
	 * Finds the sign of the input value for every element in the 2d array,
	 * where the sign is -1, 0, or 1.
	 * Returns 0 if the input is neither > 0 or < 0.
	 * @param valArray The input values
	 * @return A new 2d array where the sign function has been applied to every value in
	 * the input array
	 */
	static int[][] sign(double[][] valArray)
	{
		int[][] signArray = new int[valArray.length][];
		for(int nindex = 0; nindex < valArray.length; nindex++)
		{
			signArray[nindex] = sign(valArray[nindex]);
		}
		
		return signArray;
	}


	/**
	 * A recording contain an array of interaction values and annother array containing the indicies they
	 * correspond to
	 */
	static class BindingGeneRec
	{
		double[] bindingValGene;
		int[] bindingGeneindex;
	}

	/**
	 * Assume geneindex1 and geneindex2 are sorted arrays
	 * Returns a BindingGeneRec where bindingGeneindex has the union of geneindex1 and geneindex2
	 * and bindingValGene has the corresponding values using the maximum to resolve disagreements
	 */
	private static BindingGeneRec mergeArrays(double[] geneval1, int[] geneindex1, double[] geneval2, int[] geneindex2)
	{

		int nindex1 = 0;
		int nindex2 = 0;
		int nmatch = 0;

		while ((nindex1 < geneindex1.length)&&(nindex2 < geneindex2.length))
		{
			if (geneindex1[nindex1] == geneindex2[nindex2])
			{
				nmatch++;
				nindex1++;
				nindex2++;
			}
			else if (geneindex1[nindex1] < geneindex2[nindex2])
			{
				nindex1++;
			}
			else if (geneindex2[nindex2] < geneindex1[nindex1])
			{
				nindex2++;
			}
		}

		int nmergecount = geneindex1.length + geneindex2.length - nmatch;

		BindingGeneRec theBindingGeneRec = new BindingGeneRec();
		theBindingGeneRec.bindingValGene = new double[nmergecount];
		theBindingGeneRec.bindingGeneindex = new int[nmergecount];

		int nmergeindex = 0;
		nindex1 = 0;
		nindex2 = 0;
		while ((nindex1 < geneindex1.length)&&(nindex2 < geneindex2.length))
		{
			if (geneindex1[nindex1] == geneindex2[nindex2])
			{
				theBindingGeneRec.bindingValGene[nmergeindex] = Math.max(geneval1[nindex1], geneval2[nindex2]);
				theBindingGeneRec.bindingGeneindex[nmergeindex] = geneindex2[nindex2];
				nindex1++;
				nindex2++;
			}
			else if (geneindex1[nindex1] < geneindex2[nindex2])
			{
				theBindingGeneRec.bindingValGene[nmergeindex] = geneval1[nindex1];
				theBindingGeneRec.bindingGeneindex[nmergeindex] = geneindex1[nindex1];
				nindex1++;
			}
			else if (geneindex2[nindex2] < geneindex1[nindex1])
			{
				theBindingGeneRec.bindingValGene[nmergeindex] = geneval2[nindex2];
				theBindingGeneRec.bindingGeneindex[nmergeindex] = geneindex2[nindex2];
				nindex2++;
			}
			nmergeindex++;
		}

		while (nindex2 < geneindex2.length)
		{
			theBindingGeneRec.bindingValGene[nmergeindex] = geneval2[nindex2];          
			theBindingGeneRec.bindingGeneindex[nmergeindex] = geneindex2[nindex2];
			nmergeindex++;
			nindex2++;
		}

		while (nindex1 < geneindex1.length)
		{
			theBindingGeneRec.bindingValGene[nmergeindex] = geneval1[nindex1];
			theBindingGeneRec.bindingGeneindex[nmergeindex] = geneindex1[nindex1];
			nmergeindex++;
			nindex1++;
		}

		return theBindingGeneRec;
	} 

	/**
	 * Stores into htbinding the mapping of the gene szid to its non-zero TF-gene interactions given by
	 * nonzeroindex and nonzerovals 
	 */
	public void loadBinding(String szid, int[] nonzeroindex, double[] nonzerovals,  int numnonzero)
	{
		String szfull = szid.toUpperCase(Locale.ENGLISH);
		StringTokenizer stIDs = new StringTokenizer(szfull,SZDELIM);

		BindingGeneRec rec = new BindingGeneRec();
		rec.bindingGeneindex = new int[numnonzero];
		rec.bindingValGene = new double[numnonzero];
		for (int i = 0; i < numnonzero; i++)
		{
			rec.bindingGeneindex[i] = nonzeroindex[i];
			rec.bindingValGene[i] = nonzerovals[i];
		}

		while (stIDs.hasMoreTokens())
		{
			//union if gene matches multiple hits in binding file, takes last one seen if multiple non-zero
			String sztoken =stIDs.nextToken();
			BindingGeneRec currec = (BindingGeneRec) htBinding.get(sztoken);

			if (currec != null)
			{
				rec = mergeArrays(rec.bindingValGene, rec.bindingGeneindex,
						currec.bindingValGene, currec.bindingGeneindex);
			}
			htBinding.put(sztoken, rec);

			StringTokenizer stu = new StringTokenizer(sztoken,"_");
			if (stu.countTokens() > 1)
			{
				String szfirsttoken = stu.nextToken();
				currec =(BindingGeneRec) htBinding.get(szfirsttoken);
				if (currec != null)
				{
					rec = mergeArrays(rec.bindingValGene, rec.bindingGeneindex,
							currec.bindingValGene, currec.bindingGeneindex);
				}
				htBinding.put(szfirsttoken, rec);
			}
		}
	}

	/**
	 * Returns an entry of htBinding corresponding to szprobename or szgenename
	 */
	Object getBindingObject(String szgenename, String szprobename)
	{

		StringTokenizer st = new StringTokenizer(szprobename,SZDELIM); 
		String sztoken;
		Object obj = null;

		while ((st.hasMoreTokens()) && (obj == null))
		{
			sztoken = st.nextToken();
			obj = htBinding.get(sztoken);
		}

		if (obj == null)
		{
			st = new StringTokenizer(szgenename,SZDELIM); 
			while ((st.hasMoreTokens()) && (obj == null))
			{
				sztoken = st.nextToken();
				obj = htBinding.get(sztoken);
			}
		}

		return obj;
	}


	//////////////////////////////////////////////////

	/**
	 * Computes the average and standard deviation expression level at each time point
	 */    
	public void computeDataStats(double[][] data, int[][] pma, double[] dsigmas,double[] dmeans)
	{
		double dsum ;
		double dsumsq ;
		int npoints;
		double dsigmaavg = 0;
		for (int ncol = 0; ncol < data[0].length; ncol++)
		{
			dsum = 0;
			dsumsq = 0;
			npoints = 0;
			//need to handle 0 or 1 point error
			for(int nrow = 0; nrow < data.length; nrow++)
			{
				if (pma[nrow][ncol] != 0)
				{
					dsum += data[nrow][ncol];
					dsumsq += Math.pow(data[nrow][ncol],2);
					npoints++;
				}
			}
			dmeans[ncol]= dsum/npoints;
			dsigmas[ncol] = Math.sqrt((dsumsq- Math.pow(dsum,2)/npoints)/(npoints-1));
			if (BDEBUG)
			{
				System.out.print(dsigmas[ncol]+"\t");
			}
			dsigmaavg+= dsigmas[ncol];
			if (BDEBUG)
			{
				System.out.println("###"+dsigmas[ncol]);
			}
		}
		DEFAULTSIGMA = dsigmaavg/data[0].length;
		if (BDEBUG)
		{
			System.out.println("[[[]]]"+DEFAULTSIGMA);
			System.out.println();
		}
	}

	///////////////////////////////////////////////////
	/**
	 * Uses DREM_NaiveBayes to build a classifier which predicts given the set of transcription
	 * factors predicted to regulate a gene whether it would appear filtered or not
	 */    
	public void buildFilteredClassifier()
	{
		Iterator itrgenes = theDataSet.htFiltered.keySet().iterator();
		int numfiltered = theDataSet.htFiltered.size();

		int[][] filteredinput = new int[numfiltered][];
		int[][] filteredinputIndex = new int[numfiltered][];

		int nfilteredhits = 0;
		int[] ALLZEROES = new int[0];

		int nindex = 0;
		while (itrgenes.hasNext())
		{
			//going through filtered genes
			Object geneobj = itrgenes.next();
			String szprobe = (String) theDataSet.htFiltered.get(geneobj);
			String szgeneobj = (String) geneobj;
			//getBindingObject finds the first match for szgeneobj and szprobe
			BindingGeneRec theBindingGeneRec = (BindingGeneRec) getBindingObject(szgeneobj, szprobe);
			if (theBindingGeneRec != null)
			{
				filteredinput[nindex] = sign(theBindingGeneRec.bindingValGene);
				filteredinputIndex[nindex] = theBindingGeneRec.bindingGeneindex;
				nfilteredhits++;
			}
			else
			{
				if (!bfilterbinding)
				{
					//not filtering binding using all zeros insted
					filteredinput[nindex] = ALLZEROES;
					filteredinputIndex[nindex] = ALLZEROES; //really empty
				}
				else
				{
					filteredinput[nindex] = null;
					filteredinputIndex[nindex] = null;
				}
			}
			nindex++;
		}

		int nfinalfilter;
		if (bfilterbinding)
		{
			nfinalfilter = nfilteredhits;
		}
		else
		{
			nfinalfilter = numfiltered;
		}
		ntotalcombined = bindingSignGene.length + nfinalfilter;

		if (BDEBUG)
		{
			System.out.println(bindingSignGene.length +" $$$$$$$$$$ "+nfinalfilter);
		}
		//making new array with filtered and non-filtered
		//need to do this TF wise
		int[][] combinedbinding = new int[ntotalcombined][];
		int[][] combinedbindingIndex = new int[ntotalcombined][];
		int[] filteredlabel = new int[ntotalcombined];
		double[] trainweight = new double[ntotalcombined];
		for (int ni = 0; ni < bindingSignGene.length; ni++)
		{
			filteredlabel[ni] = 1;
			combinedbinding[ni] = bindingSignGene[ni];
			combinedbindingIndex[ni] = bindingGeneindex[ni];	     
		}

		for (int ni = bindingSignGene.length; ni < filteredlabel.length; ni++)
		{
			filteredlabel[ni] = 0;
		}

		int nfilteredindex = 0;
		int ntotalindex = bindingSignGene.length;
		while ((nfilteredindex < filteredinput.length)&&(ntotalindex < combinedbinding.length))
		{
			if (filteredinput[nfilteredindex] != null)
			{
				combinedbinding[ntotalindex] = filteredinput[nfilteredindex];
				combinedbindingIndex[ntotalindex] = filteredinputIndex[nfilteredindex];	      
				ntotalindex++;
			}
			nfilteredindex++;
		}

		for (int ni = 0; ni < trainweight.length; ni++)
		{
			trainweight[ni] = 1;
		}

		//combinedbinding input values
		//filtered label y-labels
		if (dBindingSigns.length > 0)
		{
			filteredClassifier = new DREM_NaiveBayes(combinedbinding,combinedbindingIndex,numbits,
					filteredlabel, 2);
		}
		else
		{
			filteredClassifier = null;
		}

		if (BDEBUG)
		{
			System.out.println(filteredClassifier);
		}
	}


	////////////////////////////////////////////////////
	/**
	 * Displays the current temporary DREM model, from which a search
	 * is performed for an improved model
	 */
	public void displayTempMap() throws Exception
	{
		DREM_IO.bdisplaycurrent = false;
		final DREM_Timeiohmm fthishmm = this;
		if (BDEBUG)
		{
			System.out.println("before hmmgui");
		}

		final Treenode treecopy =(Treenode) treeptr.clone();

		computeOrders(treecopy);

		viterbi(treecopy);
		if (bindingSignGene !=null)
		{
			computeStats(treecopy,treecopy);
			computeminparentlevel(treecopy);
		}

		((DREM_GoAnnotations) theDataSet.tga).buildRecDREM(treecopy,theDataSet.genenames);

		traverse(bestTree, 0,false);

		try
		{         
			javax.swing.SwingUtilities.invokeAndWait( 
					new Runnable() 
					{
						public void run() 
						{
							progressDREMGui = new DREMGui(fthishmm,treecopy, brealXaxisDEF,dYaxisDEF,dXaxisDEF,
									nKeyInputTypeDEF,dKeyInputXDEF,dpercentDEF,"(Not Final Model)",dnodekDEF);
							edu.umd.cs.piccolo.PCanvas.CURRENT_ZCANVAS = null;
							progressDREMGui.setLocation(15,40);
							progressDREMGui.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);

							progressDREMGui.setVisible(true);
							progressDREMGui.addWindowListener(new WindowAdapter() 
							{
								public void windowClosing(WindowEvent we) 
								{
									progressDREMGui.closeWindows();
								}
							});        
						}
					});
		}
		catch (InterruptedException iex)
		{
			iex.printStackTrace(System.out);
		}
		catch (java.lang.reflect.InvocationTargetException itex)
		{
			itex.printStackTrace(System.out);
		}   
		currentButton.setEnabled(true);

	}


	//////////////////////////////////////////////////
	/**
	 * Executes the first phase of the structure search
	 * Considers adding and deleting paths, does not consider 
	 * merges or delays
	 */
	public void searchstage1() throws Exception
	{
		String szDIFF = "";
		String szpercent = "";
		int numintervals = traindata[0].length-1;
		int[] path = new int[numintervals];

		traverse(treeptr, 0,false);


		if (BVITERBI)
		{
			dbestlog =trainhmmV(treeptr,true);
		}
		else
		{
			dbestlog = trainhmm(treeptr,true);
		}
		traverse(treeptr, 0,false);

		//////////////////////////////

		double dprevbestlog;
		dprevouterbestlog=Double.NEGATIVE_INFINITY;

		if (!bmodeldiff)
		{
			if (BVITERBI)
			{
				dbestlog = testhmmV(treeptr);
			}
			else
			{
				dbestlog = testhmm(treeptr);
			}
		}
		dbesttrainlike = Double.NEGATIVE_INFINITY;
		boolean bendsearchlocal; 

		do 
		{
			if (BDEBUG)
			{
				System.out.println("&&&& "+numtotalPath);
			}
			//outer loop handles the re-splitting of the data
			//and whether we should try to add and then delete nodes

			if (dprevouterbestlog != Double.NEGATIVE_INFINITY)
			{
				if (statusLabel != null)
				{
					statusLabel.setText(" Current score: "+nf2.format(dbestlog));
					statusLabel15.setText(" Current score improvement: "+szDIFF+" ("+szpercent+")");
					statusLabel3.setText(" Next score: "+nf2.format(dbestlog));
					statusLabel2.setText(" Next score improvement: 0 (0.000%)");
					statuscountLabel.setText(" Number of paths in model so far: "+numtotalPath);
				}
				String szepsilonpercent = nf3.format(100*EPSILON) + "%";
			}
			dprevouterbestlog = dbestlog;

			traverse(treeptr, 0,false);

			double dtemplog;

			bestTree = treeptr;

			if (BDEBUG)
			{
				System.out.println("trying to add");
			}
			dprevbestlog = dbestlog;

			traverse(bestTree, 0,false);
			traverseandadd(path,treeptr,treeptr);
			if (BDEBUG)
			{
				System.out.println("****** adding best log likelihood "+dbestlog);
			}
			traverse(bestTree, 0,false);

			if (BDEBUG)
			{ 
				System.out.println("after traverse");
			}
			treeptr = bestTree;

			if (BDEBUG)
			{ 
				System.out.println("$$$$$$$$$$$$$$$terminate "+dbestlog+"\t"+dprevbestlog+"\t"+
						(dbestlog-dprevbestlog)/Math.abs(dbestlog)+"\t"+(dbestlog-dprevbestlog));
			}

			if (BDEBUG)
			{
				System.out.println("YYY\tadd\t"+dbesttrainlike+"\t"+dbestlog);
			}

			if (dbestlog > dprevbestlog)
			{
				numtotalPath++;
			}

			while ((dbestlog - dprevbestlog)>0)
			{
				if (BDEBUG)
				{
					System.out.println("trying to delete");
				}

				dprevbestlog = dbestlog; 

				traverseanddelete(path,treeptr,treeptr,false,0,0);

				if (BDEBUG)
				{
					System.out.println("****** delete  best log likelihood "+dbestlog);
				}

				traverse(bestTree, 0,false);

				treeptr = bestTree;

				if (BDEBUG)
				{
					System.out.println("YYY\tdelete\t"+dbesttrainlike+"\t"+dbestlog);
				}

				if (dbestlog > dprevbestlog)
				{
					numtotalPath--;
				}
			}

			if (BDEBUG)
			{
				System.out.println("after delete "+dbestlog+"\t"+dprevouterbestlog);
				System.out.println("*************terminate "+(dbestlog-dprevouterbestlog)+"\t"+
						(dbestlog-dprevouterbestlog));
				System.out.println("ZZZ\t"+dbesttrainlike+"\t"+dbestlog);
			}

			szpercent = nf3.format(100*(dbestlog-dprevouterbestlog)/Math.abs(dbestlog)) + "%";
			szDIFF = nf3.format(dbestlog-dprevouterbestlog);

			bendsearchlocal = DREM_IO.bendsearch;

			if (DREM_IO.bdisplaycurrent)
			{
				displayTempMap();  
			}

			if (BDEBUG)
			{
				System.out.println("$"+EPSILON+" "+dbestlog+" "+dprevouterbestlog+" "+EPSILONDIFF);
			}

		}
		while ((((dbestlog-EPSILON*Math.abs(dbestlog)-dprevouterbestlog)>EPSILONDIFF)
				&&(!bendsearchlocal)&&!bmodeldiff) ||
				(bmodeldiff && (dbestlog>dprevouterbestlog)&&(!bendsearchlocal)));

		if (statusLabel15 != null)
		{
			statusLabel15.setText(" Final Main Score Improvement: "+szDIFF+" ("+szpercent +")");
			statusLabel.setText(" Final Main Score: "+nf2.format(dbestlog));
			statuscountLabel.setText(" Number of paths in model so far: "+numtotalPath);
		}
	}


	////////////////////////////////////////////////////////////////////////

	/**
	 * Assigns to the nminparentlevel in all nodes the level of the most immediate ancestor
	 * with two or more children
	 */
	public void computeminparentlevel(Treenode ptr)
	{
		if (ptr != null)
		{
			if ((ptr.parentptrA == null)||(ptr.parentptrA.length == 1))
			{ 
				if (ptr.parent == null)
				{
					ptr.nminparentlevel = ptr.ndepth;
					ptr.minparentptr = ptr.parent;
				}
				else  
				{
					if (ptr.parent.numchildren >= 2)
					{
						ptr.nminparentlevel = ptr.ndepth;
						ptr.minparentptr = ptr.parent;
					}
					else
					{
						ptr.nminparentlevel = ptr.parent.nminparentlevel;
						ptr.minparentptr = ptr.parent.minparentptr;
					}
				}
			}
			else
			{	
				Treenode[] tempptrA = new Treenode[ptr.parentptrA.length];
				int ndepth = ptr.ndepth;
				for (int nindex =0 ; nindex < tempptrA.length; nindex++)
				{
					tempptrA[nindex] = ptr.parentptrA[nindex];
				}

				boolean bfoundit = false;
				while (!bfoundit)
				{
					int nj = 1;
					boolean ballsame = true;
					while ((ballsame)&&nj<tempptrA.length)
					{
						if (tempptrA[nj] != tempptrA[nj-1])
						{
							ballsame = false;
						}
						else 
						{
							nj++;
						}
					}

					ndepth--;
					bfoundit = ballsame;
					for (int nindex =0 ; nindex < tempptrA.length; nindex++)
					{
						tempptrA[nindex] = tempptrA[nindex].parent;
					}
				}

				ptr.nminparentlevel = ndepth;
				ptr.minparentptr = tempptrA[0];
			}

			for (int nindex = 0; nindex < ptr.nextptr.length; nindex++)
			{
				computeminparentlevel(ptr.nextptr[nindex]);
			} 
		}
	}


	/**
	 * Responsible with helper functions for merging paths. Paths must share a common most
	 * recent split and no splits after a merge are allowed
	 */
	public void traverseandmerge(int[] path) throws Exception
	{

		computeNumLeaves(bestTree);
		if (BDEBUG)
		{
			System.out.println("in traverse and merge path.length = "+path.length);
		}

		traverse(bestTree,0,false);

		for (int nparentlevel = path.length-2; nparentlevel >=0; nparentlevel--)
		{
			for (int nmergingdepth = path.length; nmergingdepth > nparentlevel+1; nmergingdepth--)
			{
				//nmergingdepth is the depth of the merging
				//parent level is the common ancestor
				Treenode origbestTree;
				do
				{
					dprevouterbestlog = dbestlog; 
					dbestlog = Double.NEGATIVE_INFINITY;

					if (BDEBUG)
					{
						System.out.println("calling traverseandmergehelp "+nmergingdepth+" "+nparentlevel);
					}
					origbestTree = bestTree;
					traverseandmergehelp(path,0,nmergingdepth, nparentlevel,origbestTree,origbestTree);
					if (dbestlog == Double.NEGATIVE_INFINITY)
						dbestlog = dprevouterbestlog;
				}
				while (origbestTree != bestTree);
			}
		}

	}

	/**
	 * Helper function for traverseandmerge used to recursively search for places to merge nodes
	 */
	private void traverseandmergehelp(int[] path, int nlevel,int nmergingdepth, int nparentlevel, Treenode ptr,
			Treenode origroot) throws Exception
			{

		if (nlevel <= nmergingdepth)
		{
			for (int nindex = 0; nindex < ptr.numchildren; nindex++)
			{
				path[nlevel] = nindex;
				traverseandmergehelp(path,nlevel+1,nmergingdepth,
						nparentlevel,ptr.nextptr[nindex], origroot);
			}
		}

		if (nlevel == nmergingdepth)
		{
			if (BDEBUG)
			{
				System.out.println("nparentlevel = "+nparentlevel+" nmergindepth = "+nmergingdepth+
						" nlevel = "+nlevel+" mean = "+ptr.dmean);
			}
			mergenode(path,nmergingdepth, nparentlevel,origroot);
		} 
			}


	/**
	 * Computes the number of paths to leaves accessible from each node in the 
	 * tree pointed to by ptr
	 */
	public int computeNumLeaves(Treenode ptr)
	{
		if (ptr.numchildren == 0)
		{
			return 1;    
		}
		else
		{
			int nsum = 0;
			for (int nindex = 0; nindex < ptr.numchildren; nindex++)
			{
				nsum += computeNumLeaves(ptr.nextptr[nindex]);
			}

			if (BDEBUG)
			{
				System.out.println("ptr.dmean = "+ptr.dmean+" nsum = "+nsum);
			}

			ptr.numdownleaves = nsum;

			return nsum;
		}
	}

	//////////////////////////////////////////////////////////////////////////////////////////
	/**
	 * Returns true iff all nodes at nleveltogo only have one path to the leaves
	 * Stores in hsLevel all nodes at nleveltogo is 0
	 */
	public boolean findnodesAtLevel(Treenode treeptr, int nleveltogo, HashSet hsLevel)
	{
		if (nleveltogo == 0)
		{
			hsLevel.add(treeptr);
			if (treeptr.numdownleaves >= 2) 
				return false;
			else
				return true;
		}
		else
		{
			boolean ballsingles  = true;
			for (int nindex = 0; nindex < treeptr.numchildren; nindex++)
			{
				ballsingles= (ballsingles) && (findnodesAtLevel(treeptr.nextptr[nindex],nleveltogo-1,hsLevel));
			}
			return ballsingles;
		}
	}

	//////////////////////////////////////////////////////////////////////////////////////////////////
	/**
	 * Helper function of traverseandmergehelp responsible for finding valid merges and 
	 * evaluating them that are along path1 and at depth nmergingdepth
	 */
	private void mergenode(int[] path1, int nmergingdepth, int nparentlevel, Treenode origroot) throws Exception
	{
		double dlog;

		if (BDEBUG)
		{
			System.out.println("nmerginglevel = "+nmergingdepth);
			System.out.println("nparentlevel = "+nparentlevel);
			System.out.print("Path: ");
			for (int i = 0; i < nparentlevel; i++)
			{
				System.out.print(path1[i]+" ");
			}
			System.out.println();
		}

		Treenode tempptr = origroot;
		for (int nindex  = 0; nindex <nparentlevel; nindex++)
		{
			tempptr = tempptr.nextptr[path1[nindex]];
		}

		if ((tempptr.numchildren >= 2))
		{
			HashSet hsLevel = new HashSet();
			boolean ballsingles = findnodesAtLevel(tempptr,(nmergingdepth-nparentlevel),hsLevel);
			//ballsingles - is true if no more splits among nodes being merged	   

			int numatLevel = hsLevel.size();
			//only binary for now
			if (BDEBUG)
			{
				System.out.println("numatLevel = "+numatLevel);
			}

			if ((numatLevel == tempptr.numchildren)&&(ballsingles))
			{
				//not allowing any downstream splits from parent node
				//either before the level being merged, should have been accepted before
				//among merged nodes, not allowing any resplits
				//require ballsingles - to be true
				//also cannot have any splits 
				//assume binary for now

				if (BDEBUG)
				{
					System.out.println("before clone");
				}

				traverse(origroot,0,false);
				Treenode treeroot = (Treenode) origroot.clone();
				if (BDEBUG)
				{
					System.out.println("after clone");
				}
				traverse(treeroot,0,false);

				//traversing to the parent
				tempptr = treeroot;
				for (int nindex  = 0; nindex <nparentlevel; nindex++)
				{
					tempptr = tempptr.nextptr[path1[nindex]];
				}

				Treenode parentptr = tempptr;

				//travering to the nodes for merging
				Treenode splitptr = null;
				for (int nindex = nparentlevel; nindex < nmergingdepth; nindex++) 
				{
					tempptr = tempptr.nextptr[path1[nindex]];	     
				}

				if (BDEBUG)
				{
					System.out.println(nmergingdepth+" parent mean "+parentptr.dmean+" tempptr.dmean "
							+tempptr.dmean+" split mean ");
				}

				boolean bvalidmerge = false;
				ArrayList newparentptrA = null;
				for (int notherchild = 0; notherchild < parentptr.numchildren; notherchild++)
				{
					if (BDEBUG)
					{
						System.out.println(notherchild+"\t"+bvalidmerge);
					}

					if (notherchild != path1[nparentlevel])
					{  
						splitptr = parentptr.nextptr[notherchild];

						for (int nindex = nparentlevel+1; nindex < nmergingdepth; nindex++) 
						{
							splitptr = splitptr.nextptr[0];
						}

						if (splitptr != tempptr)
						{
							//its possible two nodes have already been merged on this level
							//but there are two paths two them, we do not need to adjust parent
							//pointers if just taking a different path to the main merging node 
							if (!bvalidmerge)
							{
								//we haven't tried merging yet, just do this once to load the merging node 
								//parents into the list
								newparentptrA = new ArrayList();
								if (tempptr.parentptrA == null)
								{
									//just single parent add to parent list
									newparentptrA.add(tempptr.parent);
									if (BDEBUG)
									{
										System.out.println("A "+tempptr.dmean+" now links to "+tempptr.parent.dmean);
									}
								}
								else
								{ 
									for (int i = 0; i < tempptr.parentptrA.length; i++)
									{
										//add all parent links to merging state
										newparentptrA.add(tempptr.parentptrA[i]);
										if (BDEBUG)
										{
											System.out.println("B "+tempptr.dmean+" now links to "+tempptr.parentptrA[i].dmean);
										}
									}
								}

								bvalidmerge = true;

								if (BDEBUG)
								{
									if (tempptr.nextptr[0]!=null) 
										System.out.println("mean is "+tempptr.nextptr[0].dmean);
									System.out.println("VALID MERGE");
								}
							}

							if (BDEBUG)
							{
								System.out.println(splitptr.dmean+" ");
							}

							if (splitptr.parentptrA == null)
							{
								splitptr.parent.nextptr[0] = tempptr;
								newparentptrA.add(splitptr.parent);
								if (BDEBUG)
								{
									System.out.println("C "+tempptr.dmean+" now links to "+splitptr.parent.dmean);
								}
							}
							else
							{
								for (int i = 0; i < splitptr.parentptrA.length; i++)
								{
									splitptr.parentptrA[i].nextptr[0] = tempptr;
									newparentptrA.add(splitptr.parentptrA[i]);
									if (BDEBUG)
									{
										System.out.println("D "+tempptr.dmean+" now links to "+splitptr.parentptrA[i].dmean);
									}
								}
							}

							if (tempptr.nextptr[0] != null)
							{
								tempptr.nextptr[0].parent = tempptr;
								tempptr.nextptr[0].parentptrA = null;
							}

							if (parentptr == splitptr.parent)
							{
								parentptr.numchildren--;
							}
						}
					}
				}

				if (bvalidmerge)
				{
					//copying over new parent nodes
					int nsize = newparentptrA.size();
					tempptr.parentptrA = new Treenode[nsize];
					for (int i = 0; i < nsize; i++)
					{
						tempptr.parentptrA[i] = (Treenode) newparentptrA.get(i);
					}
					if (BDEBUG)
					{
						System.out.println("MERGED TREE");
					}
					traverse(treeroot,0,false);

					if (BVITERBI)
					{
						dlog = trainhmmV(treeroot,true);
					}
					else
					{
						dlog =trainhmm(treeroot,true);
					}

					if (BDEBUG)
					{
						System.out.println("after training:");
					}

					traverse(treeroot,0,false);

					if (!bmodeldiff)
					{
						if (BVITERBI)
						{
							dlog = testhmmV(treeroot);
						}
						else
						{
							dlog =testhmm(treeroot);
						}
					}

					double dimprovemerge = (dlog-dprevouterbestlog)/Math.abs(dlog);
					if (BDEBUG)
					{		
						System.out.println("merge* "+dimprovemerge+"\t"+MERGEMIN+"\t"
								+dlog+"\t"+dprevouterbestlog+"\t"+dbestlog);
					}

					if (((!bmodeldiff)&&
							(((dlog+MERGEMIN*Math.abs(dlog)-dprevouterbestlog)>=MERGEMINDIFF)&&(dlog>dbestlog)))||
							(bmodeldiff&&(dlog>=dprevouterbestlog)&&(dlog>dbestlog)))
					{
						bestTree = treeroot;
						computeNumLeaves(bestTree);
						dbestlog = dlog;
						if (BDEBUG)
						{
							System.out.println("in merge "+dbestlog+" best tree is");
							traverse(treeroot,0,false);
						}

						double dimprove = (dbestlog-dprevouterbestlog)/Math.abs(dbestlog);
						String szimprove = nf3.format(100*dimprove) + "%";
						String szimproveDIFF = nf3.format(dbestlog-dprevouterbestlog);
						if (statusLabel2 != null)
							statusLabel2.setText(" Merge change: "+szimproveDIFF+" ("+szimprove+")"+
									" (Score "+nf2.format(dbestlog)+")");
					}
				}
			}
		}
	}

	///////////////////////////////////////////////////////////////////////
	/**
	 * A record containing a node id, mean, and standard deviation
	 * used for ordering the children node
	 */
	static class OrderRec
	{
		int nid;
		double dmean;
		double dsigma;

		OrderRec(double dmean, double dsigma, int nid)
		{
			this.nid = nid;
			this.dmean = dmean;
			this.dsigma = dsigma;
		}
	}

	/**
	 * Comparator for OrderRec sorts first on dmean, then dsigma, then id all
	 * in increasing order.
	 */
	static class OrderRecCompare implements Comparator
	{
		public int compare(Object o1, Object o2)
		{
			OrderRec or1 = (OrderRec) o1;
			OrderRec or2 = (OrderRec) o2;

			if (or1.dmean < or2.dmean)
				return -1;
			else if (or1.dmean > or2.dmean)
				return 1;
			else if (or1.dsigma < or2.dsigma)
				return -1;
			else if (or1.dsigma > or2.dsigma)
				return 1;
			else if (or1.nid < or2.nid)
				return -1;
			else if (or1.nid > or2.nid)
				return 1;
			else 
				return 0;
		}
	}

	///////////////////////////////////////////////////////////////////////
	/**
	 * Sorts the children of each node in treeptr and descendants
	 * based on OrderRecCompare
	 */
	public void computeOrders(Treenode treeptr)
	{
		if (treeptr.numchildren >=1)
		{
			if (treeptr.numchildren >= 2)
			{
				OrderRec[] theOrderRecs = new OrderRec[treeptr.numchildren];
				for (int nindex = 0; nindex < treeptr.numchildren; nindex++)
				{
					theOrderRecs[nindex] = new OrderRec(treeptr.nextptr[nindex].dmean, 
							treeptr.nextptr[nindex].dsigma, nindex);
				}
				Arrays.sort(theOrderRecs, new OrderRecCompare());
				treeptr.orderA = new int[theOrderRecs.length];
				treeptr.revOrderA = new int[theOrderRecs.length];

				for (int nindex = 0; nindex < treeptr.orderA.length; nindex++)
				{
					treeptr.orderA[nindex] = theOrderRecs[nindex].nid;
					treeptr.revOrderA[theOrderRecs[nindex].nid] = nindex;
				}
			}


			for (int nindex = 0; nindex < treeptr.numchildren; nindex++)
			{
				computeOrders(treeptr.nextptr[nindex]);
			}
		}
	}


	///////////////////////////////////////////////////////////////////////


	/**
	 * Builds the initial tree which is just a single chain with mean and standard
	 * deviation for each node being the global mean and standard deviation at 
	 * the time point
	 */
	public void buildEmptyTree(int ncurrdepth,
			Treenode currnode, double[] dsigma, double[] dmean)
	{
		currnode.dsigma = Math.max(dsigma[ncurrdepth],MINSIGMA);
		if (ncurrdepth+1 < dsigma.length)
		{
			currnode.numchildren = 1;
			currnode.nextptr[0] = new Treenode(currnode);
			currnode.nextptr[0].dmean = dmean[ncurrdepth];
			currnode.nextptr[0].parent = currnode;
			if (BREGDREM)
			{
				currnode.ptrans[0] = 1;
			}
			buildEmptyTree(ncurrdepth+1,currnode.nextptr[0], dsigma,dmean);
		}
	}

	/**
	 * Initializes the model parameters to those in szinitfileval
	 */    
	public void readInitTree(Treenode currnode) throws IOException
	{

		BufferedReader brinit = new BufferedReader(new FileReader(szinitfileval));
		String szfirstline = brinit.readLine();
		boolean bmergedtree = false;

		if (szfirstline.equalsIgnoreCase("MERGED TREE"))
		{
			szfirstline = brinit.readLine();
			bmergedtree = true;
		}

		if (bmergedtree)
		{
			if (ninitsearchval == 1)
			{
				throw new IllegalArgumentException("Saved model has merges, but option allowed to start "+
						"search from model is selected.\n"+
						"Search may not begin from model with merges.  "+
				"Consider option to use the saved model as is.");
			}

			if (!ballowmergeval)
			{
				throw new IllegalArgumentException("Saved model has merges, but merges are not selected to be allowed.");
			}
		}

		StringTokenizer stfirstline = new StringTokenizer(szfirstline,"\t");
		if (stfirstline.countTokens() < 2)
		{
			throw new IllegalArgumentException("Not a valid saved model file!");
		}

		String szfirsttoken = stfirstline.nextToken();
		int numcoeff = Integer.parseInt(stfirstline.nextToken());
		if ((numcoeff != tfNames.length)&&(!szfirsttoken.startsWith("REGDREM")))
		{
			throw new IllegalArgumentException("Number of coefficients ("+numcoeff+") in the saved model does not "+                               
					"match the number of coefficients in the static input file ("+tfNames.length+")");
		}

		if (!(szfirsttoken.startsWith("REGDREM"))&&(BREGDREM))
		{
			throw new IllegalArgumentException("'Use Transcription Factor-gene Interaction Data to Build Model' under 'Main Search Options' is unselected, "+
			"but saved model was built using the data.");
		}
		else if ((szfirsttoken.startsWith("REGDREM"))&&(!BREGDREM))
		{
			throw new IllegalArgumentException("'Use Transcription Factor-gene Interaction Data to Build Model' under 'Main Search Options' is selected, "+
			"but saved model was built without the data.");
		}

		int ndepth = readInitTreeHelp(currnode, brinit,numcoeff,bmergedtree);

		if (ndepth != numcols)
		{
			throw new IllegalArgumentException("Number of time points in the initial model ("+
					ndepth+") does not match number of time points in the data ("+numcols+")");
		}

		String szLine = brinit.readLine();
		if ((szLine != null) && (szLine.equals("COLORS")))
		{
			readColors(brinit);
		}

		if (bmergedtree)
		{
			HashMap createdNodes = new HashMap();
			mergeDuplicates(currnode,createdNodes);
		}
	}

	/**
	 * Helper function for reading the model parameters from a file
	 */
	private int readInitTreeHelp(Treenode currnode, BufferedReader brinit, 
			int numcoeff,boolean bmergedtree) throws IOException
			{
		String szLine = brinit.readLine();
		StringTokenizer st = new StringTokenizer(szLine,"\t");
		int nID;

		if (bmergedtree)
		{
			nID = Integer.parseInt(st.nextToken());
			currnode.nID = nID;
		}

		currnode.dmean = Double.parseDouble(st.nextToken());
		currnode.dsigma = Double.parseDouble(st.nextToken());
		currnode.numchildren = Integer.parseInt(st.nextToken());

		int numcoefftoread;

		if (currnode.numchildren <= 1)
		{
			numcoefftoread = 0;
			currnode.binit = false;
		}
		else if (currnode.numchildren == 2)
		{
			numcoefftoread = numcoeff + 1;
			currnode.binit = true;
			numtotalPath++;
		}
		else
		{
			if (currnode.numchildren > nmaxchild)
			{
				throw new IllegalArgumentException("Number of paths out of a node in the initial model file ("+
						currnode.numchildren+") exceeds the maximum allowed "+nmaxchild);

			}
			numtotalPath += currnode.numchildren-1;
			numcoefftoread = (numcoeff+1)*(currnode.numchildren-1);
			currnode.binit = true;
		}

		if ((numcoefftoread==0)&&(BREGDREM))
		{
			currnode.ptrans[0] = 1;
		}
		else if (numcoefftoread > 0)
		{
			if (BREGDREM)
			{
				for (int nindex = 0; nindex < currnode.numchildren; nindex++)
				{
					currnode.ptrans[nindex] = Double.parseDouble(brinit.readLine());
				}
			}
			else
			{
				double[] dcoeff = new double[numcoefftoread];
				String szcoeff;
				StringTokenizer stcoeff;
				for (int ncoeffindex = 0; ncoeffindex < numcoefftoread; ncoeffindex++)
				{
					szcoeff = brinit.readLine();
					stcoeff = new StringTokenizer(szcoeff,"\t");
					String szcoeffname = stcoeff.nextToken();
					int nmodrow = ncoeffindex % (tfNames.length + 1); 
					if ((nmodrow!= 0)&&(!szcoeffname.equals(tfNames[nmodrow-1])))
					{
						throw new IllegalArgumentException(ncoeffindex+" "+(tfNames.length+1)+" "+
								"Transcription factor "+szcoeffname +" in saved model  "+
								"does not match transcription factor "+tfNames[nmodrow-1]+".");  
					}
					dcoeff[ncoeffindex] = Double.parseDouble(stcoeff.nextToken());
				}

				currnode.tranC = new DREM_FastLogistic2(theInstancesIndex[currnode.numchildren-2],
						theInstances[currnode.numchildren-2],
						theInstancesTFIndex[currnode.numchildren-2],
						theInstancesTF[currnode.numchildren-2], 
						ylabels[currnode.numchildren-2],
						currnode.numchildren,dcoeff,numbits);
			}
		}

		if (currnode.numchildren == 0)
		{
			return 1;
		}
		else
		{
			int ndepth;
			currnode.nextptr[0] = new Treenode(currnode);
			currnode.nextptr[0].parent = currnode;
			ndepth = readInitTreeHelp(currnode.nextptr[0],brinit,numcoeff,bmergedtree);
			for (int nchild = 1; nchild < currnode.numchildren; nchild++)
			{
				currnode.nextptr[nchild] = new Treenode(currnode);
				currnode.nextptr[nchild].parent = currnode;
				if (ndepth != readInitTreeHelp(currnode.nextptr[nchild],brinit,numcoeff,bmergedtree))
				{
					throw new IllegalArgumentException("Invalid saved model file.  Two paths are not of the same length.");
				}
			}
			return (ndepth +1);
		}
			}


	/**
	 * When reading a model with merges from a file some of the states can be listed more than once
	 * this procedures makes it so that in the internal representation the state only exists once
	 */
	private void mergeDuplicates(Treenode currnode,HashMap createdNodes)
	{
		if (currnode != null)
		{
			createdNodes.put(new Integer(currnode.nID),currnode);
			for (int nchild = 0; nchild < currnode.numchildren; nchild++)
			{
				Treenode tnode = (Treenode) createdNodes.get(new Integer(currnode.nextptr[nchild].nID));
				if (tnode == null)
				{
					mergeDuplicates(currnode.nextptr[nchild],createdNodes);
				}
				else
				{
					currnode.nextptr[nchild] = tnode;
					if (tnode.parentptrA == null)
					{
						tnode.parentptrA = new Treenode[2];
						tnode.parentptrA[0] = tnode.parent;
						tnode.parentptrA[1] = currnode;
					}
					else
					{
						Treenode[] temptr = tnode.parentptrA;
						tnode.parentptrA = new Treenode[tnode.parentptrA.length];
						for (int nindex = 0; nindex < temptr.length; nindex++)
						{
							tnode.parentptrA[nindex] = temptr[nindex];
						}
						tnode.parentptrA[temptr.length-1] = currnode;
					}
				}
			}
		}
	}

	/**
	 * Loads saved colors for a DREM map
	 */
	public void readColors(BufferedReader brinit) throws IOException
	{
		savedColors = new ArrayList();
		String szLine;
		while ((szLine=brinit.readLine()) != null)
		{
			if (!szLine.trim().equals(""))
			{
				StringTokenizer st = new StringTokenizer(szLine,"\t");
				if (st.countTokens() != 4)
				{
					throw new IllegalArgumentException("A color line has "+st.countTokens()+" values expecting 4");
				}
				else
				{
					float f1,f2,f3,f4;
					f1 = Float.parseFloat(st.nextToken());
					f2 = Float.parseFloat(st.nextToken());
					f3 = Float.parseFloat(st.nextToken());
					f4 = Float.parseFloat(st.nextToken());
					savedColors.add(new Color(f1,f2,f3,f4));
				}
			}
		}
	}


	//////////////////////////////////////////////////////////////////////
	/**
	 * Helper function to traverseanddelayhelp that manipulates a split so that it occurs
	 * one time point later
	 */
	private Treenode delaysplit(int[] path, int ndesiredlevel, Treenode root,int nchildcombine)
	{
		Treenode treeroot = (Treenode) root.clone();
		Treenode ptr = treeroot;

		for(int nindex = 0; nindex < ndesiredlevel-1; nindex++)
		{
			ptr = ptr.nextptr[path[nindex]];
		}
		Treenode removeptr = ptr.nextptr[path[ndesiredlevel-1]];

		if (!removeptr.bchange)
		{
			return null;
		}
		else
		{
			for (int nj = path[ndesiredlevel-1]+1; nj < ptr.numchildren; nj++)   
			{
				ptr.nextptr[nj-1] = ptr.nextptr[nj];
			}

			ptr.numchildren--; 
			ptr.binit = false;

			if (BREGDREM)
			{
				double dsum = 0; 
				for (int nindex = 0; nindex < ptr.numchildren; nindex++)
				{
					dsum += ptr.ptrans[nindex]; 
				}

				for (int nindex = 0; nindex < ptr.numchildren; nindex++)
				{
					ptr.ptrans[nindex]/=dsum; 
				}
			}

			int ncurrnumchild = ptr.nextptr[nchildcombine].numchildren;

			Treenode combineptr = ptr.nextptr[nchildcombine];
			for (int nindex = 0; nindex < removeptr.numchildren ; nindex++)
			{
				combineptr.nextptr[nindex+ncurrnumchild]= removeptr.nextptr[nindex];
				removeptr.nextptr[nindex].parent = combineptr;
			}


			if (BREGDREM)
			{
				for (int nindex = 0; nindex < combineptr.numchildren; nindex++)
				{    
					combineptr.ptrans[nindex] = combineptr.ptrans[nindex]*0.5; 
				} 

				for (int nindex = 0; nindex < removeptr.numchildren; nindex++)
				{ 
					combineptr.ptrans[nindex+combineptr.numchildren] = removeptr.ptrans[nindex]*0.5;
				}
			}
			combineptr.numchildren += removeptr.numchildren;
			combineptr.binit = false;
		}
		return treeroot;
	}

	//////////////////////////////////////////////////////////////////////
	/**
	 * Deletes a child from the specified path on a cloned version of the root node
	 */
	public  Treenode deletepath(int[] path, int ndesiredlevel, Treenode root)
	{
		//used to be call deleteleaf
		Treenode treeroot = (Treenode) root.clone();
		Treenode ptr = treeroot;

		for(int nindex = 0; nindex < ndesiredlevel-1; nindex++)
		{
			ptr = ptr.nextptr[path[nindex]];
		}

		if (!ptr.nextptr[path[ndesiredlevel-1]].bchange)
		{
			return null;
		}
		else
		{
			for (int nj = path[ndesiredlevel-1]+1; nj < ptr.numchildren; nj++)   
			{
				ptr.nextptr[nj-1] = ptr.nextptr[nj];
			}

			ptr.numchildren--; 
			ptr.binit = false;

			if (BREGDREM)
			{
				double dsum = 0; 
				for (int nindex = 0; nindex < ptr.numchildren; nindex++)
				{
					dsum += ptr.ptrans[nindex]; 
				}

				for (int nindex = 0; nindex < ptr.numchildren; nindex++)
				{
					ptr.ptrans[nindex]/=dsum; 
				}
			}
		}
		return treeroot;
	} 

	////////////////////////////////////////////////////////////////////////
	/**
	 * With its helper functions adds the specified path to the model
	 */
	public void traverseandadd(int[] path,Treenode ptr,Treenode origroot) throws Exception
	{
		traverseandaddhelp(path,0,ptr,origroot);
	}
	//////////////////////////////////////////////////////////////////////


	/**
	 * Helper function for traverseandaddhelp recursively tries to add paths 
	 * to the model, and then evaluates the improvement
	 */
	private void traverseandaddhelp(int[] path, int nlevel,
			Treenode ptr, Treenode origroot) throws Exception
			{
		//searches for the best place to add a node at current level
		//tries adding path at current level       
		if ((nlevel < path.length)&&(!DREM_IO.bendsearch))
		{
			for (int nindex = 0; nindex < ptr.numchildren; nindex++)
			{
				path[nlevel] = nindex;
				traverseandaddhelp(path,nlevel+1,ptr.nextptr[nindex],origroot);
			}
		}

		if (((ptr!=null) &&((ptr.numchildren < nmaxchild)))&&(!DREM_IO.bendsearch))
		{
			if (DREM_IO.bdisplaycurrent)
			{
				displayTempMap();
			}

			if (nlevel < path.length)
			{
				//mark to try to split
				if (BDEBUG)
				{
					for (int ni =0; ni < nlevel; ni++)
					{
						System.out.print(path[ni]+"\t");
					}
				}

				Treenode splittree = splitnode(path,nlevel,origroot);  

				if (BDEBUG)
				{
					System.out.println(path[nlevel]);
				}

				if (splittree != null)
				{
					traverse(splittree,0,false);

					double dlog;
					if (BVITERBI)
					{
						dlog = trainhmmV(splittree,false);
					}
					else
					{
						dlog = trainhmm(splittree,false);
					}   


					if (!bmodeldiff)
					{
						if (BVITERBI)
						{
							dlog =testhmmV(splittree);
						}
						else
						{
							dlog =testhmm(splittree);
						}
					} 
					traverse(splittree,0,false);

					if (BDEBUG)
					{
						System.out.println("*test "+dlog);
					}

					if ((dlog>dbestlog)&&(bmodeldiff||(dlog-dprevouterbestlog)>EPSILONDIFF))
					{
						bestTree = splittree;
						dbestlog = dlog;
						dbesttrainlike = dtrainlike;
						double dimprove = (dbestlog-dprevouterbestlog)/Math.abs(dbestlog);
						String szimprove = nf3.format(100*dimprove) + "%";
						String szimproveDIFF = nf3.format(dbestlog-dprevouterbestlog);

						if (statusLabel3 != null)
						{
							statusLabel3.setText(" Next score: "+nf2.format(dbestlog));
						}

						if (statusLabel2 != null)
						{
							statusLabel2.setText(" Next score improvement: "+szimproveDIFF+" ("+szimprove+")");
						}
					}
				}
			}
		}
			}

	/////////////////////////////////////////////////////////////////////
	/**
	 * Helper function for adding a path in traverseandaddhelp. 
	 * Responsible for splitting one path out of a node into two. 
	 */
	private Treenode splitnode(int[] path, int ndesiredlevel, Treenode origroot) 
	{
		Treenode treeroot = (Treenode) origroot.clone();
		Treenode ptr = treeroot;

		for(int nindex = 0; nindex < ndesiredlevel; nindex++)
		{
			//advancing pointer to the node from which a new path should be added
			ptr = ptr.nextptr[path[nindex]];
		}

		if (!ptr.bchange)
		{
			return null;
		}
		else
		{
			//increasing the number of children
			ptr.numchildren++;
			ptr.binit = false;

			//readjusting the transition probabilities
			////////////////////////////////////////      

			if (BREGDREM)
			{
				for (int nindex = 0; nindex < ptr.numchildren-1; nindex++)
				{

					ptr.ptrans[nindex] = ptr.ptrans[nindex]*(ptr.numchildren-1)/(ptr.numchildren); 
				}  
				ptr.ptrans[ptr.numchildren-1] = 1.0/(ptr.numchildren);

				if (BDEBUG)
				{
					System.out.println("probs A "+ptr.ptrans[0]+"\t"+ptr.ptrans[1]);
				}
			}

			if (ndesiredlevel < path.length)
			{
				path[ndesiredlevel] = ptr.numchildren-1; 
			}
			////////////////////////////////////////////////////
			//more sophisticated initialization schemes possible here
			///////////////////////////////////////  
			for (int nlevel = ndesiredlevel; nlevel <path.length; nlevel++)
			{
				ptr.nextptr[ptr.numchildren-1]= new Treenode(ptr);
				if (ptr.numchildren > 1)
				{  
					double dk =0.02;
					double ddiff = ptr.dmean-ptr.nextptr[ptr.numchildren-2].dmean;

					if (Math.abs(ddiff) < dk*ptr.dsigma)
					{
						if (ptr.nextptr[ptr.numchildren-2].dmean > ptr.dmean)
						{
							ptr.nextptr[ptr.numchildren-1].dmean = ptr.dmean -
							(dk*ptr.dsigma - (ptr.nextptr[ptr.numchildren-2].dmean-ptr.dmean));

							if (BDEBUGMODEL)
							{
								System.out.println("A: MEAN "+ptr.nextptr[ptr.numchildren-1].dmean+" "+
										ptr.nextptr[ptr.numchildren-2].dmean+" "+ptr.dsigma);
							}
						}
						else
						{
							ptr.nextptr[ptr.numchildren-1].dmean = ptr.dmean +
							(dk*ptr.dsigma - (ptr.dmean-ptr.nextptr[ptr.numchildren-2].dmean));

							if (BDEBUGMODEL)
							{
								System.out.println("B: MEAN "+ptr.nextptr[ptr.numchildren-1].dmean+" "+
										ptr.nextptr[ptr.numchildren-2].dmean+" "+ptr.dsigma);
							}
						}
					}
					else
					{
						ptr.nextptr[ptr.numchildren-1].dmean = ptr.dmean;

						if (BDEBUGMODEL)
						{
							System.out.println("C: MEAN "+ptr.nextptr[ptr.numchildren-1].dmean+" "+
									ptr.nextptr[ptr.numchildren-2].dmean+" "+ptr.dsigma);
						}
					}
				}
				else
				{
					ptr.nextptr[ptr.numchildren-1].dmean = ptr.dmean;
				}

				ptr.nextptr[ptr.numchildren-1].dsigma = ptr.dsigma;


				ptr = ptr.nextptr[ptr.numchildren-1];
				ptr.numchildren = 1;
				if (BREGDREM)
				{
					ptr.ptrans[0] = 1;
				}
			}
			ptr.numchildren = 0;
		}
		return treeroot;     
	}


	//////////////////////////////////////////////////////////////////////
	/**
	 * With its helper function finds the path with the fewest genes assigned 
	 * this will become a candidate for deletion
	 */
	public MinPathRec traverseanddeleteMinPath(int[] path, int[] bestpath, Treenode ptr) throws Exception
	{
		return traverseanddeletehelpMinPath(path,0,ptr);
	}

	//////////////////////////////////////////////////////////////////////
	/**
	 * A record used for storing the path with the fewest genes assigned
	 */
	static class MinPathRec
	{
		int nval;
		int nlevel;
		int[] bestpath;

		MinPathRec(int nval, int nlevel, int[] bestpath)
		{
			this.nval = nval;
			this.nlevel = nlevel;
			this.bestpath = bestpath;
		}
	}

	/**
	 * Deletes the specificed path from the model starting from ndesiredlevel
	 */
	public void deleteMinPath(int[] path, int ndesiredlevel, Treenode root)
	{
		Treenode ptr = root;

		for(int nindex = 0; nindex < ndesiredlevel-1; nindex++)
		{
			ptr = ptr.nextptr[path[nindex]];
		}

		for (int nj = path[ndesiredlevel-1]+1; nj < ptr.numchildren; nj++)   
		{
			ptr.nextptr[nj-1] = ptr.nextptr[nj];
		}

		ptr.numchildren--; 
		ptr.binit = false;

		if (BREGDREM)
		{
			double dsum = 0; 
			for (int nindex = 0; nindex < ptr.numchildren; nindex++)
			{
				dsum += ptr.ptrans[nindex]; 
			}

			for (int nindex = 0; nindex < ptr.numchildren; nindex++)
			{
				if (dsum > 0)
				{
					ptr.ptrans[nindex]/=dsum; 
				}
				else
				{
					ptr.ptrans[nindex] = 0;
				}
			}
		}	
	} 

	/**
	 * Recursively searches for the path with the fewest genes assigned this will become a candidate
	 * for deletion
	 */
	private MinPathRec traverseanddeletehelpMinPath(int[] path, int nlevel,
			Treenode ptr) throws Exception
			{
		int nmin = Integer.MAX_VALUE;
		if (nlevel < path.length)
		{
			MinPathRec theminMinPathRec=null;
			for (int nindex = 0; nindex < ptr.numchildren; nindex++)
			{
				path[nlevel] = nindex;
				MinPathRec thetempMinPathRec
				= traverseanddeletehelpMinPath(path,nlevel+1,ptr.nextptr[nindex]);

				if (thetempMinPathRec.nval < nmin)
				{
					nmin = thetempMinPathRec.nval;
					theminMinPathRec = thetempMinPathRec;
					theminMinPathRec.bestpath[nlevel] = nindex;
					if (thetempMinPathRec.nval == ptr.numPath)
					{
						thetempMinPathRec.nlevel = nlevel;
					}
				}
			}
			return theminMinPathRec;
		}
		else
		{
			MinPathRec theMinPathRec = new MinPathRec(ptr.numPath,nlevel,new int[path.length]);
			return theMinPathRec;
		}
			}


	////////////////////////////////////////////////////////////////////////////
	/**
	 * With its helper method traverseanddelayhelp
	 * searches for splits to delay one time point
	 * Returns true if an eligible split to delay was found
	 */
	public boolean traverseanddelay(int[] path, Treenode ptr,int ndesiredlevel, 
			Treenode origroot,boolean bresplit) throws Exception
			{
		if (BDEBUG)
		{
			System.out.println("entering traverse and delay");
		}
		return traverseanddelayhelp(path,0,ndesiredlevel,ptr,origroot,bresplit,0);
			}


	/////////////////////////////////////////////////////////////////////////////
	/**
	 * Recursively searches for splits to delay one time point
	 */
	private boolean traverseanddelayhelp(int[] path, int nlevel,int ndesiredlevel,
			Treenode ptr, Treenode origroot,boolean bresplit,int nchild) throws Exception
			{

		if (BDEBUG)
		{
			System.out.println("traverseanddelayhelp "+ndesiredlevel+" "+ndesiredlevel);
		}

		boolean breturnval = false;
		if ((nlevel < path.length-1)&&(nlevel<ndesiredlevel))
		{
			for (int nindex = 0; nindex < ptr.numchildren; nindex++)
			{
				if (DREM_IO.bdisplaycurrent)
				{
					displayTempMap();
				}
				path[nlevel] = nindex;
				boolean breturnvaltemp;
				breturnvaltemp = traverseanddelayhelp(path,nlevel+1,ndesiredlevel,
						ptr.nextptr[nindex],origroot,bresplit,nindex);
				breturnval = breturnval || breturnvaltemp;
			}
		}

		//assuming ptr.parent also not null
		if ((ptr !=null)&&(nchild >= 1)&&(nlevel==ndesiredlevel))
		{
			for (int nchildcombine =0; nchildcombine < nchild; nchildcombine++)
			{
				if (ptr.numchildren + ptr.parent.nextptr[nchildcombine].numchildren <= nmaxchild)
				{
					double dabsmeandiff = Math.abs(ptr.dmean-ptr.parent.nextptr[nchildcombine].dmean);
					NumberFormat nf4 = NumberFormat.getInstance();
					nf4.setMinimumFractionDigits(4);
					nf4.setMaximumFractionDigits(4);
					if (BDEBUG)
					{
						System.out.println(nlevel+" "+nf4.format(ptr.dmean)+" "
								+nf4.format(ptr.parent.nextptr[nchildcombine].dmean)+" "
								+nf4.format(dabsmeandiff)+" "+nf4.format(ptr.dsigma)+" "+
								nf4.format(ptr.parent.nextptr[nchildcombine].dsigma));
						System.out.println("-----------");
					}

					//mark to try to delete
					if (BDEBUG)
					{
						System.out.print("dy ");
						for (int ni =0; ni <nlevel; ni++)
						{
							System.out.print(path[ni]+" ");
						}
						System.out.println();
						System.out.println("before");
					}

					traverse(origroot,0, false);
					Treenode delaytree = delaysplit(path,nlevel,origroot,nchildcombine);  

					if (delaytree != null)
					{
						double dlog;	      
						if (BVITERBI)
						{
							dlog = trainhmmV(delaytree,true);
						}
						else
						{
							dlog = trainhmm(delaytree,true);
						}

						if (!bmodeldiff)
						{
							if (BVITERBI)
							{ 
								dlog = testhmmV(delaytree);
							}
							else
							{
								dlog = testhmm(delaytree);
							}
						}

						if (BDEBUG)
						{
							System.out.println("after "+dlog);
						}
						traverse(delaytree,0, false);

						if (BDEBUG)
						{
							System.out.println("test "+dlog);
						}

						double ddelayimprove = (dlog-dprevouterbestlog)/Math.abs(dlog);

						if (BDEBUG)
						{
							System.out.println("DELAY\t"+dlog+"\t"+dbestlog+"\t"+
									ddelayimprove+"\t"+(dlog-dprevouterbestlog)/Math.abs(dlog)+"\t"+DELAYMIN);
						}

						if (BDEBUG)
						{
							System.out.println("delay = "+ddelayimprove+"\t"+DELAYMIN);
						}

						if ((((dlog+DELAYMIN*Math.abs(dlog)-dprevouterbestlog)>=DELAYMINDIFF)
								&&(dlog>dbestlog)&&!bmodeldiff)||
								((bmodeldiff)&&(dlog >= dprevouterbestlog)&&(dlog>dbestlog)))
						{
							bestTree = delaytree;
							dbestlog = dlog;
							if (BDEBUG)
							{
								System.out.println("in delay "+dbestlog);
							}

							double dimprove = (dbestlog-dprevouterbestlog)/Math.abs(dbestlog);
							String szimprove = nf3.format(100*dimprove) + "%";
							String szimproveDIFF = nf3.format(dbestlog-dprevouterbestlog);
							if (statusLabel2 != null)
							{
								statusLabel2.setText(" Delay next change: "+szimproveDIFF+" ("+szimprove+")"+
										" (Score "+nf2.format(dbestlog)+")");

							}

							dbesttrainlike = dtrainlike;
							breturnval = true;
						}
					}
				}
			}
		}
		return breturnval;
			}


	//////////////////////////////////////////////////////////////////////
	/**
	 * With its helper function searches for the best path to delete
	 * Returns true if an eligible path to delete was found
	 */
	public boolean traverseanddelete(int[] path, Treenode ptr, 
			Treenode origroot,boolean bresplit,
			double dimprovemin,double dimprovemindiff) throws Exception
			{
		return traverseanddeletehelp(path,0,ptr,origroot,
				bresplit,dimprovemin,dimprovemindiff);
			}

	//////////////////////////////////////////////////////////////////////
	/**
	 * Helper function that searches for the best path to delete
	 */
	private boolean traverseanddeletehelp(int[] path, int nlevel,
			Treenode ptr, Treenode origroot,boolean bresplit,
			double dimprovemin,double dimprovemindiff) throws Exception
			{

		boolean breturnval = false;
		if (nlevel < path.length)
		{
			for (int nindex = 0; nindex < ptr.numchildren; nindex++)
			{
				if (DREM_IO.bdisplaycurrent)
				{
					displayTempMap();
				}

				path[nlevel] = nindex;
				boolean breturnvaltemp = 
					traverseanddeletehelp(path,nlevel+1,ptr.nextptr[nindex],
							origroot,bresplit,dimprovemin,dimprovemindiff);
				breturnval = (breturnvaltemp || breturnval);
			}
		}

		if ((ptr !=null)&&(ptr.parent!=null)&&(ptr.parent.numchildren > 1))
		{
			if (BDEBUG)
			{
				//mark to try to delete
				for (int ni =0; ni <nlevel; ni++)
				{
					System.out.print(path[ni]+" ");
				}
				System.out.println();
			}

			Treenode deletetree = deletepath(path,nlevel,origroot);  

			if (deletetree != null)
			{
				double dlog;

				if (BVITERBI)
				{
					dlog = trainhmmV(deletetree,false);
				}
				else
				{
					dlog = trainhmm(deletetree,false);
				} 

				if (!bmodeldiff)
				{
					if (BVITERBI)
					{
						dlog = testhmmV(deletetree);
					}
					else
					{
						dlog =testhmm(deletetree);
					}
				}

				if (BDEBUG)
				{
					System.out.println("test "+dlog);
				}

				double dimprovedelete = (dlog-dprevouterbestlog)/Math.abs(dlog);
				if (BDEBUG)
				{
					System.out.println("del* "+dimprovedelete+"\t"+dimprovemin+"\t"+(dlog-dprevouterbestlog));
					System.out.println("~~"+dlog+"\t"+dbestlog+"\t"+dprevouterbestlog);

					System.out.println("DELETE\t"+dlog+"\t"+dbestlog+"\t"+dimprovedelete+
							"\t"+(dlog-dprevouterbestlog)/Math.abs(dlog)+"\t"+dimprovemin);
				}

				if ((((dlog+dimprovemin*Math.abs(dlog)-dprevouterbestlog)>=dimprovemindiff)
						&&(dlog>dbestlog)&&!bmodeldiff)||
						(bmodeldiff && (dimprovemin>=dprevouterbestlog)&&(dlog>dbestlog)))

				{
					bestTree = deletetree;
					dbestlog = dlog;
					if (BDEBUG)
					{
						System.out.println("in delete "+dbestlog);
					}

					if (!bresplit)
					{
						double dimprove = (dbestlog-dprevouterbestlog)/Math.abs(dbestlog);
						String szimprove = nf3.format(100*dimprove) + "%";
						String szimproveDIFF = nf3.format((dbestlog-dprevouterbestlog));
						if (statusLabel3 != null)
						{
							statusLabel3.setText(" Next score: "+nf2.format(dbestlog));
						}

						if (statusLabel2 != null)
						{
							statusLabel2.setText(" Next score improvement: "+szimproveDIFF+" ("+szimprove+")");
						}
					}
					else
					{
						double dimprove = (dbestlog-dprevouterbestlog)/Math.abs(dbestlog);
						String szimprove = nf3.format(100*dimprove) + "%";
						String szimproveDIFF = nf3.format((dbestlog-dprevouterbestlog));

						if (statusLabel2 != null)
						{
							statusLabel2.setText(" Delete path change : "+szimproveDIFF+" ("+szimprove+")");
						}
					}
					dbesttrainlike = dtrainlike;
					breturnval = true;
				}
			}
		}

		return breturnval;
			}



	///////////////////////////////////////////////////////////////////////
	/**
	 * Prints the most significant TF-association scores for each TF
	 */
	public void printMinPvals()
	{
		double[] bestpsplit = new double[tfNames.length];
		double[] bestpfull = new double[tfNames.length];

		int[] bestsplitdepth = new int[tfNames.length];
		int[] bestfulldepth = new int[tfNames.length];

		for (int ni = 0; ni < bestpsplit.length; ni++)
		{
			bestpsplit[ni] = 1;
		}

		for (int ni = 0; ni < bestpsplit.length; ni++)
		{
			bestpfull[ni] = 1;
		}

		getMinPvals(treeptr, bestpsplit, bestpfull,bestsplitdepth,bestfulldepth,1);

		try
		{
			PrintWriter pw = new PrintWriter(new FileOutputStream(SCORESDIR+"/"+nrandomseed+".txt"));
			for (int ni = 0; ni < tfNames.length; ni++)
			{
				pw.println(tfNames[ni]+"\t"+bestpfull[ni]+"\t"+bestpsplit[ni]
				                                                          +"\t"+bestfulldepth[ni]+"\t"+bestsplitdepth[ni]);
			}
			pw.close();
		}
		catch (IOException ioex)
		{
			ioex.printStackTrace(System.out);
		}
	}


	/**
	 * Helper procedure for getting the most significant TF-association scores for each TF
	 */
	private void getMinPvals(Treenode ptr, double[] bestpsplit, double[] bestpfull,
			int[] bestsplitdepth, int[] bestfulldepth, int ndepth)
	{
		if (ptr.numchildren >=1)
		{
			for (int nchild = 0; nchild < ptr.numchildren; nchild++)
			{
				if (ptr.dpvalEdgeSplit != null)
				{
					for (int ni = 0; ni < bestpsplit.length; ni++)
					{
						if (ptr.dpvalEdgeSplit[nchild][ni] < bestpsplit[ni])
						{
							bestpsplit[ni] = ptr.dpvalEdgeSplit[nchild][ni];
							bestsplitdepth[ni] = ndepth;
						}
					}
				}

				if (ptr.dpvalEdgeFull.length >= 1)
				{
					for (int ni = 0; ni < bestpfull.length; ni++)
					{
						if (ptr.dpvalEdgeFull[nchild][ni] < bestpfull[ni])
						{
							bestpfull[ni] = ptr.dpvalEdgeFull[nchild][ni];
							bestfulldepth[ni] = ndepth;
						}
					}
				}

				getMinPvals(ptr.nextptr[nchild], bestpsplit, bestpfull,
						bestsplitdepth, bestfulldepth, ndepth+1);
			}
		}
	}


	//////////////////////////////////////////////////////////////////////
	/**
	 * A Treenode corresponds to a state in the model 
	 */
	class Treenode
	{
		int ncurrtime;//use to count visited
		int nID;
		int niteration;
		int ninstanceiteration;
		int nrec;
		int ninstancerec;
		int nminparentlevel;
		Treenode minparentptr;

		Treenode[] parentptrA;
		/**number of leaves descendants of this node, only relevant if split*/
		int numdownleaves;
		double dens;
		/**expected sum of the emissions from this node */
		double dEsum; 
		/**expected sum of the square of emissions from this node*/
		double dEsumsq;
		double dPsum;
		/**emission mean*/
		double dmean;
		/**emission std*/
		double dsigma;
		double dpredictweight;
		String szgolabel="";
		String szgenesetlabel = "";
		String sztfsetlabel = "";
		/**non zero recs*/
		int nprime;
		/**number of children*/
		int numchildren;
		/**can change state's values*/
		boolean bchange;
		int ndepth;
		int nparentcolorindex;
		boolean[] bvector;

		/**The logistic regression classifier*/
		DREM_FastLogistic2 tranC;
		double[] recweight;
		double[] dTotalScore;
		double[][] dpvals;
		/** A sorted list of the Treenode's children such that orderA[0] gives the
		 * index of the 0th (lowest) child node */
		int[] orderA;
		/** The reverse mapping of orderA such that revOrderA[i] gives the 0-indexed
		 * order (lowest to highest) of child i.
		 */
		int[] revOrderA;
		/**Stores the number of genes regulated by a TF on a given path with the given
		 * interaction value.  Dimensions:
		 * <br>[num TFs][num paths out of split][num distinct binding interaction values]*/
		int[][][] ncountvals;
		/**Stores the number of genes regulated by a TF on the paths other than
		 * the given path with the given interaction value.  Dimensions:
		 * <br>[num TFs][num paths out of split][num distinct binding interaction values]*/
		int[][][] nothersum;
		double[][] dscoreother;
		double[][] dscore;
		double[][] ddiff;
		/** A TreeSet used to display important TFs based on activity score,
		 * which isn't a measure of significance in the statistical sense. */
		TreeSet tsSigTFActivity;
		TreeSet tsSigTF;
		TreeSet[] tsSigTFEdgeSplit;
		TreeSet[] tsSigTFFull;
		Hashtable[] htScore;
		/**Stores the number of genes regulated by a TF on all paths out of a split
		 * with the given interaction value.  Dimensions:
		 * <br>[num TFs][num distinct binding interaction values]*/
		int[][] ncountTotals;
		double[][] dpvalEdgeSplit;
		double[][] ddiffEdgeSplit;
		double[][] dpvalEdgeFull;
		double[][] ddiffEdgeFull;
		double[][] dexpectEdgeFull;
		double[][] dexpectEdgeSplit;
		double[][] dratioSplit;
		/**An array of activity scores for all TFs based on how the
		 * genes the TF regulates transition to the next states compared
		 * to how all genes into the split transistion.*/
		double[] dTFActivityScore;


		boolean[] bInNode;
		/**Number of genes going through this state???*/
		int numPath;
		/**probability of transitioning to each of the next states*/
		double[] ptrans;
		/**expected number of transitions given the training sequences*/
		double[] aprobs;
		/**pointer to each of the next states*/
		Treenode[] nextptr;
		/**pointer back to the parent*/
		Treenode parent;

		PText thepredictText;
		PText goText;
		PText genesetText;
		PText tfsetText;
		/**forward variable for current variables
		 * <br>p(x_1,...,x_i,pi_i=k)*/
		double df;
		/**backward variable for current variables
		 * <br>p(x_(i+1),...,x_L|pi_i=k)*/
		double db; 
		boolean binit; 

		/**
		 * Calls inittreenode with a null parent
		 */
		Treenode()
		{
			inittreenode(null);
		}


		Treenode(Treenode parent)
		{
			inittreenode(parent);
		}

		/**
		 * Calls inittreenode with parent
		 */
		void inittreenode(Treenode parent)
		{
			niteration = 0;
			ninstanceiteration = 0;
			nrec = -1;
			dmean = 0;
			dsigma = .5;
			this.parent = parent;   
			dEsum = 0;
			dEsumsq = 0;
			dPsum =0;
			nprime = 0;
			bchange = true;
			binit = false;
			if (parent == null)
			{
				ndepth = 0;
			}
			else
			{
				ndepth = parent.ndepth +1;
			}
			nextptr = new Treenode[nmaxchild];


			if (BREGDREM)
			{
				ptrans   = new double[nmaxchild];
				aprobs = new double[nmaxchild];
			}

			bvector = new boolean[numbits];

			if (bindingSignGene != null)
			{
				dTotalScore = new double[numbits];
			}

			// Not used???
			//int numPath;

			for (int nindex = 0; nindex < numbits; nindex++)
			{
				bvector[nindex] = false;
			}

			numchildren = 0; 
			for (int nindex = 0; nindex < nextptr.length; nindex++)
			{
				nextptr[nindex] = null;
			}
		} 

		/**
		 * For making copy of nodes
		 */
		public Object clone()
		{
			if (parent == null)
			{
				if (BDEBUG)
				{
					System.out.println("in clone........");
				}
				htBackPts = new Hashtable();
			}

			Treenode tnode = new Treenode();
			tnode.nID = nID;
			tnode.dmean = dmean;
			tnode.dsigma = dsigma;

			tnode.ndepth = ndepth;
			tnode.dEsum = dEsum;
			tnode.dEsumsq = dEsumsq;
			tnode.dPsum = dPsum;
			tnode.nprime = nprime;
			tnode.numchildren = numchildren;

			tnode.thepredictText = thepredictText;
			tnode.parent = null;  
			tnode.bchange = bchange;
			tnode.df =df;  
			tnode.db = db;
			tnode.binit = binit;
			tnode.nextptr = new Treenode[nmaxchild];

			if (BREGDREM)
			{
				tnode.ptrans   = new double[nmaxchild];
				tnode.aprobs = new double[nmaxchild];
			}
			else if (tranC != null)
			{
				tnode.tranC =  (DREM_FastLogistic2) tranC.clone();  //do we need a clone yes
			}

			tnode.recweight = new double[recweight.length];


			for (int nindex = 0; nindex < nextptr.length; nindex++)
			{
				if (nextptr[nindex] == null)
				{
					tnode.nextptr[nindex] = null;
				}
				else
				{
					BackPtrRec childNodeRec = (BackPtrRec) htBackPts.get(nextptr[nindex]);

					if (childNodeRec != null)
					{
						//we've already cloned the child node, just have it link to it
						tnode.nextptr[nindex] = childNodeRec.childNode;
						//have the child link back to this node
						if (BDEBUG)
						{
							System.out.println("[["+
									childNodeRec+" "+childNodeRec.childNode+" "+childNodeRec.childNode.parentptrA);
						}

						int nfindparentindex = 0;
						while (nextptr[nindex].parentptrA[nfindparentindex] != this)
						{
							nfindparentindex++;
						}
						if (BDEBUG)
						{
							System.out.println("nfindparentindex = "+nfindparentindex);
						}
						childNodeRec.childNode.parentptrA[nfindparentindex] = tnode;
					}
					else
					{
						//cloning a new node
						tnode.nextptr[nindex] = (Treenode) nextptr[nindex].clone();
						tnode.nextptr[nindex].parent = tnode;

						if (nextptr[nindex].parentptrA != null)
						{
							tnode.nextptr[nindex].parentptrA = new Treenode[nextptr[nindex].parentptrA.length]; 
							for (int nparentindex =0; nparentindex < 
							nextptr[nindex].parentptrA.length; nparentindex++)
							{
								if (nextptr[nindex].parentptrA[nparentindex] == this)
								{
									if (BDEBUG)
									{
										System.out.println("~~``"+tnode.nextptr[nindex].dmean+" "+
												nparentindex);
									}
									tnode.nextptr[nindex].parentptrA[nparentindex] = tnode;
								}
							}

							if (nextptr[nindex].parentptrA.length >= 2)
							{
								if (BDEBUG)
								{
									System.out.println(dmean+" "+nextptr[nindex].dmean);
								}
								htBackPts.put(nextptr[nindex], new BackPtrRec(tnode.nextptr[nindex]));
							}
						}
					}
				}

				if (BREGDREM)
				{
					tnode.ptrans[nindex] = ptrans[nindex];
					tnode.aprobs[nindex] =aprobs[nindex];
				}
			}

			tnode.bvector = new boolean[numbits];
			for (int nindex = 0; nindex < numbits; nindex++)
			{
				tnode.bvector[nindex] = bvector[nindex];
			}

			return tnode;
		}
		
		/**
		 * 
		 * @return a map from TF names (excluding the " [#]" used to indicate
		 * the primary path out of the split for the TF) to the best activity
		 * score seen at any of this node's ancestors.
		 */
		public HashMap getAncestorActivityScores()
		{
			HashMap bestSeenScores = new HashMap();
			// TODO need to handle the case where there are multiple parents?
			if(parent != null)
			{
				return parent.getAncestorActivityScores(bestSeenScores);
			}
			else
			{
				return bestSeenScores;
			}
		}
		
		/**
		 * Helper function for the other version of getAncestorActivityScores.
		 * @param bestSeenScores The best activity scores seen
		 * in ancestors of the original node before the search reached
		 * the current node.
		 * @return a map from TF names (excluding the " [#]" used to indicate
		 * the primary path out of the split for the TF) to the best activity
		 * score seen at any of this nodes ancestors or in the bestSeenScores.
		 */
		private HashMap getAncestorActivityScores(HashMap bestSeenScores)
		{
			// Enumerate all TFs in the tree set and add their activty
			// score to the map if there was no score for that TF or if
			// the score at the current node is better (i.e. lower, because
			// the activity scores have been transformed to resemble p-values)
			Iterator itr = tsSigTFActivity.iterator();
			DREM_Timeiohmm.SigTFRecv2 theTFRec;
			while (itr.hasNext())
			{
				theTFRec =  (DREM_Timeiohmm.SigTFRecv2) itr.next();

				String tfName = theTFRec.szname;
				// Check if the path out of the split has been added as " [#]"
				// to the TF name, and if so remove it
				// TODO store the path index separately so this does not
				// need to be done
				if(tfName.endsWith("]"))
				{
					tfName = tfName.substring(0,tfName.lastIndexOf(" ["));
				}
				
				// Add the TF and score to the map if they are not already present
				if(!bestSeenScores.containsKey(tfName))
				{
					// The "pval" is actually the transformed activity score
					bestSeenScores.put(tfName, Double.valueOf(theTFRec.dpval));
				}
				// Otherwise add them if the current score is lower
				// than the previously seen best lower score
				else
				{
					double bestScore = ((Double) bestSeenScores.get(tfName)).doubleValue();
					
					if(theTFRec.dpval < bestScore)
					{
						bestSeenScores.put(tfName, Double.valueOf(theTFRec.dpval));
					}					
				}
			} // end iteration over all TF activity scores at this node
		
			
			// TODO need to handle the case where there are multiple parents?
			if(parent != null)
			{
				bestSeenScores = parent.getAncestorActivityScores(bestSeenScores);
			}
			
			return bestSeenScores;
		}
	}

	/**
	 * A record with information about a state a back pointer points to
	 */
	static class BackPtrRec
	{
		BackPtrRec(Treenode childNode)
		{
			this.childNode = childNode;
		}
		Treenode childNode;
	}

	/////////////////////////////////////////////////////////////////////
	/**
	 * Returns true iff treeptr or a descendent has two or more parents
	 */
	public boolean hasMerge(Treenode treeptr)
	{
		boolean bhasmerge = ((treeptr.parentptrA != null) && (treeptr.parentptrA.length >= 2));
		for (int nchild = 0; (nchild < treeptr.numchildren)&&(!bhasmerge); nchild++)
		{
			bhasmerge = hasMerge(treeptr.nextptr[nchild]);
		}

		return bhasmerge;
	}

	//////////////////////////////////////////////////////////////////////
	/**
	 * Returns a string with the model parameters
	 */
	public String saveString(Treenode treecopy)
	{
		StringBuffer szbuf = new StringBuffer("");;

		bhasmerge = hasMerge(treecopy);

		if (bhasmerge)
		{
			szbuf.append("MERGED TREE\n");
		}

		if (BREGDREM)
		{
			szbuf.append("REGDREM Num. Coefficients\t"+tfNames.length+"\n");
		}
		else
		{
			szbuf.append("Num. Coefficients\t"+tfNames.length+"\n");
		}

		HashMap ptrMap = new HashMap();
		ntotalID = -1;
		szbuf.append(saveStringHelp(treecopy,ptrMap));
		return szbuf.toString();
	}

	/**
	 * Helper function for generating a string with model parameters
	 * that traverses all the states
	 */
	private StringBuffer saveStringHelp(Treenode ptr,HashMap ptrMap)
	{
		StringBuffer szbuf = new StringBuffer();
		if (ptr != null)
		{
			if (bhasmerge)
			{
				Integer obj = (Integer) ptrMap.get(ptr);
				int nID;
				if (obj != null)
				{
					nID = obj.intValue();
				}
				else
				{
					ntotalID++;
					nID = ntotalID;
					ptrMap.put(ptr,new Integer(nID));
				}
				szbuf.append(nID+"\t");
			}

			szbuf.append(ptr.dmean+"\t"+ptr.dsigma+"\t"+ptr.numchildren+"\n");
			if (ptr.numchildren > 1)
			{
				if (BREGDREM)
				{
					for (int nchild = 0;nchild < ptr.numchildren; nchild++)
					{
						szbuf.append(ptr.ptrans[nchild]+"\n");
					}  
				}
				else
				{
					szbuf.append(ptr.tranC.saveClassifier(tfNames));
				}

			}

			for (int nindex = 0; nindex < ptr.numchildren; nindex++)
			{
				if (ptr.nextptr[nindex] != null)
				{
					szbuf.append(saveStringHelp(ptr.nextptr[nindex],ptrMap));
				}
			}
		}

		return szbuf;
	}

	
	/**
	 * Save the activity scores at this node and its children
	 * @return
	 */
	public StringBuffer saveActivityScores(Treenode node)
	{
		StringBuffer buf = new StringBuffer();
		
		if(node != null && node.numchildren >= 1)
		{
			if(node.numchildren >= 2)
			{
				for(int t = 0; t < tfNames.length; t++)
				{
					buf.append(tfNames[t]).append("\t");
					buf.append(node.dTFActivityScore[t]).append("\t");
					buf.append(node.ndepth).append("\n");
				}
			}
			
			for (int nindex = 0; nindex < node.numchildren; nindex++)
			{
				if (node.nextptr[nindex] != null)
				{
					buf.append(saveActivityScores(node.nextptr[nindex]));
				}
			}
		}
		
		return buf;
	}
	

	//////////////////////////////////////////////////////////////////////
	/**
	 * Record of information about predictions
	 */
	static class PredictRec
	{
		double[] mean;
		double[] meansq;
		double[] var;
		double[] entropy;
	}


	/**
	 * Used to set probability of each state a gene with a TF-binding signature of binding will be in
	 * Also computes a record of other statistics based on these predictions
	 */
	public PredictRec predictTime(int[] binding, double dclassprob,DREM_Timeiohmm.Treenode treeptr) throws Exception
	{
		PredictRec thePredictRec = new PredictRec();
		thePredictRec.mean = new double[traindata[0].length];
		thePredictRec.var = new double[traindata[0].length];
		thePredictRec.meansq = new double[traindata[0].length];
		thePredictRec.entropy = new double[traindata[0].length];
		predictTimeHelp(treeptr,binding,dclassprob,0,
				thePredictRec.mean,thePredictRec.var, 
				thePredictRec.meansq,thePredictRec.entropy);

		return thePredictRec;

	}


	/**
	 * A recursive helper function to predictTime
	 */
	private void predictTimeHelp(Treenode node, int[] theInstance, double dweight,
			int nstep,double[] meanpredict, double[] varpredict,
			double[] meansqpredict,double[] entropypredict)  throws Exception
			{
		if (nstep < meanpredict.length)
		{
			NumberFormat nf1 = NumberFormat.getInstance(Locale.ENGLISH);

			if (dweight < 1)
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
			node.dpredictweight = dweight;

			if (node.thepredictText != null)
			{
				node.thepredictText.setText(nf1.format(dweight));
			}

			meanpredict[nstep] += dweight * node.dmean;
			varpredict[nstep] +=  Math.pow(dweight*node.dsigma, 2);
			meansqpredict[nstep] += dweight*node.dmean*node.dmean;
			entropypredict[nstep] -= dweight*Math.log(dweight)/Math.log(2);

			if (node.binit)
			{
				double[] pdist;
				if (BREGDREM)
				{
					pdist = node.ptrans;
				}
				else
				{
					pdist = node.tranC.distributionForInstance(theInstance,-1);
				}

				if (BDEBUG)
				{
					if (pdist.length==1)
						System.out.println("warning "+pdist[0]);
				}

				for (int nchild = 0; nchild < node.numchildren; nchild++)
				{
					predictTimeHelp(node.nextptr[nchild],theInstance,
							dweight*pdist[nchild],nstep+1,meanpredict,varpredict,meansqpredict,entropypredict);
				}
			}
			else
			{
				predictTimeHelp(node.nextptr[0],  theInstance,dweight,nstep+1,meanpredict,varpredict,
						meansqpredict,entropypredict);
			}
		}
			}



	///////////////////////////////////////////////////////////////////////
	/**
	 * Used for debugging to see the contents of the model
	 */
	public void traverse(Treenode root, int ndepth, boolean bprintC) 
	{
		if (BDEBUG)
		{
			for (int nindex =0; nindex < ndepth; nindex++)
			{
				System.out.print("\t");
			}

			System.out.print(root.numPath+"\t"+root.dmean+"\t"+root.dsigma+"\t"+	       
					root.numchildren+"\t"+root.binit+"\t"+
					root.dpredictweight+"\t"+root.ndepth+"\t#"+root.numdownleaves+"\t");

			if (root.parent != null)
				System.out.print("p "+root.parent.dmean);

			if (root.parentptrA != null)
			{
				for (int i = 0; i < root.parentptrA.length; i++)
					System.out.print("\t"+nf2.format(root.parentptrA[i].dmean));
			}

			if (root.ptrans != null)
				System.out.print("\t"+root.ptrans[0]+"\t"+root.ptrans[1]);

			for (int nindex =ndepth; nindex < theDataSet.data[0].length-1; nindex++)
			{
				System.out.print("\t");
			}

			if (root.nextptr[0] != null)
			{
				System.out.println();

				if ((bprintC)&&(root.numchildren > 1))
				{  
					System.out.println(root.tranC);    
				}
			}
			else
			{  
				System.out.println();
			}

			for (int nindex = 0; nindex < root.numchildren; nindex++)
			{
				if (root.nextptr[nindex] != null)
				{
					traverse(root.nextptr[nindex], ndepth+1, bprintC);
				}
			}
		}
	}

	/////////////////////////////////////////////////////////////////////////


	/**
	 * Combines the results of the forward and backward algorithm to obtain the probability of
	 * each gene going through each state when using viterbi training
	 */
	public void instanceAEV(Treenode root,double[] vals,int[] pma, int[] bpvals,int[] bestpath,
			double dpj,int nrec) throws Exception
			{
		Treenode tempptr = root;

		for (int ndepth = 0; ndepth < bestpath.length; ndepth++)
		{
			tempptr.dEsum += vals[ndepth]; 
			tempptr.dPsum += 1;
			tempptr.dEsumsq += Math.pow(vals[ndepth],2);
			tempptr.nprime++;

			if (BREGDREM)
			{
				tempptr.aprobs[bestpath[ndepth]]++;
			}

			if (tempptr.numchildren >= 1)
			{
				tempptr.recweight[tempptr.numchildren*nrec+bestpath[ndepth]] = 1;
			}

			tempptr = tempptr.nextptr[bestpath[ndepth]];
		}
			}

	///////////////////////////////////////////////////////////////////////////////////
	/**
	 * Initializes all the fields that will be used by instanceAE
	 */
	public void initAE(Treenode root)
	{
		if (root != null)
		{
			for (int nchild = 0; nchild < root.numchildren; nchild++)
			{
				if (BREGDREM)
				{
					root.aprobs[nchild] = 0;
				}
				initAE(root.nextptr[nchild]);             
			}

			if ((root.recweight==null)|| (root.recweight.length != root.numchildren*traindata.length))
			{
				root.recweight = new double[root.numchildren*traindata.length];
			}

			root.dEsum = 0;
			root.dEsumsq = 0;
			root.dPsum = 0;
			root.nprime = 0;
			root.niteration = 0;
			root.ninstancerec = 0;
			root.ninstanceiteration = 0;
			root.nrec = -1;
		}
	}


	/**
	 * Combines the results of the forward and backward algorithm to obtain the probability of
	 * each gene going through each state when using baum-welch training
	 */
	public void instanceAE(Treenode root,int ndepth,double[] vals,int[] pma, int[] bpvals,
			double dpj,int nrec) throws Exception
			{

		if ((root.ninstancerec != nrec)||(root.ninstanceiteration != nglobaliteration))
		{
			root.ninstanceiteration = nglobaliteration;
			root.ninstancerec = nrec;

			if (root.nextptr[0] != null)
			{
				for (int nchild = 0; nchild < root.numchildren; nchild++)
				{
					double dmultval = root.nextptr[nchild].df *
					root.nextptr[nchild].db/dpj;

					//P(A,B)P(C|B)=P(A,B)P(C|A,B)=P(A,B,C)

					//double df;  //p(x_1,...,x_i,pi_i=k)
					//double db; //p(x_(i+1),...,x_L|pi_i=k)
					//p(x_1,...,x_L,pi_i=k)
					//=p(x_1,...,x_i,pi_i=k)*p(x_(i+1),...,x_L|x_1,...,x_i,pi_i=k)
					//=p(x_1,...,x_i,pi_i=k)*p(x_(i+1),...,x_L|pi_i=k)
					//p(x_1,...,x_L) = p(pi_i=k|x_1,...x_L)
					root.recweight[root.numchildren*nrec+nchild] = dmultval;

					if (BREGDREM)
					{
						root.aprobs[nchild] += dmultval;
					}

					instanceAE(root.nextptr[nchild],ndepth+1,vals,pma,bpvals,dpj,nrec);             
				}
			}

			double dweight =root.df*root.db/dpj;
			if ((dweight > 0)&& (pma[ndepth] != 0)&&(!Double.isInfinite(dweight)))
			{
				// Not used???
				//double dtemp = 0;

				root.dEsum += vals[ndepth]*dweight;
				root.dPsum += dweight;
				root.dEsumsq +=vals[ndepth]*vals[ndepth]*dweight;
				root.nprime++;
			}
		}
			}

	/////////////////////////////////////////////////////////////////////////
	/**
	 * Initializes all the fields that will be used by instanceAEV
	 */
	public void initAEV(Treenode root)
	{
		if (root != null)
		{
			if (BREGDREM)
			{
				for (int nchild = 0; nchild < root.aprobs.length; nchild++)
				{
					root.aprobs[nchild] = 0;
				}
			}

			for (int nchild = 0; nchild < root.numchildren; nchild++)
			{
				if (BREGDREM)
				{
					root.aprobs[nchild] = 0;
				}
				initAEV(root.nextptr[nchild]);             
			}


			if ((root.recweight==null)|| (root.recweight.length != root.numchildren*traindata.length))
			{
				root.recweight = new double[root.numchildren*traindata.length];
			}
			else
			{
				for (int i=0; i < root.recweight.length; i++)
				{
					root.recweight[i] = 0;
				}
			}

			root.dEsum = 0;
			root.dEsumsq = 0;
			root.dPsum = 0;
			root.nprime = 0;
		}
	}

	////////////////////////////////////////////////////////////////////////
	/**
	 * Forces children of a parent to both have the average of their standard
	 * deviations
	 */
	public void averageChildrenSigmas(Treenode root)
	{
		if (root.numchildren >=2)
		{
			double dsigmasum = 0;
			double dpsum = 0;
			for (int nchild = 0; nchild< root.numchildren; nchild++)
			{
				dsigmasum +=  root.nextptr[nchild].dPsum*root.nextptr[nchild].dsigma;
				dpsum += root.nextptr[nchild].dPsum;
				if (BDEBUG)
				{
					System.out.println("##"+root.nextptr[nchild].dsigma+"\t"+root.nextptr[nchild].dPsum);
				}
			}
			double dsigmaavg = dsigmasum/dpsum;
			if (BDEBUG)
			{
				System.out.println("&&&"+dpsum+"\t"+root.dPsum+"\t"+dsigmasum+"\t"+dsigmaavg);
			}

			for (int nchild = 0; nchild< root.numchildren; nchild++)
			{
				root.nextptr[nchild].dsigma = dsigmaavg;
			}
		}

		for (int nchild = 0; nchild< root.numchildren; nchild++)
		{
			averageChildrenSigmas(root.nextptr[nchild]);
		}   
	}

	/////////////////////////////////////////////////////////////////////////

	/**
	 * Updates the mean, standard deviation, and classifier parameters
	 * based on the current expectation of which genes will go through which paths
	 */
	public void updateParams(Treenode root, int nchildnum, double[][] data,int[][] pma, int ndepth) 
	throws Exception
	{
		if (!root.bchange)
		{
			//not allowed to change states params
			//maybe we can get rid of this below
			//since if there is no change in root can assume all children not changeable
			if (root.nextptr[0] != null)
			{
				for (int nchild = 0; nchild< root.numchildren; nchild++)
				{
					updateParams(root.nextptr[nchild],nchild,data,pma, ndepth+1);
				}
			}     
		}
		else
		{

			double doldmean = root.dmean;

			if (root.dPsum > 0)
			{
				root.dmean = root.dEsum/root.dPsum;
			}
			else
			{
				root.dmean = 0;
			}	       

			if (Double.isNaN(root.dmean))
			{
				System.out.println("not a valid mean "+root.dEsum+"\t"+root.dPsum+"\t"+doldmean);
				System.out.println(root.tranC);
				System.out.println("--------------------");
				throw new IllegalArgumentException("not a valid mean");
			}
			if (root.nprime <= 1)
			{
				root.dsigma = DEFAULTSIGMA;
			}
			else
			{
				double dval = (double) root.nprime/(double)(root.nprime - 1)*
				(root.dEsumsq/root.dPsum - Math.pow(root.dEsum/root.dPsum,2));
				if (dval <= 0)
				{
					root.dsigma = DEFAULTSIGMA;
				}
				else
				{
					root.dsigma =Math.sqrt(dval);
				}
			}		


			if (!(root.dsigma > 0))
			{
				System.out.println("## "+root.dsigma+" "+root.nprime+" "+root.dEsumsq +" "+root.dEsum+" "+root.dPsum);
			}

			if (root.nextptr[0] != null)
			{
				if (root.numchildren > 1)
				{
					if (BREGDREM)
					{
						root.binit = true;
					}
					else
					{
						if (root.tranC == null)
						{
							root.tranC = 
								new DREM_FastLogistic2(theInstancesIndex[root.numchildren-2],
										theInstances[root.numchildren-2], theInstancesTFIndex[root.numchildren-2],theInstancesTF[root.numchildren-2],
										ylabels[root.numchildren-2],root.recweight,root.numchildren,numbits);
							root.tranC.setRidge(RIDGE);
						}
						else
						{
							root.tranC.reinit(theInstancesIndex[root.numchildren-2],theInstances[root.numchildren-2],
									theInstancesTFIndex[root.numchildren-2],theInstancesTF[root.numchildren-2],
									ylabels[root.numchildren-2],root.recweight,root.numchildren);
						}

						root.tranC.train();

						root.binit = true;
					}
				}

				double dsum = 0;

				if (BREGDREM)
				{ 
					for (int nchild = 0; nchild< root.numchildren; nchild++)
					{
						if (root.nextptr[nchild].bchange)
						{
							dsum += root.aprobs[nchild];
						}
					}
				}


				for (int nchild = 0; nchild< root.numchildren; nchild++)
				{
					if (BREGDREM)
					{
						if (dsum >0)
						{
							root.ptrans[nchild] = root.aprobs[nchild]/dsum;
							if (BDEBUG)
							{
								System.out.println("probs: "+root.aprobs[nchild]+" "+dsum+" "+root.ptrans[nchild]);
							}
						}
						else
						{
							root.ptrans[nchild] = 0;
						}

						if (Double.isNaN(root.ptrans[nchild]))
							System.out.println("%%%%%% "+root.ptrans[nchild]+"\t"+root.aprobs[nchild]+"\t"+dsum);
					}
					updateParams(root.nextptr[nchild],nchild,data,pma,ndepth+1);
				}
			}
		}
	}

	//////////////////////////////////////////////////////////////////////////////
	/**
	 * Executes the forward algorithm portion of the baum-welch algorithm
	 */
	public void forwardalg(double[] vals,int[] pma, int[] theInstanceIndex,
			double[] theInstance, int ni,double da, Treenode node, 
			int nrec,Treenode parentptr) throws Exception
			{
		//da is trans prob
		//f_k(i) = P(x_1,...,x_i,pi_i=k)
		//f_0(0) = 1, f_k(0)=0 for k>0

		if (node != null)
		{
			//sum_m  P(x_1,...,x_(i-1)|pi_(i-1)=m)a_(mk)*prob(vals[ni])
			//but were only assuming there is one path to each state

			// Not used???
			//double dens;

			if (pma[ni] != 0)
			{
				double dxmudiff = (vals[ni]-node.dmean);
				node.dens = Math.exp(-dxmudiff*dxmudiff/(2*node.dsigma*node.dsigma))/(node.dsigma*StatUtil.TWOPISQRT);
			}
			else
			{
				node.dens = 1;
			}

			if (node.parentptrA == null)
			{
				node.df = node.parent.df*da*node.dens;
			}
			else
			{
				if ((node.nrec != nrec)||(node.niteration != nglobaliteration))
				{
					node.df = parentptr.df*da*node.dens;
					if (node.niteration != nglobaliteration)
						node.niteration = nglobaliteration;
					if (node.nrec != nrec)
						node.nrec = nrec;
				}
				else
				{
					node.df += parentptr.df*da*node.dens;
					node.niteration++;
				}
			}


			double[] ptranlogit;

			if (BREGDREM)
			{
				ptranlogit = node.ptrans;
			}
			else if ((node.numchildren >1)&&(node.binit))
			{
				ptranlogit = node.tranC.distributionForInstance(theInstanceIndex,theInstance,nrec);
			}
			else
			{
				ptranlogit = CONSTANTA[node.numchildren];
			}

			for (int nchild = 0; nchild < node.numchildren; nchild++)
			{
				try
				{	 
					forwardalg(vals,pma, theInstanceIndex,theInstance,ni+1,ptranlogit[nchild],
							node.nextptr[nchild],nrec,node);
				}
				catch (Exception ex)
				{
					ex.printStackTrace();
					System.out.println(ptranlogit.length+" ------abc--"+node.nextptr[0]+"------------- "
							+node.nextptr.length+" "+nchild+"\t**"+node.numchildren+"\t"+node.tranC.numclasses
							+"\t"+node.binit+"\t"+node.dmean+"\t"+node.dsigma);
					System.out.println(node.tranC);
				}
			}
		}
			}

	/////////////////////////////////////////////////////////////////////////////
	/**
	 * Executes the backward portion of the baum-welch algorithm
	 */
	public void backalg(double[] vals,int[] pma, int[] theInstanceIndex,
			double[] theInstance,int ni, Treenode node, int nrec,boolean bforward) throws Exception
			{
		//b_k(i) = P(x_(i+1)...x_L|pi=k)
		//P(x_1,...,x_i,pi_i=k)
		if (node.nextptr[0] == null)
		{
			node.db = 1;
		}
		else
		{ 
			node.db = 0;
			double[] ptranlogit;

			if (BREGDREM)
			{
				ptranlogit = node.ptrans;
			}
			else if ((node.numchildren >1)&&(node.binit))
			{
				ptranlogit = node.tranC.distributionForInstance(theInstanceIndex, theInstance,nrec);
			}
			else
			{
				ptranlogit = CONSTANTA[node.numchildren];
			}

			for (int nchild = 0; nchild < node.numchildren; nchild++)
			{
				Treenode nextnode = node.nextptr[nchild];
				backalg(vals, pma, theInstanceIndex,theInstance,ni+1, nextnode,nrec,bforward);

				double dens;

				if (bforward)
				{
					dens = nextnode.dens;
				}
				else if (pma[ni+1] != 0)
				{
					double dxmudiff = (vals[ni+1]-nextnode.dmean);
					dens = Math.exp(-dxmudiff*dxmudiff/(2*nextnode.dsigma*nextnode.dsigma))/
					(nextnode.dsigma*StatUtil.TWOPISQRT);
				}
				else
				{
					dens = 1;
				}

				if (BDEBUG)
				{
					if (nchild >= ptranlogit.length)
					{
						System.out.println("**********"+nchild+"\t"+ptranlogit.length);
					}
				}

				node.db += ptranlogit[nchild]* nextnode.db *dens;

				if (node.db < 0)
				{
					System.out.println("##"+node.db+"\t"+nextnode.dens);
					for (int i = 0; i < ptranlogit.length; i++)
						System.out.print(ptranlogit[i]+"\t");
					System.out.println();

				}
			}
		}
			}

	////////////////////////////////////////////////////////////////
	/**
	 * Clears any prior assignments of genes to paths through the model
	 */
	public void clearCounts(Treenode ptr)
	{
		if (ptr != null)
		{
			ptr.bInNode = new boolean[theDataSet.data.length];
			ptr.numPath = 0;


			if (bindingSignGene != null)
			{
				ptr.htScore = new Hashtable[numbits];
				ptr.dTotalScore = new double[numbits];


				for (int nbindingindex = 0; nbindingindex < ptr.htScore.length; nbindingindex++)
				{
					ptr.htScore[nbindingindex] = new Hashtable();
				}
			}

			for (int nchild = 0; nchild < ptr.numchildren; nchild++)
			{
				clearCounts(ptr.nextptr[nchild]);
			}
		}
	}

	//////////////////////////////////////////////////////////////////////////////
	/**
	 * Determines the most likely path of each gene through the model
	 */
	public void viterbi(Treenode treeptr) throws Exception
	{
		int[] bestpath;
		bestpath = new int[theDataSet.data[0].length];
		storedbestpath = new int[theDataSet.data.length];
		clearCounts(treeptr);
		treeptr.numPath = theDataSet.data.length;
		if (BDEBUG)
		{
			System.out.print("!!!!!!"+theDataSet.data.length);
			if (bindingValGene !=null)
				System.out.println(" "+bindingValGene.length);
		}

		// Loop through all genes to find the best path for each one
		// one at a time
		for (int nrow = 0; nrow < theDataSet.data.length; nrow++)
		{
			if (bindingValGene != null)
			{
				computevlogit(theDataSet.data[nrow],theDataSet.pmavalues[nrow], bindingValGene[nrow],
						bindingGeneindex[nrow],0,treeptr, bestpath);
			}
			else
			{
				computevlogit(theDataSet.data[nrow],theDataSet.pmavalues[nrow], null,null,
						0,treeptr, bestpath);
			}

			int nsum = 0;
			Treenode tempptr = treeptr;

			tempptr.bInNode[nrow] = true;

			for (int nindex =0; nindex < bestpath.length-1; nindex++)
			{
				tempptr = tempptr.nextptr[bestpath[nindex]];
				tempptr.bInNode[nrow] = true;
				tempptr.numPath++;

				if (bindingSignGene != null)
				{
					int ntfindex = 0;
					for (int nbit = 0; nbit < numbits; nbit++)
					{
						while ((ntfindex < bindingGeneindex[nrow].length)&&
								(nbit > bindingGeneindex[nrow][ntfindex]))
						{
							ntfindex++;
						}

						Double objMap;
						if ((ntfindex < bindingGeneindex[nrow].length)&&
								(nbit == bindingGeneindex[nrow][ntfindex]))
						{
							// Binding signs are used here instead of binding values
							// because dTotalScore and htScore are used to
							// count genes along a path ???
							tempptr.dTotalScore[nbit] += bindingSignGene[nrow][ntfindex];
							objMap =new Double(bindingSignGene[nrow][ntfindex]);
						}
						else
						{
							objMap = new Double(0.0);
						}

						Integer objCount = (Integer) tempptr.htScore[nbit].get(objMap);
						if (objCount != null)
						{
							tempptr.htScore[nbit].put(objMap, new Integer(objCount.intValue()+1));
						}
						else
						{
							tempptr.htScore[nbit].put(objMap, new Integer(1));
						}
					}
				}
				nsum += (int) Math.pow(nmaxchild,bestpath.length-nindex-1)*bestpath[nindex];
			}
			storedbestpath[nrow] = nsum;
		}
	}

	//////////////////////////////////////////////////////////////////////////////////
	/**
	 * Computes association scores for the transcription factors overall on paths and
	 * on paths out of splits conditional on the set of genes going into the split
	 */
	public void computeStats(Treenode treeptr, Treenode rootptr)
	{
		if (treeptr != null)
		{
			treeptr.dscore = new double[numbits][treeptr.numchildren];
			treeptr.ddiff = new double[numbits][treeptr.numchildren];
			treeptr.dpvals = new double[numbits][treeptr.numchildren];

			treeptr.ncountvals = new int[numbits][treeptr.numchildren][dBindingSigns.length];

			if (treeptr.numchildren >= 3)
			{
				treeptr.nothersum = new int[numbits][treeptr.numchildren][dBindingSigns.length];
				treeptr.dscoreother = new double[numbits][treeptr.numchildren];
			}

			int[][] ncountgrid = new int[2][dBindingSigns.length];
			treeptr.ncountTotals = new int[numbits][dBindingSigns.length];

			treeptr.dTFActivityScore = new double[numbits];

			// Create a TreeSet that stores important TFs at this node
			// based on activity scores instead of pvalues
			treeptr.tsSigTFActivity = new TreeSet(new SigTFRecv2Compare());
			treeptr.tsSigTF = new TreeSet(new SigTFRecv2Compare());
			for (int ntf = 0; ntf < numbits; ntf++)
			{
				double dnumpaths = 0;
				double dscoretotal = 0;
				for (int nindex = 0; nindex < treeptr.ncountTotals[ntf].length; nindex++)
				{
					treeptr.ncountTotals[ntf][nindex] = 0;
				}

				for (int nchild = 0; nchild < treeptr.numchildren; nchild++)
				{
					if (treeptr.nextptr[nchild].htScore == null)
					{
						for (int nvalindex = 0; nvalindex < treeptr.ncountvals[ntf][nchild].length; nvalindex++)
						{
							treeptr.ncountvals[ntf][nchild][nvalindex] = 0;
						}
					}
					else
					{
						for (int nel = 0; nel < dBindingSigns.length; nel++)
						{
							Object key = new Double(dBindingSigns[nel]);
							Object objcount = treeptr.nextptr[nchild].htScore[ntf].get(key);

							if (objcount != null)
							{  
								int ntempval =((Integer) objcount).intValue();
								treeptr.ncountvals[ntf][nchild][nel] = ntempval;
								treeptr.ncountTotals[ntf][nel] += ntempval;  
							}
							else
							{
								treeptr.ncountvals[ntf][nchild][nel] = 0;
							}
						}
					}

					treeptr.dscore[ntf][nchild] = 
						treeptr.nextptr[nchild].dTotalScore[ntf]/treeptr.nextptr[nchild].numPath;
					dnumpaths += treeptr.nextptr[nchild].numPath;
					dscoretotal += treeptr.nextptr[nchild].dTotalScore[ntf];
				}

				if (treeptr.numchildren == 2)
				{
					treeptr.ddiff[ntf][0] =  treeptr.dscore[ntf][0]- treeptr.dscore[ntf][1];
					treeptr.dpvals[ntf][0] = computepval( treeptr.ncountvals[ntf],  Math.abs(treeptr.ddiff[ntf][0]));
					treeptr.tsSigTF.add(new SigTFRecv2(tfNames[ntf],treeptr.dpvals[ntf][0]));

				}
				else  if (treeptr.numchildren >= 3)
				{
					for (int nchild = 0; nchild < treeptr.numchildren; nchild++)
					{
						//divide by 0 need to check about that
						double ddiff = (dnumpaths-treeptr.nextptr[nchild].numPath);
						double davgother;

						if (ddiff > 0)
						{
							davgother =  (dscoretotal- treeptr.nextptr[nchild].dTotalScore[ntf])/ddiff;
						}
						else
						{
							davgother = 0;
						}
						treeptr.dscoreother[ntf][nchild] = davgother;

						treeptr.ddiff[ntf][nchild] = treeptr.dscore[ntf][nchild]-davgother;

						for (int nindex = 0; nindex < treeptr.ncountTotals[ntf].length; nindex++)
						{
							ncountgrid[0][nindex] = treeptr.ncountvals[ntf][nchild][nindex];
							ncountgrid[1][nindex] = treeptr.ncountTotals[ntf][nindex] - ncountgrid[0][nindex];
							treeptr.nothersum[ntf][nchild][nindex] =  ncountgrid[1][nindex];
						}

						treeptr.dpvals[ntf][nchild] = computepval(ncountgrid, Math.abs(treeptr.ddiff[ntf][nchild]));
					}

					for (int nchild = 0; nchild < treeptr.numchildren; nchild++)
					{
						treeptr.tsSigTF.add(
								new SigTFRecv2(tfNames[ntf],treeptr.dpvals[ntf][treeptr.orderA[nchild]]));						 
					}
				}

				
				// Activity scores are only calculated if there is a split at this node
				if(treeptr.numchildren > 1)
				{
					/**The number of genes on a particular path out of a split regulated by a TF.
					 * Regulated means the interaction value is not 0.  Dimension:<br>
					 * [num paths out of split]*/
					int[] nRegThisPath = new int[treeptr.numchildren];
					/**The number of genes on all paths except the given path
					 * out of a split that are regulated by a TF.
					 * Regulated means the interaction value is not 0.  Dimension:<br>
					 * [num paths out of split]*/
					int[] nRegOtherPaths = new int[treeptr.numchildren];
					/**The number of genes on all paths
					 * out of a split that are regulated by a TF.
					 * Regulated means the interaction value is not 0.*/
					int nRegAllPaths = 0;

					// Initialize nRegThisPath
					for (int nchild = 0; nchild < treeptr.numchildren; nchild++)
					{
						nRegThisPath[nchild] = 0;

					}

					// Populate nRegThisPath and nRegAllPaths
					for (int nchild = 0; nchild < treeptr.numchildren; nchild++)
					{
						for (int nel = 0; nel < dBindingSigns.length; nel++)
						{
							// Only count genes the TF regulates, which are those
							// genes with an interaction value != 0
							if(dBindingSigns[nel] != 0)
							{
								nRegThisPath[nchild] += treeptr.ncountvals[ntf][nchild][nel];
								nRegAllPaths += treeptr.ncountvals[ntf][nchild][nel];
							}
						}
					}

					// Populate nRegOtherPaths
					for (int nchild = 0; nchild < treeptr.numchildren; nchild++)
					{
						nRegOtherPaths[nchild] = nRegAllPaths - nRegThisPath[nchild];
					}

					// Now determine which path is the primary path out of the split.
					// The primary path is the path with the most regulated
					// genes.  In the case of the tie, the path with the fewest genes
					// overall out of the split is the primary path.
					int nMaxPathVal = -1;
					for (int nchild = 0; nchild < treeptr.numchildren; nchild++)
					{
						nMaxPathVal = Math.max(nMaxPathVal, nRegThisPath[nchild]);
					}

					// See how many paths have this max value
					boolean[] bMaxPaths = new boolean[treeptr.numchildren];
					int nMaxPathCount = 0;
					for (int nchild = 0; nchild < treeptr.numchildren; nchild++)
					{
						if(nRegThisPath[nchild] == nMaxPathVal)
						{
							bMaxPaths[nchild] = true;
							nMaxPathCount++;
						}
						else
						{
							bMaxPaths[nchild] = false;
						}
					}

					// The index of the primay path
					int nPrimPathInd = -1;
					if(nMaxPathCount == 1)
					{
						for (int nchild = 0; nchild < treeptr.numchildren; nchild++)
						{
							if(bMaxPaths[nchild])
							{
								nPrimPathInd = nchild;
							}
						}
					}
					else
					{
						// There was a tie so look at all paths out of the split
						// there were involved in the tie
						// to find the path with the fewest genes
						// If there is another tie, it won't affect the score
						// calculation so pick the first one found
						int nMinPathVal = Integer.MAX_VALUE;

						for (int nchild = 0; nchild < treeptr.numchildren; nchild++)
						{
							if(bMaxPaths[nchild] && treeptr.nextptr[nchild].numPath < nMinPathVal)
							{
								nMinPathVal = treeptr.nextptr[nchild].numPath;
								nPrimPathInd = nchild;							
							}
						}
					}

					// Calculate the number of genes on the paths that are not
					// the primary path
					int nNumOtherPaths = 0;
					for (int nchild = 0; nchild < treeptr.numchildren; nchild++)
					{
						if(nPrimPathInd != nchild)
						{
							nNumOtherPaths += treeptr.nextptr[nchild].numPath;
						}
					}

					// The formula for a TF's activity score is:
					// score = numer / denom
					// where numer = (prob binding is functional)^(num genes TF regulates on primary path)
					//                 * (prob binding is not functional)^(num genes TF regulates on other paths)
					//                 * (TF activity prior)
					// and denom = (ratio of all genes on TFs primary path)^(num genes TF regulates on primary path)
					//               * (ratio of all genes on other paths)^(num genes TF regulates on other paths)
					//               * (1 - TF activity prior)
//					double scoreNumer = Math.pow(dProbBindingFunctional, nRegThisPath[nPrimPathInd]);
//					scoreNumer *= Math.pow(1 - dProbBindingFunctional, nRegOtherPaths[nPrimPathInd]);
//					scoreNumer *= dTFActivityPrior[ntf];
//					double scoreDenom = Math.pow(treeptr.nextptr[nPrimPathInd].numPath / ((double) treeptr.numPath), nRegThisPath[nPrimPathInd]);
//					scoreDenom *= Math.pow(nNumOtherPaths / ((double) treeptr.numPath), nRegOtherPaths[nPrimPathInd]);
//					scoreDenom *= (1 - dTFActivityPrior[ntf]);
//					treeptr.dTFActivityScore[ntf] = scoreNumer / scoreDenom;
					
					// Calculations need to be done in logspace, otherwise if there are too
					// many bound genes the denominator can go to 0 or the score will
					// be imprecise
					double scoreNumer = Math.log(dProbBindingFunctional) * nRegThisPath[nPrimPathInd];
					scoreNumer += Math.log(1 - dProbBindingFunctional) * nRegOtherPaths[nPrimPathInd];
					scoreNumer += Math.log(dTFActivityPrior[ntf]);
					double scoreDenom = Math.log(treeptr.nextptr[nPrimPathInd].numPath / ((double) treeptr.numPath)) * nRegThisPath[nPrimPathInd];
					scoreDenom += Math.log(nNumOtherPaths / ((double) treeptr.numPath)) * nRegOtherPaths[nPrimPathInd];
					scoreDenom += Math.log(1 - dTFActivityPrior[ntf]);
					// Take out of logspace for display and output
					treeptr.dTFActivityScore[ntf] = Math.exp(scoreNumer - scoreDenom);
					
					// Max activity score is either the new TF activity score at this node
					// or the pre-existing max activity score
					dMaxTFActivityScore[ntf] = Math.max(treeptr.dTFActivityScore[ntf], dMaxTFActivityScore[ntf]);
					
					// Add the activity score at this node to the TreeSet
					// Instead of storing the score directly, store it in a form that resembles
					// a pvalue so that the existing signifance score threshold and sorting
					// code work
					double pseudoPvalue = Double.POSITIVE_INFINITY;
					if (treeptr.dTFActivityScore[ntf] != 0)
					{
						pseudoPvalue = 1 / treeptr.dTFActivityScore[ntf];
					}
					// Also with the TF name, store the index of the primary path for that
					// TF coming out of the split
					// TODO store the path index separately
					// Use the reverse order map to go from the actual index of the child node (path)
					// to the sorted order of that node based on the mean expression levels
					treeptr.tsSigTFActivity.add(new SigTFRecv2(tfNames[ntf] + " [" + (treeptr.revOrderA[nPrimPathInd] + 1) + "]",pseudoPvalue));
				}
				else // treeptr.numchildren is 0 or 1 so there is no split to calculate the score
				{
					treeptr.dTFActivityScore[ntf] = 0;
					// No need to add a SigTFRecv2 object to the tsSigTFActivity TreeSet
					// if we are not at a split
				}

			} // end loop through all TFs


			int numrowstable;
			int nsize = hsUniqueInput.size()-1;

			numrowstable = tfNames.length * nsize;

			boolean bsplit;
			if (treeptr.numchildren >= 2)
			{
				treeptr.dpvalEdgeSplit= new double[treeptr.numchildren][numrowstable];
				treeptr.ddiffEdgeSplit = new double[treeptr.numchildren][numrowstable];
				treeptr.dexpectEdgeSplit = new double[treeptr.numchildren][numrowstable];
				treeptr.tsSigTFEdgeSplit = new TreeSet[treeptr.numchildren];
				treeptr.dratioSplit = new double[treeptr.numchildren][numrowstable];

				for (int nindex= 0; nindex < treeptr.numchildren; nindex++)
				{
					treeptr.tsSigTFEdgeSplit[nindex] = new TreeSet(new SigTFRecv2Compare());
				}
				bsplit = true;
			}
			else
			{
				bsplit = false;
			}

			treeptr.dpvalEdgeFull= new double[treeptr.numchildren][0];
			treeptr.ddiffEdgeFull = new double[treeptr.numchildren][0];
			treeptr.dexpectEdgeFull = new double[treeptr.numchildren][0];
			treeptr.tsSigTFFull = new TreeSet[treeptr.numchildren];
			int nsumpathfull;

			nsumpathfull = rootptr.numPath;

			if ((treeptr.numchildren ==1) && (treeptr != rootptr))
			{
				boolean bfound = false;
				int nindex = 0;
				while (!bfound)
				{
					if (treeptr.parent.nextptr[nindex] == treeptr)
					{
						bfound = true;
						if (BDEBUG)
						{
							System.out.println(treeptr+"\t"+treeptr.parent+
									"\t"+treeptr.parentptrA+"\t"+treeptr.parent.dpvalEdgeFull);
						}

						if (treeptr.parentptrA == null)
						{
							treeptr.dpvalEdgeFull[0]= treeptr.parent.dpvalEdgeFull[nindex];
							treeptr.ddiffEdgeFull[0] =  treeptr.parent.ddiffEdgeFull[nindex];
							treeptr.dexpectEdgeFull[0] =  treeptr.parent.dexpectEdgeFull[nindex];
							treeptr.tsSigTFFull[0] = treeptr.parent.tsSigTFFull[nindex];
						}
					}
					else
					{
						nindex++;
					}
				}
			}
			else
			{
				treeptr.dpvalEdgeFull= new double[treeptr.numchildren][numrowstable];
				treeptr.ddiffEdgeFull = new double[treeptr.numchildren][numrowstable];
				treeptr.dexpectEdgeFull = new double[treeptr.numchildren][numrowstable];
				for (int nindex= 0; nindex < treeptr.numchildren; nindex++)
				{
					treeptr.tsSigTFFull[nindex] = new TreeSet(new SigTFRecv2Compare());
				}

				for (int nchild = 0; nchild < treeptr.numchildren; nchild++)
				{
					int nrowindex = 0;
					for (int nrow = 0; nrow < tfNames.length; nrow++)
					{
						int nsumthispath = 0;
						int[] dvals = treeptr.ncountvals[nrow][nchild];
						nsumthispath = treeptr.nextptr[nchild].numPath;

						for (int nel = 0; nel < dBindingSigns.length; nel++)
						{  
							if (dBindingSigns[nel] != 0)
							{
								int nsumfullel = filteredClassifier.nBaseCount[nrow][filteredClassifier.getFeatureIndex(dBindingSigns[nel])];
								double dexpectfull = (nsumthispath* nsumfullel)/(double) ntotalcombined;
								treeptr.dexpectEdgeFull[nchild][nrowindex] =dexpectfull;

								double ddiffull =  treeptr.ncountvals[nrow][nchild][nel] - dexpectfull;
								treeptr.ddiffEdgeFull[nchild][nrowindex] = ddiffull;

								treeptr.dpvalEdgeFull[nchild][nrowindex] = 
									StatUtil.hypergeometrictail(treeptr.ncountvals[nrow][nchild][nel]-1, 
											nsumfullel, ntotalcombined-nsumfullel,nsumthispath);
								//assuming 0 always unique input

								if (hsUniqueInput.size()==2)
								{
									treeptr.tsSigTFFull[nchild].add(new SigTFRecv2(tfNames[nrow],
											treeptr.dpvalEdgeFull[nchild][nrowindex]));
								}
								else
								{
									treeptr.tsSigTFFull[nchild].add(
											new SigTFRecv2(tfNames[nrow]+" "+(int)dBindingSigns[nel],
													treeptr.dpvalEdgeFull[nchild][nrowindex]));
								}

								if (bsplit)
								{
									int nsumthisel = treeptr.ncountTotals[nrow][nel];

									int nsumthisparentpath = treeptr.numPath;

									double dexpect = (nsumthispath* nsumthisel)/(double) nsumthisparentpath;
									treeptr.dexpectEdgeSplit[nchild][nrowindex] = dexpect;
									double ddiff =  treeptr.ncountvals[nrow][nchild][nel] - dexpect;
									treeptr.ddiffEdgeSplit[nchild][nrowindex] = ddiff;
									treeptr.dpvalEdgeSplit[nchild][nrowindex] = 
										StatUtil.hypergeometrictail(treeptr.ncountvals[nrow][nchild][nel]-1, 
												nsumthisel, nsumthisparentpath-nsumthisel,nsumthispath);

									treeptr.dratioSplit[nchild][nrowindex] = (double)treeptr.ncountvals[nrow][nchild][nel]/treeptr.ncountTotals[nrow][nel];

									if (hsUniqueInput.size()==2)
									{
										treeptr.tsSigTFEdgeSplit[nchild].add(
												new SigTFRecv2(tfNames[nrow], treeptr.dpvalEdgeSplit[nchild][nrowindex],
														treeptr.dratioSplit[nchild][nrowindex]));
									}
									else
									{
										treeptr.tsSigTFEdgeSplit[nchild].add(
												new SigTFRecv2(tfNames[nrow]+" "+(int)dBindingSigns[nel],
														treeptr.dpvalEdgeSplit[nchild][nrowindex],
														treeptr.dratioSplit[nchild][nrowindex]));
									}
								}
								nrowindex++;
							}
						}
					}
				}
			}

			for (int nchild = 0; nchild < treeptr.numchildren; nchild++)
			{
				computeStats(treeptr.nextptr[nchild],rootptr);
			}
		}
	}

	/////////////////////////////////////////////////////////////////////////////////
	/**
	 * Record used in the ranking of transcription factors
	 */
	static class SigTFRecv2 
	{
		String szname;
		double dpval;
		double dpercent;

		SigTFRecv2(String szname, double dpval,double dpercent)
		{
			this.szname = szname;
			this.dpval = dpval;
			this.dpercent = dpercent;
		}

		SigTFRecv2(String szname, double dpval)
		{
			this.szname = szname;
			this.dpval = dpval;
			this.dpercent = -1;
		}
	}

	/**
	 * Compares TFs first based on the value first lower dpval, then great dpercent, then name
	 */
	static class SigTFRecv2Compare implements Comparator 
	{
		public int compare(Object c1, Object c2)
		{
			SigTFRecv2 cr1 = (SigTFRecv2) c1;
			SigTFRecv2 cr2 = (SigTFRecv2) c2;

			if (cr1.dpval < cr2.dpval)
				return -1;
			else if (cr1.dpval > cr2.dpval)
				return 1;
			else if (cr1.dpercent > cr2.dpercent)
				return -1;
			else if (cr1.dpercent < cr2.dpercent)
				return 1;
			else
				return cr1.szname.compareTo(cr2.szname);

		}
	}

	/////////////////////////////////////////////////////////////////////////////////////////
	/**
	 * Computes the fraction of combinations for vals that could lead to a greater difference
	 * than ddiff
	 */
	public double computepval(int[][] vals, double diff)
	{
		// Not used ???
//		double dcol0 = 0;
//		double dcol1 = 0;

		// The parameter vals stores the number of genes regulated by a
		// particular TF on a given path with the given interaction value.
		// Dimensions: [num paths out of split][num distinct binding interaction values] ???
		// Assumes 2 paths out of the split.  For more, the paths must
		// be grouped so that all but one path form a single group ???

		// vals[0].length is the number of distinct binding interaction values ???
		int[] rowSum = new int[vals[0].length];
		double dtimes = 0;
		int nsum0 = 0;
		int nsum1 = 0;

		for (int nindex = 0; nindex < vals[0].length; nindex++)
		{
			// nsum0 stores the number of genes that are on the
			// path at index 0 ???
			nsum0 += vals[0][nindex];
			// nsum1 is the same as nsum0 for the path at index 1 ???
			nsum1 += vals[1][nindex];
			// rowSum stores the number of genes that have the same binding
			// value regardless of what path they are assigned to ???
			rowSum[nindex] = vals[0][nindex] + vals[1][nindex];
		}

		double dpvalsum = 0;
		int[] ncol0 = new int[vals[0].length];

		// Two distinct binding values  (e.g. 0,1)
		if (ncol0.length == 2)
		{
			// Iterate through all possible assignments of genes with binding value
			// at rowSum[0] to path with index 0 ???
			for (ncol0[0] = 0; ncol0[0] <= Math.min(rowSum[0],nsum0);  ncol0[0]++)
			{
				double dinnersum = 0;
				// ncol0[1] is the number of genes with the other binding value
				// assigned to the path with index 0 ???
				ncol0[1] = nsum0-ncol0[0];

				// Make sure there are more genes assigned to the path than
				// there are genes to assign
				if (ncol0[1] <= rowSum[1])
				{
					double dweightsum0 = 0;
					double dweightsum1 = 0;
					for (int i = 0; i < ncol0.length; i++)
					{
						dweightsum0 += dBindingSigns[i]*ncol0[i]; 
						dweightsum1 += dBindingSigns[i]*(rowSum[i]-ncol0[i]); 
					}
					double dcurrdiff = Math.abs(dweightsum0/nsum0-dweightsum1/nsum1);
					if (dcurrdiff >= diff)
					{
						dtimes++;
						dinnersum++;
					}
				}
				dpvalsum += StatUtil.hypergeometric(ncol0[0], nsum0, nsum1, rowSum[0])* dinnersum;
			}
		}
		// Three distinct binding values (e.g. -1,0,1)
		else if (ncol0.length==3)
		{
			for (ncol0[0] = 0; ncol0[0] <= Math.min(rowSum[0],nsum0);  ncol0[0]++)
			{
				double dinnersum = 0;
				int nremainder0 =  nsum0-ncol0[0];
				int nremainder1 = nsum1 - (rowSum[0]- ncol0[0]);
				for (ncol0[1] = 0; ncol0[1] <= Math.min(rowSum[1],nremainder0); ncol0[1]++)
				{
					ncol0[2] = nsum0;
					for (int i = 0; i < ncol0.length - 1; i++)
					{
						ncol0[2] -= ncol0[i];
					}
					if (ncol0[2] <= rowSum[2])
					{
						double dweightsum0 = 0;
						double dweightsum1 = 0;
						for (int i = 0; i < ncol0.length; i++)
						{
							dweightsum0 += dBindingSigns[i]*ncol0[i]; 
							dweightsum1 += dBindingSigns[i]*(rowSum[i]-ncol0[i]); 
						}
						double dcurrdiff = Math.abs(dweightsum0/nsum0-dweightsum1/nsum1);
						if (dcurrdiff >= diff)
						{
							dtimes++;
							dinnersum += StatUtil.hypergeometric(ncol0[1], nremainder0, nremainder1,rowSum[1]);
						}
					}
				}
				dpvalsum += StatUtil.hypergeometric(ncol0[0], nsum0, nsum1, rowSum[0])* dinnersum;
			}
		}
		else if (ncol0.length == 1)
		{
			return 1;
		}
		else
		{
			throw new IllegalArgumentException("can only handle two or three unique inputs");
		}

		return dpvalsum;
	}

	/////////////////////////////////////////////////////////////////////////

	/**
	 * Recursively determines the best path through the model and its likelihood
	 * for a single gene
	 */
	public double computevlogit(double[] vals,int[] pma,
			double[] theInstance,int[] theInstanceIndex, 
			int ndepth,Treenode node,int[] bestpath) throws Exception
			{

		double dmax = Double.NEGATIVE_INFINITY;
		int bestchild = 0;
		double dlogout = 0;
		double dval = 0;

		if (node.dsigma > 0)
		{
			if (pma[ndepth] != 0)
			{

				if (ndepth == 0)
				{
					dlogout = 0;
				}
				else
				{
					dlogout = Math.log(StatUtil.normaldensity(vals[ndepth], node.dmean, node.dsigma));
				}
			}
			else
			{
				dlogout = 0;
			}
		}

		// The base case is when there are no states to move to so return
		// only the emission probability ???
		if (node.numchildren == 0)
		{
			return dlogout;
		}

		int[] currbestpath = new int[bestpath.length];

		// log (a) + c(log(b_1) + log(b_2) + ... + log(b_n))
		//= log(a*(b_1*b_2*...*b_n)^c) 

		double [] ptranlogit;

		// Determine the probability of the gene moving to each of the
		// possible subsequent states
		if (BREGDREM)
		{
			ptranlogit = node.ptrans;
		}
		else if ((node.numchildren >1)&&(node.binit))
		{
			ptranlogit = node.tranC.distributionForInstance(theInstanceIndex,theInstance,-1);
		}
		else
		{
			ptranlogit = CONSTANTA[node.numchildren];
		}

		for (int nchild = 0; nchild < node.numchildren; nchild++)
		{
			dval = dlogout;
			if (node.nextptr[nchild] != null)
			{
				double dvi = computevlogit(vals,pma,theInstance,theInstanceIndex,
						ndepth+1,node.nextptr[nchild],bestpath);
				if (ptranlogit[nchild]==0)
				{
					dval += Math.log(MINPROB);
				}
				else
				{
					dval += Math.log(ptranlogit[nchild])+dvi;
				}
			}

			if (dval > dmax)
			{
				dmax = dval;
				bestchild = nchild;
				for (int nindex = 0; nindex <bestpath.length; nindex++)
				{
					currbestpath[nindex] = bestpath[nindex];
				}
			}
		}

		for (int nindex = 0; nindex <bestpath.length; nindex++)
		{
			bestpath[nindex] = currbestpath[nindex];
		}
		bestpath[ndepth] = bestchild;

		return dmax;
			}

	//////////////////////////////////////////////////////////////////////////
	/**
	 * Evaluates the likelihood of the viterbi assignments of the holdoutdata 
	 * under the current settings of the model parameters
	 */
	public double testhmmVHoldOut(Treenode treehmm) throws Exception
	{
		double dlike = 0;
		int[] bestpath = new int[traindata[0].length];
		for (int nrow = 0; nrow <holdoutdata.length; nrow++)
		{
			double[] theInstance = holdoutVal[nrow];
			int[] theInstanceindex = holdoutIndex[nrow];

			double dbest = computevlogit(holdoutdata[nrow],holdoutpma[nrow],theInstance,
					theInstanceindex, 0,treehmm, bestpath);
			{
				dlike += dbest;
			}
		}
		return dlike;
	}

	//////////////////////////////////////////////////////////////////////////
	/**
	 * Evaluates the likelihood of the viterbi assignments of the testdata 
	 * under the current settings of the model parameters
	 */
	public double testhmmV(Treenode treehmm) throws Exception
	{
		double dlike = 0;
		int[] bestpath = new int[traindata[0].length];
		for (int nrow = 0; nrow <testdata.length; nrow++)
		{

			double[] theInstance = testVal[nrow];
			int[] theInstanceIndex = testIndex[nrow];
			double dbest = computevlogit(testdata[nrow],testpma[nrow],theInstance,
					theInstanceIndex, 0,treehmm, bestpath);
			dlike += dbest;
		}
		return dlike;
	}


	//////////////////////////////////////////////////////////////////////////
	/**
	 * Evaluates the likelihood of the holdoutdata under the current settings of the
	 * model parameters
	 */
	public double testhmmHoldOut(Treenode treehmm) throws Exception
	{
		double dlike = 0;
		for (int nrow = 0; nrow <holdoutdata.length; nrow++)
		{
			backalg(holdoutdata[nrow],holdoutpma[nrow],holdoutIndex[nrow],
					holdoutVal[nrow], 0, treehmm,nrow,false);
			if (treehmm.db <= 0)
			{
				dlike += Math.log(MINPROB);
				System.out.println("B warning row "+nrow+" is "+treehmm.db);
			}
			else
			{
				dlike += Math.log(treehmm.db);
			}
		}
		return dlike;
	}

	//////////////////////////////////////////////////////////////////////////
	/**
	 * Evaluates the likelihood of the testdata under the current settings of the
	 * model parameters
	 */
	public double testhmm(Treenode treehmm) throws Exception
	{
		double dlike = 0;
		for (int nrow = 0; nrow <testdata.length; nrow++)
		{
			backalg(testdata[nrow],testpma[nrow],testIndex[nrow],
					testVal[nrow],0, treehmm,nrow,false);

			if (treehmm.db <= 0)
			{
				dlike += Math.log(MINPROB);
			}
			else
			{
				dlike += Math.log(treehmm.db);
			}
		}

		return dlike;
	}


	/**
	 * Returns the number of nodes in the tree pointed by root for which its
	 * ncurrtime field does not equal the ncurrtime parameter
	 */
	public int countNodes(Treenode root, int ncurrtime)
	{

		if (root == null)
		{
			return 0;
		}
		else 
		{
			int nsum;
			if (root.ncurrtime != ncurrtime)
			{
				root.ncurrtime = ncurrtime;
				nsum = 1;
			}
			else
			{
				nsum = 0;
			}

			for (int nchild = 0; nchild < root.numchildren; nchild++)
			{
				nsum += countNodes(root.nextptr[nchild],ncurrtime);
			}

			return nsum;
		}
	}


	//////////////////////////////////////////////////////////////////////////
	/**
	 * Implements a viterbi training method for the parameters with the
	 * current model structure
	 */
	public double trainhmmV(Treenode treehmm, boolean bpruneexempt) throws Exception
	{

		double dlike=0;
		double doldlike;

		//first time through we will need to initialize the hmm without the classifiers
		//then we will need to train the classifiers

		initAEV(treehmm);
		int[] bestpath = new int[traindata[0].length];

		double[] ptranlogit;

		for (int nrow = 0; nrow <traindata.length; nrow++)
		{
			double[] theInstance = trainVal[nrow];
			int[] theInstanceIndex = trainIndex[nrow];

			if (BREGDREM)
			{
				ptranlogit = treehmm.ptrans;
			}
			else if ((treehmm.numchildren >1)&&(treehmm.binit))
			{
				ptranlogit = treehmm.tranC.distributionForInstance(theInstanceIndex,theInstance,-1);
			}
			else
			{
				ptranlogit = CONSTANTA[treehmm.numchildren];
			}

			double dbest = computevlogit(traindata[nrow],trainpma[nrow],theInstance,
					theInstanceIndex, 0,treehmm, bestpath);

			instanceAEV(treehmm,traindata[nrow],trainpma[nrow],trainSign[nrow],
					bestpath,Math.exp(dbest),nrow);

			dlike += dbest;

		}
		updateParams(treehmm,0,traindata,trainpma,0);

		if (BEQUALSTD)
		{
			averageChildrenSigmas(treehmm);
		}

		int ncount = 0;

		//very first time going to assign labels 
		double dpredictitr;

		do
		{
			//train while there is still a big enough improvement
			doldlike = dlike;
			initAEV(treehmm);
			dlike = 0;

			for (int nrow = 0; nrow <traindata.length; nrow++)
			{
				double[] theInstance = trainVal[nrow];
				int[] theInstanceIndex = trainIndex[nrow];

				if (BREGDREM)
				{
					ptranlogit = treehmm.ptrans;
				}
				else if ((treehmm.numchildren >1)&&(treehmm.binit))
				{
					ptranlogit = treehmm.tranC.distributionForInstance(theInstanceIndex,theInstance,nrow);
				}
				else
				{
					ptranlogit = CONSTANTA[treehmm.numchildren];
				}

				double dbest = computevlogit(traindata[nrow],trainpma[nrow],theInstance,
						theInstanceIndex, 0,treehmm, bestpath);


				instanceAEV(treehmm,traindata[nrow],trainpma[nrow],trainSign[nrow],
						bestpath,Math.exp(dbest),nrow);

				dlike += dbest;

				if (Double.isNaN(dlike))
				{
					traverse(treehmm,0,true);

					System.out.println(nrow+"\t");
					for (int na = 0; na < trainSign[nrow].length; na++)
					{
						System.out.print(trainSign[nrow][na]+"\t");
					}
					System.out.println();
					if (!BREGDREM)
					{
						for (int na = 0; na < treehmm.tranC.dcoeff.length; na++)
						{
							System.out.print(treehmm.tranC.dcoeff[na]+"\t");
						}
					}
					System.out.println();
					System.out.println("treehmm.db = "+treehmm.db);
					throw new Exception();
				}
			}

			updateParams(treehmm,0,traindata,trainpma,0);
			if (BEQUALSTD)
			{
				averageChildrenSigmas(treehmm);
			}

			ncount++;
			dpredictitr= (dbesttrainlike-dlike)/(dlike-doldlike);
		} 
		while (((dlike-doldlike)/Math.abs(dlike) > BEPSILON)
				&&((bpruneexempt)||(dpredictitr<MAXFUTUREITR)));


		if (BDEBUG)
		{
			System.out.println("Likelihood: "+dlike);
		}

		dtrainlike = dlike;   

		return (dlike);
	}

	//////////////////////////////////////////////////////////////////////////
	/**
	 * Implements a baum-welch method for training the parameters with the
	 * current model structure
	 */
	public double trainhmm(Treenode treehmm, boolean bpruneexempt) throws Exception
	{

		double dlike=0;
		double doldlike;


		//first time through we will need to initialize the hmm without the classifiers
		//then we will need to train the classifiers

		initAE(treehmm);
		nglobaliteration = 1;
		dlike = 0;
		double[] ptranlogit;

		traverse(treehmm,0,false);
		for (int nrow = 0; nrow <traindata.length; nrow++)
		{
			treehmm.df = 1;  

			double[] theInstance = trainVal[nrow];
			int[] theInstanceIndex = trainIndex[nrow];

			if (BREGDREM)
			{
				ptranlogit = treehmm.ptrans;
			}
			else if ((treehmm.numchildren >1)&&(treehmm.binit))
			{
				ptranlogit = treehmm.tranC.distributionForInstance(theInstanceIndex,theInstance,-1);
			}
			else
			{
				ptranlogit = CONSTANTA[treehmm.numchildren];
			}

			for (int nchild =0;  nchild < treehmm.numchildren; nchild++)
			{
				forwardalg(traindata[nrow],trainpma[nrow],theInstanceIndex,theInstance,
						1,ptranlogit[nchild],treehmm.nextptr[nchild],nrow,treehmm);
			}    
			backalg(traindata[nrow],trainpma[nrow],theInstanceIndex,theInstance, 0, treehmm,nrow,true);
			double dpj =  treehmm.db;
			if (dpj == 0) 
			{
				dlike += Math.log(MINPROB);
			}
			else
			{
				instanceAE(treehmm,0,traindata[nrow],trainpma[nrow],
						trainSign[nrow],dpj,nrow);

				dlike += Math.log(treehmm.db);
			}
		}

		traverse(treehmm,0,false);
		updateParams(treehmm,0,traindata,trainpma,0);

		traverse(treehmm,0,false);
		if (BEQUALSTD)
		{
			averageChildrenSigmas(treehmm);
		}

		int ncount = 0;

		//very first time going to assign labels 
		double dpredictitr;

		do
		{
			nglobaliteration++;
			if (BDEBUG)
			{
				System.out.println("Global iteration is now "+nglobaliteration);
			}

			//train while there is still a big enough improvement
			doldlike = dlike;
			initAE(treehmm);
			dlike = 0;

			for (int nrow = 0; nrow <traindata.length; nrow++)
			{
				double[] theInstance = trainVal[nrow];
				int[] theInstanceIndex = trainIndex[nrow];

				if (BREGDREM)
				{
					ptranlogit = treehmm.ptrans;
				}
				else if ((treehmm.numchildren >1)&&(treehmm.binit))
				{
					ptranlogit = treehmm.tranC.distributionForInstance(theInstanceIndex,theInstance,nrow);
				}
				else
				{
					ptranlogit = CONSTANTA[treehmm.numchildren];
				}

				treehmm.df = 1;  
				for (int nchild =0;  nchild < treehmm.numchildren; nchild++)
				{
					forwardalg(traindata[nrow],trainpma[nrow],theInstanceIndex,theInstance,
							1,ptranlogit[nchild],treehmm.nextptr[nchild],nrow,treehmm);
				}

				backalg(traindata[nrow],trainpma[nrow],theInstanceIndex,theInstance, 0, 
						treehmm,nrow,true);

				double dpj =  treehmm.db;

				if (dpj == 0) 
				{
					dlike += Math.log(MINPROB);
				}
				else
				{
					instanceAE(treehmm,0,traindata[nrow],trainpma[nrow], trainSign[nrow],dpj,nrow);
					dlike += Math.log(treehmm.db);
				}

				if (Double.isNaN(dlike))
				{
					traverse(treehmm,0,true);

					System.out.println(nrow+"\t");
					for (int na = 0; na < trainSign[nrow].length; na++)
					{
						System.out.print(trainSign[nrow][na]+"\t");
					}
					System.out.println();
					if(!BREGDREM)
					{
						for (int na = 0; na < treehmm.tranC.dcoeff.length; na++)
						{
							System.out.print(treehmm.tranC.dcoeff[na]+"\t");
						}
					}
					System.out.println();
					System.out.println("treehmm.db = "+treehmm.db);
					throw new Exception();
				}
			} 

			traverse(treehmm,0,false);
			updateParams(treehmm,0,traindata,trainpma,0);
			if (BEQUALSTD)
			{
				averageChildrenSigmas(treehmm);
			}

			traverse(treehmm,0,false); 
			ncount++;
			dpredictitr= (dbesttrainlike-dlike)/(dlike-doldlike);
			if (BDEBUGMODEL)
			{
				System.out.println(ncount +" - Train hmm Likelihood: "+dlike+
						" best "+dbesttrainlike+" "+"oldlike"+doldlike+" "+dpredictitr);
			}
		} 
		while (((dlike-doldlike)/Math.abs(dlike) > BEPSILON)
				&&((bpruneexempt)||(dpredictitr<MAXFUTUREITR)));

		if (BDEBUG)
		{
			System.out.println("Likelihood: "+dlike);
		}

		dtrainlike = dlike;   


		nglobaltime++;
		int nparams =countNodes(treehmm,nglobaltime);
		dgloballike = dlike;

		if (BDEBUG)
		{
			System.out.print(dlike+"\t"+nparams);
		}
		dlike -= nodepenalty*nparams;

		if (BDEBUG)
		{
			System.out.println("\t"+dlike);
		}

		return (dlike);
	}
}