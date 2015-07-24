package edu.cmu.cs.sb.drem;
import java.util.HashMap;

import edu.cmu.cs.sb.core.*;

/**
 * Implements a naive bayes classifier
 * This class is used for predicting whether a gene with a specified
 * set of regulators will be filtered or not. 
 */
public class DREM_NaiveBayes
{
	/** The number of TFs ??? */
	int numcols;
	int numclasses;
	/** Stores the probability that P(C=c) for all c ??? */
	double[] pClass;
	/** Used for mapping training data values to feature indices */
	HashMap featureMap;

	//p(C|F_1,...,F_n)
	//P(C,F_1,...,F_n)/P(F_1,...,F_n)
	//P(C)P(F_1|C)...P(F_n|C)/(P(C=1)P(F_1,...,F_n|C=1)+P(C=0)P(F_1,...,F_n|C=0))
	/**first class, second feature, third feature val.
	 * Stores P(F_i=f|C=c) for all values of c, i, and f ??? */
	double[][][] pNB;
	// nBaseCount is used in many classes outside of DREM_NaiveBayes
	/**first feaure, second feature val*/
	int[][] nBaseCount;

	/**
	 * Class constructor - builds the classifier
	 */
	public DREM_NaiveBayes(int[][] traindata, int[][] traindataIndex,
			int numcols, int[] y, int numclasses)
	{
		this.numclasses = numclasses;
		this.numcols = numcols;
		int numrows = traindata.length;
		
		// Create a new HashMap and populate it with a mapping from
		// each possible feature value to an index
		featureMap = new HashMap();
		int numfeaturevals = 0;
		
		// By default, 0 will always be a possible feature value
		// It must be added separately because traindata only contains
		// nonzero values, the zeros are implicit
		featureMap.put(Integer.valueOf(0),Integer.valueOf(numfeaturevals));
		numfeaturevals++;
		
		for(int nrow = 0; nrow < numrows; nrow++)
		{
			for(int nindex = 0; nindex < traindata[nrow].length; nindex++)
			{
				// If this feature value has not been observed before, add
				// it to the map
				if(!featureMap.containsKey(Integer.valueOf(traindata[nrow][nindex])))
				{
					featureMap.put(Integer.valueOf(traindata[nrow][nindex]),
							Integer.valueOf(numfeaturevals));
					numfeaturevals++;
				}
			}
		}		

		pNB = new double[numclasses][numcols][numfeaturevals];
		nBaseCount = new int[numcols][numfeaturevals];
		pClass = new double[numclasses];

		for (int nrow = 0; nrow < numrows; nrow++)
		{
			int ntfindex = 0;
			// Iterate through all TFs ???
			for (int ncol = 0; ncol < numcols; ncol++)
			{
				// Iterate through all TFs that regulate the gene at nrow
				// until finding one with index >= ncol ???
				while ((ntfindex < traindataIndex[nrow].length)&&
						(ncol > traindataIndex[nrow][ntfindex]))
				{
					ntfindex++;
				}

				int nfeatureindex;
				if ((ntfindex<traindataIndex[nrow].length)&&
						(ncol == traindataIndex[nrow][ntfindex]))
				{
					// The TF indexed by ncol regulates the gene at nrow ???
					nfeatureindex = getFeatureIndex(traindata[nrow][ntfindex]);
				}
				else
				{
					// The TF does not regulate the gene so its feature
					// value is 0 (only nonzero feature values are stored
					// in traindata) ???
					nfeatureindex = getFeatureIndex(0);
				}
				
				// Increment a counter because a gene with this class assignment
				// and feature value is bound by the TF with index ncol ???
				pNB[y[nrow]][ncol][nfeatureindex]++;
				nBaseCount[ncol][nfeatureindex]++;
			}
			// Increment the class counter for the class assigned to this gene
			pClass[y[nrow]]++;
		}

		// At this point, pNB = count(F_i=f & C=c)
		// and pClass = count(C=c)
		
		for (int nclass = 0; nclass < numclasses; nclass++)
		{
			for (int ncol =0; ncol < numcols; ncol++)
			{
				for (int nfeatureval = 0; nfeatureval < numfeaturevals; nfeatureval++)
				{
					// No need to involve the count(genes) terms because they cancel ???
					// P(F_i=f & C=c) = count(F_i=f & C=c)/count(genes)
					// P(C=c) = count(C=c)/count(genes)
					// P(F_i=f|C=c) = P(F_i=f & C=c)/P(C=c)
					//  = (count(F_i=f & C=c)/count(genes))/(count(C=c)/count(genes))
					//  = count(F_i=f & C=c)/count(C=c)
					pNB[nclass][ncol][nfeatureval] /= pClass[nclass];
				}
			}
			// P(C=c) = count(C=c)/count(genes)
			pClass[nclass] /= numrows;
		}
	}

	/**
	 * Returns the probability of each class for the feature values in theInstance
	 * under the naive bayes models.
	 */
	public double[] distributionForInstance(int[] theInstance)
	{
		double[] dist = new double[numclasses];
		double dsum = 0;
		for (int nclass = 0; nclass < numclasses; nclass++)
		{
			dist[nclass] = pClass[nclass];
			for (int nfeature =0; nfeature < numcols; nfeature++)
			{
				// theInstance[nfeature] is the binding sign of the TF
				// at index nfeature ???
				// Getting the index for the feature value assumes
				// that theInstance does not have feature values that
				// were not in the training data.				
				dist[nclass] *= pNB[nclass][nfeature][getFeatureIndex(theInstance[nfeature])];
			}
			dsum += dist[nclass];
		}

		// No need to calculate P(F_1,...,F_n|C=c) because
		// the probabilities can be normalized instead ???
		for (int nclass = 0; nclass < numclasses; nclass++)
		{
			if (dist[nclass] != 0)
			{
				dist[nclass] /= dsum;
			}
		}

		return dist;
	}
	
	/**
	 * Takes a feature value and returns the appropriate index
	 * @param nkey the feature value
	 * @return the index for that feature value or -1 if the feature value
	 * was not in the training data
	 */
	public int getFeatureIndex(int nkey)
	{
		if(featureMap.containsKey(Integer.valueOf(nkey)))
		{
			Integer tmpIndexInt = (Integer) featureMap.get(Integer.valueOf(nkey));
			return tmpIndexInt.intValue();
		}
		else
		{
			System.err.println("Naive Bayes feature map couldn't find: " + nkey);
			return -1;
		}
	}
}

