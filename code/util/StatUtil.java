package util;
// package edu.cmu.cs.stem;
import java.util.*;


/**
 *The class implements static methods to compute values related to
 *the binomial, hypergeometric, and normal distribution
 * @author Jason Ernst
 */
public class StatUtil
{
    /**
     *for low accuracy calculations any p-value below this value will be set to 0.  
     *any p-value within 1 of this value will be set to 1
     */
    static double ERROR = Math.pow(10, -12); 

    /**
     *this hashtable caches previously computed binomial coefficient values so they do not need to be recomputed
     */
    static Hashtable ht =new Hashtable();
    
    /** Cache the value of log_10(2) */
    static final double log2base = Math.log(2);

    /**
     *Returns the log of the binomial coefficient N choose ni
     */
    static double logbinomcoeff(int ni, int N)
    {
	String sz = ni+";"+N;
	Double dobj = (Double) ht.get(sz);
        double dsum;

	if (dobj != null)
       	{
	   dsum = ((Double) dobj).doubleValue();
       	}
        else
	{
           dsum = 0;
           int dmax = Math.max(ni,N-ni);
           int dmin = Math.min(ni,N-ni);

	   //the log of the part of the numerator not cancelled by the larger factorial in the denominator
           for (int nj =dmax+1; nj <=N; nj++)
           {
              dsum += Math.log(nj);
           }

	   //subtract off the log of the denominator of the smaller term
           for (int nj =2 ; nj <=dmin; nj++)
           {
              dsum -= Math.log(nj);
           }
	   //store it
           ht.put(sz,new Double(dsum));
       }   

       return dsum;
    }


    public static void main(String[] args)
    {
	//System.out.println(StatUtil.hypergeometrictail(27, 451, 13198-451, 118));
	System.out.println(StatUtil.hypergeometrictail(99, 451, 13198-451, 362));
    }

    /**
     * Returns the probability of seeing x of type A, when there nA objects of type A
     * nB objects of type B, and nm objects total drawn
     */
    static double hypergeometric(int nx, int nA, int nB, int nm)
    {
       if (nx < 0)
       {
          return 0;
       }
        
       //System.out.println(nx+" "+nA+" "+nB+" "+nm);
       double dprob = 0;
      
       int nminx = Math.max(nm-nB, 0); //how many must be of type A
       int nmaxx = Math.min(nx, nA); //the most that can be of type A
    
       if (nminx > nmaxx)
       {
          return 0;
       }

       double dsum = -logbinomcoeff(nm, nA + nB)+ logbinomcoeff(nx, nA)
	             + logbinomcoeff(nm-nx, nB);
       //old way if nm is greater than nB, then log 1
       //correct way nB, then still log 1

       dprob = Math.pow(Math.E,dsum);

       if (dprob <=0)// ERROR)
       {
          dprob = 0;
       }
       else if (dprob >= 1)//-ERROR)
       {
          dprob = 1;
       }
	
       return dprob;
    }

    // Call with nx - 1 if you want the probability of x or more objects (which
    // is usually the case)
    /**
     * Returns the probability of seeing more than x objects of type A, when there nA objects of type A
     * nB objects of type B, and nm objects total drawn
     */
    //static double hypergeometriccumulative(int nx, int nA, int nB, int nm)
    public static double hypergeometrictail(int nx, int nA, int nB, int nm)
    {
       if (nx < 0)
       {
	   //return 0;
	  return 1;
       }
        
       //System.out.println(nx+" "+nA+" "+nB+" "+nm);
       double dprob = 0;
      
       //int nminx = Math.max(nm-nB, 0); //how many must be of type A
       //int nmaxx = Math.min(nx, nA); //the most that can be of type A
   
       int nminx = Math.max(nx+1,0); //the first element to the right of x or 0
       int nmaxx = Math.min(nA,nm);  //the max of type A there can be 

       if (nminx > nmaxx)
       {
	   //if nx approaches infinity tail should be 0
	   return 0;
	   //had to have seen more of type A
	   //return 1;
       }

       double dsum = -logbinomcoeff(nm, nA + nB)+ logbinomcoeff(nminx, nA)
	             + logbinomcoeff(nm-nminx, nB);
       //old way if nm is greater than nB, then log 1
       //correct way nB, then still log 1

       //dprob = Math.pow(Math.E,dsum);
       double dlogprob=dsum;
       //System.out.println("!!!"+dsum);
       for (int ni = nminx+1; ni <= nmaxx; ni++)
       {
	   //computing the increase in probability mass
	   //numerator has nA!/(ni!(nA-ni)!) * nB!/((nm-ni)!(nB-nm+ni)!)
	   //denominator has (nA+nB)!/(nm!(nA+nB-nm)!)
	   //numerator has nA!/((ni-1)!(nA-ni+1)!) * nB!/((nm-ni+1)!(nB-nm+ni-1)!)
	   //denominator has (nA+nB)!/(nm!(nA+nB-nm)!)
	   //cancelling gives
	   //numerator has 1/(ni!(nA-ni)!) * 1/((nm-ni)!(nB-nm+ni)!)
	   //numerator has 1/((ni-1)!(nA-ni+1)!) * 1/((nm-ni+1)!(nB-nm+ni-1)!)
	   //dval += Math.log(nA-ni+1)-Math.log(nB-nm+ni)
	   //      + Math.log(nm-ni+1)-Math.log(ni);
          dsum  += Math.log(nA-ni+1)-Math.log(nB-nm+ni)
	         + Math.log(nm-ni+1)-Math.log(ni);

	  //System.out.println((nA-ni+1)+" "+(nB-nm+ni)+" "+(nm-ni+1)+" "+ni+" "+dsum+" "+dlogprob);
	  //a+b+c+d+e
	  //log(a+b+c+d+e)
	  //log(a+b+c+d)+log(e)

	  //log(e) + log(a+b+c+d+e) - log(e)
	  //log(e) + log((a+b+c+d+e)/e)
	  //log(e) + log(1+(a+b+c+d)/e)
	  //log(e) + log(1+Math.exp(log(a+b+c+d)-log(e)))
	  //log(e) + log(1+exp(log(a+b+c+d)/log(e)))
	  //log(a+b+c+d+e)
	  //log(a)  + log(b) --> log(a+b)
	  //log(a*b    
	  if (dsum >= dlogprob)
	  {   
	     dlogprob = dsum + Math.log(1+Math.pow(Math.E,dlogprob-dsum));     
	  }
	  else
	  {
	     dlogprob = dlogprob + Math.log(1+Math.pow(Math.E,dsum-dlogprob)); 
	  }
       }

       dprob = Math.pow(Math.E,dlogprob);

       if (dprob <= 0)
       {
	   return 0;
       }
       else if (dprob >= 1)
       {
	   return 1;
       }
       else
       {
	   //System.out.println(dprob);
	   return dprob;
       }
    }


    /**
     * Returns the probability of seeing x or few objects of type A, when there nA objects of type A
     * nB objects of type B, and nm objects total drawn
     */
    static double hypergeometriccumulative(int nx, int nA, int nB, int nm)
    {
       if (nx < 0)
       {
          return 0;
       }
        
       //System.out.println(nx+" "+nA+" "+nB+" "+nm);
       double dprob = 0;
      
       int nminx = Math.max(nm-nB, 0); //how many must be of type A
       int nmaxx = Math.min(nx, nA); //the most that can be of type A
   
       if (nminx > nmaxx)
       {	   
	   //if nx less < 0 cumulative should be 0
          return 0;
       }

       double dsum = -logbinomcoeff(nm, nA + nB)+ logbinomcoeff(nminx, nA)
	             + logbinomcoeff(nm-nminx, nB);
       //old way if nm is greater than nB, then log 1
       //correct way nB, then still log 1

       dprob = Math.pow(Math.E,dsum);
       for (int ni = nminx+1; ni <= nmaxx; ni++)
       {
	   //computing the increase in probability mass
	   //numerator has nA!/(ni!(nA-ni)!) * nB!/((nm-ni)!(nB-nm+ni)!)
	   //denominator has (nA+nB)!/(nm!(nA+nB-nm)!)
	   //numerator has nA!/((ni-1)!(nA-ni+1)!) * nB!/((nm-ni+1)!(nB-nm+ni-1)!)
	   //denominator has (nA+nB)!/(nm!(nA+nB-nm)!)
	   //cancelling gives
	   //numerator has 1/(ni!(nA-ni)!) * 1/((nm-ni)!(nB-nm+ni)!)
	   //numerator has 1/((ni-1)!(nA-ni+1)!) * 1/((nm-ni+1)!(nB-nm+ni-1)!)
          dsum  += Math.log(nA-ni+1)-Math.log(nB-nm+ni)
	         + Math.log(nm-ni+1)-Math.log(ni);
          
          dprob += Math.pow(Math.E,dsum);
       }

       if (dprob <= ERROR)
       {
          dprob = 0;
       }
       else if (dprob >= 1-ERROR)
       {
          dprob = 1;
       }
	
       return dprob;
    }

    /**
     *Computes the probability of seeing x or smaller successes in dN trials 
     *where the probability of a success is dp.
     */
    static double binomialcumulative(double x, double dN, double dp)
    {
	//System.out.println(x+" "+N+" "+dp);
	//returns the amount of probability less than or equal to x
            
	if (x > dN)
	{
	    return 1;
	}
	else if ((x < 0)||(dp<=0)||(dp>=1))
	{
	    return 0;
	}

        int N = (int) Math.ceil(dN);
        double dterm = logbinomcoeff(0, N);
        double dpv1 = Math.log(dp);
        double dpv2 = Math.log(1-dp);

        dterm += N*dpv2;
        double dprob = Math.pow(Math.E,dterm);
	double dpdiff =  dpv1-dpv2;
        for (int ni = 1; ni <= x; ni++)
	{
	    //N!/(ni!(N-ni)!)
	    //N!/((ni-1)!(N-ni+1)!)
	   dterm += Math.log(N-ni+1) - Math.log(ni);
	   dterm += dpdiff;
           dprob += Math.pow(Math.E,dterm);
	}
        if (dprob <= ERROR)
	{
           dprob =0;
	}
        else if (dprob >= 1-ERROR)
	{
           dprob =1;
	}
	
        //if (!((dprob<=1) &&(dprob >=0)))
	//   System.out.println(dp+" "+dterm);

        return dprob;
    }



    static double binomialtail(int x,int N, double dp)
    {
	//System.out.println(x+" "+N+" "+dp);
	//returns the amount of probability less than or equal to x
         
	if (x > N)
	{
	    return 0;
	}  
	else if ((x < 0)||(dp<=0)||(dp>=1))
	{
	    return 1;
	}
	//remove dN
        //int N = (int) Math.ceil(dN);
	x++;
        double dterm = logbinomcoeff(x, N);
        double dpv1 = Math.log(dp);
        double dpv2 = Math.log(1-dp);

        dterm += x*dpv1+(N-x)*dpv2;
        //double dprob = Math.pow(Math.E,dterm);
	double dlogprob = dterm;
	double dpdiff =  dpv1-dpv2;
	double dprob;
        for (int ni = x+1; ni <= N; ni++)
	{
	    //N!/(ni!(N-ni)!)
	    //N!/((ni-1)!(N-ni+1)!)
	   dterm += Math.log(N-ni+1) - Math.log(ni) + dpdiff;
	   //dterm += dpdiff;

	  if (dterm >= dlogprob)
	  {   
	     dlogprob = dterm + Math.log(1+Math.pow(Math.E,dlogprob-dterm));     
	  }
	  else
	  {
	     dlogprob = dlogprob + Math.log(1+Math.pow(Math.E,dterm-dlogprob)); 
	  }
	}
        dprob = Math.pow(Math.E,dlogprob);

       if (dprob <= 0)
       {
	   return 0;
       }
       else if (dprob >= 1)
       {
	   return 1;
       }
       else
       {
	   return dprob;
       }
        //if (!((dprob<=1) &&(dprob >=0)))
	//   System.out.println(dp+" "+dterm);

        //return dprob;
    }

    /**
     *computes the value of f(x) where f is a density for a normal distribution with
     *mean dmu and standard deviation dsigma.  Not used by STEM.
     */
    static double TWOPISQRT = Math.sqrt(2*Math.PI);
    //static double LOGTWOPISQRT = Math.log(Math.sqrt(2*Math.PI));
    //static HashSet calls = new HashSet();
    //static int ncalls = 0;
    static double normaldensity(double x, double dmu, double dsigma)
    {
	//calls.add(x+";"+dmu+";"+dsigma);
	//ncalls++;
	//if (ncalls % 1000 == 0)
	//    System.out.println(ncalls);//+"\t"+calls.size());

	double dens;
	double dxmudiff = (x-dmu);
	dens = Math.exp(-dxmudiff*dxmudiff/(2*dsigma*dsigma))/(dsigma*TWOPISQRT);
	//System.out.println("@@@\t"+x+"\t"+dmu+"\t"+dsigma+"\t"+dens);
	return dens;
    }


	/**
	 * Use Fisher's Exact test to obtain a p-value for the overlap
	 * of the bound and differentially expressed genes.  Case-sensitve.
	 * @param overlap - genes bound and expressed
	 * @param boundGenes - genes bound by a target
	 * @param expressedGenes - genes differentially expressed
	 * @param possibleBoundGenes - all genes considered in the binding experiment
	 * @param possibleExpressedGenes - all genes considered in the expression experiment
	 * @return
	 */
	public static double overlapSignificance(Set<String> overlap,
			Set<String> boundGenes,
			Set<String> expressedGenes, Set<String> possibleBoundGenes,
			Set<String> possibleExpressedGenes)
	{
		return StatUtil.overlapSignificance(overlap, boundGenes, expressedGenes,
				MapUtil.union(possibleBoundGenes, possibleExpressedGenes));
	}


	/**
	 * Use Fisher's Exact test to obtain a p-value for the overlap
	 * between two sets of genes.  Case-sensitive.
	 * @param overlap - genes in set1 and set2
	 * @param geneSet1
	 * @param geneSet2
	 * @param allGenes - all genes that could have possibly been in either
	 * set 1 or set 2
	 * @return
	 */
	public static double overlapSignificance(Set<String> overlap,
			Set<String> geneSet1,
			Set<String> geneSet2,
			Set<String> allGenes)
	{
		return hypergeometrictail(overlap.size() - 1,
				geneSet1.size(),
				allGenes.size() - geneSet1.size(),
				geneSet2.size());
	}
	
	/**
	 * Use Fisher's Exact test to obtain a p-value for the overlap
	 * between two sets of objects
	 * @param overlap number of objects in both sets (first row, first col
	 * in 2x2 contingency table)
	 * @param set1 number of objects in set 1 (sum of first row in contingency table)
	 * @param set2 number of objects in set 2 (sum of first col in contingency table)
	 * @param allPossible total number of possible objects (sum of all entries in
	 * contingency table)
	 * @return
	 */
	public static double overlapSignificance(int overlap,
			int set1, int set2, int allPossible)
	{
		// Subtract 1 because we want the probability of seeing the overlap
		// or more than the overlap at random (tail only considers more than
		// the overlap)
		return hypergeometrictail(overlap - 1,
				set1,
				allPossible - set1, // sum of second col in contingency table
				set2);
	}
	
	/**
	 * Calculate log_2(x)
	 * @param x
	 * @return
	 */
	public static double log2(double x)
	{
		return Math.log(x)/log2base;
	}

}