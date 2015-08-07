***********************************************
Version 1.2.1
***********************************************

The SDREM software is described in:

Linking the signaling cascades and dynamic regulatory networks controlling stress responses. 
Anthony Gitter, Miri Carmi, Naama Barkai, Ziv Bar-Joseph. 
Genome Research. 23:2, 2013.
http://www.genome.org/cgi/doi/10.1101/gr.138628.112

Identifying proteins controlling key disease signaling pathways. 
Anthony Gitter, Ziv Bar-Joseph. 
Bioinformatics. 29:13, 2013.
http://bioinformatics.oxfordjournals.org/content/29/13/i227.full

The SDREM method for reconstructing signaling and regulatory response networks: Applications for studying disease progression.
Anthony Gitter, Ziv Bar-Joseph.
Systems Biology of Alzheimer's Disease.  Methods in Molecular Biology. Volume 1303, 2015.
http://link.springer.com/protocol/10.1007/978-1-4939-2627-5_30

Please cite one or more of those manuscripts if using the software in published work.  The 2015 book chapter provides a step-by-step protocol for using the software.  Contact gitter@biostat.wisc.edu with any questions or for a copy of the book chapter PDF if you cannot access it.

Use of this code implies the user has accepted the terms of license.txt.  The code requires Java 5.0 or above.

Properties files are used to specify the parameters and input data.  See the DREM manual (http://sb.cs.cmu.edu/drem/DREMmanual.pdf) for details regarding the gene expression and protein-DNA binding file formats.  The TF-gene priors have the same format as the DREM protein-DNA binding data file, but have the value 0.5 instead of 1 for all non-zero entries (assuming uniform initial priors).

Examples of how to run SDREM on a PBS cluster are included.  The maximum Java heap size must be increased when working with large networks.

Many intermediate output files are generated during each iteration of SDREM.  The SDREM executables, code, input files, output produced, and example data are described below:


***********************************************
Precomputing and storing paths (optional)
***********************************************
This is a recommended  optional step that can be performed before running SDREM.

StorePaths.jar - The executable for preprocessing the network data.  It searches for paths from the sources to each TF and writes them all to disk, optionally filtering them to keep only the highest confidence paths.  The stored files are large and can require many GBs of disk space in total.  Without this preprocessing, SDREM has to search for the paths many times as it iterates, which is wasteful and makes it very, very, very slow for large networks.

allPathsEgfPriors.props - An example properties file for StorePaths.jar.  store/filter means the user wants to enumerate paths and remove low-confidence paths (recommended for large PPI networks to speed network orientation).  The next lines define the sources, targets, and networks.  Node priors are used to give weights to vertices in the network.  If the user stores and filters paths, then he/she can delete the intermediate output in stored.paths.dir once the filtering step is done (give the final directory with the filtered paths to SDREM as input).

allPathsEgfPriors.qub - An example showing how to call StorePaths.jar.  This is a submission script for a PBS cluster and the last line shows how to call the jar.


***********************************************
SDREM
***********************************************
The SDREM algorithm executable.

sdrem.jar - The SDREM executable.

sdremEgfPriors.props - A sample SDREM properties file.  model.dir is where the output will be written and where some of the input files must be located.  stored.paths.dir is the location of the filtered paths from StorePaths.jar if it was used.  It's important to use the same settings (sources, targets, node priors, path length, etc.) for the StorePaths and SDREM properties.  predict.knockdown is used to enable the knockdown effect prediction described in the SDREM Bioinformatics paper and can be commented out if this feature is not desired.  The rest of the parameters can usually be left at the default values and are described in the SDREM Genome Research paper supplement.

DREM_defaults.txt - A separate DREM properties file has to be configured.  It's the same format that the original DREM software uses (as described in http://sb.cs.cmu.edu/drem/DREMmanual.pdf) except the TF-gene_Interactions_File will be generated dynamically and Active_TF_influence is a new SDREM parameter.  This file must be located in the model.dir.

sdremEgfPriors.qub - An example showing how to call sdrem.jar on a PBS cluster.


***********************************************
Modified DREM
***********************************************
The DREM software (http://sb.cs.cmu.edu/drem/) that SDREM was built upon allows visualization of the
active TFs and gene expression profiles after SDREM is run.  There are several differences between the version distributed here and DREM 2.0:
None of the new features described in the DREM 2.0 manuscript are present yet (e.g. support for motif finding).
To view the final output (assuming 10 iterations), load 10.model as the saved model and tfActivityPriors_round9.txt as the TF-gene interactions.  Use these priors instead of tfActivityPriors_round10.txt because tfActivityPriors_round10.txt are the updated priors after running the 10th round of network orientation as opposed to the file that was given to DREM as input at the start of the 10th SDREM iteration.
The split table will show the activity score for each TF at that split and the max activity score across all splits.
Key TF Labels includes options to display TFs at each split based on activity score.  If the user uses activity scores to choose which TFs to show, the slider will be used to calculate 10^X (instead of 10^-X) and all TFs with activity scores greater than 10^X at a split are shown.  
The output file 10.targetsStd can be used to help choose the activity score threshold.  The last column in this file gives the max activity score for each active TF.  Therefore, the user can find the minimum of these values and use it (or a value slightly less than it) as the threshold for display purposes.
There is also an option to show the TF only at the earliest split at which it exceeds this threshold on each path.
When using activity score to display TFs, there will be [1] or [2] after the TF id if binary splits were used.  [1] means the TF primarily controls the lower path out of the split and [2] is for the higher path.  For higher order splits (3+), [1] is the lowest path out of the split and [3] is the highest.

drem.jar - The DREM executable, which can be run without command line parameters.  Use the GUI to load the model file from SDREM and the data.  Optionally provide the DREM_defaults.txt file at the command line to load the same settings that SDREM used.  See http://sb.cs.cmu.edu/drem/DREMmanual.pdf for details.


***********************************************
Source code
***********************************************
The executable jar files correspond to following classes:
scripts.StorePaths
alg.SDREM
edu.cmu.cs.sb.drem.DREM_IO


***********************************************
Output files
***********************************************
Intermediate output files are generated at every iteration of SDREM, which can be used for debugging and analysis.  At each iteration 'N', SDREM writes the following files (some filenames may vary slightly based on the parameters used):

/N - A subdirectory that contains output from the DREM runs that use randomized TF binding data.  These can be ignored or deleted.

N.model - The DREM model file for iteration N, which can be loaded by drem.jar.
N.model.activities - Activity information for the TFs used to determine which TFs should be targets.
N.model.activitiesDynamic - Activity information for the TFs used to determine which TFs should be targets.
N.model.activitiesStd - Activity information for the TFs used to determine which TFs should be targets.

N.targets -The TFs that were selected as targets at iteration N and their target weights.  These are the proteins that were used for the network orientation.
N.targetsStd - Like N.targets but contains information about the distribution of random TF activity scores.

itrN.out - A log file.  It also prints the version of SDREM that was run.

conflictOrientations_itrN.txt - A code that represents how each PPI was oriented.  See pathEdges_iterN.txt for a human-readable version.

nodeScores_iterN_Pathweight_10_1000.txt - A summary of the oriented network at iteration N.  It gives the sources, target TFs (same ones as N.targets) and various measures of how many oriented paths use a particular protein.  The '% top 1000 paths through node' column was used to select nodes in the SDREM Genome Research paper.  From this file, the user can extract the sources, targets, and all other proteins that have a value >= 0.01 in this column (called the 'internal' or 'signaling' proteins), and they are also written as a Cytoscape-formatted file at the final iteration.  That score means that at least 1% of the highest confidence oriented paths go through that particular protein.

pathEdges_iterN.txt - The predicted PPI orientation for all edges that were used on at least one source-target path (i.e. some edges are excluded because they are not "between" the sources and targets).

scores_itrN_r1.0_PathWeight_10_1000.txt - Connectivity of the true target TFs in the network orientation versus randomly selected targets, which is used to determine which TFs' priors should be increased or decreased at the next iteration of SDREM.

tfActivityPriors_roundN.txt - The new TF-gene binding file that will be used as input at iteration N+1.

statisfiedPaths_itrN.txt.gz - Lists every oriented, satisfied path that connects a source and target and has less than 5 edges.  It doesn't contain any information that couldn't be reconstructed from pathEdges_itrN.txt but is more convenient for downstream analyses that examine individual paths.


After the final iteration, iteration M, SDREM writes:
postProcessing.out - A log file.

topPathEdges_itrM.sif - A file that can be loaded into Cytoscape v2.8 to visualize the high-confidence paths SDREM inferred.  Only edges between sources, targets, or internal nodes are shown.
topPathNodes_itrM.noa - A file that can be loaded into Cytoscape v2.8 to annotate the nodes on the high-confidence path with their role in the network (Source, Target, or Internal).

droppedTargets.txt - Targets that were present at iteration N-1 but not iteration N.
targetsByIteration.txt - All targets at each iteration.
newTargets.txt - Targets that were not present at iteration N-1 but are at iteration N.

singleKnockdown_itrM.txt - File that is generated only if single or double knockdown effects were requested.  The SDREM Bioinformatics paper defines the scoring metrics.
doubleKnockdown_itrM.txt - File that is generated only if double knockdown effects were requested.  The SDREM Bioinformatics paper defines the scoring metrics.


***********************************************
Sample data
***********************************************
The example property files refer to these data files, which demonstrate the expected file formats.  These are human datasets, and all ids are NCBI Gene ids (http://www.ncbi.nlm.nih.gov/gene).  DREM always requires that all proteins and genes use the same type of identifier, and it is especially important that the same types of ids are used in the PPI network and the TF columns of the TF-gene interaction files.  For example, do not use gene symbols in one and ids in the other.

Yarden_MCF10A_expr.txt - Sample expression data in the DREM format (see http://sb.cs.cmu.edu/drem/DREMmanual.pdf).  The data are from MCF10A cells stimulated with EGF.  Please cite PMID 17322878 if using this data.

ppi_ptm_pd_edges.txt - The human interaction network, which uses PPI from BioGRID and HPRD (pp), predicted protein-DNA binding edges (pd), and post-translational modifications from HPRD (ptm).  Please cite the SDREM Bioinformatics paper and the original data sources if using this data.

tfList.txt - Gene ids for the TFs that have protein-DNA interactions.  This is used for precomputing paths with StorePaths.jar, and these proteins are also the set of possible random targets when SDREM builds a background distribution of TF network connectivity scores (if using precomputed paths).

tfActivityPriors_round0.txt  - A DREM-style protein-DNA binding grid (see http://sb.cs.cmu.edu/drem/DREMmanual.pdf) but with 0.5 as a prior for all interactions.  This will be updated at each iteration as the network is used to refine the priors.  Please cite PMIDs 22897824 and 20219943 if using this data.

sources.txt - The ids of the proteins that are used as the sources for the network orientation in all of the iterations.  These two proteins were selected from the EGFR pathway.

egfPriors.txt - A tab-separated file that lists proteins and the prior probability of their involvement in the signaling pathway.