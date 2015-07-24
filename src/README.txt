Use of this code implies you have accepted the terms of orientation-license.txt
as well as the DREM and Piccolo licenses in the drem subdirectory.
The code requires Java 5.0 or above.


A properties file is used to specify the paramaters and data to use.
See sample.props for an example and sampleEdges.txt and sampleSources.txt
for examples of the formatting required for the network and source files.
See the DREM manual for details regarding the gene expression and
protein-DNA binding data files.  The TF-gene priors have the same
format as the DREM protein-DNA binding data file, but should have the value 0.5
instead of 1 for all non-zero entries (assuming uniform initial priors).

The files SGD_standardToOrf.txt and SGD_orfToStandard.txt are used to
resolve yeast gene synonyms and must be placed in the working directory
in which the code is run.


The usage of the orientation algorithms is:
java SDREM <properties_file>

It is recommended that the maximum heap size is increased when working with large networks:
java -Xmx<memory> SDREM <properties_file> 

where <memory> includes both the heap size and units, e.g. 1500m


Many intermediate output files are generated during each iteration of SDREM.
The *.model file produced during the last iteration can be loaded into DREM to
view the regulatory paths.  A pathEdges_itr*.txt file is produced at each
iteration, which provides the orientation of all network edges that are on
at least one source-target path.