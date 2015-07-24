Version 1.0.6b of the Dynamic Regulatory Events Miner (DREM)
-------------------------------------------------------------------------------
The DREM software was developed by Jason Ernst in collaboration with Ziv Bar-Joseph. 
Any questions about the software or bugs found should be emailed to jernst@cs.cmu.edu.
The latest version of DREM will be available at www.sb.cs.cmu.edu/drem.
More information about the DREM software can be found in the file DREMmanual.pdf.  The
method and application to yeast is presented in:

Ernst J., Vainas O., Harbison C.T., Simon I., and Bar-Joseph, Z. 
Reconstructing Dynamic Regulatory Maps.  Molecular Systems Biology, 3:74, 2007 

-------------------------------------------------------------------------------

This readme file as well as the DREM manual and documentation apply to the original version
of DREM only and not the SDREM extensions.  Updated documentation is forthcoming.
Please contact Anthony Gitter with any questions or comments specific to SDREM at
agitter@cs.cmu.edu

-------------------------------------------------------------------------------
DREM requires Java 1.4 or later to be installed.  Java can be downloaded from
http://www.java.com.

-------------------------------------------------------------------------------
Once Java is installed, to start DREM in Microsoft Windows double click on drem.cmd.  
Otherwise to start DREM, from a command line while in the drem directory enter the 
command

java -mx1024M  -jar drem.jar 

Append "-d defaults.txt" to the above command to have DREM start with its initial
settings be those specified in the defaults.txt file.  

-------------------------------------------------------------------------------
To start DREM with default settings to analyze the sample Heat shock data included in the 
drem directory, if in Microsoft Windows double click on the file 
dremHeatSample.cmd, otherwise from a command line while in the drem
directory type the command:

java -mx1024M -jar drem.jar -d defaultsHeatSample.txt

The sample expression data, expressionHeat.txt, is based on time series expression 
data from Gasch et al, 2000.  The transcription factor-gene association file,
tfinput_heat.txt, is based on data from Harbison et al, 2004 and Hahn et al, 2004.
model_heat.txt contains a saved DREM model for this condition inferred based on
these input data sources.
-------------------------------------------------------------------------------
Source code and Javadoc api for the source code are included in the sourcedoe 
directory.

-------------------------------------------------------------------------------
The TFInput directory contains transcription factor-gene associaton files.  
The ones for yeast that are included with the DREM download are based on Harbison et al, 2004 
and MacIsaac et al., 2006, and are described in the Appendix of the manual.
For E. coli there are curated interactions with direct evidence
based on EcoCyc version 11.5 (Keseler et al, 2005) 
and this input input extended with computational predictions as described
in (Ernst et al, 2008).
If a user adds their own files to this directory they will automatically
appear on the menu of available transcription factor-gene association files.  

-------------------------------------------------------------------------------
DREM is distributed to academic and non-profit users under a non-commercial 
research use license.  

DREM makes use of the Piccolo toolkit which is distributed under a BSD license.  
More information about Piccolo can be found at www.cs.umd.edu/hcil/piccolo

DREM makes use of the Gene Ontology and gene annotations provided by Gene Ontology
Consortium members.  More information about the Gene Ontology can be found at
www.geneontology.org.

-------------------------------------------------------------------------------
References
Ernst J., Beg Q.K., Kay K.A., Balázsi G., Oltvai Z.N., Bar-Joseph Z.
A Semi-Supervised Method for Predicting Transcription Factor-Gene Interactions in Escherichia coli.
PLoS Computational Biology, in press.

Gasch A.P., Spellman P.T., Kao C.M., Carmel-Harel O., Eisen M.B., Storz G., Botstein D., Brown P.O.
Genomic expression programs in the response of yeast cells to environmental changes. 
Mol. Biol. Cell 11, 4241-4257, 2000.

Hahn JS, Hu Z, Thiele DJ, Iyer VR.  Genome-wide analysis of the biology of stress responses through 
heat shock transcription factor.  Mol Cell Biol. 2004 Jun;24(12):5249-56.

Harbison CT, Gordon DB, Lee TI, Rinaldi NJ, Macisaac KD, Danford TW, Hannett NM, Tagne JB,
Reynolds DB, Yoo J, Jennings EG, Zeitlinger J, Pokholok DK, Kellis M, Rolfe PA, Takusagawa KT,
Lander ES, Gifford DK, Fraenkel E, Young RA Transcriptional regulatory code of a eukaryotic genome.
Nature 431: 99104, 2004.

Keseler IM, Collado-Vides J, Gama-Castro S, Ingraham J, Paley S, et al. (2005) 
EcoCyc: A comprehensive database resource for Escherichia coli. Nucleic Acids Res 33: D334-337.

MacIsaac KD, Wang T, Gordon DB, Gifford DK, Stormo GD, Fraenkel E. An improved map of conserved
regulatory sites for Saccharomyces cerevisiae. BMC Bioinformatics 7:113, 2007.