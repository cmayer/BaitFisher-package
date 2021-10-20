BaitFisher-package: Design DNA hybrid enrichment baits
======================================================

Table of contents:
- [About the BaitFisher package](#about-the-baitfisher-package)
- [Documentation](#documentation)
- [Compiling and installing BaitFisher and BaitFilter](#compiling-and-installing-baitfisher-and-baitfilter)
  * [System requirements:](#system-requirements)
- [Quickstart](#quickstart)
- [Frequently aksed questions](#Frequently-aksed-questions)



About the BaitFisher package  <a id="about-the-baitfisher-package"></a>
----------------------------
The BaitFisher package consists of two programs: BaitFisher and BaitFilter. 

BaitFisher was been designed to construct hybrid enrichment baits from multiple sequence alignments (MSAs) or annotated features in MSAs. The main goal of BaitFisher is to avoid redundancy in the construction of baits by designing fewer baits in conserved regions of the MSAs and designing more baits in variable regions. This makes use of the fact that hybrid enrichment baits can differ to some extends from the target region, which they should capture in the enrichment procedure. By specifying the allowed distance between baits and the sequences in the MSAs the user can control the allowed bait-to-target distance and the degree of reduction in the number of baits that are designed. See the BaitFisher paper for details.

BaitFilter was designed (i) to determine whether baits bind unspecifically to a reference genome, (ii) to filter baits that only have partial length matches to a reference genome, (iii) to determine the optimal bait region in a MSA and to convert baits to a format that can be uploaded at a bait constructing company. The optimal bait region can be the most conserved region in the MSA or the region with the highest number of sequences without gaps or ambiguous nucleotides.

Since performance was one of the major design goals, the BaitFisher and BaitFilter programs are both written in C++. 

**Software development:**   
Christoph Mayer, 
Forschungsmuseum Alexander Koenig, 
53113 Bonn, 
Germany

**Main contributors:**
The main contributors to the development of the BaitFisher software are:  
Christoph Mayer (\*), Oliver Niehuis (\*\*,\*), Manuela Sann (\*\*,\*)  
(\*): Forschungsmuseum Alexander Koenig, 53113 Bonn, Germany  
(\*\*): Institut für Biologie I (Zoologie), Albert-Ludwigs-Universität Freiburg, Germany 


**When using BaitFisher please cite:**   
Mayer, C., Sann, M., Donath, A., Meixner, M., Podsiadlowski, L., Peters, R.S., Petersen, M., Meusemann, K., Liere, K., Wägele, J.-W., Misof, B., Bleidorn, C., Ohl, M., Niehuis, O., 2016. BaitFisher: A Software Package for Multispecies Target DNA Enrichment Probe Design. Mol. Biol. Evol. doi:10.1093/molbev/msw056

Documentation  <a id="documentation"></a>
-------------

A full manual for the BaitFisher package as PDF can be found in the "Documentation" folder.
Furthermore, example analyses can be found in the Example folder.

Online help for the BaitFisher and BaitFilter programs can be obtained by called he programs
as follows:

PATH/BaitFisher -h  
PATH/BaitFileter -h



Compiling and installing BaitFisher and BaitFilter <a id="compiling-and-installing-baitfisher-and-baitfilter"></a>
--------------------------------------------------
### System requirements:  <a id="system-requirements"></a>

BaitFisher can be compiled on all platforms, which provide a C++ compiler that includes the following system header files: unistd.h, sys/stat.h, sys/types.h, dirent.h. These header files are standard header files on all unix systems such as Linux and Mac Os X. I also managed to compile BaitFisher and BaitFilter on Windows using the Mingw compiler (www.mingw.org).

### Compiling BaitFisher and BaitFilter:

Please, unpack the BaitFisher archive you downloaded. 

On the command line enter the BaitFisher-package-master directory. Now type "make" on the command line. This compiles the BaitFisher as well as the BaitFilter program. The executables of these to programs care called: BaitFisher-v1.2.7 and BaitFilter-v1.0.5. They can be copied on your computer to any place you like. To use these executables they have to be addressed with its full path, or they have to be placed somewhere where in your system path.


Quickstart  <a id="quickstart"></a>
----------
More content will me added soon.
For details please see the manual in the "Documentation" folder.

Frequently asked questions  <a id="Frequently-aksed-questions"></a>
--------------------------
1. **BaitFisher does not design baits for all input files or all excised features. What is the reason for this?**

   BaitFisher does not design baits for input files or excised features, which cannot host a tiling design of the required length.


2. **How can I determine the input files or excised features for which no baits have been designed?**

   There are different reasons and therefore different messages why no baits have been designed for a given input file or excised feature.
   A full list of the messages written to the screen (i.e. to standard output) and the underlying problem are listed here:

>  **Message:**    
   *NOTICE: Skipping alignment since it is too short. Alignment length: the-number*   
   **Reason:**    
   In this case the alignment is shorter than the requested tiling design.   
   **Solution:**   
   If the size of the requested tiling design is decreased, it might be possible to design baits for this alignment.
   Be sure that decreasing the size of the tiling design makes sense.

>  **Message 1:**    
   *ALERT: Skipping feature-name feature since its overlap with the transcript is too short. Length of overlap: a-number*   
   **Message 2:**    
   *ALERT: Skipping feature-name feature since the alignment is too short after removing gap only positions. Length of overlap: a-number*   
   **Reason:**     
   The feature, the transcript or only the overlap of the two are too short to host a complete tiling design.   
   **Solution:**   
   If the size of the requested tiling design is decreased, it might be possible to design baits for this feature.   

>  **Message:**   
   *ALERT: Alignment contains #number features but no single valid bait region. (Details: see above.)*   
   **Reason:**   
   This message hints at the avbove messages.

>  **Message:**   
   *ALERT: No valid bait region found for alignment:*   
   **Reason:**    
   There are two possibilities: If too many gaps appear in effect in every sequence, BaitFisher might not find a single sequence without gaps
   in length of the requested tiling design.
   If required taxa are specified, no region might exist in which the required taxa are present in the required length without gaps.   
   **Solution:**
   If the size of the requested tiling design is decreased, it might be possible to design baits for this alignment.
   If required taxa are removed,  it might be possible to design baits for this alignment of it.

>  **Message:**   
   *ALERT: No valid bait regions found for feature: name-of-feature*   
   **Reason:**    
   There are two possibilities: If too many gaps appear in effect in every sequence of an excised feature, BaitFisher might not find a single sequence without gaps
   in length of the requested tiling design.
   If required taxa are specified, no region might exist in which the required taxa are present in the required length without gaps.   
   **Solution:**
   If the size of the requested tiling design is decreased, it might be possible to design baits for this feature.
   If required taxa are removed,  it might be possible to design baits for this features.

>  **Message:**   
   *WARNING: For gene " gene-name " no entries where found in the gff file.   
             While this is possible, in most cases the reason is either a wrong attribute name for the   
	     gene-name-field in the gff-file or the gene name specified in the sequence   
	     name of the transcript alignment is incompatible to the gene names in the gff file.   
	     It is recommended to check whether gene name " gene-name " occurs in the   
	     gff file under the attribute name " the-gff-record-attribute-name-for-geneID ".   
	     Note that the gene name as specified in the sequence name and the gene name in   
	     the gff file have to match exactly. The attribute name specified in the parameter    
	     file must be an exact (case insensitive) match to the attribute name in the gff file.*   
    **Reason:**    
    This message can only occur if features are excised from genes/alignments. In the configuration file the
    user specifies which attribute contains the gene name. A gff file contains lines of features. Each line has
    nine tab delimited columns. The last column contains the so called attributes. E.g. "Parent=FBtr0081133".
    If BaitFisher found a gene name in the sequence name in the alignment file, it will search for all features
    that belong to this gene. The features are identified by the gene name saved in an attribute.
    If for a given gene name no entry is found in the gff file, this is in most cases due to a nameing problem,
    or incompatibility of the gene names and the gff file. If no features with the given gene name are found,
    no baits will be designed for this gene.
   **Solution:**   
   Search for the gene-name given in the error message und verify that it can be found in the gff file.
   Make sure it can be found under the  the-gff-record-attribute-name-for-geneID specified in the configuration file.
   In most cases this name is "parent".

>  **Message:**   
   *WARNING: For gene " gene-name " no entries where found in the gff file that matched the specified feature name. 
	     Since entries for this gene name have been found, but with different feature names, the problem could be the feature name 
	     or the fact that features exist for this gene.*
    **Reason:**    
    Here BaitFisher found entries with gene-name for the given attribute name, but they could not be found for the feaute
    that was specified.   
   **Solution:**   
   Have a closer look at the lines containing the features you want to excise. If you want to excise CDS features
   look at lines in the gff file containing this feature and check that the gene name is found under the correct attribute name.
   Note that different features often have different attribute names for the name of the gene.
   E.g. the gene or mRNA feature usually does not mention the gene name under the attribute name parent, but under the attribute name ID
   since these features have no parent features. In contrast, the exon, CDS or intron features are usually
   assigned to a parent gene and the attribute containing the gene name is often the parent attribute.
   Check that the given gene-name is found for at least one feature of the given type (e.g. CDS) under
   the specified attribute. If the gene name does not occur at all double check that the gff file and
   the genome file are compatible and that the gene name has been determined correctly.


2. **How can I determine the input files or excised features for which baits have been designed?**

   The bait-loci.txt file list all bait regions for which baits have been designed. The first column in this file is the alignment file name.
   The second column is the gene name if applicalble, or again the alignment file name if no features have been excised.

   On Unix system such as Linux or Mac, the unique filenames, gene names, or features can be easily determined on the command line:
   In order to determine the alignment files for which baits have been designed you can use:   
   tail -n +2  loci_baits.txt | cut -f 1 | sort -u

   In order to determine the genes for which baits have been designed you can use:   
   tail -n +2  loci_baits.txt | cut -f 2 | sort -u

   In order to determine the feautes for which baits have been designed you can :   
   tail -n +2  loci_baits.txt | cut -f 2,3 | sort -u

&nbsp;


In case of questions not answered here, please contact the author Christoph Mayer.

>Last updated: September, 2017 by Christoph Mayer

