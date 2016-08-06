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
Christoph Mayer (\*), Oliver Niehuis (\*), Manuela Sann (\*,\*\*)  
(\*): Forschungsmuseum Alexander Koenig, 53113 Bonn, Germany  
(\*\*): Museum für Naturkunde, Berlin, Germnay.  


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

On the command line enter the BaitFisher-package-master directory. Now type make on the command line. This compiles the BaitFisher as well as the BaitFilter program. The executables of these to programs care called: BaitFisher-v1.2.7 and BaitFilter-v1.0.5. They can be copied on your computer to any place you like. To use these executables they have to be addressed with its full path, or they have to be placed somewhere where in your system path.


Quickstart  <a id="quickstart"></a>
----------
More content will me added soon.
For details please see the manual in the "Documentation" folder.

Frequently asked questions  <a id="Frequently-aksed-questions"></a>
--------------------------
1. **BaitFisher does not design baits for all input files or all excised features. What is the reason for this?**

   BaitFisher does not design baits for input files which cannot host a tiling design of the required length.


2. **How can I determine the input files or excised features for which no baits have been designed?**

   There are multple reasons and therefore multiple messages why no baits have been designed for a given input file or excised feature.
   A full list of the messages written to the screen (i.e. to standard output) and the underlying problem are listed here:

>  **Message:**    
   *NOTICE: Skipping alignment since it is too short. Alignment length: the-number*   
   **Reason:**    
   In this case the alignment is shorter than the requested tiling design.   
   **Solution:**   
   If the size of the requested tiling design is decreased, it might be possible to design baits for this alignment.
   Be sure that decreasing the size of the tiling design makes sense.

>  **Message 1:**    
   *ALERT: Skipping <feature name> feature since its overlap with the transcript is too short. Length of overlap: a-number*   
   **Message 2:**    
   *ALERT: Skipping <feature name> feature since the alignment is tool short after removing gap only positions. Length of overlap: a-number*   
   **Reason:**     
   The feature, the transcript or only the overlap of the two are too short to host a complete tiling design.   
   **Solution:**   
   If the size of the requested tiling design is decreased, it might be possible to design baits for this feature.   

>  **Message:**   
   *ALERT: Alignment contains ... features but no single valid bait region. (Details: see above.)*   
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



2. **How can I determine the input files or excised features for which baits have been designed?**

   The bait-loci.txt file list all bait regions for which baits have been designed. The first column in this file is the alignment file name.
   The second column is the gene name if applicalble, or again the alignment file name if no features have been excised.

   On Unix system such as Linux or Mac, the number of unique filenames, gene names, or features can be easily determined on the command line:
   In order to determine the alignment files for which baits have been designed you can use:   
   tail -n +2  loci_baits.txt | cut -f 1 | sort -u

   In order to determine the genes for which baits have been designed you can use:   
   tail -n +2  loci_baits.txt | cut -f 2 | sort -u

   In order to determine the feautes for which baits have been designed you can :   
   tail -n +2  loci_baits.txt | cut -f 2,3 | sort -u

&nbsp;


In case of questions not answered here, please contact the author Christoph Mayer.

>Last updated: August 5th, 2016 by Christoph Mayer

