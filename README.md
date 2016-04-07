BaitFisher-package: Design DNA hybrid enrichment baits
======================================================


About the BaitFisher package
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
Christoph Mayer (*), Oliver Niehuis (*), Manuala Sann (*,**)  
(*): Forschungsmuseum Alexander Koenig, 53113 Bonn, Germany  
(**): Museum für Naturkunde, Berlin, Germnay.  


**When using BaitFisher please cite:**   
Mayer, C., Sann, M., Donath, A., Meixner, M., Podsiadlowski, L., Peters, R.S., Petersen, M., Meusemann, K., Liere, K., Wägele, J.-W., Misof, B., Bleidorn, C., Ohl, M., Niehuis, O., 2016. BaitFisher: A Software Package for Multispecies Target DNA Enrichment Probe Design. Mol. Biol. Evol. doi:10.1093/molbev/msw056

Documentation
-------------

A full manual for the BaitFisher package can be found in the "Documentation" folder.
Furthermore, example analyses can be found in the Example folder.

Online help for the BaitFisher and BaitFilter programs can be obtained by called he programs
as follows:

PATH/BaitFisher -h  
PATH/BaitFileter -h



Compiling and installing BaitFisher and BaitFilter
--------------------------------------------------
### System requirements:

BaitFisher can be compiled on all platforms, which provide a C++ compiler that includes the following system header files: unistd.h, sys/stat.h, sys/types.h, dirent.h. These header files are standard header files on all unix systems such as Linux and Mac Os X. I also managed to compile BaitFisher and BaitFilter on Windows using the Mingw compiler (www.mingw.org).

### Compiling BaitFisher and BaitFilter:

Please, unpack the BaitFisher archive you downloaded. 

On the command line enter the BaitFisher-package-master directory. Now type make on the command line. This compiles the BaitFisher as well as the BaitFilter program. The executables of these to programs care called: BaitFisher-v1.2.7 and BaitFilter-v1.0.5. They can be copied on your computer to any place you like. To use these executables they have to be addressed with its full path, or they have to be placed somewhere where in your system path.


Quickstart
----------
More content will me added soon.



>Last updated: April 4th, 2016 by Christoph Mayer

