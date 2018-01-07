Version history of the BaitFisher and BaitFilter programs
=========================================================

BaitFisher:
-----------

BaitFisher version 1.2.7:
*   First public release

BaitFisher version 1.2.8:	
*   New g++ (GNU c++) compilers are more restrictive than older versions. The code has been changed to work with newer compilers.
*   Typos have been removed in the screen output as well is in the manual.


BaitFilter:
-----------

BaitFilter version 1.0.5:
*   First public release

BaitFilter version 1.0.6:
* Modified comparison functions when sorting baits. The problem was that different implementations of the sort function lead to different orders of sorted bait regions if the characteristics by which bait regions are sorted are different. Therefore bait regions selected by BaitFilter in the modes -m as, -m ab, etc. could differ for different implementations of the library functions, in particular on different platforms. The reason is that often loci are equivalent according to the optimality criterion. Selecting one or the other loci is possible, but it would be better if this would happen in a consistent and predictable way. This has been changed so that loci at the beginning of the sequences are preferred. -> Results can differ in comparison to version 1.0.5. The selected bait regions  can differ, but will be equally good.
* Corrected typos in output. Improved output.
* Log output and verbosity flag threshold values have changed.
* The B flag (and the --blast-executable) did not work. This has been fixed.
* The default value for the evalue passed to blastn as been changed. The default now is twice the second hit value. (Was first hit value, which was error.) If a coverage filter is specified, the value is increased to 0.001 if it is smaller. -> This has an effect the number of hits one obtains and thus the filtering result.
* The stderr log file now also contains the information how BaitFilter was called. I.e. the exact command line is written to the log file.
* Added the feature that a blast hit coverage filter can be conducted independently of other blast filters.
* Now uses multimap instead of map in CBlastParser.h. Should not make a difference in our context. But makes the CBlastParser more usable.
* Cleaned and improved stats output. Stats were correct, but not easy to read and partially redundant. (Number of bait regions appeared twice. For blast and coverage filter, one result was missing.)
* Added a new thinning mode. This allows us to have multiple bait regions at a specific distance in an alignment or feature.
* Added option: --blast-result-file. A filtering step with BaitFilter can take a long time if a blast analysis is required. Trying different parameters then wastes a lot of time if one and the same blast is done multiple times. With the option a blast result file from a previous run can be reused in the current run. TODO: Do some tests whether the blast file indeed corresponds to the loci-baits.txt file.
