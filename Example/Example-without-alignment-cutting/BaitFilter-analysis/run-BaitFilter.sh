BaitFilter-v1.0.6  -m ab -i ../BaitFisher-results/loci_baits.txt  -o baits_filtered_best-in-alignment_minNumBaits.txt  1> console_o_filter_bestBaits 2> console_e_filter_bestBaits
BaitFilter-v1.0.6  -m as -i ../BaitFisher-results/loci_baits.txt  -o baits_filtered_best-in-alignment_maxNumSeq.txt    1> console_o_filter_bestSeqs  2> console_e_filter_bestSeqs

## Produce bait file that can be uploaded at a bait producing company. So war all companies accepted the file produced with the following command.
BaitFilter-v1.0.6 -c Agilent-Germany -i baits_filtered_best-in-alignment_maxNumSeq.txt   -o Final_baits_filtered_best-in-alignment_maxNumSeq.txt   --ID-prefix My-project-name- 1> conversion_maxNumSeq.o   2> conversion_maxNumSeq.e
BaitFilter-v1.0.6 -c Agilent-Germany -i baits_filtered_best-in-alignment_minNumBaits.txt -o Final_baits_filtered_best-in-alignment_minNumBaits.txt --ID-prefix My-project-name- 1> conversion_minNumBaits.o 2> conversion_minNumBaits.e