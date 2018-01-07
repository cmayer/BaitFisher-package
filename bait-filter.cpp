/*  BaitFisher (version 1.2.8) a program for designing DNA target enrichment baits
 *  BaitFilter (version 1.0.6) a program for selecting optimal bait regions
 *  Copyright 2013-2017 by Christoph Mayer
 *
 *  This source file is part of the BaitFisher-package.
 * 
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with BaitFisher.  If not, see <http://www.gnu.org/licenses/>.
 *
 *
 *  For any enquiries send an Email to Christoph Mayer
 *  c.mayer.zfmk@uni-bonn.de
 *
 *  When publishing work that is based on the results please cite:
 *  Mayer et al. 2016: BaitFisher: A software package for multi-species target DNA enrichment probe design
 *  
 */
#include <iostream>
#include <fstream>
#include <ctime>

#include "global-types-and-parameters.h"
#include "CBaitRecord.h"
#include "faststring2.h"
#include "CBlastParser.h"
#include "print_container.h"
#include "Ctriple.h"
#include "climits"

using namespace std;

bool     skip_blast = false;
unsigned start_time = 0;
unsigned blast_time = 0;


// Verbosity rules:                                                                                                                  
//      verbosity 0 -> only error messages that lead to exit()                                                                       
//      verbosity 1 -> only warnings
//      verbosity 2 -> information on baits, time used
//      verbosity 3 -> progress                                                                                                      
//      verbosity 4 -> more progress                                                                                                 
//      verbosity

void welcome(std::ostream &os)
{
  os << endl << endl;
  os << "Welcome to " << PROGNAME << ", version " << VERSION << ".\n\n";
}

void good_bye(std::ostream &os)
{
  os << "Your BaitFilter analysis finished successfully.\n"; 
}


void get_hours_minutes_seconds(unsigned secs, unsigned &h, unsigned &m, unsigned &s)
{
  h     =  secs/3600;
  secs  %= 3600;
  m     =  secs/60;
  s     =  secs%60;
}

void print_time(unsigned secs, ostream &os)
{
  char prev = os.fill('0');

  unsigned h, m, s;
  get_hours_minutes_seconds(secs, h, m, s);

  os << h << ":";
  os.width(2);
  os << m << ":";
  os.width(2);
  os << s;

  os.fill(prev);
}

void print_stats(const CBaitLoci_collection &blc, const char * header= NULL)
{
  if (global_verbosity >= 2)
  {
    CBaitRegionSelectorSets brss;
    if (header)
      cout << header << "\n";
    blc.fill_locus_selector(brss);
    brss.print_num_genes_num_features(cout, "  ", blc.has_feature_information());
    
    unsigned NumSeqs, NumBaits;
    blc.get_seqs_baits_count(NumSeqs, NumBaits);
    //   cout << "  Number of bait regions:                      " << blc.get_number_of_loci() << "\n";
    cout << "  Total number of baits:                                     " << NumBaits << "\n";
    cout << "  Total number of bait windows fully covered by sequences    " << NumSeqs  << "\n";
    cout << "  Proportion of baits saved (1-#baits/#sequences):           " << (1.0 -(double)NumBaits/NumSeqs)*100 << "%\n";

    cout.precision(6);
    cout << "  Mean maximum distance of baits to MSAs:                    " << blc.get_max_dist_mean() << "\n";
    cout << "  Overall maximum distance of baits to MSA:                  " << blc.get_overall_max_dist() << "\n";
    cout.precision(2);
    cout << "  Mean CG content of baits:                                  " << blc.get_CG_mean() << "\n";
  }

}


void extract_loci_information(const faststring &header, faststring &alignment_name, faststring &gene_name, unsigned &feature_num, unsigned &start)
{
  faststring rest1;
  faststring rest2;
  faststring tmp;

  header.divide_at('|', alignment_name, rest1);

  rest1.divide_at('|',  gene_name,      rest2);

  rest2.divide_at('|',  tmp,            rest1);

  if (!tmp.isAnUnsigned())
  {
    cerr << "Error: Header does not obey formatting conventions: " << header << "\n";
    cerr << "The third field (separator \"|\") should be an unsigned integer number, but found \""  << tmp << "\" instead.\n";
    exit(-3);
  }
  feature_num = tmp.ToUnsigned();

  rest1.divide_at('|',  tmp,            rest2);

  if (!tmp.isAnUnsigned())
  {
    cerr << "Error: Header does not follow formatting conventions: " << header << "\n";
    cerr << "The fourth field (separator \"|\") should be an unsigned integer number, but found \"" << tmp << "\" instead.\n";
    exit(-3);
  }
  start = tmp.ToUnsigned();

//   rest2.divide_at('|', tmp, rest1);

//   if (!tmp.isAnUnsigned())
//   {
//     cerr << "Error: Header does not follow formatting conventions: " << header << "\n";
//     cerr << "The fifth field (separator \"|\") should be an unsigned integer number, but found \""  << tmp << "\" instead.\n";
//     exit(-3);
//   }

//   if (!rest1.isAnUnsigned())
//   {
//     cerr << "Error: Header does not follow formatting conventions: " << header << "\n";
//     cerr << "The sixth field (separator \"|\") should be an unsigned integer number, but found \""  << rest1 << "\" instead.\n";
//     exit(-3);
//   }

//   tiling_index = tmp.ToUnsigned();
//   counter      = rest1.ToUnsigned();

}


int main(int argc, char ** argv)
{

  start_time = time(NULL);

  cout.setf(ios::fixed);
  cout.precision(2);

  welcome(cout);
  read_and_init_parameters(argc, argv, cerr);

#ifdef SKIPBLAST
  skip_blast = true;
  global_blast_result_file = "blast_result.txt";
#endif

  if (!global_blast_result_file.empty())
  {
    skip_blast = true;
  }


  if (global_verbosity >= 2)
  {
    cerr << '\n';
    print_parameters(cerr, "");
    print_calling_command_line(cerr, argc, argv);
    cerr << '\n';
  }

  // Parameter tests:
  if (0)
  {
    faststring cmd = global_blast_exe + " -version";
    cout << "Blast exe test: " << cmd << "\n";
    int err = system(cmd.c_str());

    cout << err << "\n";

    exit(0);
  }


  ifstream is;
  ofstream os;

  is.open(global_bait_filename.c_str());

  if (is.fail())
  {
    cerr << "Error: Input file: " << global_bait_filename << " could not be opened. I guess it does not exist.\n";
    exit(-1);
  }

  if (global_mode != 'S')
  {
    os.open(global_output_filename.c_str());
    if (os.fail())
    {
      cerr << "\n";
      cerr << "Error: Failed to open or create output file. " << global_output_filename << "\n";
      cerr << "Exiting.\n";
      exit(-1);
    }
  }

  if (global_verbosity >= 3)
  { 
    cout << "PROGRESS: Parsing bait file: ... " << flush;
  }
  CBaitLoci_collection blc(is);
  is.close();

  if (global_verbosity >= 3)
  {
    cout << "Done. Time used since start in seconds: " << (time(NULL)-start_time) << "\n";
  }

  // Note: Bait loci are sorted in the constructor. Usually, they should be sorted in the bait file,
  //       if the user has not edited the file.

  //  DEBUG code:
  //  blc.print_stats(cout);
  //  blc.print_all(cout, 's');
  //  blc.print_all(cout, 'f');
  //  exit(0);


  if (global_mode == 'S')
  {
    faststring msg;

    msg = "Characteristics of the set of bait regions in input file: ";

    if (global_verbosity < 2) // Can be the case if (global_mode == 'S')
      global_verbosity = 2; // Minimum value required by print_stats() to print something. Could be handled with a parameter!
    print_stats(blc, msg.c_str());
    exit(0);
  }
  else if (global_verbosity >= 2)
  {
    faststring msg;
    msg = "Characteristics of the set of bait regions in input file (before filtering): ";
    print_stats(blc, msg.c_str());
  }



  if (global_mode == 'a')
  {
    if (global_verbosity >= 3)
    {
      cout << "\n\n";
      cout << "PROGRESS: Determining: Best bait region in each alignment/gene with criterion minimum number of baits\n";
    }
    CBaitLoci_collection blc_best_baitRegions_in_genes_crit_numBaits(blc, 'g', 'b');
    blc_best_baitRegions_in_genes_crit_numBaits.print_all(os, 's');
    if (global_verbosity >= 2)
    {
      cout << "\n";
      cout << "Characteristics of the set of bait regions after determining the best bait region per alignment/gene under the criterion to minimize the number of baits we have\n";
      print_stats(blc_best_baitRegions_in_genes_crit_numBaits);
      cout << "\n";
    }
  }
  else if (global_mode == 'A')
  {
    if (global_verbosity >= 3)
    {
      cerr << "\n\n";
      cerr << "PROGRESS: Determining: Best bait region in each gene with criterion maximum number of sequences.\n";
    }
    CBaitLoci_collection blc_best_baitRegions_in_genes_crit_numSeqs(blc, 'g', 's');
    blc_best_baitRegions_in_genes_crit_numSeqs.print_all(os, 's');
    if (global_verbosity >= 2)
    {
      cout << "\n";
      cout << "Characteristics of the set of bait regions after determining the best bait region per alignment/gene under the criterion to maximize the number of sequences we have\n";
      print_stats(blc_best_baitRegions_in_genes_crit_numSeqs);
      cout << "\n";
    }
  }
  else if (global_mode == 'f')
  {
    if (global_verbosity >= 3)
    {
      cerr << "\n\n";
      cerr << "PROGRESS: Determining: Best bait region in each feature with criterion minimum number of baits\n";
    }
    CBaitLoci_collection blc_best_baitRegions_in_features_crit_numBaits(blc, 'f', 'b');
    blc_best_baitRegions_in_features_crit_numBaits.print_all(os, 's');
    if (global_verbosity >= 2)
    {
      cout << "\n";
      cout << "After determining the best bait region per feature under the criterion to minimize the number of baits we have\n";
      print_stats(blc_best_baitRegions_in_features_crit_numBaits);
      cout << "\n";
    }
  }
  else if (global_mode == 'F')
  {
    if (global_verbosity >= 3)
    {
      cerr << "\n\n";
      cerr << "PROGRESS: Determining: Best bait region in each feature with criterion maximum number of sequences.\n";
    }
    CBaitLoci_collection blc_best_baitRegions_in_features_crit_numSeqs(blc, 'f', 's');
    blc_best_baitRegions_in_features_crit_numSeqs.print_all(os, 's');
    if (global_verbosity >= 2)
    {
      cout << "\n";
      cout << "Characteristics of the set of bait regions after determining the best bait region per feature under the criterion to maximize the number of sequences we have\n";
      print_stats(blc_best_baitRegions_in_features_crit_numSeqs);
      cout << "\n";
    }
  }
  else if (global_mode == 't' || global_mode == 'T')
  {
    unsigned  sum_seqs, sum_baits;
    
    unsigned  start_count = 0;
    unsigned  best_index  = 1;
    unsigned  best_num    = 0; // This initialisation is only done to silence the compiler warning for using an uninitialised variable.

    // This is essential output: Restrict it by verbosity??
    if (global_verbosity >= 2)
    {
      cout << "Running thinning mode:\n";
      cout << "\nThe columns are: starting offset, sum of sequences the result is based on, number of required baits.\n";
    }
    // Run thinning mode:
    while (start_count < global_thinning_step_width)
    {
      blc.count_seqs_baits_for_start_and_offset(start_count, global_thinning_step_width, sum_seqs, sum_baits);
      ++start_count;

      if (start_count == 1)  // In this block we will surely initialize best_num:
      {
	if (global_mode == 't')
	  best_num = sum_baits;
	else // global_mode == 'T'
	  best_num = sum_seqs;
      }

      if (global_mode == 't')
      {
	if (best_num > sum_baits)
	{
	  best_index = start_count;
	  best_num = sum_baits;
	}
      }
      else // global_mode == 'T'
      {
	if (best_num < sum_seqs)
	{
	  best_index = start_count;
	  best_num = sum_seqs;
	}
      }
      if (global_verbosity >= 2)
      {
	cout << start_count << " " << sum_seqs << " " << sum_baits << "\n";
      }
    }
    
    cout << "Best starting offset: " << best_index << "\n";
    cout << "A bait file that has been thinned out is written to the output file: " << global_output_filename << "\n";
    cout << "This file contains each " << global_thinning_step_width << "th bait region starting with bait-region " << best_index << " of the input bait file.\n"; 
    cout << "Numbers correspond to the line numbers in the input file +1, since all bait files have a headline.\n";

    blc.print_with_start_stepwidth(os, 's', best_index, global_thinning_step_width);
  }
  else if (global_mode == 'h' || global_mode == 'H')
  {
    unsigned  sum_seqs, sum_baits;
    
    unsigned  start_count;
    unsigned  best_index  = UINT_MAX;
    unsigned  best_num    = 0; // This initialisation is only done to silence the compiler warning for using an uninitialised variable.
    unsigned  best_penalty;
    unsigned  penalty_diff;

    std::vector<CBaitLocus*>::iterator it_bait_vector, it_next_alignment_start, it_this_alignment_start;

    // Loop over all alignments or features stored in this collection.
    // Bait loci of all alignments or features are stored this one file.
    // We will handle one alignment/feature a time.
    // The current alignment or feature ends at it_tmp.
    // The file ends at it_end

    unsigned alignment_counter = 0;
    bool     finished          = false;

    // For each alignment
    while (!finished) // The variable finished is being set inside blc.count_seqs_baits_for_start_and_offset_multi_gene_feature_aware
    {
      ++ alignment_counter; // First value in loop is 1.

      // For alignment_counter==1 the variables it_this_alignment_start and it_next_alignment_start are not initialised, but this is intended.
      // They will be initialised with the first call to blc.count_seqs_baits_for_start_and_offset_multi_gene_feature_aware.
      // If alignment_counter > 1, we start the next iteration.
      it_this_alignment_start = it_next_alignment_start;

      if (global_verbosity >= 2)
      {
	faststring original_file_name;
	faststring gene_name;

	blc.get_original_file_name_gene_name(alignment_counter, it_this_alignment_start, original_file_name, gene_name);

	cout << "\n";
	cout << "Running thinning mode for next alignment-file/feature: " << original_file_name << "/" << gene_name << "\n";
	//	cout << "\nThe columns are: starting offset, sum of sequences the result is based on, number of required baits.\n";
      }
      // Run thinning mode:

      start_count  = 0;
      penalty_diff = 0;
      best_penalty = UINT_MAX;

      unsigned best_sum_seqs, best_sum_baits;

      // Test all starting points in this alignment
      while (start_count < global_thinning_step_width)
      {
	// In the following call, it_this_alignment_start is only changed if alignment_counter==1.
	// In this case it is changed to the "beginning of the data set" and therefore to the beginning of the first alignment.
	// Otherwise it remains unchanged.
	// it_next_alignment_start is set to the start of the next alignment.
	blc.count_seqs_baits_for_start_and_offset_multi_gene_feature_aware(alignment_counter,
									   it_this_alignment_start,
									   it_next_alignment_start,
									   finished,
									   start_count,
									   global_thinning_step_width,
									   sum_seqs, sum_baits, penalty_diff);

	if (start_count == 0)  // In this block we will initialise best_num:
	{
	  if (global_mode == 'h')
	  {
	    best_num = sum_baits;
	  }
	  else // global_mode == 'H'
	  {
	    best_num = sum_seqs;
	  }
	  best_penalty = penalty_diff;
	  best_index   = 0;
	  best_sum_seqs  = sum_seqs;
	  best_sum_baits = sum_baits; 
	}
	else
	{
	  if (global_mode == 'h')
	  {
	    // Penalty is the primary criterion to choose the best starting position 
	    if (best_penalty > penalty_diff || (best_penalty == penalty_diff && best_num > sum_baits) )
	    { 
	      best_num = sum_baits;
	      best_index = start_count;
	      best_penalty = penalty_diff;

	      best_sum_seqs  = sum_seqs;
	      best_sum_baits = sum_baits; 
	    }
	  }
	  else // global_mode == 'H'
	  {
	    // Penalty is the primary criterion to choose the best starting position 
	    if (best_penalty > penalty_diff || (best_penalty == penalty_diff && best_num < sum_seqs) )
	    {
	      best_num = sum_seqs;
	      best_index = start_count;
	      best_penalty = penalty_diff;

	      best_sum_seqs  = sum_seqs;
	      best_sum_baits = sum_baits; 
	    }
	  }
	}
	if (global_verbosity >= 5)
	{
	  cout << ": " << start_count << " " << sum_seqs << " " << sum_baits << " " << penalty_diff << "\n";
	  if (global_verbosity >= 7)
	  {
	    cout << "best index, penalty, score: " << best_index << " " << best_penalty << " " << best_num << "\n";
	  }
	}
	++start_count;
      } // END while (start_count < global_thinning_step_width) // Test all starting points in this alignment
    
      cout << "Best starting offset (add 1 to get the index): " << best_index << "\n";
      cout << "Sequences & baits for best starting offset: " << best_sum_seqs << " " << best_sum_baits << "\n";
      cout << "Thinned bait loci have been appended to the output file: " << global_output_filename << "\n";
      //      cout << "This file contains each " << global_thinning_step_width << "th bait region starting with bait-region " << best_index << " of the input bait file." << "\n"; 
      //      cout << "Numbers correspond to the line numbers in the input file +1, since all bait files have a headline." << "\n";

      // Write result for this alignment file/feature, before we go to the next a/f or finish.
      blc.print_with_start_stepwidth_multi_gene_feature_aware(os, 's', best_index, global_thinning_step_width, it_this_alignment_start);

    } // END while (!finished)
  } // END else if (global_mode == 'h' || global_mode == 'H')
  else if (global_mode == 'b' || global_mode == 'B' || global_mode == 'x' || global_mode == 'C')
  {
    faststring blastfilename = "blast_result.txt";

    if (skip_blast)
    {
      blastfilename = global_blast_result_file;
    }

    // Conduct the blast search.
    if (!skip_blast)
    {
      // Write baits to fasta file:
      ofstream os_tmp("tmp_baits.txt");
      blc.print_all(os_tmp, 'f');
      os_tmp.close();

      // Conduct a blast of the baits sequences against a reference genome:

      // In the blast call we use an E-value threshold of:
      //     global_blast_max_second_hit_evalue*2
      //     This will report all relevant hits.

      faststring e_value_command;

      if (!global_blast_evalue_commandline.empty())
      {
	double tmp = global_blast_evalue_commandline.ToDouble();
	if (tmp <= 10 && tmp > 0)
	{
	  e_value_command = global_blast_evalue_commandline;
	}
	else
	{
	  if (global_blast_min_hit_coverage_of_bait > 0)   // For the coverage filter to work, the E-value should not be too small.
	  {
	    double temp = global_blast_max_second_hit_evalue*2;
	    if (temp < 0.001)
	      temp = 0.001;

	    e_value_command = faststring(temp, 15, true);
	  }
	  else
	  {
	    e_value_command = faststring(global_blast_max_second_hit_evalue*2, 15, true);
	  }
	  cerr << "Warning: The Evalue that can be specified with the blast-evalue-cutoff option has to be greater than 0 and smaller or equal to 10."
	       << " The evalue that has been specified is not in this range."
	       << " Therefore we are falling back to the default value: " << e_value_command << "\n";
	}
      }
      else
      {
	if (global_blast_min_hit_coverage_of_bait > 0)   // For the coverage filter to work, the E-value should not be too small.
	{
	  double temp = global_blast_max_second_hit_evalue*2;
	  if (temp < 0.001)
	    temp = 0.001;
	  
	  e_value_command = faststring(temp, 15, true);
	}
	else
	{
	  e_value_command = faststring(global_blast_max_second_hit_evalue*2, 15, true); 
	}
      }
      faststring blastcommand  =   global_blast_exe + " " 
	                         + global_blast_extra_commandline
                                 + " -outfmt 6 -task blastn -query tmp_baits.txt -db " + global_blast_db 
	                         + " -out " + blastfilename + " -evalue " + e_value_command;

      if (global_verbosity >= 2)
      {
	cout << "Running blast with command:\n";
	cout << blastcommand << "\n";
      }

      blast_time = time(NULL);

      //    cerr << "Please press return to continue.\n";
      //  cin.get();

      int blast_command_error;
      blast_command_error = system(blastcommand.c_str());
      
      if (blast_command_error != 0)
      {
	cerr << "ERROR: The blast command returned with non-zero error status. This indicates that blast encountered an error. The status blast returned via the system command was: " << blast_command_error << "\n";
	cerr << "I suggest to run blast command manually:\n";
	cerr << blastcommand << "\n";
	cerr << "Exiting.\n";
	exit(-1);
      }

      blast_time = time(NULL)-blast_time;
    } // END if (!skip_blast)

    if (global_verbosity >= 50)
      cout << "PROGRESS: Reading the blast result file: " << blastfilename << " ..." << flush;
    CBlast_parser bp(blastfilename.c_str());
    if (global_verbosity >= 50)
      cout << " Done.\n";

    // How to use the CBalstParser to filter hits:
    // Mode 0:
    // Filter gene if first hit is below global_blast_max_first_hit_evalue and if second hit is below best_e_value_this_query*factor_max_increase
    // Mode 1:
    // Filter gene if first hit is below global_blast_max_first_hit_evalue and if second hit is below global_blast_max_second_hit_evalue 

    // double_hits contain only one of the 2 good hits that have been found.
    // In case more than 2 good hits are found, it contains all but the first good hit, which is more than what we need.
    // But this redundancy is removed when we add the bait information to the locus_selector data structure.
    // 
    CBlast_parser double_hits; 
    bp.detect_2_good_hits_per_query(global_blast_max_first_hit_evalue, global_blast_max_second_hit_evalue, double_hits, 1);

    set<faststring> set_of_names;
    double_hits.set_of_query_names(set_of_names);

    // The locus_selector contains information about 
    CBaitRegionSelectorSets locus_selector;

    set<faststring>::iterator it, it_end;
    it      = set_of_names.begin();
    it_end  = set_of_names.end();

    //    cout << "Start Time: " << time(NULL) << "\n";

    // For all baits
    while (it != it_end)
    {
      faststring alignment_name;
      faststring gene_name;
      unsigned   feature_num;
      unsigned   start;
//       unsigned   tiling_index;  // Just dummy, will not be used.
//       unsigned   counter;       // Just dummy, will not be used.

// Extract locus information from bait sequence name:
      extract_loci_information(*it, alignment_name, gene_name, feature_num, start);

      if (0)
      {
	if (!locus_selector.exists(alignment_name) )
	  cout << "New alignmentname:                        " << alignment_name << "\n";
	else
	  cout << "Existing alignmentname:                   " << alignment_name << "\n";

	if (!locus_selector.exists(alignment_name, feature_num) )
	  cout << "New alignmentname feature_num:            " << alignment_name << " " << feature_num << "\n";
	else
	  cout << "Existing alignmentname feature_num:       " << alignment_name << " " << feature_num << "\n";

	if (!locus_selector.exists(alignment_name, feature_num, start) )
	  cout << "New alignmentname feature_num start:      " << alignment_name << " " << feature_num << " " << start << "\n";
	else
	  cout << "Existing alignmentname feature_num start: " << alignment_name << " " << feature_num << " " << start << "\n";
      }

      locus_selector.add(alignment_name, feature_num, start);
      ++it;
    } // END For all baits: while (it != it_end)

    if (0)
    {
      cout << "\n";
      cout << "Print locus selectors:\n";
      locus_selector.print_all(cout);
      cout << "\n";
    }

    //    cout << "End Time: " << time(NULL) << "\n";

    // Filter bait loci according to double_hits:
    
    if (global_verbosity >= 100)
    {
      cout << "START: List of double hits:\n";
      double_hits.print_all(cout);
      cout << "END: List of double hits:\n\n";
    }

    // The following code was moved to the top.
//     if (global_verbosity >= 2)
//     {
//               //      cout << "Number of bait regions before searching/filtering for baits with multiple blast hits: " << blc.get_number_of_loci() << "\n";
//       cout << "Characteristics of input file (before filtering): ";
//       print_stats(blc);
//       cout << "\n";
//     }

    // Allowed modes in this constructor:
    // Filter modes: b, B, x:
    // No filter: C. Means that filtered_blc obtains a copy of blc.

    CBaitLoci_collection filtered_blc(blc, locus_selector, global_mode);

    // Print stats with multiple hits results: (Do not print stats in coverage filter only mode. This case is handled below.)
    if (global_verbosity >= 2  && global_mode != 'C')
    {
      //       cout << "Number of bait regions after searching/filtering for baits with multiple blast hits:  " << filtered_blc.get_number_of_loci() << "\n";
      cout << "\n";
      cout << "Characteristics of set of bait regions after filtering for baits with multiple blast hits to reference genome:\n";
      print_stats(filtered_blc);
      cout << "\n";
    }

    if (global_blast_min_hit_coverage_of_bait > 0)
    {
      CBaitLoci_collection coverage_filter(filtered_blc, bp, global_blast_min_hit_coverage_of_bait);
//       if (global_verbosity >= 2)
//       {
// 	cout << "Number of bait regions after hit coverage filtering:  " << coverage_filter.get_number_of_loci() << "\n";
//       }
      coverage_filter.print_all(os, 's');

      if (global_verbosity >= 2)
      {
	if (global_blast_min_hit_coverage_of_bait > 0)
	{
	  cout << "\n";
	  cout << "Characteristics of set of bait regions after filtering for bait regions without sufficient hit coverage values in blast against a reference genome:\n";
	  print_stats(coverage_filter);
	  cout << "\n";
	}
      }
    }
    else
    {
      // Report what we have found:
      filtered_blc.print_all(os, 's');
    }
  } // END if (global_mode == 'b' || global_mode == 'B' || global_mode == 'x' || global_mode == 'C')
  else if (global_mode == 'c')
  {
    if (global_conversion_mode == 1)
    {
      // Write baits to file in Agilent-Germany format
      blc.print_all(os, 'A', global_ProbeID_prefix.c_str());

      if (global_verbosity >= 2)
      {
	cout << "Bait statistics for the converted bait file:\n";
	print_stats(blc);
	cout << "\n";
      }
    }
    else
    {
      cerr << "INTERNAL ERROR: Unknown conversion mode : " << global_conversion_mode << "\n";
      cerr << "Please report this problem to the developer.\n";
      exit(-13);     
    }
  } // END if (global_mode == 'c')
  else
  {
    cerr << "INTERNAL ERROR: Unknown internal mode: " << global_mode << "\n";
    cerr << "Please report this problem to the developer.\n";
    exit(-13);
  }

  unsigned total_seconds = time(NULL)-start_time;
  if (global_verbosity >= 2)
  {
    cout << "Total time used: ";
    print_time(total_seconds, cout);
    cout << "\n";
  }

  if (global_mode == 'b' || global_mode == 'B' || global_mode == 'x' || global_mode == 'C')
  {
    if (global_verbosity >= 2)
    {
      cout << "Time used excluding the Blast search: ";
      print_time(total_seconds-blast_time, cout);
      cout << "\n";
    }
  }

  os.close();
  good_bye(cout);
}
