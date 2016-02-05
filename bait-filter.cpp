/*  BaitFisher (version 1.2.7) a program for designing DNA target enrichment baits
 *  Copyright 2013-2016 by Christoph Mayer
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
#include "global-types-and-parameters.h"
#include "CBaitRecord.h"
#include "faststring2.h"
#include "CBlastParser.h"
#include "print_container.h"
#include "Ctriple.h"

#include <ctime>

using namespace std;


bool skip_blast = false;
unsigned start_time=0;
unsigned blast_time=0;


// Verbosity rules:                                                                                                                  
//      verbosity 0 -> only error messages that lead to exit()                                                                       
//      verbosity 1 -> only warnings
//      verbosity 2 -> information on baits, time used
//      verbosity 3 -> progress                                                                                                      
//      verbosity 4 -> more progress                                                                                                 


// The very first version of this program was called filter-double-hits.cpp.
//  double max_e               = 1.0/1000000.0;
//  double factor_max_increase = 1000;

//  faststring blast_db_name = "/Volumes/Mac2/Daten2/Bait-Fisher-Tests/05-2015/Baitfishing_Vespidae/Apis_mellifera_genome/Amel_4.5_scaffolds.fa";

void welcome(std::ostream &os)
{
  os << endl << endl;
  os << "Welcome to " << PROGNAME << ", version " << VERSION << "."
     << endl << endl << endl;
}

void good_bye(std::ostream &os)
{
  os << "Your BaitFilter analysis finished successfully." << endl; 
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
  unsigned h, m, s;
  get_hours_minutes_seconds(secs, h, m, s);
  os << h << ":" << m << ":" << s;
}

void print_stats(const CBaitLoci_collection &blc, const char * header= NULL)
{
  if (global_verbosity >= 2)
  {
    CBaitRegionSelectorSets brss;
    if (header)
      cout << header << endl;
    blc.fill_locus_selector(brss);
    brss.print_num_genes_num_features(cout, "  ");
    
    unsigned NumSeqs, NumBaits;
    blc.get_seqs_baits_count(NumSeqs, NumBaits);

    cout << endl;
    cout << "  Total number of baits:             " << NumBaits << endl;
    cout << "  Total number of sequences at loci: " << NumSeqs  << endl;
    cout << "  Proportion of baits saved:         " << (1.0 -(double)NumBaits/NumSeqs)*100 << "%" << endl;

    cout.precision(6);
    cout << "  Mean maximum distance of baits to MSAs:    " << blc.get_max_dist_mean() << endl;
    cout << "  Overall maximum distance of bait to MSA:   " << blc.get_overall_max_dist() << endl;
    cout.precision(2);
    cout << "  Mean CG content of baits:               " << blc.get_CG_mean() << endl;
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
    cerr << "Error: Header does not follow formating conventions: " << header << endl;
    cerr << "The third field (separator \"|\") should be an unsigned integer number, but found \""  << tmp << "\" instead." << endl;
    exit(-3);
  }
  feature_num = tmp.ToUnsigned();

  rest1.divide_at('|',  tmp,            rest2);

  if (!tmp.isAnUnsigned())
  {
    cerr << "Error: Header does not follow formating conventions: " << header << endl;
    cerr << "The fourth field (separator \"|\") should be an unsigned integer number, but found \"" << tmp << "\" instead." << endl;
    exit(-3);
  }
  start = tmp.ToUnsigned();

//   rest2.divide_at('|', tmp, rest1);

//   if (!tmp.isAnUnsigned())
//   {
//     cerr << "Error: Header does not follow formating conventions: " << header << endl;
//     cerr << "The fifth field (separator \"|\") should be an unsigned integer number, but found \""  << tmp << "\" instead." << endl;
//     exit(-3);
//   }

//   if (!rest1.isAnUnsigned())
//   {
//     cerr << "Error: Header does not follow formating conventions: " << header << endl;
//     cerr << "The sixth field (separator \"|\") should be an unsigned integer number, but found \""  << rest1 << "\" instead." << endl;
//     exit(-3);
//   }

//   tiling_index = tmp.ToUnsigned();
//   counter      = rest1.ToUnsigned();

}




int main(int argc, char ** argv)
{

  start_time = time(NULL);

#ifdef SKIPBLAST
skip_blast = true;
#endif

  cout.setf(ios::fixed);
  cout.precision(2);

  if (global_verbosity >= 2)
  {
    welcome(cout);
  }
  read_and_init_parameters(argc, argv);

  if (global_verbosity >= 2)
  {
    putc('\n', stderr);
    print_parameters(stderr,"");
    putc('\n', stderr);
  }

  // Parameter tests:
  if (0)
  {
    faststring cmd = global_blast_exe + " -version";
    cout << "Blast exe test: " << cmd << endl;
    int err = system(cmd.c_str());

    cout << err << endl;

    exit(0);
  }


  ifstream is;
  ofstream os;

  is.open(global_bait_filename.c_str());

  if (is.fail())
  {
    cerr << "Error: Input file: " << global_bait_filename << " could not be opened. I guess it does not exist." << endl;
    exit(-1);
  }

  if (global_mode != 'S')
  {
    os.open(global_output_filename.c_str());
    if (os.fail())
    {
      cerr << endl;
      cerr << "Error: Failed to open or create output file. " << global_output_filename << endl;
      cerr << "Exiting." << endl;
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
    cout << "Done. Time used since start in seconds: " << (time(NULL)-start_time) << endl;
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
    faststring msg = "Stats for input file";


    if (global_mode < 3)
      global_verbosity = 3; // Minimum value required by print_stats() to print something. Could be handled with a parameter!
    print_stats(blc, msg.c_str());
    exit(0);
  }


  if (global_mode == 'a')
  {
    if (global_verbosity >= 3)
    {
      cout << endl << endl;
      cout << "PROGRESS: Determining: Best bait region in each alignment/gene with criterion minimum number of baits" << endl;
    }
    CBaitLoci_collection blc_best_baitRegions_in_genes_crit_numBaits(blc, 'g', 'b');
    blc_best_baitRegions_in_genes_crit_numBaits.print_all(os, 's');
    if (global_verbosity >= 2)
    {
      cout << "After determining the best bait region per alignment/gene under the criterion to minimize the number of baits we have" 
	   << endl;
      print_stats(blc_best_baitRegions_in_genes_crit_numBaits);
      cout << endl;
    }
  }
  else if (global_mode == 'A')
  {
    if (global_verbosity >= 3)
    {
      cerr << endl << endl;
      cerr << "PROGRESS: Determining: Best bait region in each gene with criterion maximum number of sequences." << endl;
    }
    CBaitLoci_collection blc_best_baitRegions_in_genes_crit_numSeqs(blc, 'g', 's');
    blc_best_baitRegions_in_genes_crit_numSeqs.print_all(os, 's');
    if (global_verbosity >= 2)
    {
      cout << "After determining the best bait region per alignment/gene under the criterion to maximize the number of sequences we have" 
	   << endl;
      print_stats(blc_best_baitRegions_in_genes_crit_numSeqs);
      cout << endl;
    }
  }
  else if (global_mode == 'f')
  {
    if (global_verbosity >= 3)
    {
      cerr << endl << endl;
      cerr << "PROGRESS: Determining: Best bait region in each feature with criterion minimum number of baits" << endl;
    }
    CBaitLoci_collection blc_best_baitRegions_in_features_crit_numBaits(blc, 'f', 'b');
    blc_best_baitRegions_in_features_crit_numBaits.print_all(os, 's');
    if (global_verbosity >= 2)
    {
      cout << "After determining the best bait region per feature under the criterion to minimize the number of baits we have" 
	   << endl;
      print_stats(blc_best_baitRegions_in_features_crit_numBaits);
      cout << endl;
    }
  }
  else if (global_mode == 'F')
  {
    if (global_verbosity >= 3)
    {
      cerr << endl << endl;
      cerr << "PROGRESS: Determining: Best bait region in each feature with criterion maximum number of sequences." << endl;
    }
    CBaitLoci_collection blc_best_baitRegions_in_features_crit_numSeqs(blc, 'f', 's');
    blc_best_baitRegions_in_features_crit_numSeqs.print_all(os, 's');
    if (global_verbosity >= 2)
    {
      cout << "After determining the best bait region per feature under the criterion to maximize the number of sequences we have"
	   << endl;
      print_stats(blc_best_baitRegions_in_features_crit_numSeqs);
      cout << endl;
    }
  }
  else if (global_mode == 't' || global_mode == 'T')
  {
    unsigned  sum_seqs, sum_baits;
    
    unsigned  start_count = 0;
    unsigned  best_index  = 1;
    unsigned  best_num    = 0; // This initialistion is only done to silence the compiler warning for using an unitialized variable.

    // This is essential output: Restrict it by verbostiy??
    if (global_verbosity >= 2)
    {
      cout << "Running thinning mode:" << endl;
      cout << "\nThe columns are: starting offset, sum of sequences the result is baised on, number of required baits." << endl;
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
	else
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
      else
      {
	if (best_num < sum_seqs)
	{
	  best_index = start_count;
	  best_num = sum_seqs;
	}
      }
      if (global_verbosity >= 2)
      {
	cout << start_count << " " << sum_seqs << " " << sum_baits << endl;
      }
    }
    
    cout << "Best starting offset: " << best_index << endl;
    cout << "A bait file that has been thinned out is written to the output file: " << global_output_filename << endl;
    cout << "This file contains each " << global_thinning_step_width << "th bait region starting with bait-region " << best_index << " of the input bait file." << endl; 
    cout << "Numbers correspond to the line numbers in the input file +1, since all bait files have a headline." << endl;

    blc.print_with_start_stepwidth(os, 's', best_index, global_thinning_step_width);
  }
  else if (global_mode == 'b' || global_mode == 'B' || global_mode == 'x' )
  {
    faststring blastfilename = "blast_result.txt";

    if (!skip_blast)
    {
      // Write baits to fasta file:
      ofstream os_tmp("tmp_baits.txt");
      blc.print_all(os_tmp, 'f');
      os_tmp.close();

      // Conduct a blast of the baits sequences against a reference genome:

      // In the blast call we use an evalue threshold of:
      //     global_blast_max_second_hit_evalue*2
      //     This will report all relevant hits.

      faststring e_value_command;

      if (!global_blast_evalue_commandline.empty())
      {
	double tmp = global_blast_evalue_commandline.ToDouble();
	if (tmp <= 10 && tmp > 0)
	{
	  e_value_command = faststring(tmp, 15, true);
	}
	else
	{
	  e_value_command = faststring(global_blast_max_second_hit_evalue*2, 15, true); 
	  cerr << "Warning: The Evalue that can be specified with the blast-evalue-cutoff option has to be greater than 0 and smaller or equal to 10."
	       << " The evalue that has been specified is not in this range."
	       << " Therefore we are falling back to the default value: " << e_value_command << endl;
	}
      }
      else
      {
	  e_value_command = faststring(global_blast_max_second_hit_evalue*2, 15, true); 
      }
      faststring blastcommand  = "blastn -outfmt 6 -task blastn -query tmp_baits.txt -db " + global_blast_db + " -out " + blastfilename + " -evalue " + e_value_command;

      if (!global_blast_extra_commandline.empty())
      {
	blastcommand += " " + global_blast_extra_commandline;
      }

      if (global_verbosity > 50)
      {
	cout << "Running blast with command: " << endl;
	cout << blastcommand << endl;
      }

      blast_time = time(NULL);

      //    cerr << "Please press return to continue." << endl;
      //  cin.get();

      int blast_command_error;
      blast_command_error = system(blastcommand.c_str());
      
      if (blast_command_error != 0)
      {
	cerr << "ERROR: The blast command returned with non-zero error status. This indicates that blast encountered an error. The status blast returned via the system command was: " << blast_command_error << endl;
	cerr << "I suggest to run blast command manually: " << endl;
	cerr << blastcommand << endl;
	cerr << "Exiting." << endl;
	exit(-1);
      }

      blast_time = time(NULL)-blast_time;
    }

    if (global_verbosity > 50)
      cout << "PROGRESS: Reading the blast result file: " << blastfilename << " ..." << flush;
    CBlast_parser bp(blastfilename.c_str());
    if (global_verbosity > 50)
      cout << " Done." << endl;

    // How to use the CBalstParser to filter hits:
    // Mode 0:
    // Filter gene if first hit is below global_blast_max_first_hit_evalue and if second hit is below best_e_value_this_query*factor_max_increase
    // Mode 1:
    // Filter gene if first hit is below global_blast_max_first_hit_evalue and if second hit is below global_blast_max_second_hit_evalue 

    // double_hits contain only one of the 2 good hits that have been found.
    // In case more than 2 good hits are found, it contains all but the first good hit, which is more than what we need.
    // But this redudancy is removed when we add the bait information to the locus_selector data structure.
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

    //    cout << "Start Time: " << time(NULL) << endl;

    while (it != it_end)
    {
      faststring alignment_name;
      faststring gene_name;
      unsigned   feature_num;
      unsigned   start;
//       unsigned   tiling_index;  // Just dummy, will not be used.
//       unsigned   counter;       // Just dummy, will not be used.

      extract_loci_information(*it, alignment_name, gene_name, feature_num, start);

      if (0)
      {
	if (!locus_selector.exists(alignment_name) )
	  cout << "New alignmentname:                        " << alignment_name << endl;
	else
	  cout << "Existing alignmentname:                   " << alignment_name << endl;

	if (!locus_selector.exists(alignment_name, feature_num) )
	  cout << "New alignmentname feature_num:            " << alignment_name << " " << feature_num << endl;
	else
	  cout << "Existing alignmentname feature_num:       " << alignment_name << " " << feature_num << endl;

	if (!locus_selector.exists(alignment_name, feature_num, start) )
	  cout << "New alignmentname feature_num start:      " << alignment_name << " " << feature_num << " " << start << endl;
	else
	  cout << "Existing alignmentname feature_num start: " << alignment_name << " " << feature_num << " " << start << endl;
      }

      locus_selector.add(alignment_name, feature_num, start);
      ++it;
    }

    if (0)
    {
      cout << endl;
      cout << "Print locus selectors:" << endl;
      locus_selector.print_all(cout);
      cout << endl;
    }

    //    cout << "End Time: " << time(NULL) << endl;

    // Filter bait loci according to double_hits:
    
    if (global_verbosity >= 100)
    {
      cout << "START: List of double hits:" << endl;
      double_hits.print_all(cout);
      cout << "END: List of double hits:" << endl;
      cout << endl;
    }

    if (global_verbosity >= 2)
    {
      cout << "Number of bait regions before searching/filtering for baits with multiple blast hits: " << blc.get_number_of_loci() << endl;
      cout << "Before filtering: ";
      print_stats(blc);
      cout << endl;
    }

    CBaitLoci_collection filtered_blc(blc, locus_selector, global_mode);

    if (global_verbosity >= 2)
    {
      cout << "Number of bait regions after searching/filtering for baits with multiple blast hits:  " << filtered_blc.get_number_of_loci() << endl;
    }

    if (global_blast_min_hit_coverage_of_bait > 0)
    {
    // Now filter bait regions to ensure minimum coverage of hits:
      CBaitLoci_collection coverage_filter(filtered_blc, bp, global_blast_min_hit_coverage_of_bait);
      if (global_verbosity >= 2)
      {
	cout << "Number of bait regions after hit coverage filtering:  " << coverage_filter.get_number_of_loci() << endl;
      }
      coverage_filter.print_all(os, 's');
    }
    else
    {
      // Report what we have found:
      filtered_blc.print_all(os, 's');
    }
    
    // Print collectively for all modes (global_mode == 'b' || global_mode == 'B' || global_mode == 'x'  ):
    if (global_verbosity >= 2)
    {
      cout << "Bait statistics after filtering: ";
      print_stats(filtered_blc);
      cout << endl;
    }

  } // END if (global_mode == 'b' || global_mode == 'B' || global_mode == 'x'  )
  else if (global_mode == 'c')
  {
    if (global_conversion_mode == 1)
    {
      // Write baits to file in Agilent-Germany format
      blc.print_all(os, 'A', global_ProbeID_prefix.c_str());

      if (global_verbosity >= 2)
      {
	cout << "Bait statistics for the converted bait file:" << endl;
	print_stats(blc);
	cout << endl;
      }
    }
    else
    {
      cerr << "INTERNAL ERROR: Unknown converson mode : " << global_conversion_mode << endl;
      cerr << "Please report this problem to the developer." << endl;
      exit(-13);     
    }
  } // END if (global_mode == 'c')
  else
  {
    cerr << "INTERNAL ERROR: Unknown internal mode: " << global_mode << endl;
    cerr << "Please report this problem to the developer." << endl;
    exit(-13);
  }

  unsigned total_seconds = time(NULL)-start_time;
  if (global_verbosity >= 2)
  {
    cout << "Total time used: ";
    print_time(total_seconds, cout);
    cout << endl;
  }

  if (global_mode == 'b' || global_mode == 'B' || global_mode == 'x' )
  {
    if (global_verbosity >= 2)
    {
      cout << "Time used excluding the Blast search: ";
      print_time(total_seconds-blast_time, cout);
      cout << endl;
    }
  }

  os.close();
  good_bye(cout);
}
