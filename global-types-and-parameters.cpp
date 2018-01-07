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
#include "global-types-and-parameters.h"
#include "tclap/CmdLine.h"
#include <string>
#include "faststring2.h"

using namespace TCLAP;
using namespace std;

faststring                       global_bait_filename;
char                             global_mode;
faststring                       global_blast_exe;
faststring                       global_output_filename;

double                           global_blast_min_hit_coverage_of_bait;
double                           global_blast_max_first_hit_evalue;
double                           global_blast_max_second_hit_evalue;
faststring                       global_blast_db;

faststring                       global_blast_extra_commandline;
faststring                       global_blast_evalue_commandline;

unsigned                         global_thinning_step_width;

bool                             global_use_GUI;
bool                             global_stats_mode;

unsigned                         global_conversion_mode;
faststring                       global_ProbeID_prefix;
unsigned                         global_verbosity;
faststring                       global_blast_result_file;


bool file_exists(const char* fn)
{
  ifstream is;
  is.open(fn);
  if (is.fail())
  {
    return false;
  }
  is.close();
  return true;
}


bool blast_exe_exists()
{
  faststring cmd = global_blast_exe + " -version 1> /dev/null 2> /dev/null";
  int err = system(cmd.c_str());

  if (err == 0)
    return true;
  else
    return false;
}

bool blast_db_files_exist(faststring &blast_db_base_name)
{
  faststring tmp;

  tmp = blast_db_base_name + ".nhr";
  if (!file_exists(tmp.c_str()))
    return false;
  tmp = blast_db_base_name + ".nin";
  if (!file_exists(tmp.c_str()))
    return false;
  tmp = blast_db_base_name + ".nsq";
  if (!file_exists(tmp.c_str()))
    return false;
  return true;
}


void good_bye_and_exit(FILE *of, int error)
{
  if (error == 0)
    fprintf(of, "#\n#\n## Finished successfully. ##\n");
  fprintf(of, "\n" PROGNAME " says goodbye.\n\n");
  exit(error);
}

void good_bye_and_exit(ostream &os, int error)
{
  if (error == 0)
    os << "#\n#\n## Finished successfully. ##\n";
  os << "\n" PROGNAME " says goodbye.\n\n";
  exit(error);
}



void init_param()
{
  global_blast_exe                      = "blastn";
  global_blast_min_hit_coverage_of_bait = 0;
  global_blast_max_first_hit_evalue     = 0.000001;
  global_blast_max_second_hit_evalue    = 0.000001;
  global_thinning_step_width            = 0;

  global_use_GUI         = false;
  global_stats_mode      = false;

  global_mode            = 0; // 0: Undefined, 'c': conversion, ....
  global_conversion_mode = 0; // 0: No conversion, 1: Agilent-Germany

  global_verbosity       = 2;
}



void read_and_init_parameters(int argc, char** argv, ostream &logerr)
{
  init_param();

  try
  {
    CmdLine cmd("The " PROGNAME  " program has been designed to post process the output of the BaitFisher program in order select appropriate bait regions and to create the final bait set. "
                "BaitFilter offers several filtering and conversion modes. If multiple filtering steps and a final conversion are required, BaitFilter will have to be started multiple times and the output "
		"of the different runs are used as input in the next step.\n"
                "The BaitFisher program designs baits for every locus for which a bait design is possible for a full bait region. A bait region can start at every nucleotide as long as the remaining sequence is long enough. "
		"This output has to be reduced and the purpose of BaitFilter is to find for each feature, gene or alignment the optimal locus or the optimal loci for the bait regions. "
		"Before determining the locus with the fewest number of baits or the largest sequence coverage, one might want to determine which baits are expected to bind specifically in a given reference genome. "
                "This is achieved by conducting a Blast search of the baits against a genome. Baits which are highly similar to at least two loci of the genome can be determined and their bait regions can be removed. "
                "The blast search result can also be used to specify a minimum hit coverage of the baits in a bait region against the reference genome. "
		"After removing bait regions at inferior loci, the optimal bait region starting locus (start coordinate) can be inferred with the aid of different criteria in a subsequent run of BaitFilter. "
                "As input, BaitFilter requires a bait file generated by the BaitFisher program or a BaitFile generated by a previous filtering run of BaitFilter. "
		"This bait file is specified with the -i command line parameter (see below). Furthermore, the user has to specify an output file name with the -o parameter and a filter mode with "
		"the -m parameter.\n"
		"To convert a file to final and uploadable output format, see the -c option below.\n"
		"To compute a bait file statistics of an input file, see the -S option below.\n" 
		"The different filter modes provided by BaitFilter are the following:\n" 
		"1a) Retain only the best bait locus per alignment file. Criterion: Minimize number of required baits.\n"
		"1b) Retain only the best bait locus per alignment file. Criterion: Maximize number of sequences.\n"
		"2a) Retain only best bait locus per feature (requires that features were selected in BaitFisher). Criterion: Minimize number of required baits.\n"
		"2b) Retain only best bait locus per feature (requires that features were selected in BaitFisher). Criterion: Maximize number of sequences.\n"
		"3) Use a blast search of the bait sequences against a reference genome to detect putative non-unique target loci. Non unique target sites will have multiple good hits against the reference genome. "
		"  Furthermore, a minimum coverage of the best blast hit of bait sequence against the genome can be specified. Note that all blast modes require additional command line parameters!"
		"  These modes remove bait regions for which multiple good blast hits where found or for which baits have insufficiently long hits. Different versions of this mode are available:\n" 
		" 3a) If a single bait is not unique, remove all bait regions from the current gene.\n" 
		" 3b) If a single bait is not unique, remove all bait regions from the current feature (if applicable).\n"
		" 3c) If a single bait is not unique, remove only the bait region that contains this bait.\n"
		" 4) Thin out the given bait file: Retain only every Nth bait region, where N has to be specified by the user. Two submodes are available:\n" 
		" 4a) Thin out bait regions by retaining only every Nth bait region in a bait file. The starting offset will by chosen such that the number of required baits is minimized.\n"
		" 4b) Thin out bait regions by retaining only every Nth bait region in a bait file. The starting offset will by chosen such that the number of sequences the result is based on is maximized.\n",
		' ', VERSION);

    
    ValueArg<string> blast_result_Arg("b", "blast-result-file",
       "Conducting a blast analysis of all baits against a reference genome can take a long time. If different filtering parameters, "
       "e.g. different coverage thresholds are to be compared, the same blast has to be done multiple times. "
       "With this argument, the blast will be skipped and the specified blast result file will be used. "
       "This option has to be used with caution! No checks are done (so far) to ensure that the blast result file "
       "corresponds to the specified bait file. "
       "If a BaitFilter run was conducted which did a blast search, BaitFilter will not delete the blast result file "
       "after the run was completed. "
       "The result file with the name blast_result.txt will remain in the working directory. "
       "It can be moved or renamed and with this option it can be specified as the input file for further BaitFilter runs. "
       "If you have the slightest doubt whether you are using the correct blast result file, you should not use this option. "
       "This option is only allowed in modes that would normally do a blast search. This option cannot be specified together "
       "with the blast-executable, blast-evalue-cutoff, blast-extra-commandline, ref-blast-db options, "
       "since these are options specific to runs in which a blast search is conducted.",
	false, "", "string");
    cmd.add( blast_result_Arg );

    ValueArg<unsigned> verbosity_Arg("", "verbosity",
				     "The verbosity option controls the amount of information " PROGNAME 
				     " writes to the console while running. 0: Print only welcome message and essential error messages that lead to exiting the program. "
				     "1: report also warnings, 2: report also progress, 3: report more detailed progress, >10: debug output. Maximum 10000: write all possible diagnostic output."
				     " A value of 2 is required if startup parameters should be reported.\n",
				     false, global_verbosity, "unsigned integer");
    cmd.add( verbosity_Arg );

    SwitchArg global_stats_Arg("S", "stats", "Compute bait file characteristics for the input file and report these. This mode is automatically used for all modes specified with -m option or the conversion mode "
			       "specified with -c option. "
			       "The purpose of the -S option is to compute stats without having to filter or convert the input file. In particular, the -S mode does not require specifying an output file.\n"
			       "This option has no effect if combined with the -m or -c modes.\n", false);
    cmd.add( global_stats_Arg );


    ValueArg<string> conversion_IDprefix_Arg("", "ID-prefix",
        "In the conversion mode to the four-column-upload file format, each converted file should get a unique ProbeID prefix, since even among multiple files, "
        "ProbeIDs are not allowed to be identical. With this option the user is able to specify a prefix string to all probe IDs in the four-column-upload file created by BaitFilter.",
	false, "", "string");
    cmd.add( conversion_IDprefix_Arg );


    ValueArg<unsigned> thinning_step_width_Arg("t", "thinning-step-width",
					      "Thin out the bait file by retaining only every Nth bait region. The integer after the option specifies the step width N. If one of the modes thin-b (thin-b-old), or" 
					      " thin-s (thin-s-old) is active, this option is required, otherwise it is not allowed to set this parameter.",
					      false, global_thinning_step_width, "positive integer");
    cmd.add( thinning_step_width_Arg );

    ValueArg<string> blast_exe_Arg("B", "blast-executable",
       "Name of or path+name to the blast executable. Default: blastn. Minimum blast version number: Blast+ 2.2.x. Default: blastn. Cannot be specified together with the blast-result-file option. ",
       false, global_blast_exe.c_str(), "string");
    cmd.add( blast_exe_Arg );

    ValueArg<string> blast_evalue_extra_Arg("", "blast-evalue-cutoff",
       "When conducting a blast search, a maximum E-value can be specified when calling the blast program. "
       "The effect is that hits with a higher E-value are not reported. "
       "BaitFilter always specifies such an E-value when calling the blast program. The default E-value passed by BaitFilter "
       "to the blast program is twice the --blast-second-hit-evalue. "
       "If a coverage filter is requested the default value is set to 0.001 if twice the value of "
       "--blast-second-hit-evalue is smaller than 0.001. "
       "This should guarantee that all hits necessary for the blast and/or coverage filter are found. "
       "If the user wants to set a different E-value threshold, this can be specified with this option. "
       "With version 1.0.6 of this program, the value is automatically changed to be larger or equal to 0.001 if the coverage filter is used. This makes the usage of this option unnecessary in most cases.",
       false, "", "floating point number");
    cmd.add( blast_evalue_extra_Arg );

    ValueArg<string> blast_extra_commandline_Arg("", "blast-extra-commandline",
       "When invoking the blast command, extra command line parameters can be passed to the blast program with the aid of this option. As an example, "
       "this option allows to specify the number of threads the blast program should use. Example: --blast-extra-commandline \"-num_threads 20\" sets the number of threads to 20. ",
       false, "", "string");
    cmd.add( blast_extra_commandline_Arg );

    ValueArg<string> blast_db_Arg("", "ref-blast-db",
       "Base name to a blast data base file. This name is passed to the blast command. This is the name of the fasta file of your reference genome."
       " IMPORTANT: The makeblastdb program has to be called before starting the " PROGNAME " program. makeblastdb takes the fasta file and "
       " creates data base files out of it. Cannot be specified together with the blast-result-file option.",
	false, "", "string");
    cmd.add( blast_db_Arg );

    ValueArg<double> blast_hit_coverage_Arg("", "blast-min-hit-coverage-of-baits-in-tiling-stack",
       "Can be specified together with the following modes (-m option): blast-a, blast-f, blast-l, blast-c. In all these modes, a blast analysis of all baits against a reference genome is conducted. "
       "This option specifies a minimum query hit coverage which at least one bait has to have in each tiling stack (i.e. the column in the tiling design). Otherwise the bait region is discarded. "
       "If not specified, no hit coverage is checked. The coverage is determined for each bait by dividing the length of the best hit of this bait against the specified genome "
       "by the length of this bait. Then the highest coverage is determined for each bait stack of the tiling design. "

       " If this option is used together with another filter, it is important to know the order in which the two are applied, since the order matters for the final result:"
					    //       "For the mode options: ab, as, fb, fs the coverage is checked before "
					    //       "determining the optimal bait region. "
       "For the mode options: blast-a, blast-f, blast-l the hit coverage is checked after filtering for baits with multiple "
       "good hits to the reference genome.", false, global_blast_min_hit_coverage_of_bait, "floating point number");
    cmd.add( blast_hit_coverage_Arg );

    ValueArg<double> blast_primary_hit_evalue_Arg("", "blast-first-hit-evalue",
						  "Maximum E-value for the first or best hit of the bait against the genome. A bait is characterized to bind ambiguously, if we have at least two good hits to different loci of the genome. "
						  "This option is the E-value threshold for the first/best hit. Default: 0.000001",
	false, global_blast_max_first_hit_evalue, "floating point number");
    cmd.add( blast_primary_hit_evalue_Arg );

    ValueArg<double> blast_secondary_hit_evalue_Arg("", "blast-second-hit-evalue",
						    "Maximum E-value for the second or second best hit. A bait is characterised to bind ambiguously, if we have at least two good hits. This option is the E-value threshold for the second best hit to different loci of the genome."
						  "This option is the E-value threshold for the second best hit. Default: 0.000001",
	false, global_blast_max_second_hit_evalue, "floating point number");
    cmd.add( blast_secondary_hit_evalue_Arg );

 
    // The most important options must be mentioned first. Not they are specified in (reverse) order:

    // (#4) in option list 
    ValueArg<string> bait_filter_mode_Arg("m","mode",
				      "Apart from the input file option, the mode option is the most important option. This option specifies which filter mode BaitFilter uses. (See the user manual for more details):\n"
				      "\"ab\":      Retain only the best bait locus for each alignment file when using the optimality criterion"
				      "             to minimize the total number of required baits.\n"
				      "\"as\":      Retain only the best bait locus for each alignment file when using the optimality criterion"
				      "             to maximize the number of sequences the result is based on.\n"
				      "\"fb\":      Retain only the best bait locus for each feature (e.g. CDS) when using the optimality criterion"
				      "             to minimize the total number of required baits. Only applicable if alignment cutting has been used in BaitFisher.\n"
				      "\"fs\":      Retain only the best bait locus for each feature (e.g. CDS) when using the optimality criterion"
				      "             to maximize the number of sequences the result is based on. Only applicable if alignment cutting has been used in BaitFisher. \n"
				      "\"blast-a\": Remove all bait regions of all ALIGNMENTs for which one or more baits have at least two good hits to a reference genome. (Not recommended.)\n"
				      "\"blast-f\": Remove all bait regions of all FEATUREs for which one or more baits have at least two good hits to a reference genome. (Not recommended.)\n"
				      "\"blast-l\": Remove only the bait REGIONs that contain a bait that has multiple good hits to a reference genome. (Recommended over blast-f and blast-a.)\n"
				      "\"blast-c\": Conduct a coverage filter run without a search for multiple hits. Requires the blast-min-hit-coverage-of-baits-in-tiling-stack option to be specified.\n"
                                      "\"thin-b\":  Thin out a bait file to every Nth bait region, by finding the start position that minimizes the number of baits.\n"
                                      "\"thin-s\":  Thin out a bait file to every Nth bait region, by finding the start position that maximizes the number of sequences.\n"
                                      "\"thin-b-old\":  Similar to thin-b, but treats all loci as if they come from one alignment file. Identical to behaviour of thin-b in version 1.0.5 or earlier.\n"
				      "\"thin-s-old\":  Similar to thin-s, but treats all loci as if they come from one alignment file. Identical to behaviour of thin-b in version 1.0.5 or earlier.\n",
					  false, "","string");
    cmd.add(bait_filter_mode_Arg);    

    // (#3) in option list 
    ValueArg<string> conversion_mode_Arg("c", "convert",
       "Allows the user to produce the final output file which can be uploaded at a bait producing company. "
       "In this mode, BaitFilter reads the input bait file and instead of doing a filtering step, it produces a "
       "custom bait file that can be uploaded at the baits producing company. In order to avoid confusion, "
       "a filtering step cannot be done in the same run as the conversion. If you want to filter a bait file and "
       "convert the output, you will need to call this program more than once, first to do the filtering and second to do "
       "the conversion. Allowed conversion parameters currently are: \"four-column-upload\".\n"
       "New output formats can be added upon request. Please contact the author: Christoph Mayer, Email: Mayer Christoph <c.mayer.zfmk@uni-bonn.de>",
	false, "", "string");
    cmd.add( conversion_mode_Arg );

    // (#2) in option list 
    ValueArg<string> out_filename_Arg("o", "output-bait-file-name",
       "Name of the output bait file. All modes, except the conversion mode, produce files in the BaitFisher format.",
	false, "", "string");
    cmd.add( out_filename_Arg );

    // (#1) in option list 
    ValueArg<string> bait_filename_Arg("i", "input-bait-file-name",
       "Name of the input bait locus file. This is the bait file obtained from the BaitFisher program or from a previous filter run with BaitFilter.",
	true, "", "string");
    cmd.add( bait_filename_Arg );
    



//************************************************
    cmd.parse( argc, argv );
//************************************************

    faststring tmp, tmp1;

    // Assigning parameters to variables:
    global_bait_filename                  = bait_filename_Arg.getValue().c_str();
    tmp                                   = bait_filter_mode_Arg.getValue().c_str();
    global_output_filename                = out_filename_Arg.getValue().c_str();
    global_blast_exe                      = blast_exe_Arg.getValue().c_str();
    global_blast_db                       = blast_db_Arg.getValue().c_str();
    global_blast_min_hit_coverage_of_bait = blast_hit_coverage_Arg.getValue();
    global_blast_max_first_hit_evalue     = blast_primary_hit_evalue_Arg.getValue();
    global_blast_max_second_hit_evalue    = blast_secondary_hit_evalue_Arg.getValue();

    global_blast_evalue_commandline       = blast_evalue_extra_Arg.getValue().c_str();
    global_blast_extra_commandline        = blast_extra_commandline_Arg.getValue().c_str();
    global_thinning_step_width            = thinning_step_width_Arg.getValue();
    tmp1                                  = conversion_mode_Arg.getValue().c_str();
    global_ProbeID_prefix                 = conversion_IDprefix_Arg.getValue().c_str();
    global_verbosity                      = verbosity_Arg.getValue();

    global_stats_mode                     = global_stats_Arg.getValue();
    global_blast_result_file              = blast_result_Arg.getValue().c_str();

    if (!global_blast_result_file.empty())
    {
      skip_blast = true;
      if ( !file_exists(global_blast_result_file.c_str()) )
      {
	logerr << "Error: You specified a blast result file which does not exist or cannot be opened. Exiting.\n";
	good_bye_and_exit(stdout, -2);
      }
    }

    if (global_output_filename.empty() && (!tmp.empty() || !tmp1.empty()) )
    {
      logerr << "Error: The -o parameter, used to specify an output file, is a required parameter for all modes specified with the -m and -c parameter.\n"; 
      good_bye_and_exit(stdout, -2);
    }

    if (tmp.empty() && tmp1.empty() && !global_stats_mode )
    {
      logerr << "Error: When using the BaitFilter program you have to specify either the -m (--mode) option, the -c (--convert) option or the -S option.\n";
      logerr << "These are the three main running modes and one of them has to be chosen.\n";
      logerr << "See the BaitFilter manual for help. ";
      logerr << "Type " << PROGNAME << " -h for a quick online help.\n";
      good_bye_and_exit(stdout, -1);
    }

    if (!tmp.empty() && !tmp1.empty() )
    {
      logerr << "Error: You specified a -m (--mode) option together with the -c (--convert) option.\n";
      logerr << "Filtering runs and conversion runs cannot be combined in one run, but have to be done consecutively, by invoking " << PROGNAME << " twice.\n";
      logerr << "See the BaitFilter manual for help.\n";
      logerr << "Type " << PROGNAME << " -h for a quick online help.\n";
      good_bye_and_exit(stdout, -1);
    }

    if (global_blast_max_first_hit_evalue > global_blast_max_second_hit_evalue)
    {
      logerr << "Error: You specified that the E-value threshold of the first (best) hit is larger than that of the second hit for a bait against a genome.\n"
	     << "This does not make sense. Equal E-value thresholds for the first and the second hit or a smaller E-value for the first hit are allowed.\n";
      good_bye_and_exit(stdout, -1);
    }

    if (!global_blast_result_file.empty() && global_blast_exe != "blastn")
    {
      logerr << "Error: You specified a blast result file and a blast executable. However, specifying a blast result file has the effect that no blast analysis is done\n";
      logerr << "Therefore the two parameters cannot be specified in the same run.\n";
      good_bye_and_exit(stdout, -1);
    }

    if (!global_blast_result_file.empty() && !global_blast_db.empty())
    {
      logerr << "Error: You specified a blast result file and a blast data base file. However, specifying a blast result file has the effect that no blast analysis is done\n";
      logerr << "Therefore the two parameters cannot be specified in the same run.\n";
      good_bye_and_exit(stdout, -1);
    }

    if (!global_blast_result_file.empty() && !global_blast_extra_commandline.empty() )
    {
      logerr << "Error: You specified a blast result file and the blast-extra-commandline option. Since no blast analysis is done when specifying a blast result file, "
	        "no blast-extra-commandline option is allowed if a blast result file is specified.\n";
      good_bye_and_exit(stdout, -1);
    }

    if (!global_blast_result_file.empty() && !global_blast_evalue_commandline.empty() )
    {
      logerr << "Error: You specified a blast result file and a blast-evalue-cutoff. Since no blast analysis is done when specifying a blast result file, "
	        "no blast-evalue-cutoff option is allowed if a blast result file is specified.\n";
      good_bye_and_exit(stdout, -1);
    }

    if (global_stats_mode)
    {
      global_mode = 'S';
    }
    else if (!tmp.empty())
    {
      global_conversion_mode = 0;

      if (tmp == "ab")      // Best per alignment, number of baits
      {
	global_mode = 'a';
      }
      else if (tmp == "as")
      {
	global_mode = 'A';
      }
      else if (tmp == "fb")
      {
	global_mode = 'f';
      }
      else if (tmp == "fs")
      {
	global_mode = 'F';
      }
      else if (tmp == "blast-a")
      {
	global_mode = 'b';
      }
      else if (tmp == "blast-f")
      {
	global_mode = 'B';
      }
      else if (tmp == "blast-l")
      {
	global_mode = 'x';
      }
      else if (tmp == "blast-c")
      {
	global_mode = 'C';
      }
      else if (tmp == "thin-b-old")
      {
	global_mode = 't';
	if (global_thinning_step_width == 0)
	{
	  logerr << "Error: In \"thin-b\" mode, the --thinning-step-width is a required parameter and the specified value must be greater than 0,\n\n";
	  good_bye_and_exit(stdout, -11);
	}
      }
      else if (tmp == "thin-s-old")
      {
	global_mode = 'T';
	if (global_thinning_step_width == 0)
	{
	  logerr << "Error: In \"thin-s-old\" mode, the --thinning-step-width is a required parameter and the specified value must be greater than 0,\n\n";
       	  good_bye_and_exit(stdout, -11);
	}
      }
      // BEGIN: Non permanent test code
      else if (tmp == "thin-b")
      {
	global_mode = 'h';
	if (global_thinning_step_width == 0)
	{
	  logerr << "Error: In \"thin-b\" mode, the --thinning-step-width is a required parameter and the specified value must be greater than 0,\n\n";
	  good_bye_and_exit(stdout, -11);
	}
      }
      else if (tmp == "thin-s")
      {
	global_mode = 'H';
	if (global_thinning_step_width == 0)
	{
	  logerr << "Error: In \"thin-s\" mode, the --thinning-step-width is a required parameter and the specified value must be greater than 0,\n\n";
	  good_bye_and_exit(stdout, -11);
	}
      }
      // END: Non permanent test code
      else
      {
	logerr << "Unknown option for parameter -m <string> : " << tmp << " encountered.\n";
	logerr << "Allowed modes can be viewed with: " << PROGNAME << " -h\n";
	good_bye_and_exit(stdout, -3);
      }
    } // END if (!tmp.empty())
    else // if ( !tmp1.empty() )  // The if is not necessary, since either tmp or tmp1 are ensured to be non-empty.
    {
      global_mode = 'c';

      tmp1.ToLower(); // This argument shall be case insensitive. It is still sufficiently sensitive. 

      if (tmp1 == "agilent-germany" || tmp1 == "four-column-upload")
      {
	global_conversion_mode = 1;
      }
      else 
      {
	logerr << "Unknown option for parameter -c <string> : " << tmp1 << " encountered.\n";
	logerr << "Allowed modes can be viewed with: " << PROGNAME << " -h\n";
	good_bye_and_exit(stdout, -3);
      }
    }
  }
  catch (ArgException &e)
  {
    logerr << "Error: " << e.error() << " for arg " << e.argId() << "\n";
    good_bye_and_exit(stdout, -1);
  }

  // Do some checks here:

  // Check whether the blast executable can be found and can be executed.
  // Check whether the data base file exists.

  if (global_thinning_step_width > 0 && (tolower(global_mode) != 't' &&  tolower(global_mode) != 'h' ) ) 
  {
    logerr << "Error: The --thinning-step-width option only has an effect in the \"-m thin-b\", \"-m thin-s\", \"-m thin-b-old\" or \"-m thin-s-old\" modes. Specifying this value in other modes is not allowed.\n";
    good_bye_and_exit(stdout, -12);
  }

  if (global_mode == 'C' && global_blast_min_hit_coverage_of_bait == 0)
  {
    logerr << "Error: You specified the mode -m blast-c but did not specify the blast-min-hit-coverage-of-baits-in-tiling-stack option which is required in this mode. Exiting.\n";
    good_bye_and_exit(stdout, -14);
  }

  // In case we have to conduct a blast search:
  if (global_mode == 'b' || global_mode == 'B' || global_mode == 'x' ||  global_mode == 'C' )
  {
    //    #ifndef SKIPBLAST
    if (!skip_blast)
    {
      bool good = blast_exe_exists();

      if (!good)
      {
	logerr << "Error: You try to run a blast filter, but the blast executable does not seem to work.\n";
	logerr << "       I was not able to run the blast program with the following command line option:\n";
	logerr << "       " << global_blast_exe << " -version\n";
	logerr << "       This command should not give an error message, but display the version number of blast program.\n";
	logerr << "       Exiting.\n";
	good_bye_and_exit(stdout, -1);
      }

      if (global_blast_db.empty())
      {
	logerr << "Error: You try to run a blast filter, but the index files of your reference genome is not specified.\n";
	logerr << "       Use the --ref-blast_db option to specify a blast data base index file.\n";
	good_bye_and_exit(stdout, -1);
      }

      good = blast_db_files_exist(global_blast_db);

      if (!good)
      {
	logerr << "Error: You try to run a blast filter, but the index files of your reference genome does not exist.\n";
	logerr << "       If you only have a fasta file at the moment, you have to use the makeblastdb program which comes with your blast installation.\n";
	logerr << "       The creation of the blast index was successful if you where able to created file with the extension: nin, nsq and nhr.\n";
	logerr << "       Exiting.\n";
	good_bye_and_exit(stdout, -1);
      }
    }
  }
  // #endif

  if (!(global_mode == 'b' || global_mode == 'B' || global_mode == 'x' ||  global_mode == 'C' )) // No blast mode
  {
    if (!global_blast_db.empty())
    {
      logerr << "Warning: You specified a blast database file, but you run " << PROGNAME
	     << " in a non blast mode. Thus, the specified data base file will not be used.\n";
    }

    if (!global_blast_result_file.empty())
    {
      logerr << "Error: You specified a blast result file together with a filtering mode that does require a blast search. Exiting.\n";
      good_bye_and_exit(stdout, -1);
    }
  }


  if ( !global_ProbeID_prefix.empty() && global_mode != 'c')
  {
    logerr << "Warning: You specified a ProbeID prefix string, which has an effect only in the conversion mode.\n"
	   << "         In non conversion mode this prefix will be ignored.\n";
  }
}



void print_parameters(FILE *of, const char *s)
{
  fprintf(of,   "%sBaitFilter has been started with the following command line parameters:\n",                      s);
  fprintf(of,   "%sInput bait file name:                             %s\n",   s, global_bait_filename.c_str());
  fprintf(of,   "%sOutput bait file name:                            %s\n",   s, global_output_filename.c_str());
  if (global_mode == 'a')
  {
    fprintf(of, "%sFilter mode:                                      Retain only the best locus per alignment file; criterion: minimize number of baits.\n", s);
  }
  else if (global_mode == 'A')
  {
    fprintf(of, "%sFilter mode:                                      Retain only the best locus per alignment file; criterion: maximize sequence coverage in bait region.\n", s);
  }
  else if (global_mode == 'f')
  {
    fprintf(of, "%sFilter mode:                                      Retain only the best locus per feature; criterion: minimize number of baits in bait region.\n", s);
  }
  else if (global_mode == 'F')
  {
    fprintf(of, "%sFilter mode:                                      Retain only the best locus per feature; criterion: maximize sequence coverage in bait region.x\n", s);
  }
  else if (global_mode == 'b')
  {
    fprintf(of, "%sFilter mode:                                      Blast filter: Remove all bait regions from output that belong to the ALIGNMENT in which one bait has at least two good hits to the reference genome.\n", s);
  }
  else if (global_mode == 'B')
  {
    fprintf(of, "%sFilter mode:                                      Blast filter: Remove all bait regions from output that belong to the FEATUE in which one bait has at least two good hits to the reference genome.\n", s);
  }
  else if (global_mode == 'x')
  {
    fprintf(of, "%sFilter mode:                                      Blast filter: Remove all bait REGIONs in which one bait has at least two good hits to the reference genome.\n", s);
  }
  else if (global_mode == 'C')
  {
    fprintf(of, "%sFilter mode:                                      Blast filter: Conduct coverage filter without any other blast filter.\n", s);
  }
  else if (global_mode == 't')
  {
    fprintf(of, "%sFilter mode:                                      Thinning mode: Retain only every %u th bait region. Choose the starting offset such that the number of baits is minimized.\n", s, global_thinning_step_width);    
  }
  else if (global_mode == 'T')
  {
    fprintf(of, "%sFilter mode:                                      Thinning mode: Retain only every %u th bait region. Choose the starting offset such that the number of sequences the result is based on is maximized.\n", s, global_thinning_step_width);    
  }
  else if (global_mode == 'h')
  {
    fprintf(of, "%sFilter mode:                                      Thinning mode new algorithm: Retain only every %u th bait region. Choose the starting offset such that the number of baits is minimized.\n", s, global_thinning_step_width);    
  }
  else if (global_mode == 'H')
  {
    fprintf(of, "%sFilter mode:                                      Thinning mode new algorithm: Retain only every %u th bait region. Choose the starting offset such that the number of sequences the result is based on is maximized.\n", s, global_thinning_step_width);    
  }
  else if (global_mode == 'c')
  {
    if (global_conversion_mode == 1)
    {
      fprintf(of, "%sConversion mode:                                Convert bait file to four-column-upload format.\n", s);
    }
    else
    {
      fprintf(of, "%sERROR: Invalid internal conversion mode: %c.\nPlease report this bug.\n", s, global_mode);
      good_bye_and_exit(stdout, -99);
    }
  }
  else if (global_mode == 'S')
  {
    fprintf(of, "%sMode:                                             Prints stats mode.\n", s); 
  }
  else
  {
    fprintf(of, "%sERROR: Invalid internal mode:                     %c.\n", s, global_mode);
    good_bye_and_exit(stdout, -99);
  }

  if (global_mode == 'c' && !global_ProbeID_prefix.empty())
  {
    fprintf(of, "%sProbeID suffix string in conversion mode:         %s\n", s, global_ProbeID_prefix.c_str());
  }

  fprintf(of, "%sVerbosity:                                        %d\n", s, global_verbosity);


  if (global_mode == 'b' || global_mode == 'B' || global_mode == 'x' ||  global_mode == 'C')
  {
    faststring str;

    fprintf(of, "%sBlast search settings and parameters:\n",   s);

    if (!skip_blast)
      fprintf(of, "%sBlast executable:                                 %s\n",  s, global_blast_exe.c_str());

    fprintf(of, "%sBlast data base file name:                        %s\n",  s, global_blast_db.c_str());

    fprintf(of, "%sBlast minimum hit coverage of bait in first hit:  %lf\n", s, global_blast_min_hit_coverage_of_bait);

    fprintf(of, "%sBlast maximum E-value of first/best hit:          %lg\n", s, global_blast_max_first_hit_evalue);
    fprintf(of, "%sBlast maximum E-value of second best hit:         %lg\n", s, global_blast_max_second_hit_evalue);

    if (!global_blast_extra_commandline.empty())
    {
      fprintf(of, "%sBlast extra command line argument:                %s\n", s, global_blast_extra_commandline.c_str());
    }
    else
    {
      fprintf(of, "%sBlast extra command line argument:                <not specified>\n", s);
    }

    if (!global_blast_evalue_commandline.empty())
    {
      fprintf(of, "%sBlast E-value cut-off:                            %s\n", s, global_blast_evalue_commandline.c_str());
    }
    else
    {
      fprintf(of, "%sBlast E-value cut-off:                            <not specified>\n", s);
    }

    if (!global_blast_result_file.empty())
    {
      fprintf(of, "%sDo not conduct a new blast search, but use the    \n"
	          "%sfollowing existing blast result file:             %s\n",  s,s, global_blast_result_file.c_str() );
    }


  } // END if (global_mode == 'b' || global_mode == 'B' || global_mode == 'x' ||  global_mode == 'C')
}


void print_parameters(ostream &os, const char *s)
{
  os << s <<   "BaitFilter has been started with the following command line parameters:\n";
  os << s <<   "Input bait file name:                             " <<  global_bait_filename.c_str() << '\n';
  os << s <<   "Output bait file name:                            " <<  global_output_filename.c_str() << '\n';
  if (global_mode == 'a')
  {
    os << s << "Filter mode:                                      Retain only the best locus per alignment file; criterion: minimize number of baits.\n";
  }
  else if (global_mode == 'A')
  {
    os << s << "Filter mode:                                      Retain only the best locus per alignment file; criterion: maximize sequence coverage in bait region.\n";
  }
  else if (global_mode == 'f')
  {
    os << s << "Filter mode:                                      Retain only the best locus per feature; criterion: minimize number of baits in bait region.\n";
  }
  else if (global_mode == 'F')
  {
    os << s << "Filter mode:                                      Retain only the best locus per feature; criterion: maximize sequence coverage in bait region.x\n";
  }
  else if (global_mode == 'b')
  {
    os << s << "Filter mode:                                      Blast filter: Remove all bait regions from output that belong to the ALIGNMENT in which one bait has at least two good hits to the reference genome.\n";
  }
  else if (global_mode == 'B')
  {
    os << s << "Filter mode:                                      Blast filter: Remove all bait regions from output that belong to the FEATUE in which one bait has at least two good hits to the reference genome.\n";
  }
  else if (global_mode == 'x')
  {
    os << s << "Filter mode:                                      Blast filter: Remove all bait REGIONs in which one bait has at least two good hits to the reference genome.\n";
  }
  else if (global_mode == 'C')
  {
    os << s << "Filter mode:                                      Blast filter: Conduct coverage filter without any other blast filter.\n";
  }
  else if (global_mode == 't')
  {
    os << s << "Filter mode:                                      Thinning mode: Retain only every " << global_thinning_step_width << " th bait region. Choose the starting offset such that the number of baits is minimized.\n";    
  }
  else if (global_mode == 'T')
  {
    os << s << "Filter mode:                                      Thinning mode: Retain only every " << global_thinning_step_width << " th bait region. Choose the starting offset such that the number of sequences the result is based on is maximized.\n";
  }
  else if (global_mode == 'h')
  {
    os << s << "Filter mode:                                      Thinning mode new algorithm: Retain only every " << global_thinning_step_width  << " th bait region. Choose the starting offset such that the number of baits is minimized.\n";
  }
  else if (global_mode == 'H')
  {
    os << s << "Filter mode:                                      Thinning mode new algorithm: Retain only every " << global_thinning_step_width  << " th bait region. Choose the starting offset such that the number of sequences the result is based on is maximized.\n";
  }
  else if (global_mode == 'c')
  {
    if (global_conversion_mode == 1)
    {
      os << s << "Conversion mode:                                Convert bait file to four-column-upload format.\n";
    }
    else
    {
      os << s << "ERROR: Invalid internal conversion mode: " << global_conversion_mode << " .\nPlease report this bug.\n";
      good_bye_and_exit(cout, -99);
    }
  }
  else if (global_mode == 'S')
  {
    os << s << "Mode:                                             Prints stats mode.\n";
  }
  else
  {
    os << s << "ERROR: Invalid internal mode:                     " << global_mode << '\n';
    good_bye_and_exit(cout, -99);
  }

  if (global_mode == 'c' && !global_ProbeID_prefix.empty())
  {
    os << s << "ProbeID suffix string in conversion mode:         " << global_ProbeID_prefix.c_str() << '\n';
  }

  os << s << "Verbosity:                                        " << global_verbosity << '\n';


  if (global_mode == 'b' || global_mode == 'B' || global_mode == 'x' ||  global_mode == 'C')
  {
    faststring str;

    os << s << "Blast search settings and parameters:\n";
    os << s << "Blast executable:                                 " << global_blast_exe.c_str() << '\n';
    os << s << "Blast data base file name:                        " << global_blast_db.c_str() << '\n';

    os << s << "Blast minimum hit coverage of bait in first hit:  " << global_blast_min_hit_coverage_of_bait << '\n';

    os << s << "Blast maximum E-value of first/best hit:           " << global_blast_max_first_hit_evalue << '\n';
    os << s << "Blast maximum E-value of second best hit:          " << global_blast_max_second_hit_evalue << '\n';

    if (!global_blast_extra_commandline.empty())
    {
      os << s << "Blast extra command line argument:               " << global_blast_extra_commandline.c_str() << '\n';
    }
    else
    {
      os << s << "Blast extra command line argument:               <not specified>\n";
    }

    if (!global_blast_evalue_commandline.empty())
    {
      os << s << "Blast E-value cut-off:                           " << global_blast_evalue_commandline.c_str() << '\n';
    }
    else
    {
      os << s << "Blast E-value cut-off:                           <not specified>\n";
    }

    if (!global_blast_result_file.empty())
    {
      os << s << "Do not conduct a new blast search, but use the    " << '\n'
	 << s << "following existing blast result file:             " << global_blast_result_file.c_str()  << '\n';
    }
  } // END if (global_mode == 'b' || global_mode == 'B' || global_mode == 'x' ||  global_mode == 'C')
}



void print_calling_command_line(FILE *of, unsigned argc, char ** argv)
{
  faststring cmd;
  
  for (unsigned i=0; i<argc; ++i)
  {
    cmd += argv[i];
    cmd += " ";
  }  

  fprintf(of, "\nThe BaitFilter program was called with the following command line:\n");
  fputs(cmd.c_str(), of);
  fputc('\n',of);  
}


void print_calling_command_line(ostream &os, unsigned argc, char ** argv)
{
  faststring cmd;
  
  for (unsigned i=0; i<argc; ++i)
  {
    cmd += argv[i];
    cmd += " ";
  }  

  os << "\nThe BaitFilter program was called with the following command line:\n";
  os << cmd.c_str() << '\n';  
}



void print_output_creation_message(FILE *of, char *s)
{
  fprintf(of, "%sOutput created by " PROGNAME ", version " VERSION "\n", s);
  fprintf(of, "%s\n", s);
  print_parameters(of, s);
  fprintf(of, "%s\n", s);
}

void print_output_creation_message(ostream &os, char *s)
{
  os << s << "Output created by " PROGNAME ", version " VERSION "\n";
  os << s << '\n';
  print_parameters(os, s);
  os << s << '\n';
}
