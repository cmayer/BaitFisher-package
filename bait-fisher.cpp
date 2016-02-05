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
#define  NODEBUG
#include "faststring2.h"
#include <iostream>
#include "GFF-class.h"
#include "GFF-collection.h"
#include "mydir-unix.h"
#include <list>
#include "CFile/CFile2_1.h"
#include "CSequences2.h"
#include "CSequence_Mol2_1.h"
#include "CDnaString2.h"
#include <climits>
#include <ctime>

#include "CAligner/CAligner.h"
#include "scoring-matrices/CScoreMatrix.h"
#include "basic-DNA-RNA-AA-routines.h"

#include "CDistance_matrix.h"
#include "CSeqNameList.h"
#include "CTaxonNamesDictionary.h"
#include "CRequiredTaxon.h"

//#include "pipeline-step3.h"
#include "Csequence_cluster_and_center_sequence.h"
#include "DEBUG_STUFF.h"


#define PROGNAME  "BaitFisher"
#define VERSION   "1.2.7"


using namespace std;

#include "bait-fisher-helper.h"

// Verbosity rules:
//      verbosity 0 -> only error messages that lead to exit() (goes to cerr)
//      verbosity 1 -> warnings, ALERT (goes to cerr)
//      verboisty 2 -> more warinings, e.g. skipping features
//      verbosity 3 -> progress (goes to cout)
//      verbosity 4 -> more progress (goes to cout)
//      verbosity 6 -> per alignment/feautre success NOTES
//      verbosity 8 -> per coordinate success NOTES
//      verbosity >=20: DEGUB
// Default: 5

unsigned global_VERBOSITY = 0; // Will be initialized in parse_parameter_file. This value will never be used.




// Goal of this program:
//


int DEBUG_COUNT_NEW_CALLS_FEATURE_ALIGNMENT = 0;
int DEBUG_COUNT_NEW_CALLS_CBAIT_REGION_INFORMATION = 0;
int DEBUG_COUNT_NEW_CALLS_CSEQUENCE_LOCI_CLUSTER_COLLECTION = 0;
int DEBUG_COUNT_NEW_CALLS_MSA = 0;
int DEBUG_COUNT_NEW_ALIGNMENT_FROM_FILE = 0;

int DEBUG_COUNT_DELETE_CALLS_FEATURE_ALIGNMENT = 0;
int DEBUG_COUNT_DELETE_CALLS_CBAIT_REGION_INFORMATION = 0;
int DEBUG_COUNT_DELETE_CALLS_CSEQUENCE_LOCI_CLUSTER_COLLECTION = 0;
int DEBUG_COUNT_DELETE_CALLS_MSA = 0;
int DEBUG_COUNT_DELETE_ALIGNMENT_FROM_FILE = 0;

void print_DEBUG_COUNT (ostream &os, const char * comment)
{
  os << "DEBUG: NEW/DELETE COUNTs:              " << comment << endl;
  os << "DEBUG: NUMBER OF ALIGNMENTS CREATED/DELETED:  " << DEBUG_COUNT_NEW_ALIGNMENT_FROM_FILE                     << "/" <<  DEBUG_COUNT_DELETE_ALIGNMENT_FROM_FILE << endl;
  os << "DEBUG: NUMBER OF FEATUREALIGNMENT C/D:        " << DEBUG_COUNT_NEW_CALLS_FEATURE_ALIGNMENT                 << "/" <<  DEBUG_COUNT_DELETE_CALLS_FEATURE_ALIGNMENT << endl;
  os << "DEBUG: MSA  CREATED/DELETED                   " << DEBUG_COUNT_NEW_CALLS_MSA                               << "/" <<  DEBUG_COUNT_DELETE_CALLS_MSA << endl;
  os << "DEBUG: CSEQUENCE_LOCI_CLUSTER_COLLECTION C/D: " << DEBUG_COUNT_NEW_CALLS_CSEQUENCE_LOCI_CLUSTER_COLLECTION << "/" <<  DEBUG_COUNT_DELETE_CALLS_CSEQUENCE_LOCI_CLUSTER_COLLECTION << endl;
  os << "DEBUG: CBAIT_REGION_INFORMATION C/D:          " << DEBUG_COUNT_NEW_CALLS_CBAIT_REGION_INFORMATION          << "/" <<  DEBUG_COUNT_DELETE_CALLS_CBAIT_REGION_INFORMATION << endl;
}

void get_time_string(faststring &str)
{
  time_t      rawtime;
  struct tm   *timeinfo;
  static char buffer[80];

  std::time (&rawtime);
  timeinfo = std::localtime(&rawtime);
  std::strftime (buffer,80,"%F %T",timeinfo);
  str = faststring(buffer);
}

template<class T>
int delete_pointers_in_vector(vector<T> &v)
{
  int i, N, count=0;
  N = v.size();

  for (i=0; i<N; ++i)
  {
    if (v[i])
    {
      ++count;
      delete v[i];
    }
  }
  v.clear();
  return count;
}



// Add more checks!! Currently we do not check whether parameters occur twice, etc.
void parse_parameter_file(faststring      & fname_parameter_file,
			  faststring      & fname_gff,
			  faststring      & fname_genome,
			  faststring      & string_required_taxa,
			  faststring      & dname_alignments,
			  faststring      & genome_name,
			  unsigned short  & bait_length,
			  double          & cluster_threshold,
			  faststring      & output_directory,
			  unsigned short  & bait_offset,
			  unsigned short  & bait_number,
			  bool            & remove_reference_species_from_alignments,
			  bool            & alignment_cutting,
			  faststring      & gff_feature_name,
			  faststring      & gff_record_attribute_name_for_geneID,
			  char            & sequence_name_field_delimiter,
			  unsigned short  & sequence_name_taxon_field_number,
			  unsigned short  & sequence_name_geneID_field_number,
			  bool            & sequence_name_tokenizer_overwritten,
			  char            & center_computation_mode,
			  unsigned        & verbosity
			 )
{
  fname_gff.clear();
  fname_genome.clear();
  string_required_taxa.clear();
  dname_alignments.clear();
  genome_name.clear();
  bait_length                              = 0;
  cluster_threshold                        = -1;
  output_directory.clear();
  bait_offset                              = USHRT_MAX;
  bait_number                              = USHRT_MAX;
  remove_reference_species_from_alignments = false;
  alignment_cutting                        = false;

  sequence_name_field_delimiter            = '|';
  sequence_name_taxon_field_number         = 2;
  sequence_name_geneID_field_number        = 3;
  sequence_name_tokenizer_overwritten      = false;

  center_computation_mode                = 1; // 0:exhaustive, 1:heuristic, 3: special comparison
  verbosity                              = 5;

  faststring line, pre, post;

  ifstream is;
  is.open(fname_parameter_file.c_str());
  if (is.fail() )
  {
    cerr << "Error: Parameter file \""  << fname_parameter_file << "\" could not be opened. I guess it does not exist." << endl;
    exit(-1);
  }

  getline(is, line);

  while (!is.fail() )
  {
    line.shorten_to_first_occurrence_of('#');
    line.removeSpacesBack();
    line.removeSpacesFront();

    if ( line.empty() )
    {
      getline(is, line);
      continue;
    }

    if(line.divide_at_first_of(" \t",pre, post))
    {
      pre.removeSpacesBack();
      post.removeSpacesFront();
    }
    else
    {
      pre = line;
      post.clear();
    }

    if (pre == "gff-file-reference-genome")
    {
      fname_gff = post;
    }
    else if (pre == "reference-genome-file")
    {
      fname_genome = post;
    }
    else if (pre == "required-taxonomic-groups-string")
    {
      string_required_taxa = post;
    }
    else if (pre == "directory-transcript-files")
    {
      dname_alignments = post;
    }
    else if (pre == "reference-genome-name")
    {
      genome_name = post;
    }
    else if (pre == "bait-length")
    {
      bait_length = post.ToUnsigned();
    }
    else if (pre == "bait-offset")
    {
      bait_offset = post.ToUnsigned();
    }
    else if (pre == "bait-number")
    {
      bait_number = post.ToUnsigned();
    }
    else if (pre == "cluster-threshold")
    {
      cluster_threshold = post.ToDouble();
    }
    else if (pre == "output-directory")
    {
      output_directory = post;
    }
    else if (pre == "remove-reference-sequence")
    {
      post.ToLower();
      if (post == "yes")
	remove_reference_species_from_alignments = true;
      else if (post == "no")
	remove_reference_species_from_alignments = false;
      else
      {
	cerr << "Error: Unknown option " << post << " for parameter " << pre << endl;
	cerr << "The " << pre << " parameter will be ignored." << endl;
	exit(-7);
      }
    }
    else if (pre == "alignment-cutting")
    {
      post.ToLower();
      if (post == "yes")
	alignment_cutting = true;
      else if (post == "no")
	alignment_cutting = false;
      else
      {
	cerr << "Error: Unknown option " << post << " for parameter " << pre << endl;
	cerr << "The " << pre << " parameter will be ignored." << endl;
	exit(-7);
      }
    }
    else if (pre == "gff-feature-name")
    {
      gff_feature_name = post;
    }
    else if (pre == "gff-record-attribute-name-for-geneID")
    {
      gff_record_attribute_name_for_geneID = post;
    }
    else if (pre == "sequence-name-field-delimiter")
    {
      if (post.size()!=1)
      {
	cerr << "ERROR when parsing the parameter file: sequence-name-field-delimiter must be a single character, but I found "
	     << sequence_name_field_delimiter << "." << endl;
	exit(-1);
      }
      sequence_name_field_delimiter = post[0];
    }
    else if (pre == "sequence-name-taxon-field_number")
    {
      sequence_name_taxon_field_number = post.ToUnsigned();
    }
    else if (pre == "sequence-name-geneID-field_number")
    {
      sequence_name_geneID_field_number = post.ToUnsigned();
    }
    else if (pre == "verbosity")
    {
      verbosity = post.ToUnsigned();
    }
    else if (pre == "Hamming-center-sequence-computation-mode" || pre == "centroid-computation-mode")
    {
      if (post == "exhaustive")
      {
	center_computation_mode = 0;
      }
      else if (post == "heuristic")
      {
	center_computation_mode = 1;
      }
      else if (post == "comparison")
      {
	center_computation_mode = 3;
      }
      else
      {
	cerr << "Unknown option for parameter \"center-computation-mode\" in parameter file. Allowed options are: \"exhaustive\" and \"heuristic\"." << endl;
	exit(-1);
      }
    }
    else
    {
      cerr << "Error: Unknown parameter in parameter file: \"" << pre << "\"."<< endl;
      exit(-1);
    }
    getline(is, line);
  } // END while (!is.fail() )

  // Basic consitency checking of parameters.
  // Are all required parameters supplied?

  if (alignment_cutting) // Then we require some parameters:
  {
    if (fname_gff.empty())
    {
      cerr << "Error: The parameter \"gff-file\" is required if alignment_cutting is set to yes. This parameter must be specified in the parameter file." << endl;
      exit(-4);
    }

    if (fname_genome.empty() )
    {
      cerr << "Error: The parameter \"genome-file\" is required if alignment_cutting is set to yes. This parameter must be specified in the parameter file." << endl;
      exit(-4);
    }

    if (genome_name.empty() )
    {
      cerr << "Error: The parameter \"genome-name\" is required if alignment_cutting is set to yes. This parameter must be specified in the parameter file." << endl;
      exit(-4);
    }
    // Parameters with default values:
    if (gff_feature_name.empty() )
      gff_feature_name                         = "CDS";
    if (gff_record_attribute_name_for_geneID.empty() )
      gff_record_attribute_name_for_geneID     = "Parent";
  }
  else // No alignment cutting
  {
    if (!fname_gff.empty())
    {
      cerr << "Error: The parameter \"gff-file\" is specified even though alignment_cutting is set to no. Please remove this parameter from the parameter file." << endl;
      exit(-4);
    }

    if (!fname_genome.empty() )
    {
      cerr << "Error: The parameter \"genome-file\" is specified even though alignment_cutting is set to no. Please remove this parameter from the parameter file." << endl;
      exit(-4);
    }

    if (!genome_name.empty() )
    {
      cerr << "Error: The parameter \"genome-name\" is specified even though alignment_cutting is set to no. Please remove this parameter from the parameter file." << endl;
      exit(-4);
    }
    if (!gff_feature_name.empty() )
    {
      cerr << "Error: The parameter \"gff-feature-name\" is specified even though alignment_cutting is set to no. Please remove this parameter from the parameter file." << endl;
      exit(-4);
    }
    if (!gff_record_attribute_name_for_geneID.empty() )
    {
      cerr << "Error: The parameter \"gff-record-attribute-name-for-geneID\" is specified even though alignment_cutting is set to no. Please remove this parameter from the parameter file." << endl;
      exit(-4);
    }
  }

//   if (string_required_taxa.empty() )
//   {
//     cerr << "Error: The required parameter \"required-taxonomic-groups-string\" is not specified in the parameter file." << endl;
//     exit(-4);
//   }
  if (dname_alignments.empty() )
  {
    cerr << "Error: The required parameter \"directory-transcript-files\" is not specified in the parameter file." << endl;
    exit(-4);
  }
  if (bait_length == 0)
  {
    cerr << "Error: Required parameter \"bait-length\" is not specified in the parameter file." << endl;
    exit(-4);
  }
  if (cluster_threshold == -1)
  {
    cerr << "Error: Required parameter \"cluster-threshold\" is not specified in the parameter file." << endl;
    exit(-4);
  }
  if (  output_directory.empty() )
  {
    cerr << "Error: Required parameter \"output-directory\" is not specified in the parameter file." << endl;
    exit(-4);
  }
  if (bait_number == SHRT_MAX)
  {
    cerr << "Error: Required parameter \"bait-number\" is not specified in the parameter file." << endl;
    exit(-4);
  }
  if (bait_number == 0)
  {
    cerr << "Error: Parameter \"bait-number\" has been set to 0 in the parameter file, which is not a valid value." << endl;
    exit(-4);
  }
  if (bait_offset == USHRT_MAX)
  {
    cerr << "Error: Required parameter \"bait-offset\" is not specified in the parameter file." << endl;
    exit(-4);
  }
  if (bait_offset == 0 && bait_number > 1)
  {
    cerr << "Error: The parameter \"bait-offset\" has been set to 0 in the parameter file, while the bait number is larger than 1. This does not make sense." << endl;
    cerr << "Please specify a valid tiling design with bait number equal to 1 and bait offset equal to 0, or bait number larger than 1 and offset larger than 0." << endl;
    exit(-4);
  }

  if (!alignment_cutting && string_required_taxa.empty() ) // In this case, sequence names do not need to contain the taxon name.
  {
    sequence_name_field_delimiter    = '\0';
    sequence_name_taxon_field_number = 1;
    sequence_name_geneID_field_number = 1;
    sequence_name_tokenizer_overwritten = true;
  }

  // The field numbers are stored 0-based. In the user interaction, numbers are 1-based.
  --sequence_name_geneID_field_number;
  --sequence_name_taxon_field_number;
}


void print_parameters(ostream         & os,
		      faststring      & fname_gff,
		      faststring      & fname_genome,
		      faststring      & string_required_taxa,
		      faststring      & dname_alignments,
		      faststring      & genome_name,
		      unsigned short    bait_length,
		      double            cluster_threshold,
		      faststring      & output_directory,
		      unsigned short    bait_offset,
		      unsigned short    bait_number,
		      bool              remove_reference_species_from_alignments,
		      bool              alignment_cutting,
		      faststring      & gff_feature_name,
		      faststring      & gff_record_attribute_name_for_geneID,
		      char            sequence_name_field_delimiter,
		      unsigned short  sequence_name_taxon_field_number,
		      unsigned short  sequence_name_geneID_field_number,
		      bool            sequence_name_tokenizer_overwritten,
		      char            center_computation_mode,
		      unsigned        verbosity
		     )
{
  os << "Welcome to the " << PROGNAME << " program, version " << VERSION << endl << endl
     << "This software and any of its components are copyright protected by the authors" << endl
     << "Christoph Mayer and Oliver Niehuis, Zoologisches Forschungsmuseum Alexander" << endl
     << "Koenig, Bonn, Germany. Copyright 2014." << endl
     << "This software is distributed as working material within a workshop given by" << endl
     << "Manuela Sann and Oliver Niehuis at the Australian National University, Canberra." << endl
     << "As long as it is not published, the software may only be used for the purposes of" << endl
     << "the workshop or with explicit permission of one of the authors." << endl;
  os << "Redistribution as source code in part or in whole or in binary form is not permitted." << endl;
  os << "The software comes without warranty of any kind." << endl << endl;

  os << "PROGRAM has been started with the following parameters:"                << endl;
  os << "Alignment cutting                             " << (alignment_cutting ? "yes" : "no");
  os << endl;

  if (alignment_cutting)
  {
    os << "Name of gff file:                             " << fname_gff             << endl;
    os << "Name of genome file:                          " << fname_genome          << endl;
    os << "Name of genome:                               " << genome_name           << endl;
    os << "GFF feature to be excised from sequences:     " << gff_feature_name      << endl;
    os << "GFF record attribute name for geneID:         " << gff_record_attribute_name_for_geneID << endl;
  }

  if (string_required_taxa.empty())
    os << "String with required taxonomic groups:        " << "no taxa specified" << endl;
  else
    os << "String with required taxonomic groups:        " << string_required_taxa  << endl;

  os << "Name of directory with transcriptome files:   " << dname_alignments      << endl;
  os << "Cluster threshold:                            " << cluster_threshold     << endl;
  os << "Bait length:                                  " << bait_length           << endl;
  os << "Bait offset (for multiple baits in region):   " << bait_offset           << endl;
  os << "Bait number (for multiple baits in region):   " << bait_number           << endl;
  os << "Output directory:                             " << output_directory      << endl;
  os << "Remove refernce species from alignment:       " << (remove_reference_species_from_alignments ? "yes": "no") << endl;
  if (sequence_name_tokenizer_overwritten)
  {
    os << "sequence_name_field_delimiter (default \"|\") " << "Not used since no required taxa specified and no aligment cutting." << endl;
    // The field numbers are stored 0-based. In the user interaction, numbers are 1-based.
    os << "Number of taxon field in sequence name:       " << "Not used since no required taxa specified and no aligment cutting." << endl;
    os << "Number of geneID field in sequence name:      " << "Not used since no required taxa specified and no aligment cutting." << endl;
  }
  else
  {
    os << "sequence_name_field_delimiter (default \"|\") " << sequence_name_field_delimiter << endl;
    // The field numbers are stored 0-based. In the user interaction, numbers are 1-based.
    os << "Number of taxon field in sequence name:       " << sequence_name_taxon_field_number+1  << endl;
    os << "Number of geneID field in sequence name:      " << sequence_name_geneID_field_number+1 << endl;
  }
  if (center_computation_mode == 0)
    os << "Center computation mode:                      " << "exhaustive" << endl;
  else if (center_computation_mode == 1)
    os << "Center computation mode:                      " << "heuristic" << endl;
 else
    os << "Center computation mode:                      " << "comparison mode" << endl;
  os <<   "Verbosity:                                    " << verbosity << endl;
  os << endl;
}





// Computes alignment and determines some characteristics:
//  - alignment_score, alignment_length,
//  - s1ins2s: starting pos of s1 in s2. If negative it is outside of s2 by the negative number (0 based)
//    s1ins2e: position after last symbol of s1 in s2. If negative it is outside of s2 by the negative number (0 based)
//    s2ins1s: similar to s1ins2s
//    s2ins1e: similar to s2ins2e
inline void respective_sequence_positions(CScoreMatrix &DNAmat, faststring &s1, faststring &s2,
					  int &s1ins2s,
					  int &s1ins2e,
					  int &s2ins1s,
					  int &s2ins1e,
					  int &alignment_score,
					  int &alignment_length,
					  int verbosity)
{
  //  CScoreMatrix DNAmat("simpledna.mat"); // It is not efficient to open this file for many million short alignments.
  short mode = 1; // local - free end gaps

  std::ofstream log;
  log.open("CAligner.log");

  CAligner CA(s1.c_str(), s1.length(), s2.c_str(), s2.length(), DNAmat, mode, log);

  faststring al;
  align_coord a, e;
  faststring l1, l2, l3;

  CA.align();
  alignment_score = CA.traceback(al, a, e);

  CA.alignment_strings(l1, l2, l3, al);

  if (verbosity >= 50)
  {
    cout << "Alignment score: " << alignment_score << endl;
    cout << "Alignment:" << endl;
    cout << l1 << endl;
    cout << l3 << endl;
    cout << l2 << endl;
  }

  CA.get_respective_start_end_coordinates(s2ins1s,
					  s2ins1e,
					  s1ins2s,
					  s1ins2e,
					  a,e,
					  s1.length(),
					  s2.length(),
					  al);

  if (verbosity > 50)
  {
    std::cout << al << std::endl;
    std::cout << "Score: " << alignment_score << std::endl;
    std::cout << "Internal representation of start: "; a.print(std::cout); std::cout << std::endl;
    std::cout << "Internal representation of end:   "; e.print(std::cout); std::cout << std::endl;

    //    pretty_print_alignment_base0(stdout, l1, l2, l3, 100);
    //    pretty_print_alignment(stdout, l1, l2, l3, 100);

    //    std::cout << l1 << std::endl;
    //    std::cout << l3 << std::endl;
    //    std::cout << l2 << std::endl;

    //   std::cout << "Start: "; a.print(std::cout); std::cout << std::endl;
    //   std::cout << "End:   "; e.print(std::cout); std::cout << std::endl;

    std::cout << "Respective start/end (internal representation): "
	 << "s2ins1s " << s2ins1s << " "
	 << "s2ins1e " << s2ins1e << " "
	 << "s1ins2s " << s1ins2s << " "
	 << "s1ins2e " << s1ins2e << std::endl;

    std::cout << "len s1 " << s1.length() << std::endl;
    std::cout << "len s2 " << s2.length() << std::endl;
    std::cout << "len al " << al.length() << std::endl;
  }


  alignment_length = al.length();
}


void find_coords(CScoreMatrix &DNAmat, CDnaString &s1, CDnaString &s2,
		 int &s1ins2s_best,
		 int &s1ins2e_best,
		 int &s2ins1s_best,
		 int &s2ins1e_best,
		 int &alignment_score_best,
		 int &alignment_length_best,
		 int verbosity)
{
  // Coordinates in alignment, in px_y: x is first and second part, y is linker index.

  char  better;

  int s1ins2s_tmp1;
  int s1ins2e_tmp1;
  int s2ins1s_tmp1;
  int s2ins1e_tmp1;
  int alignment_score_tmp1;
  int alignment_length_tmp1;

  int s1ins2s_tmp2;
  int s1ins2e_tmp2;
  int s2ins1s_tmp2;
  int s2ins1e_tmp2;
  int alignment_score_tmp2;
  int alignment_length_tmp2;

  // Try both: feature sequence and its reverse complement
  respective_sequence_positions(DNAmat, s1, s2,
				s1ins2s_tmp1, s1ins2e_tmp1, s2ins1s_tmp1, s2ins1e_tmp1, alignment_score_tmp1, alignment_length_tmp1,
				verbosity);

  CDnaString revcomp_s2;

  revcomp_s2.setToReverseComplementOf(s2);

  respective_sequence_positions(DNAmat, s1, revcomp_s2,
				s1ins2s_tmp2, s1ins2e_tmp2, s2ins1s_tmp2, s2ins1e_tmp2, alignment_score_tmp2, alignment_length_tmp2,
				verbosity);

  if (alignment_score_tmp1 > alignment_score_tmp2)
  {
    s1ins2s_best          = s1ins2s_tmp1;
    s1ins2e_best          = s1ins2e_tmp1;
    s2ins1s_best          = s2ins1s_tmp1;
    s2ins1e_best          = s2ins1e_tmp1;
    alignment_score_best  = alignment_score_tmp1;
    alignment_length_best = alignment_length_tmp1;
    better = 1;
  }
  else
  {
    s1ins2s_best          = s1ins2s_tmp2;
    s1ins2e_best          = s1ins2e_tmp2;
    s2ins1s_best          = s2ins1s_tmp2;
    s2ins1e_best          = s2ins1e_tmp2;
    alignment_score_best  = alignment_score_tmp2;
    alignment_length_best = alignment_length_tmp2;
    better = 2;
  }

  if (verbosity >5)
  {
    cout << "Best alignment:" << endl;
    cout << "Score:     " <<  alignment_score_best << endl;
    cout << "Direction: " <<  (better==1 ? "Same" : "RevComp") << endl;
  }

}


bool check_full_coverage_and_required_taxa_in_range(CSequences2 *seqs, unsigned pos1, unsigned pos2,
						    const vector<unsigned> &unique_taxon_ids,
						    const CRequiredTaxa    &required_groups_of_taxa)
{
  unsigned N = seqs->GetTaxaNum();
  unsigned i;

  CSequence_Mol*      seq;
  set<unsigned>       set_of_taxon_ids;
  faststring          taxon_name;

  // Determine set of complete sequences (taxa) of this range.

  // For each sequence in this range we check that it does not contains gaps or Ns.
  // In the set_of_taxon_ids we store the unique taxon ID, not the sequence number.
  for (i=0; i<N; ++i)
  {
    seq = seqs->get_seq_by_index(i);

    // Has no gaps, Ns in range?
    if (! seq->range_contains_gaps_or_Ns(pos1, pos2) )
    {
      // *seq is a good sequence in this range. (No gaps or Ns.)
      set_of_taxon_ids.insert(unique_taxon_ids[i]);
    }
  }

  // Checks whether there is any sequence left. In case of an empty required_taxa string
  // this is an important check. Otherwise gap only regions could pass this test.
  if (set_of_taxon_ids.empty())
  {
    if (global_VERBOSITY >= VERBOSITY_COORDINATES)
    {
      cout << "NOTE: Not a good bait window. No complete sequences in this range: " << pos1+1 << "-" << pos2
	   << " Num seqs: " << set_of_taxon_ids.size() << endl;
    }
    return false;
  }

  if (required_groups_of_taxa.is_empty() )
  {
    if (global_VERBOSITY >= VERBOSITY_COORDINATES)
    {
      cout << "NOTE: Good bait window. Valid sequence range (no required taxa specified): " << pos1+1 << "-" << pos2 << " Num seqs: " << set_of_taxon_ids.size() << endl;
    }
    return true;
  }

  // Check that all required taxonomic groups are present in this feature:
  if (required_groups_of_taxa.check_precense_of_all_groups_of_taxa(set_of_taxon_ids) )
  {
    if (global_VERBOSITY >= VERBOSITY_COORDINATES)
    {
      cout << "NOTE: Good bait window. All required taxonomic groups present in this range: " << pos1+1 << "-" << pos2  << " Num seqs: " << set_of_taxon_ids.size() << endl;
    }
    return true;
  }
  else
  {
    if (global_VERBOSITY >= VERBOSITY_COORDINATES)
    {
      cout << "NOTE: Not a good bait windows. NOT all required taxonomic groups present in this range." << pos1+1 << "-" << pos2  << " Num seqs: " << set_of_taxon_ids.size() << endl;
    }
    return false;
  }
} // END check_full_coverage_and_required_taxa_in_range(...)


// Return true if at least one good bait region was found, otherwise it returns false.
bool handle_this_alignment(CSequences2           *alignment,
			   CTaxonNamesDictionary &taxonNamesDict,
			   CRequiredTaxa         &required_groups_of_taxa,
			   unsigned short        bait_length,
			   unsigned short        number_of_baits_in_region,
			   unsigned short        bait_offset,
			   faststring            name,     // gene name
			   unsigned short        running_number,
			   unsigned short        running_number_maximum,
			   double                cluster_threshold,
			   ostream               &os_result_file,
			   char                  delim,
			   unsigned short        field_num,
			   char                  center_computation_mode,
			   const faststring      &filename)
{
  unsigned alignment_length = alignment->GetPosNum();
  int      region_length    = (number_of_baits_in_region-1)*bait_offset+bait_length;

  vector<unsigned>      seqNum_to_unique_taxonID_this_alignment; // Used to translate the sequnce number in the current file to the unique taxon number.
  vector<faststring>    full_names_this_alignment;
  vector<faststring>    taxon_names_this_alignment;
  set<faststring>       taxon_set_this_alignment;

  //DEBUG BREAK POINT - allows us to catch certain gene names for debugging purposes.
//   {
//     //    os_result_file << "name: " << name << " " << running_number << endl;
//
//     if (name == "NV18435" && running_number == 4)
//     {
//       cerr << "Found name: " << name << " " << running_number << endl;
//     }
//
//   }


  //************************************************
  // The following four arrays belong together.
  //************************************************
  // For each site in this msa is_valid_bait_window stores the information
  // whether the site represents a good starting point and whether
  // for the window starting here we need to do the clustering and determine
  // the center sequence.
  // For each site in the msa is_valid_bait_region stores the information
  // whether we have a good bait region starting at this point.
  // For each site cluster_loci_this_alignment stores the information about
  // clusters and center sequences for the bait window starting at this point.
  // If the window does not need to be computed, the pointer in this array is NULL.
  faststring is_valid_bait_window;
  faststring is_valid_bait_region;
  vector<Csequence_loci_cluster_collection *> cluster_loci_this_alignment(alignment_length, NULL);
  vector<CBait_region_information *> bait_region_info(alignment_length, NULL);

  // TODO: currently, the is_valid_bait_window and is_valid_bait_region are filled up to the end.
  //       Last valid/possible index for is_valid_bait_window: alignment_length - bait_length
  //       Last valid/possible index for is_valid_bait_region: alignment_length - region_length

  // Initialize is_valid_bait_window with fail character '-':
  is_valid_bait_window.clear();
  is_valid_bait_window.assign(alignment_length, '-');

  is_valid_bait_region.clear();
  is_valid_bait_region.assign(alignment_length, '-');

  unsigned i;
  for (i=alignment_length - bait_length+1; i< alignment_length; ++i)
    is_valid_bait_window[i] = 'e';  // Indicates that this region extends beyond the end.

  for (i=alignment_length - region_length+1; i< alignment_length; ++i)
    is_valid_bait_region[i] = 'e';     // Indicates that this region extends beyond the end.

//   cout << is_valid_bait_window.size() << endl;
//   cout << is_valid_bait_window << endl;

//   cout << is_valid_bait_region.size() << endl;
//   cout << is_valid_bait_region << endl;

  // Determine translation table: seqNum_to_unique_taxonID_this_alignment
  // Clear the data for this alignment:
  seqNum_to_unique_taxonID_this_alignment.clear();
  full_names_this_alignment.clear();
  //	      taxon_names_this_alignment.clear();
  //	      taxon_set_this_alignment.clear();

  alignment->get_sequence_names(full_names_this_alignment);

  unsigned j;
  faststring tmp_name;

  // Determine vector that translates sequence number to unique_taxonIDs
  for (j=0; j < full_names_this_alignment.size(); ++j)
  {
    tmp_name = extract_taxon_name(full_names_this_alignment[j], delim, field_num);
    //		taxon_names_this_alignment.push_back(sequence_to_taxon_name(full_names_this_alignment[j]));
    //		taxon_set_this_alignment.insert(taxon_names_this_alignment[j]);
    //		cout << taxon_names_this_alignment[j] << endl;

    seqNum_to_unique_taxonID_this_alignment.push_back(taxonNamesDict.dictionary_get_index_of_taxonname(tmp_name));
    //		cout << j << "-->" << tmp_name << "->" << taxonNamesDict.dictionary_get_index_of_taxonname(tmp_name) << endl;
  }

  // Move through alignment sequence alignment and determine those parts which have a sufficient data coverage:
  {
    unsigned pos1, pos2;
    bool     is_valid_range;

    pos1 = 0; // 0 based
    pos2 = pos1 + bait_length;

    // pos2 is the index after the last nucleotide in the bait.
    // Thus pos2 is allowed to be equal to alignment_length, which specifies the index after
    // the last nucleotide in the sequence.
    while (pos2 <= alignment_length)
    {
      is_valid_range = check_full_coverage_and_required_taxa_in_range(alignment,
								      pos1, pos2,
								      seqNum_to_unique_taxonID_this_alignment,
								      required_groups_of_taxa);
      if (is_valid_range)
	is_valid_bait_window[pos1] = 'x';  // Indicates a valid bait window.

      ++pos1;
      ++pos2;
    }
    //    cout << "Valid bait windows in alignment: " << endl << is_valid_bait_window << endl;
  } // END Move through alignment sequence alignment

//   cout << "Valid bait windows: " << endl;
//   cout << is_valid_bait_window << endl;

  unsigned num_good_regions=0;

  // Determine good bait regions, i.e. with requested number of baits and starting points separated by an offset of bait_offset
  {
    // for all starting point of a bait region:
    unsigned i, j, pos;
    bool     is_good_region;

    // i is the start index of the current range.
    for (i=0; i < alignment_length; ++i)
    {
      if (is_valid_bait_region[i] != 'e')
      {
	//	cout << "Starting good: " << i << endl;
	is_good_region = true;
	// For all number_of_baits_in_region baits in one region with offset bait_offset
	for (j=0, pos=i; j < number_of_baits_in_region; ++j, pos += bait_offset)
	{
	  // Allowed values for  is_valid_bait_window[pos] are 'c' and 'x'. Not allowed are: '-' and 'e'
	  // Some checks are redundant here: E.g. pos >= alignment_length and is_valid_bait_window[pos]=='e'
	  // should not be possilbe, since we check is_valid_bait_region[i] != 'e' before we get here.
	  // Doing these checks is more secure.
	  // TODO: In future versions check whether pos >= alignment_length and is_valid_bait_window[pos]=='e'
	  //       are necessary.
	  if (pos >= alignment_length || is_valid_bait_window[pos] == '-' || is_valid_bait_window[pos]=='e')
	  {
	    //	    cout << "Not a good bait region: " << i << " " << pos << " " <<   is_valid_bait_window[pos] << endl;
	    is_good_region = false;
	    break; // This region is at best incomplete, so we mark it as bad.
	  }
	}
	if (is_good_region)
	{
	  //	  cout << "Mark as good: " << i << endl;
	  ++num_good_regions;
	  // Mark good region.
	  is_valid_bait_region[i] = 'x';

	  //	  cout << "Interm. " << i << is_valid_bait_region << endl;
	  
	  // For all good bait regions we also allocate memory for the final colletion of the data:
	  bait_region_info[i] = new CBait_region_information(name, running_number, running_number_maximum, i, // bait_length, bait_offset,
							     number_of_baits_in_region);
	  //		    cout << "Intermediate is_valid_bait_region: " << endl << is_valid_bait_region << endl;
	  // For this valid bait region, we mark all windows as "need to be computed".

	  ++DEBUG_COUNT_NEW_CALLS_CBAIT_REGION_INFORMATION;
	  
	  for (j=0, pos=i; j < number_of_baits_in_region; ++j, pos += bait_offset)
	  {
	    if (pos >= alignment_length)
	      break;
	    is_valid_bait_window[pos] = 'c';
	  }
	  //		    cout << "Intermediate is_valid_bait_window: " << endl << is_valid_bait_window << endl;
	  if (global_VERBOSITY >= VERBOSITY_COORDINATES)
	  {
	    cout << "NOTE: Good bait region: " << i+1 << " " << i+1+region_length << endl;
	  }
	}
	else
	{
	  if (global_VERBOSITY >= VERBOSITY_COORDINATES)
	  {
	    cout << "NOTE: Not a good bait region: " << i+1  << " " << i+1+region_length << endl;
	  }
	}
      } // END if (is_valid_bait_region[i] != 'e')
    } // END for (i=0; i < alignment_length; ++i)
    //          cout << "Valid bait regions in feature: " << endl << is_valid_bait_region << endl;
    //          cout << "Bait windows that need to be computed" << endl << is_valid_bait_window << endl;
  }   // END Determine baits with requested number_of_baits_in_region good baits separated by an offset of bait_offset

  // If we come here we have the following information:
  //	      is_valid_bait_window contains a 'c' in all bait window starting indices that need to be looked at.
  //              is_valid_bait_region contains a 'x' in all bait region starting indices that we will explore,
  //                                    i.e. for which we have the necesarry successive bait windows.

  // If this file has one or more regions that should be looked at:
  // Cluster loci which need to be clustered.
  // Compute center sequences.

  if (num_good_regions > 0)
  {
    // OK - Where are we?
    // We have extracted an feature. We determined the loci of possible
    // bait regions in this feature. This number > 0.
    // For all bait windows in this feature that belong to good bait regions,
    // we cluster the sequences and compute an artificial center sequence for 
    // each cluster.

    unsigned i;
    unsigned num_taxa    = alignment->GetTaxaNum();

    // We need a recoded msa in order to compute artificial center sequences.
    // We could do this for each window, but since the windows overlap, we
    // do a lot of redudant work. Doing this per window we would only be on the positive
    // side if there would be only very few windows in a long alignment.
    // I do not see this as the general case.
    // The recoded msa is used for the complete file and only needs to be created once.

    // The vector of faststrings allocates the memory for msa:
    std::vector<faststring>  msa_faststring;
    char** msa;
    msa = new char* [num_taxa];  // remember to delete this dynamic array
    msa_faststring.clear();
    ++DEBUG_COUNT_NEW_CALLS_MSA;

    // Copy the sequence data to the msa_faststring
    // The data is stored in the alignment (type CSequences2) object.
    alignment->get_alignment_as_copy(msa_faststring);

    for (i=0; i < num_taxa; ++i)
    {
      // Important: Even though they are initially 0 terminated cstrings,
      //            this does not provide much advantage since the sequences
      //            will be recoded anyway. In the recoded form, the bit code 0
      //            is used. So recoded sequences cannot be interpreted as cstrings.
      msa[i]  = msa_faststring[i].data();
    }

    // Recode the msa:
    // IMPORTANT: Recoding is possible also for ambiguous characters.
    //            Gaps are also recoded.
    recode_MSA(msa, num_taxa, alignment_length);

    unsigned num_valid_bait_windows   = 0;
    unsigned count_valid_bait_windows = 0;
    unsigned report_interval          = 0;

    for (i=0; i < alignment_length; ++i)
    {
      if (is_valid_bait_window[i] == 'c')
      {
	++num_valid_bait_windows;
      }
    }

    report_interval = num_valid_bait_windows/10;
    if (report_interval == 0)
      report_interval = 1;

    // For all loci:
    for (i=0; i < alignment_length; ++i)
    {
      // Do we have to compute cluster and center sequence for this loci?
      if (is_valid_bait_window[i] == 'c')
      {
	++count_valid_bait_windows;

	if (global_VERBOSITY > VERBOSITY_MOREPROGRESS)  // Report progress
	{
	  if (count_valid_bait_windows % report_interval == 0)
	  {
	    int p = count_valid_bait_windows*100 /  num_valid_bait_windows;
	    faststring t;
	    get_time_string(t);
	    cout << "PROGRESS: " << t << " Percent done this alignment or gene: " << p << "%."
		 << endl << flush;
	  }
	}

	// This constructor again determins the valid sequences and their number in this range
	Csequence_loci_cluster_collection* tmp_collection = new Csequence_loci_cluster_collection(i, i+bait_length,
												  alignment,
												  // taxonNamesDict,
												  cluster_threshold, msa );
	++DEBUG_COUNT_NEW_CALLS_CSEQUENCE_LOCI_CLUSTER_COLLECTION;
	cluster_loci_this_alignment[i] = tmp_collection;
	// Compute what we can compute now for this locus.
	tmp_collection->cluster_locus(name, filename);
	tmp_collection->compute_center_for_all_clusters(center_computation_mode);

	// Debug output:
// 	cout << "*******************" << endl;
// 	tmp_collection->print_msa_backrecoded_all(cout, 1);
// 	cout << "*******************" << endl;
      }
    } // END for all loci

    // Now we can transfer the clustering and center sequence data to the bait_region_info vector.

    unsigned short tmp_num_sequences;
    unsigned short tmp_num_clusters;

    // For all possible starting positions of bait regions.
    for (i=0; i < alignment_length; ++i)
    {
      if (is_valid_bait_region[i] == 'x')
      {
	// Collect the informaiton for all number_of_baits_in_region baits in this bait region.
	for (j=0; j < number_of_baits_in_region; ++j)
	{
	  if (cluster_loci_this_alignment[i+j*bait_offset] != NULL)
	  {
	    if (bait_region_info[i] == NULL)
	    {
	      cerr << "INTERNAL ERROR: bait_region_info["<< i << "]==NULL" << endl;
	      exit(-4);
	    }
	    cluster_loci_this_alignment[i+j*bait_offset]->get_numbers_of_valid_sequences_number_of_clusters(tmp_num_sequences, tmp_num_clusters);
	    bait_region_info[i]->add_info__cluster_locus(cluster_loci_this_alignment[i+j*bait_offset],
							 tmp_num_sequences, tmp_num_clusters);
	  }
	  else
	  {
	    cerr << "INTERNAL ERROR: The cluster_loci_this_alignment is unexpectedly NULL at position: "
		 << i+j*bait_offset << endl;
	    exit(-5);
	  }
	}
      }
    }

    for (i=0; i< alignment_length; ++i)
    {
      if (is_valid_bait_region[i]=='x')
      {
	bait_region_info[i]->print_baits_result(os_result_file, filename);
      }
    }


    // The msa array is not needed any more, so we delete it.
    // The acutual sequences are stored in the msa_faststring
    // vector and they will be deleted automatically when this block will be left.
    delete [] msa;
    ++DEBUG_COUNT_DELETE_CALLS_MSA;
  } // END if (num_good_regions > 0) // Compute cluster and center sequence for loci.
  else
  {
    if (global_VERBOSITY >= VERBOSITY_WARNINGS)
    {
      if (running_number > 0) // In non cutting mode this 0, otherwise > 0
      {
	cerr << "ALERT: No valid bait regions found for feature: ";
	cerr << name << " " << running_number << "/" << running_number_maximum;
	cerr << endl;
      }
      else
      {
	cerr << "ALERT: No valid bait region found for alignment: " << name << endl;
      }
    }
  }


  // Remember to release all memory allocated in this block:

  DEBUG_COUNT_DELETE_CALLS_CBAIT_REGION_INFORMATION          += delete_pointers_in_vector(bait_region_info);
  DEBUG_COUNT_DELETE_CALLS_CSEQUENCE_LOCI_CLUSTER_COLLECTION += delete_pointers_in_vector(cluster_loci_this_alignment);

  if (num_good_regions > 0)
    return true;
  else
    return false;

} // END handle_this_alignment()

// Only used in special mode to test gff file integrety:
void gff_file_test(faststring gff_file_name)
{
  GFF_collection gff_col;

  unsigned good_lines;
  unsigned bad_lines;
  bool     file_error;

  if (global_VERBOSITY >= VERBOSITY_MOREPROGRESS)
  {
    cout << "Parsing the gff file ..." << flush;
  }

  gff_col.readGFFfile(gff_file_name.c_str(), good_lines, bad_lines, file_error);

  if (global_VERBOSITY >= VERBOSITY_DEBUG)
  {
    cout << " Done." << endl;
    cout << "Good lines: " << good_lines << endl;
    cout << "Bad lines:  " << bad_lines  << endl;
    cout << "File error: " << file_error << endl;
    cout << endl;
    cout << endl;
  }

  if (file_error)
  {
    cerr << "ERROR: The specified gff-file could not be opened. It seems it does not exist." << endl;
    cerr << "Exiting." << endl;
    exit(-3);
  }

  gff_col.print(cout, 1);
}


int main(int argc, char **argv)
{
  faststring      fname_parameter_file;

  faststring      fname_gff;
  faststring      fname_genome;
  faststring      string_required_taxa;
  faststring      dname_alignments;
  faststring      genome_name;

  unsigned short  bait_length;
  double          cluster_threshold;
  faststring      output_directory;

  unsigned  short bait_offset;
  unsigned  short number_of_baits_in_region;

  bool            remove_reference_species_from_alignments;
  bool            alignment_cutting;

  faststring      gff_feature_name;
  faststring      gff_record_attribute_name_for_geneID;

  char            sequence_name_field_delimiter;
  unsigned short  sequence_name_taxon_field_number;    // Used to find reference species and test for required taxa.
  unsigned short  sequence_name_geneID_field_number;   // Only used if extracting with the aid of a gff file.
  bool            sequence_name_tokenizer_overwritten;


  char            center_computation_mode;
  int             special_mode = 0; // 0: normal, 1: gff-test


  if (argc == 2)
  {
    fname_parameter_file = argv[1];
    parse_parameter_file(fname_parameter_file,
			 fname_gff,
			 fname_genome,
			 string_required_taxa,
			 dname_alignments,
			 genome_name,
			 bait_length,
			 cluster_threshold,
			 output_directory,
			 bait_offset,
			 number_of_baits_in_region,
			 remove_reference_species_from_alignments,
			 alignment_cutting,
			 gff_feature_name,
			 gff_record_attribute_name_for_geneID,
			 sequence_name_field_delimiter,
			 sequence_name_taxon_field_number,
			 sequence_name_geneID_field_number,
			 sequence_name_tokenizer_overwritten,
			 center_computation_mode,
			 global_VERBOSITY
			 );
  }
  else if (argc == 3)
  {
    if (faststring(argv[2])=="gffFileTest")
    {
      special_mode = 1;
    }
    else
    {
      cerr << "Welcome to the " << PROGNAME << " program, version " << VERSION << endl << endl;
      cerr << "Usage: " << argv[0] << " parameter-file [gffFileTest]" << endl;
      //    cerr << "OR" << endl;
      //    cerr << "Usage: " << argv[0] << " gff-file genome-file transcriptome-folder required_taxa_file probe-length cluster-threshold output-directory" << endl;
      exit(-3);
    }

    fname_parameter_file = argv[1];
    parse_parameter_file(fname_parameter_file,
			 fname_gff,
			 fname_genome,
			 string_required_taxa,
			 dname_alignments,
			 genome_name,
			 bait_length,
			 cluster_threshold,
			 output_directory,
			 bait_offset,
			 number_of_baits_in_region,
			 remove_reference_species_from_alignments,
			 alignment_cutting,
			 gff_feature_name,
			 gff_record_attribute_name_for_geneID,
			 sequence_name_field_delimiter,
			 sequence_name_taxon_field_number,
			 sequence_name_geneID_field_number,
			 sequence_name_tokenizer_overwritten,
			 center_computation_mode,
			 global_VERBOSITY
			 );

  }
//   else if (argc == 8)
//   {
//     fname_gff           = argv[1];
//     fname_genome        = argv[2];
//     fname_required_taxa = argv[3];
//     dname_alignments = argv[4];
//     genome_name         = "Nasonia_vitripennis";

//     bait_length        = faststring(argv[5]).ToInt();
//     cluster_threshold   = faststring(argv[6]).ToDouble();
//     output_directory    = argv[7];
//   }
  else
  {
    cerr << "Welcome to the " << PROGNAME << " program, version " << VERSION << endl << endl;
    cerr << "Usage: " << argv[0] << " parameter-file [gffFileTest]" << endl;
    //    cerr << "OR" << endl;
    //    cerr << "Usage: " << argv[0] << " gff-file genome-file transcriptome-folder required_taxa_file probe-length cluster-threshold output-directory" << endl;
    exit(-3);
  }


  print_parameters(cout,
		   fname_gff,
		   fname_genome,
		   string_required_taxa,
		   dname_alignments,
		   genome_name,
		   bait_length,
		   cluster_threshold,
		   output_directory,
		   bait_offset,
		   number_of_baits_in_region,
		   remove_reference_species_from_alignments,
		   alignment_cutting,
		   gff_feature_name,
		   gff_record_attribute_name_for_geneID,
		   sequence_name_field_delimiter,
		   sequence_name_taxon_field_number,
		   sequence_name_geneID_field_number,
		   sequence_name_tokenizer_overwritten,
		   center_computation_mode,
		   global_VERBOSITY
		  );


  if (special_mode == 1)
  {
    gff_file_test(fname_gff);
    exit(0);
  }
  //  exit(-111);


  unsigned short   min_alignment_length     = bait_length + (number_of_baits_in_region-1)*bait_offset;
  CScoreMatrix DNAmat("simpledna_N0.mat");

  // Do some initial checks.


//   if (!file_exists(fname_required_taxa.c_str()))
//   {
//     cerr << "I could not open required taxa file: " << fname_required_taxa << endl;
//     cerr << "I guess it does not exist. Exiting." << endl;
//     exit(-1);
//   }


  if (cluster_threshold <= 0 || cluster_threshold >= 1)
  {
    cerr << "The cluster threshold must be greater than 0 and smaller than 1." << endl;
    exit(-1);
  }

  if (bait_length < 1)
  {
    cerr << "Probe length cannot be 0 or smaller." << endl;
    exit(-1);
  }

  if (!file_or_dir_exists(dname_alignments.c_str()))
  {
    cerr << "I could not read the transcriptome directory: " << dname_alignments << endl;
    cerr << "I guess it does not exist. Exiting." << endl;
    exit(-1);
  }

  if (!file_or_dir_exists(output_directory.c_str()))
  {
    cerr << "The output-directory " << output_directory << endl;
    cerr << "could not be opened. I guess it does not exist. Please make sure this directory exists. Exiting." << endl;
    exit(-1);
  }

  // Open file for output:
  faststring result_filename = output_directory + "/loci_baits.txt";
  ofstream   os_result_file(result_filename.c_str());

  if (alignment_cutting)
    os_result_file << "# Corresponding sequence name\tGene-name\tfeature-number\tnumber-of-features-of-specified-type-for-this-gene-in-gff-file\tstart-coordinate-of-bait-region-in-alignment\tnumber-of-sequences-in-region-the-baits-are-based-on\tnumber-of-baits-constructed-for-this-bait-region\tList of baits - each bait is followed by CG content and the maximum distance to the sequences in the alignment." << endl;
  else
    os_result_file << "# Corresponding sequence name\tGene-name\tn.a.\tn.a.\tstart-coordinate-of-bait-region-in-alignment\t number-of-sequences-in-region-the-baits-are-based-on\tnumber-of-baits-constructed-for-this-bait-region\tList of baits - each bait is followed by CG content and the maximum distance to the sequences in the alignment." << endl;

  // Read names of alignment files:
  list<faststring> fasta_names;
  list<faststring> filter_names;

  // If we use the Qt library:
  //   filter_names.push_back("*.fas");
  //   filter_names.push_back("*.fa");
  //   filter_names.push_back("*.fasta");

  filter_names.push_back("fas");
  filter_names.push_back("fa");
  filter_names.push_back("fasta");

  get_file_names_in_dir(fasta_names, dname_alignments, filter_names);



  //****************************************************************
  // BEGIN Alignment cutting mode
  //****************************************************************
  if (alignment_cutting)
  {
    if (!file_exists(fname_gff.c_str()))
    {
      cerr << "I could not open the gff file: " << fname_gff << endl;
      cerr << "I think it does not exist. Exiting." << endl;
      exit(-1);
    }

    if (!file_exists(fname_genome.c_str()))
    {
      cerr << "I could not open the genome file: " << fname_genome << endl;
      cerr << "I think it does not exist. Exiting." << endl;
      exit(-1);
    }


    // Read gff file:

    unsigned good_lines;
    unsigned bad_lines;
    bool     file_error;

    GFF_collection gff_col;

    if (global_VERBOSITY >= VERBOSITY_PROGRESS)
    {
      cout << "Parsing the gff file ... " << flush;
    }
    gff_col.readGFFfile(fname_gff.c_str(), good_lines, bad_lines, file_error);
    if (global_VERBOSITY >= VERBOSITY_PROGRESS)
    {
      cout << "Done." << endl;
    }

    if (good_lines == 0 || bad_lines > 0 || file_error)
    {
      cerr << "A problem occured while reading the first gff file." << endl;
      cerr << "Good lines: " << good_lines << endl;
      cerr << "Bad lines:  " << bad_lines  << endl;
      cerr << "File error: " << file_error << endl;

      if (file_error)
      {
	cerr << "The specified gff-file could not be opened. It seems it does not exist." << endl;
      }
      cerr << "Exiting." << endl;
      exit(-3);
    }

    //  gff_col.print(cout, 0);

    // Read genome file:

    if (global_VERBOSITY >= VERBOSITY_PROGRESS)
      cout << "Reading the genome file ..." << flush;

    CSequences2 genome_seqs(CSequence_Mol::dna);
    CSequence_Mol::processing_flag pflag = CSequence_Mol::processing_flag(0);

    pflag = CSequence_Mol::processing_flag(pflag | CSequence_Mol::convert_toupper);
    //  pflag = CSequence_Mol::processing_flag(pflag | CSequence_Mol::convert_ambig2N);

    if (global_VERBOSITY >= 50)
    {
      cout << "with processing flag: " << pflag << "..." << endl;
    }

    CFile file;
    file.ffopen(fname_genome.c_str());
    if (genome_seqs.read_from_Fasta_File(file, pflag, 0, -1, false) < 0)
    {
      cerr << "ERROR while reading the fasta file: " << fname_genome << endl;
      cerr << "Please check the line indicated above since this line prevented me" << endl
	   << "from parsing the fasta file." << endl;
      exit(-2);
    }
    file.ffclose();

    if (global_VERBOSITY >= VERBOSITY_PROGRESS)
    {
      cout << "Done." << endl;

      if (global_VERBOSITY >= VERBOSITY_FILEINFO)
      {
	cout << "NOTE: Found this number of sequences in genome file: " << genome_seqs.GetTaxaNum()  << endl;
	cout << "NOTE: Found this number of positions in genome file: " << genome_seqs.GetPosNum()   << endl;
      }
      faststring t;
      get_time_string(t);
      cout << "PROGRESS: " << t << " Genome file has been read." << endl << flush; 
    }

    // Read fasta sequnces one by one:
    list<faststring>::iterator it, it_end;
    it     = fasta_names.begin();
    it_end = fasta_names.end();

    unsigned        Num_transcripts   = fasta_names.size();
    unsigned        count_transcripts = 0;
    unsigned        count_transcripts_no_gff_entries = 0;
    unsigned        count_transcripts_no_gff_entries_for_feature = 0;

    CSequence_Mol   *transcript_seq;
    unsigned        transcript_seq_index;
    faststring      full_alignment_file_name;
    faststring      sequence_header_target;
    faststring      gene_name;
    GFF_record_list gene_record_list;
    faststring      genomic_seq_name;
    //  unsigned        overall_counter  = 0;
    //  double          mean_match_score = DNAmat.get_mean_matchscore_of_first_N_symbols(4);

    CSequences2     *feature_transcript_alignment = NULL, *transcript_alignment = NULL;
    CSequences2     *feature_transcript_alignment_rmgo = NULL;

    unsigned        feature_number;
    unsigned        features_this_gene;
    bool            no_gff_entries_warning;


    //  cout << "Mean match score: " << DNAmat.get_mean_matchscore_of_first_N_symbols(4) << endl;

    // We have to keep a list of all taxon names that we deal with in this analysis.
    // We give the constructor the name of the directory that contains all fasta files.
    // It will read all fasta files, determine the taxon names from the sequence names and add them
    // to the list of taxon names. Furthermore, this object will keep a map that maps taxon names
    // to a unique taxon index in this analysis. Internally it is important to have such a unique index.

    // The unique taxon indeces are used:
    // -  when checking for all required taxa.
    //    The required taxa are stored as sets of unique taxon ids.
    //    To check that there is at least one representative for each required set,
    //    we determine the unique IDs for the taxa that are present.
    //    Comparing intersections of sets then tells us whether all required taxonomic groups are present.

    CTaxonNamesDictionary taxonNamesDict(dname_alignments, fasta_names, sequence_name_taxon_field_number, sequence_name_field_delimiter);
    //  Print dictionary:
    //  taxonNamesDict.print(os_result_file);

    // Read required groups of taxa:
    int error_required_taxa;
    CRequiredTaxa required_groups_of_taxa(string_required_taxa.c_str(), taxonNamesDict, error_required_taxa);

    { // Debug output:
      if (global_VERBOSITY > 1000)
      {
	cout << endl;
	taxonNamesDict.print_taxon_names_to_index(cout);
	cout << "Required groups: From each of the following sets of taxa, we require at least one to be present in the data set." << endl;
	required_groups_of_taxa.print(cout, 1);
	cout << endl;
      }
    }

    // For all transcript files:
    // We extract the gene name and use the gene name and feature name to extract all gff records corresponding to these two.
    // For each gff record we cut out the region from the alignment file that corresponds to the genomic region in the genome
    // file depicted by the gff-record. In this sub alignment we search for baits.

    // For all transcript files:
    while (it != it_end)
    {
      no_gff_entries_warning = false;
      ++ count_transcripts;
      full_alignment_file_name = dname_alignments + "/" + *it;

      if (get_file_size(full_alignment_file_name) > 0)
      {
	// Here is the main action: We handle this gene:

	// For this gene file, read the transcript for the genome we are concerned with

	if (global_VERBOSITY >= VERBOSITY_PROGRESS) // Report progress
	{
	  faststring t;
	  get_time_string(t);
	  cout << "PROGRESS: " << t << " Reading alignemnt file ("
	       << count_transcripts << "/" << Num_transcripts << "): "
	       << full_alignment_file_name << endl << flush; 
	}
	  
	read_alignment_sequence_of_species2(full_alignment_file_name,
					    genome_name,
					    transcript_seq,
					    transcript_alignment,
					    transcript_seq_index,
					    sequence_name_field_delimiter,
					    sequence_name_taxon_field_number);

	if (transcript_alignment == NULL)
	{
	  cerr << "ERROR: Could not read transcript file: "
	       << full_alignment_file_name << endl;
	  exit(-4);
	}

	if (transcript_seq == NULL) // This is rather common. But in this case we have to skip this alignment.
	{
	  if (global_VERBOSITY >= VERBOSITY_WARNINGS)
	  {
	    cerr << "WARNING: Requested species " << genome_name << " not found in sequence: "
		 << full_alignment_file_name << endl;
	    cerr << "         This means that we have to skip this alignment, since it cannot be handled in this mode." << endl;
	  }
	}
	else
	{
	  sequence_header_target = transcript_seq->getFullName();
	  gene_name = extract_gene_name(sequence_header_target, sequence_name_field_delimiter, sequence_name_geneID_field_number);

	  if (global_VERBOSITY >= VERBOSITY_DEBUG)
	  {
	    cout << "INFO: Name of sequence containing the reference species: " << sequence_header_target << endl;
	    cout << "INFO: Gene name:                                         " << gene_name << endl;

	    transcript_seq->writeSequence(cout, 100);
	  }

	  gff_record_attribute_name_for_geneID.ToLower();

	  // Determine coordinates of this gene in genome file using the gff file:
	  gff_col.get_List_of_gff_records_for_attribute(gff_record_attribute_name_for_geneID.c_str(), gene_name.c_str(), gene_record_list);

	  if (gene_record_list.size() == 0 && global_VERBOSITY >= VERBOSITY_WARNINGS)
	  {
	    cerr << "WARNING: For gene \"" << gene_name << "\" no entries where found in the gff file.\n"
	            "         While this is possible, in most cases the reason is either a wrong attribute name for the\n"
	            "         gene-name-field in the gff-file or the gene name specified in the sequence\n"
	            "         name of the transcript alignment is incompatible to the gene names in the gff file.\n"
		 << "         It is recommended to check whether gene name \"" << gene_name << "\" occurs in the\n"
	            "         gff file under the attribute name \"" << gff_record_attribute_name_for_geneID << "\".\n";
	    cerr << "         Note that the gene name as specified in the sequence name and the gene name in\n"
	            "         the gff file have to match exactly. The attribute name specified in the parameter\n"
	            "         file must be an exact (case insensitive) match to the attribute name in the gff file." << endl << endl;
	    ++count_transcripts_no_gff_entries;
	    no_gff_entries_warning = true;
	  }

	  easystring filter_feature = gff_feature_name.c_str();
	  gene_record_list.filter_by_feature(filter_feature);

	  // TODO: Check empty list.
	  //	  if (gene_record_list)


	  if (global_VERBOSITY > VERBOSITY_FILEINFO)
	  {
	    // Print gff information:
	    cout << "Features with this gene name: " << endl;
	    cout << "----" << endl;
	    gene_record_list.print(cout);
	    cout << "----" << endl;
	  }

	  // Foreach feature in this gene:
	  list<GFF_record*>::iterator gl_it, gl_it_end;

	  gl_it     = gene_record_list.begin();
	  gl_it_end = gene_record_list.end();

	  // NOTE: We count all genes, not only the once that have a sufficient length.
	  // TODO: Check that we really want to count all genes, not only the once that are long enought.

	  features_this_gene = gene_record_list.size();

	  if (features_this_gene == 0)
	  {
	    if (no_gff_entries_warning == false && global_VERBOSITY >= VERBOSITY_WARNINGS)
	    {
	      cerr << "WARNING: For gene \"" << gene_name << "\" no entries where found in the gff file that matched the specified feature name.\n"
	           << "         Since entries for this gene name have been found, but with different feature names, the problem could be the feature name\n"
		   << "         or the fact that features exist for this gene." << endl; 
	    }
	    ++count_transcripts_no_gff_entries_for_feature;
	  }

	  // For the gene in the transcript file we are working on we found exi
	  // in the gff file. I.e. the gene is annotated in the genome.
	  // The genomic sequences of the exi will be extracted and aligned
	  // to the sequence in the transcript file. This determines the coordinates
	  // of the exi in the transcript file.

	  CDnaString     feature_sequence;
	  CDnaString     transcript_sequence;

	  CSequence_Mol* genome_seq_with_feature;

	  faststring     seqid_feature_in_genome;
	  unsigned       start_of_feature_in_genome;
	  unsigned       end_of_feature_in_genome;    // Remember: This is an inclusive end

	  transcript_sequence = transcript_seq->getSeqStr();

	  // The transcript sequence can contain gaps. These cause problems when aligning features against the transcripts.
	  // Removing gaps is not an option, since we need to preserver the original coordiantes.
	  // Converting gaps to Ns will preserve coordinates. It will do no harm, since Ns will nicely align to the features and the sequence will
	  // not be used for anything except aligning features to it.

	  transcript_sequence.replace_char('-','N');

	  /*
	    vector<unsigned>      seqNum_to_unique_taxonID_this_feature; // Used to translate the sequnce number in the current file to the unique taxon number.
	    vector<faststring>    full_names_this_feature;
	    vector<faststring>    taxon_names_this_feature;
	    set<faststring>       taxon_set_this_feature;

	    unsigned   feature_transcript_alignment_length;
	  */

	  unsigned features_with_at_least_one_bait_regions = 0;
	  feature_number = 0;
	  // For each feature of the gene under consideration:  #foreach feature
	  while (gl_it != gl_it_end)
	  {
	    ++feature_number;

	    // The following coordinates are 0 based and the End lies outside of the range.

	    // Indices are stored 0 based and non inclusive of the end.
	    seqid_feature_in_genome    = (**gl_it).getSeqid().c_str();
	    start_of_feature_in_genome = (**gl_it).getStartIndex(); // 0-based.
	    end_of_feature_in_genome   = (**gl_it).getEndIndex();   // 0-based, exclusive end

	    if (global_VERBOSITY >= VERBOSITY_PROGRESS)  // Report progress
	    {
	      faststring t;
	      get_time_string(t);
	      cout << "PROGRESS: "
		   << " Working with gene: " << gene_name  << " " 
		                             << feature_number << "/" << features_this_gene
		   << " found in genome at "
		   << seqid_feature_in_genome << " "
		// This function gives the 0 based index. We want to print the 1 based index.
		   << start_of_feature_in_genome+1 << " " << end_of_feature_in_genome 
		   << endl << flush;
	    }

	    // Find sequence with short_name equal to seqid_feature_in_genome.
	    genome_seq_with_feature = genome_seqs.get_seq_by_name(seqid_feature_in_genome);
	    if (genome_seq_with_feature == NULL)
	    {
	      cerr << "ERROR: gff and genome file are not compatible: Sequence name: \""
		   << seqid_feature_in_genome
		   << "\" does not exist in the genome file."  << endl;

	      cerr << "It could be that some pipeline renamed the sequences in your genome file." << endl;

	      exit(-2);
	    }
	    else
	    {
	      // getPartialSeq: Takes 0 based start and length:
	      feature_sequence = genome_seq_with_feature->getPartialSeq(start_of_feature_in_genome, end_of_feature_in_genome - start_of_feature_in_genome);

	      if (global_VERBOSITY >= VERBOSITY_DEBUG)
	      {
		cout << "Pair:" << endl;
		cout << ">T" << endl;
		cout << transcript_sequence << endl;
		cout << ">GE" << endl;
		cout << feature_sequence << endl;
	      }
	    }

	    if ((int)feature_sequence.length() < min_alignment_length)
	    {
	      faststring tmp;

	      // Notice:
	      if (alignment_cutting)
		tmp = gff_feature_name + " feature alignment";
	      else
		tmp = "alignment";

	      if (global_VERBOSITY >= VERBOSITY_WARNINGS)
	      {
		cerr << "NOTICE: Skipping " << tmp << " in " << gene_name  << " with number " << feature_number << "/" << features_this_gene << " in " << gene_name
		     << " since it is too short. Alignment length: "  << feature_sequence.length() << ". Minimum length: " << min_alignment_length << endl;
	      }
	    }
	    else
	    {
	      // Determine feature range in transcriptome MSA:

	      // 	  transcript_sequence = "AACGTGTCGTAA";
	      // 	  feature_sequence       =   "CGTCGTCGT";
	      // 	  transcript_sequence =   "CGTCGTCGT";
	      // 	  feature_sequence       = "AACGTGTCGTAA";

	      int s1ins2s_best;
	      int s1ins2e_best;
	      int s2ins1s_best;
	      int s2ins1e_best;
	      int alignment_score_best;
	      int alignment_length_best;

	      find_coords(DNAmat, transcript_sequence, feature_sequence,
			  s1ins2s_best, s1ins2e_best, s2ins1s_best, s2ins1e_best,
			  alignment_score_best, alignment_length_best, global_VERBOSITY);

	      //**********************************************************
	      // Coordinates are 0 based, end is non inclusive.
	      //**********************************************************

	      // Start coordinate in transcriptome alignment in which the feature starts (0 based).
	      // If the number is negative, then the feature starts before the transcript starts.
	      // Basically, part of the transcript is missing. The effective starting point in our
	      // transcript is: 0.  +++
	      int in_transcript_start = macromax2(s2ins1s_best, 0);

	      // End coordinate in transcript alignment in which the feature starts. End is non inclusive.
	      int in_transcript_end   = s2ins1e_best;
	      if (in_transcript_end <= 0)
	      {
		//
		in_transcript_end = (int) transcript_sequence.size();
	      }

	      int feature_in_transcript_length = in_transcript_end - in_transcript_start;

	      if (global_VERBOSITY >= VERBOSITY_DEBUG)
	      {
		cout << "in_transcript_start:               " << in_transcript_start+1 << endl;
		cout << "in_transcript_end:                 " << in_transcript_end     << endl;
		cout << "feature length:                    " <<  feature_sequence.size()      << endl;
		cout << "feature_in_transcript_length:      " <<  feature_in_transcript_length << endl;
		cout << "Score:                             " <<  alignment_score_best      << endl;
		cout << "Score/feature-length               " << (double)alignment_score_best/feature_sequence.size()      << endl;
		cout << "Score/feature_in_transcript_length " << (double)alignment_score_best/feature_in_transcript_length << endl;
	      }

	      if ( (int)feature_sequence.size() > feature_in_transcript_length)
	      {
		if (global_VERBOSITY >= VERBOSITY_WARNINGS)
		  cerr << " ALERT: Length of excised genomic feature is longer than the transcipt, i.e. it looks as if the transcipt is incomplete.\n"
		          " This can happen. Gene name:      " << gene_name
		       << " Feature number: "                  << feature_number   
		       << " Lengths are:    "                  << feature_in_transcript_length << "/" << feature_sequence.size()
		       << endl;
	      }

	      if ( feature_in_transcript_length < min_alignment_length)
	      {
		// Notice:
		if (global_VERBOSITY >= VERBOSITY_WARNINGS)
		{
		  cerr << "ALERT: Skipping " << gff_feature_name << " feature since its overlap with the transcript is too short. Length of overlap: "
		       << feature_in_transcript_length << " Minimum: " << min_alignment_length << endl;
		}
	      }
	      else // We will now deal with this feature here:  ####
		// This is the block in which a lot of memory is allocated (with new) that needs to be freed when leaving the block.
	      {
		// We extract the sequence alignment in the range: in_transcript_start, in_transcript_end
		feature_transcript_alignment = new CSequences2(*transcript_alignment,
							       in_transcript_start,
							       in_transcript_end );
		if (transcript_seq_index != -1u && remove_reference_species_from_alignments)
		  feature_transcript_alignment->remove(transcript_seq_index);

		//	      feature_transcript_alignment_length = feature_transcript_alignment->GetPosNum();
		++DEBUG_COUNT_NEW_CALLS_FEATURE_ALIGNMENT;

		// Remove gap only positions:
		feature_transcript_alignment_rmgo = new CSequences2(CSequence_Mol::dna, feature_transcript_alignment->get_ambiguity_character());
		feature_transcript_alignment->alignment_without_gap_only_positions(*feature_transcript_alignment_rmgo);
		feature_transcript_alignment_rmgo->remove_gap_or_ambig_only_sequences();

		int feature_in_transcript_length_rmgo = feature_transcript_alignment_rmgo->GetPosNum();

		if (feature_in_transcript_length_rmgo < min_alignment_length)
		{
		  // Notice:
		  if (global_VERBOSITY >= VERBOSITY_WARNINGS)
		  {
		    cerr << "ALERT: Skipping " << gff_feature_name << " feature since the alignment is tool short after removing gap only positions. Length of overlap: "
			 << feature_in_transcript_length_rmgo << " Minimum: " << min_alignment_length << endl;
		  }
		}
		else
		{
		  ///////////////
		  // Exporting transcript alignment before removing lines containig gaps and Ns.
		  ///////////////
		  if (true) // if export is requested
		  {
		    faststring output_file_name = output_directory + "/" + gene_name + "_" + faststring(feature_number) + "_" + faststring(features_this_gene) + "_all.fas";

		    ofstream os(output_file_name.c_str());
		    if (os.fail())
		    {
		      cerr << "ERROR: Could not open file " << output_file_name << " for writing." << endl;
		      exit(-9);
		    }
		    feature_transcript_alignment_rmgo->ExportSequences(os, 'f', 100000000);
		    os.close();
		  }

		  // Debug output:
		  //	      feature_transcript_alignment_rmgo->ExportSequences(cout, 'f', 10000);

		  /// We have an alignment which we can work with:
		
		  bool found_bait_region;
		
		  found_bait_region = handle_this_alignment(feature_transcript_alignment_rmgo,
							    taxonNamesDict,
							    required_groups_of_taxa,
							    bait_length,
							    number_of_baits_in_region,
							    bait_offset,
							    gene_name,
							    feature_number,
							    features_this_gene,
							    cluster_threshold,
							    os_result_file,
							    sequence_name_field_delimiter,
							    sequence_name_taxon_field_number,
							    center_computation_mode, *it);

		  if (found_bait_region)
		    ++features_with_at_least_one_bait_regions;

		  if (feature_transcript_alignment)
		  {
		    delete feature_transcript_alignment;
		    ++DEBUG_COUNT_DELETE_CALLS_FEATURE_ALIGNMENT;
		    feature_transcript_alignment = NULL;
		  }
		  if (feature_transcript_alignment_rmgo)
		  {
		    delete feature_transcript_alignment_rmgo;
		    feature_transcript_alignment_rmgo = NULL;
		  }
		} // END ELSE if (feature_in_transcript_length_rmgo < min_alignment_length)
	      } // END ELSE if ( feature_in_transcript_length < min_alignment_length)
	    } // END ELSE: if ((int)feature_sequence.length() < min_alignment_length)    
	    ++gl_it;
	  }   // END while (gl_it != gl_it_end)  - For each feature in gene

	  // Check whether we have found something:
	  if (features_with_at_least_one_bait_regions)
	  {
	    if (global_VERBOSITY >= VERBOSITY_MOREPROGRESS)
	    {
	      cout << "NOTE: Alignment contains " << features_this_gene << " features, of which " << features_with_at_least_one_bait_regions
		   << " have at least one valid bait region." << endl;
	    }
	  }
	  else
	  {
	    if (global_VERBOSITY >= VERBOSITY_WARNINGS)
	    {
	      cerr << "ALERT: Alignment contains " << features_this_gene << " features but no single valid bait region. (Details: see above.)" 
		   << endl;
	    }
	  }

	  // Do not delete the transcript_seq pointer, since it is part of transcript_alignment.
	} // ELSE of if (transcript_seq == NULL)
      } // END  if (get_file_size(full_alignment_file_name) > 0)
      else
      {
	if (global_VERBOSITY >= VERBOSITY_WARNINGS)
	  cerr << "WARNING: Skipping file due to zero file size: " << full_alignment_file_name << endl;
      }

      // Delete alignment file:
      delete transcript_alignment;
      transcript_alignment = NULL;

      ++DEBUG_COUNT_DELETE_ALIGNMENT_FROM_FILE;
      if (global_VERBOSITY >= 1000)
      {
	print_DEBUG_COUNT(cout, "After calling all deletes of this transcript file:");
      }
      ++it;
    } // END while (it != it_end) -- For all transcript files

    if (global_VERBOSITY >= VERBOSITY_PROGRESS)
      cout << "PROGRESS: Bait-Fisher handled " << count_transcripts << " alignment files." << endl;

    if (global_VERBOSITY >= VERBOSITY_WARNINGS)
    {
      if (count_transcripts_no_gff_entries > 0)
      {
	cerr << "WARNING: For some alignment files the gene name could not be found in the gff file. These alignments had to be skipped.\n";
	cerr << "         Number of alignments affected: " << count_transcripts_no_gff_entries << ". See warnings above."  << endl;
      }
      if ((count_transcripts_no_gff_entries_for_feature-count_transcripts_no_gff_entries) > 0)
      {
	cerr << "WARNING: For some alignment files we found entries for the specified gene name in the gff file. But after filtering for the specified\n"
	        "         feature name no features where left. These alignment files had to be skipped.\n";
	cerr << "         Number of alignments affected: " << count_transcripts_no_gff_entries_for_feature-count_transcripts_no_gff_entries << ". See warnings above."  << endl;
      }
    }
  } // END if (alignment_cutting)

    //****************************************************************
    // END Alignment cutting mode
    //****************************************************************

  else // No alignment cutting
  {
    // Read fasta sequnces one by one:
    list<faststring>::iterator it, it_end;
    it     = fasta_names.begin();
    it_end = fasta_names.end();

    unsigned        Num_alignments    = fasta_names.size();
    unsigned        count_alignments = 0;

    CSequence_Mol   *alignment_seq_dummy;      // Only used as parameter. Not used in this case.
    unsigned        alignment_seq_index_dummy; // Only used as parameter. Not used in this case.
    faststring      full_alignment_file_name;
    faststring      sequence_header_target;
    //    faststring      alignment_name;

    faststring      genomic_seq_name;

    //  unsigned        overall_counter  = 0;
    //  double          mean_match_score = DNAmat.get_mean_matchscore_of_first_N_symbols(4);

    CSequences2     *alignment_from_file = NULL;

    CTaxonNamesDictionary taxonNamesDict(dname_alignments, fasta_names, sequence_name_taxon_field_number, sequence_name_field_delimiter);

    int error_required_taxa;
    CRequiredTaxa required_groups_of_taxa(string_required_taxa.c_str(), taxonNamesDict, error_required_taxa);

    if (error_required_taxa > 0)
    {
      if (error_required_taxa == 3)
      {
	cerr << "ERROR: When parsing the required taxa string. At least one of the taxon names did not occur in any of the sequence files that are analysed. Please look for earlier warning messages." << endl;
	cerr << "Exiting." << endl;
	exit(-3);
      }
    }

    { // Debug output:
      if (global_VERBOSITY > 1000)
      {
	cout << endl;

	taxonNamesDict.print_taxon_names_to_index(cout);

	cout << "Required groups: From each of the following sets of taxa, we require at least one to be present in the data set." << endl;
	required_groups_of_taxa.print(cout, 1);
	cout << endl;
      }
    }

    // For all alignment files:
    while (it != it_end)
    {
      ++count_alignments;
      full_alignment_file_name = dname_alignments + "/" + *it;

      if (get_file_size(full_alignment_file_name) > 0)
      {
	// Here is the main action: We handle this gene:

	// For this gene file, read the alignment for the genome we are concerned with
	if (global_VERBOSITY >= VERBOSITY_PROGRESS)
	{
	  faststring t;
	  get_time_string(t);
	  cout << "PROGRESS: " << t << " Reading alignment file: ("
	       << count_alignments << "/" << Num_alignments
	       << ") " << full_alignment_file_name << endl;
	}

	read_alignment_sequence_of_species2(full_alignment_file_name,
					    faststring(),
					    alignment_seq_dummy,
					    alignment_from_file,
					    alignment_seq_index_dummy,
					    sequence_name_field_delimiter,
					    sequence_name_taxon_field_number);

	if (alignment_from_file == NULL)
	{
	  cerr << "ERROR: Could not read alignment file: "
	       << full_alignment_file_name << endl;
	  exit(-4);
	}

	int alignment_length = alignment_from_file->GetPosNum();

	if (alignment_length < min_alignment_length)
	{
	  // Notice:
	  if (global_VERBOSITY >= VERBOSITY_WARNINGS)
	  {
	    cerr << "NOTICE: Skipping alignment since it is too short. Alignment length: "
		 << alignment_length << " Minimum: " << min_alignment_length << endl;
	  }
	}
	else
	{
	  // We will now deal with this alignment here:  ####
	  // This is the block in which a lot of memory is allocated (with new) that needs to be freed when leaving the block.
	    /// We have an alignment which we can work with:

	  bool found_bait_region;

	  found_bait_region = handle_this_alignment(alignment_from_file,
				  taxonNamesDict,
				  required_groups_of_taxa,
				  bait_length,
				  number_of_baits_in_region,
				  bait_offset,
				  *it,  // alignment name
				  0,
				  0,
				  cluster_threshold,
				  os_result_file,
				  sequence_name_field_delimiter,
				  sequence_name_taxon_field_number,
				  center_computation_mode, *it);


	} // END ELSE if ( alignment_length < min_alignment_length)

	// Do not delete the alignment_seq pointer, since it is part of alignment_from_file.
      } // END (get_file_size(full_alignment_file_name) > 0)
      else
      {
	if (global_VERBOSITY >= VERBOSITY_WARNINGS) 
	  cerr << "WARNING: Skipping file due to zero file size: " << full_alignment_file_name << endl;
	  //	cout << "WARNING: Skipping file due to zero file size: " << full_alignment_file_name << endl;
      }

      // Delete alignment file:
      delete alignment_from_file;
      alignment_from_file = NULL;

      ++DEBUG_COUNT_DELETE_ALIGNMENT_FROM_FILE;

      if (global_VERBOSITY > 1000)
      {
	print_DEBUG_COUNT(cout, "After calling all deletes of this alignment file.:");
      }
      ++it;
    } // END while (it != it_end) -- For all alignment files

    if (global_VERBOSITY >= VERBOSITY_PROGRESS)
      cout << "Bait-Fisher handled " << count_alignments << " alignment files." << endl;
  } // END else // No alignment cutting

  os_result_file.close();

  cout << "Bait-Fisher finished without errors." << endl;

} // END main(int argc, char **argv)


// Memory leek detection:
//
// The function read_transcript_sequence_of_species2 reserves memory for:
// - transcript alignment. This is variable "transcript_alignment" in the main function.
//   This memory is released at the end of the while loop over all transcripts that are read.
//
// In the loop over all features we reserve memory for the
// - "feature_transcript_alignment"
//   This memory is released at the the end of the block that handles features that are sufficienlty long
//   for designing baits.
//
//
//
//
//
//
//
//
//

