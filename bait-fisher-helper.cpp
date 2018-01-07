/*  BaitFisher (version 1.2.8) a program for designing DNA target enrichment baits
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
#include "faststring2.h"
#include <iostream>
#include <list>

#include "CFile/CFile2_1.h"
#include "CSequences2.h"
#include "CSequence_Mol2_1.h"

#include "bait-fisher-helper.h"

using namespace std;

// Extract species (taxon) name from transcript name:
faststring extract_taxon_name(faststring sequence_name, char delim, unsigned short field_num)
{
  list<faststring> container;  char             delim_cstring[2];
  
  delim_cstring[0] = delim;
  delim_cstring[1] = '\0';

  split(container, sequence_name, delim_cstring);

  // If the field number (0 based) is not smaller the number of fields, we have a problem.
  if (!(field_num < container.size() ) )
  {
    cerr << "ERROR: In extract_taxon_name: Found unexpected number of fields in a sequence name: " << sequence_name << endl;
    cerr << "I tried to extract field number " << field_num << " as specified in the parameter file or given as default value." << endl;
    cerr << "But this field does not exist." << endl;
    exit(-3);
  }

  list<faststring>::iterator it = container.begin();
  unsigned short i;
  
  for (i=0; i<field_num; ++i)
    ++it;

  return *it;
}


// Extract gene name from transcript name:
faststring extract_gene_name(faststring sequence_name, char delim, unsigned short field_num)
{
  list<faststring> container;
  char             delim_cstring[2];
  
  

  delim_cstring[0] = delim;
  delim_cstring[1] = '\0';

  split(container, sequence_name, delim_cstring);

  // If the field number (0 based) is not smaller the number of fields, we have a problem.
  if (!(field_num < container.size() ) )
  {
    cerr << "ERROR: In extract_gene_name: Found unexpected number of fields in a sequence name: " << sequence_name << endl;
    cerr << "I tried to extract field number " << field_num << " as specified in the parameter file or given as default value." << endl;
    cerr << "But this field does not exist." << endl;
    exit(-3);
  }

  list<faststring>::iterator it = container.begin();
  unsigned short i;
  
  for (i=0; i<field_num; ++i)
    ++it;
  
#ifdef SHORTEN_GENENAME_AFTER_MINUS
    faststring pre, post;
    it->backwards_search_divide_at_break_at('-', delim, pre, post);
    return pre;
#else
    return *it;
#endif
}

// (1) Read alignment.
// (2) Provide pointer to specific sequence (the one for which the name is specified, see (3)) as well as to the full alignment.
// (3) If species_name_query is not empty, find the sequence with this name.
//     If species_name_query is an empty string, assign NULL to pseq.
// (4) Important: The CSequence_Mol and CSequences2 objects need to be deleted by the user of this function when not used any more.
// Version 2: Read full alignment and return both: full alignment and requested sequence.
void read_alignment_sequence_of_species2(faststring       alignment_file,
					 faststring       species_name_query,
					 CSequence_Mol*&  pseq,
					 CSequences2*&    pseqs,
					 unsigned&        pseq_index,
					 char             delim, 
					 unsigned short   field_num)
{
  CFile           infile;
  //  unsigned        count_seq=0;
  CSequences2*    tmp_pseqs;


//   fputs("Reading alignment file: ", stderr);
//   fputs(alignment_file.c_str(), stderr);
//   fputs(".\n", stderr);

  infile.ffopen(alignment_file.c_str());
  if ( infile.fail() )
  {
    fprintf(stderr, "Sorry! The input file \n   \"%s\"\ncould not be opened. I guess it does not exist.\n",
	    alignment_file.c_str());
    exit(-1);
  }

  CSequence_Mol::processing_flag pflag = CSequence_Mol::processing_flag(0);
  // Convert all nucleotides to upper case letters and ambig to Ns.
  pflag = CSequence_Mol::processing_flag(pflag | CSequence_Mol::convert_toupper | CSequence_Mol::convert_ambig2N);
  tmp_pseqs = new CSequences2(CSequence_Mol::dna);

  ++DEBUG_COUNT_NEW_ALIGNMENT_FROM_FILE;

  tmp_pseqs->read_from_Fasta_File(infile, pflag, 0, -1, false);
  infile.ffclose();


//   fputsf("Alignment file has been read successfully.\n", stderr);

  pseqs = tmp_pseqs;


  if (species_name_query.empty() )
  {
    pseq = NULL;
    pseq_index = -1u;    
  }
  else
  {
    // Find sequence with taxon name equal to species_name_query
    unsigned i, N;
    N = pseqs->GetTaxaNum();
  
    for (i=0; i<N; ++i)
    {
      const char*  seq_name = pseqs->get_Seq_Name(i);
    
      if (seq_name == NULL)
      {
	cerr << "INTERNAL ERROR: Could not retrieve sequence name for sequence of index: " 
	     << i << " in read_transcript_sequence_of_species2." << endl;
	exit(-5);
      }

      faststring species_name = extract_taxon_name(seq_name, delim, field_num);

      //     cout << "Index: " << i << " " << N << endl;
      //     cout << "Comparing: " << species_name << " " << species_name_query << endl;

      if (species_name == species_name_query)
	break;
    }
    if (i == N) // There is no sequence with species name equal to the one we search for.
    {
      pseq = NULL;
      pseq_index = -1u;
    }
    else
    {
      pseq = pseqs->get_seq_by_index(i);
      pseq_index = i;
    }
  }
}
