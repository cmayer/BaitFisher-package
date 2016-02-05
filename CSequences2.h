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
// TODO: Check all new whether they have a corresponding delete



//////////////////////////////////////////////////////////////////////
// CSequences.h: interface for the CSequences class.
//
//////////////////////////////////////////////////////////////////////

#ifndef  CSEQUENCES2_H
#define  CSEQUENCES2_H

#include <vector>
#include <map>
#include "faststring2.h"
#include "fast-realloc-vector.h"
#include "CSplit2.h"
#include "CSequence_Mol2_1.h"
#include <cassert>
#include <cmath>
#include <utility>
#include "basic-DNA-RNA-AA-routines.h"

#include <algorithm>

//#include "basic-DNA-RNA-AA-routines.h"

#define PrintMessage_cerr(s1)                fputs(s1, stderr)
#define ErrorMessage(s1)                     fputs(s1, stderr)
#define flush_cerr()                         fflush(stderr)

////////////////////////////////////////////////////////////////////
// Head file for a DNA or PROTEIN sequence collection. 
//
// Characteristics:
// - most member functions assume that sequences have identical lengths, though this is not true for all.
//   If sequences are not of equal length, it is recommended to revise the code in this head file first.
//
////////////////////////////////////////////////////////////////////



// Global Functions
//////////////////////////////////

// Criteria for reference taxa:
// - Original reference taxa have 2 pipes (3 sections), other sequences have 3 pipes (4 sections) in sequence names.
// - Hamstered reference taxa are those for which section 2 and 3 are identical in the first 4 characters.


/* inline bool is_1kite_hamstered_core_seq_id(faststring seq_id) */
/* { */
/* #ifdef DEBUG */
/*   std::cerr << "Called: CSequences2.h:is_1kite_hamstered_core_seq_id with parameter " << seq_id << std::endl; */
/* #endif */

/*   //  std::cerr << "Calling outdated function: is_1kite_hamstered_core_seq_id " << std::endl; */
/*   //  exit(-1); */

/*   std::vector<faststring> l; */

/*   split(l, seq_id, "|"); */

/*   unsigned num_sections = l.size(); */

/*   if (num_sections > 4 || num_sections < 3) */
/*   { */
/*     std::cerr << "Error: Non standard 1kite seq_id encountered: " << seq_id << std::endl; */
/*     exit(0); */
/*   } */

/* #ifdef DEBUG */
/*   std::cerr << "num_sections: " << num_sections << std::endl; */
/* #endif */

/*   if (num_sections == 3) */
/*     return true; */
/*   else */
/*     return false; */
/* } */



inline bool is_1kite_original_reference_taxon_seq_id(faststring seq_id)
{
#ifdef DEBUG
  std::cerr << "Called: CSequences2.h:is_1kite_original_reference_taxon_seq_id with parameter " << seq_id << std::endl;
#endif

  std::vector<faststring> l;

  split(l, seq_id, "|");

  unsigned num_sections = l.size();

  if (num_sections > 4 || num_sections < 3)
  {
    std::cerr << "Error: Non standard 1kite seq_id encountered: " << seq_id << std::endl;
    exit(0);
  }

#ifdef DEBUG
  std::cerr << "num_sections: " << num_sections << std::endl;
#endif

  if (num_sections == 3)
    return true;
  else
    return false;
}


inline bool is_1kite_hamstered_reference_seq_id(faststring seq_id)
{
  std::vector<faststring> l;

  split(l, seq_id, "|");

#ifdef DEBUG
  std::cerr << "Called: CPfamScanParser.h:is_1kite_hamstered_reference_seq_id with parameter " << seq_id << std::endl;
#endif

  if (l.size() != 4)
  {
    std::cerr << "Error: Non standard 1kite seq_id encountered: " << seq_id << std::endl;
    exit(0);
  }
  faststring part2 = l[1];

  part2.shorten_to_first_occurrence_of('_');
  faststring part3 = l[2];

  part2.shorten(4);
  part3.shorten(4);

#ifdef DEBUG
  std::cerr << "num_sections: " << l.size() << " " << part2 << " " << part3 << std::endl;
#endif

  return (part2 == part3);
}




//template<class T>
inline void add_or_count(std::map<faststring, unsigned> &m, faststring &x)
{
  std::map<faststring,unsigned>::iterator it;

  it = m.find(x);
  if (it == m.end() )
  {
    m[x] = 1;
  }
  else
  {
    ++it->second;
  }
}


//////////////////////////////////
// Data types
//////////////////////////////////

//////////////////////////////////
// Class CSequences
//////////////////////////////////

// Class for aligned sequences

class CSequences2
{
public:
  typedef char chartype;

private:
  CSequence_Mol::DataTypesEnum	datatype;            // DNA, Protein, or other
  unsigned                  taxaNum;
  char                      ambig_char;
  bool                      originalPosNumbers_supplied;

  unsigned	            posNum;

  std::vector<unsigned>         originalPosNumbers;  // In case we excluded positions from the alignment, we want to remember the original positions.
  std::vector<CSequence_Mol*>   seqData;             // The vector of pointers to sequences


  std::map<faststring, CSequence_Mol*> sn_map;       // Obtain pointer to the sequence by sequence name.

  void add_seq(CSequence_Mol* seq)
  {
    //    std::cerr << "Adding: " << seq->getName() << std::endl;

    seqData.push_back(seq);
    sn_map[seq->getName()] = seq;
  }

  bool equal_length_of_all_sequences()
  {
    size_t len;
    unsigned i;

    if (taxaNum==0)
      return true;

    len = seqData[0]->length();

    for (i=1; i<taxaNum; ++i)
    {
      if (len != seqData[i]->length())
	return false;
    }
    return true;
  }

  void determine_map_of_sequence_names()
  {
    if (sn_map.size() == 0)
    {
      int i, N=seqData.size();
      for(i=0; i<N; ++i)
      {
	sn_map[seqData[i]->getName()] = seqData[i];
      }
    }
  }


public:
  // Minimal constructor: Empty sequences object
 CSequences2(CSequence_Mol::DataTypesEnum Dt, char set_ambig_char='?'):
  datatype(Dt), taxaNum(0), ambig_char(set_ambig_char), originalPosNumbers_supplied(false),posNum(0)

  {}

  // Constructor for a set of empty sequences with names and length.
 CSequences2(CSequence_Mol::DataTypesEnum Dt, std::vector<faststring> names, unsigned len):
  datatype(Dt), taxaNum(names.size()), ambig_char('?'), originalPosNumbers_supplied(false),posNum(0)

  {
    seqData.reserve(taxaNum);
    unsigned        i;
    CSequence_Mol   *seq;

    for (i=0; i<taxaNum; ++i)
    {
      seq = new CSequence_Mol (CSequence_Mol::dna, names[i], len, ambig_char);
      add_seq(seq);
    }
  }


  ~CSequences2()
  {
    int i, n=taxaNum;

    for (i=0; i<n; ++i)
    {
      delete seqData[i];
    }

  }


  // This constructor can be used to
  //  - general copy constructor
  //  - extract a range of sites
  // Coordinates:
  //    pos1 must be the first column starting with 0.
  //    pos2 must be the index after the last column.
  //    pos2-pos1 must be the number of bases that are copied to this sequence.

 CSequences2(const CSequences2 &s, unsigned pos1 = 0, unsigned pos2 = -1u):
  datatype(s.datatype), taxaNum(s.taxaNum), ambig_char(s.ambig_char),
    originalPosNumbers_supplied(false), posNum(0)
  {
    unsigned        i, n=taxaNum;
    CSequence_Mol   *seq;

    for (i=0; i<n; ++i)
    {
      seq = new CSequence_Mol ( *(s.seqData[i]), pos1, pos2);
      add_seq(seq);
    }
    if (i != n)
    {
      std::cerr << "Critical error in CSequences2 constructor: taxaNum and number of sequences found disagree." << std::endl;
    }
    
    if (n > 0)
      posNum = seq->length();
  }
 
  CSequence_Mol* get_seq_by_name(faststring name)
  {
    // Test code:
    {
      //      print_DEBUG(cerr, 0);
    }

    std::map<faststring, CSequence_Mol*>::iterator find_it = sn_map.find(name);
    if (find_it != sn_map.end() )
      return find_it->second;
    else
      return NULL;
  }

  CSequence_Mol* get_seq_by_index(unsigned id)
  {
    if (id >= seqData.size())
    {
      return NULL;
    }
    else
    {
      return seqData[id];
    }
  }

  // Depreciated:
  CSequence_Mol* get_seq(unsigned id)
  {
    return get_seq_by_index(id);
  }

  const char* get_Seq_Data(unsigned id)
  {
    if (id >= seqData.size())
    {
      return NULL;
    }
    else
    {
      return seqData[id]->getSeqStr();
    }
  }

 const char* get_Seq_Name(unsigned id)
  {
    if (id >= seqData.size())
    {    
      return NULL;
    }
    else
    {
      return seqData[id]->getName();
    }
  }
  
  unsigned         GetTaxaNum() { return taxaNum;}
  unsigned         GetPosNum()  { return posNum;}

  CSequence_Mol::DataTypesEnum    get_datatype()
  {
    return datatype;
  }

  chartype         GetChar(unsigned TaxaIndex,
			   unsigned PosIndex) const
  {
    assert(TaxaIndex < taxaNum);
    assert(PosIndex  < posNum);
    return seqData[TaxaIndex]->get_pos(PosIndex);
  }
  
   void print_DEBUG(std::ostream &os, unsigned flag=1)
  {
    if (flag & 1) // scalars:
    {
      os << "Debug output CSequences2 object, flag==0" << std::endl;
      os << "Data type:  " << (int)datatype<< std::endl;
      os << "taxaNum:    " << taxaNum      << std::endl;
      os << "posNum:     " << posNum       << std::endl;
      os << "ambig_char: " << ambig_char   << std::endl;
      os << "originalPosNumbers_supplied " << (int)originalPosNumbers_supplied << std::endl;

      os << "size of originalPosNumbers vector: " << originalPosNumbers.size() << std::endl;
      os << "size of seqData                  : " << seqData.size()            << std::endl;
      os << "size of sn_map                   : " << sn_map.size()             << std::endl;
    }

    if (flag & 2)
    {
      int i, n=seqData.size();

      os << "Data types of sequences:\n";
      for (i=0; i<n; ++i)
      {
	  os << i << ": " << seqData[i]->get_datatype() << " " << seqData[i]->type_as_string() << std::endl;
      }
    }
  }



/*    void             SetChar(unsigned TaxaIndex, unsigned PosIndex, chartype mychar) */
/*    { */
/*      assert(TaxaIndex < taxaNum); */
/*      assert(PosIndex < posNum); */
/*      seqData[TaxaIndex]->set_pos(PosIndex, mychar); */
/*    } */

   void SetOriginalPosNumber(unsigned index, unsigned theOriginalPosNumber)
   {
     originalPosNumbers_supplied = true;
     assert(index < posNum);
     originalPosNumbers[index] = theOriginalPosNumber;
   }

   unsigned GetOriginalPosNumber(unsigned index) const
   {
     assert(index < posNum);
     if (originalPosNumbers_supplied)
       return originalPosNumbers[index];
     else
       return -1u;
   }

   // Moves taxon with index index to top, preseving the order of all other taxa. 
   // Of course this changes all indices with the only exception that this taxon is
   // already at the top of the vector.
   // The index is of course 0 based.
   void reorder_move_seq_to_top(unsigned index)
   {
     if (index < taxaNum && index > 0)
     {
       seqData.insert(seqData.begin(), seqData[index]);
       seqData.erase(seqData.begin()+index+1); // We have to add 1 since we have one additional entry at the beginning
     }
   }

   // TODO: Not very efficient!!
   // Return true if taxon_name has been found, false otherwise
   bool reorder_move_seq_to_top(faststring &taxon_name)
   {
     unsigned i;
     for (i=0; i< taxaNum; ++i)
     {
       faststring iname = seqData[i]->getFullName();
       if (seqData[i]->getFullName() == taxon_name)
       {
	 std::cerr << "Move to top: " << i << std::endl;
	 reorder_move_seq_to_top(i);
	 return true;
       }
     }
     return false;
   }

   // Reorder the sequences such the those in the given vector are at the top
   // in the order given by the vector.
   // TODO: Not very efficient! - Very simple and thus secure implementation.
   bool reorder_sort_by(std::vector<faststring> &names)
   {
     int  i,n=names.size();
     bool success = true;

     for (i=n-1; i>=0; --i)
     {
       success = success & reorder_move_seq_to_top(names[i]);
     }
     return success;
   }

   bool get_Originalposnumbers_Supplied()
   {
     return originalPosNumbers_supplied;
   }

   char get_ambiguity_character()
   {
     return ambig_char;
   }


   // Returns the maximum sequence, length
   void get_sequence_names(std::vector<faststring> &snames)
   {
     snames.clear();
     unsigned  i;
     //     unsigned  max=0;

     for (i=0; i < taxaNum; ++i)
     {
       snames.push_back(seqData[i]->getFullName());
/*        if (snames[i].size() > max) */
/*        { */
/* 	 max = snames[i].size(); */
/*        } */
     }
     //     return max;
   }

   // Not yet implemented
/*    void unique_full_names(std::vector<faststring> &snames) */
/*    { */
     
/*    } */

   // Not yet implemented
/*    void unique_short_names(std::vector<faststring> &snames) */
/*    { */
     
/*    } */



   void unique_phylip_names(std::vector<faststring> &snames)
   {
     snames.clear();
     std::map<faststring, unsigned>  m_snames;
     std::map<faststring, unsigned>  m_shortend_names;
     std::map<faststring, unsigned>::iterator  f_it, f_it2;
     
     unsigned                             digits = (unsigned)log10(taxaNum) + 1;
     unsigned                             i;
     
     for (i=0; i < taxaNum; ++i)
     {
       // push_back first sequence phylip sequence name. This is already trimmed to 10 characters.
       snames.push_back(seqData[i]->getPhylipName());
       // Search for this name in m_snames map.
       f_it = m_snames.find(snames[i]);

       // If we find it it is not unique
       if (f_it != m_snames.end() ) // This name is not unique
       {
	 ++f_it->second; // We count this occurence
       }
       else // If we do not find it, we insert it and set the counter to 1
       {
	 m_snames.insert(std::make_pair(snames[i], 1));
       }
     }
     // We have build the map that count the number of occurences - now we need to create unique names:

     // We will shorten the names to add a number.
     // When shortening names, previously unique names might become non-unique so we have to take care of this first:

     f_it = m_snames.begin();
     while (f_it != m_snames.end() )
     {
       //       std::cout << f_it->first << " " << f_it->second << std::endl;

       if (f_it->second > 1) // If the 10 character sequence name occurs more than one, we have to number them:
       {
	 faststring temp = f_it->first;
	 temp.shorten(9-digits);
	 
	 f_it2 =  m_shortend_names.find(temp);
	 if (f_it2 == m_shortend_names.end() ) // We already have this name in the map:
	 {
	   m_shortend_names[temp] += f_it->second;
	 }
	 else // Create this entry:
	 {
	    m_shortend_names[temp] = f_it->second;
	 }
       }
       ++f_it;
     }

     // Now we have two maps:
     // The first has entries for all names
     // The second has entries for all names that need to be shortend.

     // For all names that need to be shortend:
     f_it = m_shortend_names.begin();
     while (f_it != m_shortend_names.end() )
     {
       unsigned internal_counter = 1;
       unsigned j;
       faststring temp;
       for (j=0; j<taxaNum; ++j)
       {
	 temp = snames[j];
	 if ( m_snames[temp] != 1 )
	 {
	   temp.shorten(9-digits);
	 }
	 if (temp == f_it->first)
	 {  
	   faststring nn = faststring(internal_counter, '0', digits);
	   //	   nn.replace_char(' ', '0');

	   snames[j]  = temp;
	   snames[j].append('_');

	   snames[j].append(nn);
	   ++internal_counter;
	 }
       }
       ++f_it;
     }
   }


   unsigned PairwiseSequenceSymbolMatches(unsigned taxon1, unsigned taxon2,
					  char     numDistCharacters,
					  const signed char* pos_vec,
					  signed char* match_vec) const
   {
     unsigned  count    = 0;
     unsigned  pos;
     
     for (pos = 0; pos < posNum; ++pos)
     {
       if (pos_vec[pos]  &&
	   seqData[taxon1]->get_pos(pos) == seqData[taxon2]->get_pos(pos) &&
	   seqData[taxon1]->get_pos(pos) < numDistCharacters)
       {
	 ++count;
	 match_vec[pos] = 1;
       }
       else
	 match_vec[pos] = 0;
     }
     return count;
   }
   
 
   unsigned PairwiseSequenceSymbolMatches(
					  unsigned taxon1,
					  const faststring & ref_seq,
					  char               numDistCharacters,
					  const signed char* pos_vec,
					  signed char* match_vec) const
  

   // Ist hier alles OK?????????????????????????????????????????
   {
     unsigned  count    = 0;
     unsigned  pos;
     
     for (pos = 0; pos < posNum; ++pos)
     {
       if (pos_vec[pos]  &&  seqData[taxon1]->get_pos(pos) == ref_seq[pos]
	              &&  ref_seq[pos] < numDistCharacters)
       {
	 ++count;
	 match_vec[pos] = 1;
       }
       else
	 match_vec[pos] = 0;
     }
     return count;
   }
   

   void     ConsensusSequence(faststring&           conSeq,
			      const CSplit&         setOfTaxa,
			      char                  numDistCharacters,
			      unsigned              numSymbols,
			      double                consensusThreshold )
   {
     unsigned          pos;
     unsigned          taxon, maxindex, equalindex;
     unsigned          taxaCount = 0;
     unsigned          consensusSymMinimum;
     char              i;


     std::vector<unsigned>  counterSymbolsVec(numSymbols,0);
     
     for (pos=0; pos < posNum; ++pos)
     {
       // Initialise variable
       taxaCount = 0;
       for (i=0; i < numDistCharacters; ++i)
	 counterSymbolsVec[i] = 0;
       
       // Count number of occuring symbols for this position over all taxa
       for (taxon = 0; taxon < taxaNum; ++taxon)
       {
	 if (setOfTaxa.test(taxon))
	 {
	   ++taxaCount;
	   ++counterSymbolsVec[seqData[taxon]->get_pos(pos)];
	 }
       }
       
       consensusSymMinimum = (unsigned) std::ceil(consensusThreshold * taxaCount);
       
       maxindex = 0;
       equalindex = numDistCharacters;
       for (i = 1; i < numDistCharacters; ++i)
       {
	 if (counterSymbolsVec[i] >= counterSymbolsVec[maxindex])
	 {
	   if (counterSymbolsVec[i] == counterSymbolsVec[maxindex])
	     equalindex = i;
	   maxindex = i;
	 }
       }
       if (counterSymbolsVec[maxindex] >= consensusSymMinimum &&
	   maxindex != equalindex)
	 conSeq[pos] = maxindex;
       else // Default value, in case consensus cannot be determined
	 conSeq[pos] = numDistCharacters;
     }
   }
   

   void GetSymbolFrequenciesAtPosition(unsigned pos, const CSplit &taxaSet,
				       unsigned numSymbols,
				       std::vector<unsigned>& frequencies)
   {
     unsigned i;
     //  chartype sym;
     
     //  assert(frequencies.size() <= numDistCharacters);
     
     for (i=0; i<numSymbols; ++i)
       frequencies[i]=0;
     
     for (i=0; i<taxaNum; ++i)
     {
       if (taxaSet.test(i))
       {
	 ++frequencies[seqData[i]->get_pos(pos)];
       }
     }
   }

   int read_from_Fasta_File(CFile& infile,
			    CSequence_Mol::processing_flag pflag,
			    unsigned first_seq,
			    unsigned last_seq,
			    bool report_status_to_cerr=true
			    )
   {
     CSequence_Mol   *seq;
     unsigned        count_seq = 0;
     unsigned        count_seqs_stored = 0;


     // TODO
     (void) report_status_to_cerr;

     // Remove all sequences from this object - we want to read a new data set.
     clear(datatype);

     while (!infile.eof())
     {
       seq = new CSequence_Mol (datatype);
       ++count_seq;

       seq->readFastaSequence_generic(infile, pflag);

       // PrintMessage_cerr(seq.getName());PrintMessage_cerr("-\n");PrintMessage_cerr(seq.getFullName());PrintMessage_cerr("-\n");PrintMessage_cerr(seq.getDescription());PrintMessage_cerr("-\n");
       if ( infile.fail() && !infile.eof() )
       {
	 PrintMessage_cerr("\n\n");
	 faststring errormsg;
	 errormsg   =  "An error occurred while reading the input file. It might not be a valid fasta file.\n";
	 errormsg  +=  "File position: line ";
	 errormsg  +=  faststring(infile.line());
	 ErrorMessage(errormsg.c_str());
	 PrintMessage_cerr("\n");
	 flush_cerr();
	 return -25;
       }

       if ( count_seq < first_seq || count_seq > last_seq )
       {
	 if (report_status_to_cerr)
	 {
	   PrintMessage_cerr("Skipping sequence ");
	   PrintMessage_cerr(faststring(count_seq).c_str());
	   PrintMessage_cerr(": ");
	   PrintMessage_cerr(seq->getName());
	   PrintMessage_cerr("\n");
	   flush_cerr();
	 }
	 continue;
       }

       if (report_status_to_cerr)
       {
	 PrintMessage_cerr("Processing sequence ");
	 PrintMessage_cerr(faststring(count_seq).c_str());
	 PrintMessage_cerr(": ");
	 PrintMessage_cerr(seq->getName());
	 PrintMessage_cerr("\n");
	 flush_cerr();
       }

       add_seq(seq);
       ++count_seqs_stored;
     } // End while

     // Check that all sequences have the same length - we require that they have the same length
     if (!equal_length_of_all_sequences() )
       return -1;

     // Check equal data type and ambig character:

     taxaNum = count_seqs_stored;

     if (taxaNum > 0)
     {
       posNum  = seqData[0]->length();
       ambig_char = seqData[0]->get_ambiguity_character();
     }

     return 0;
   }  // End this element function.



   // Can read sequential and interleaved phylip files
   // Cannot read multi phylip files.

   // format_flag: 0 strict phylip standard. Sequence names are exactly 10 characters long.
   //              1 relaxed sequence name length. Sequence names end with first space and can have anly length.
   //                Spaces are not allowed in sequence names.

   int read_from_Phylip_File(CFile& infile, unsigned format_flag = 1)
   {
     CSequence_Mol   *seq;
     //     unsigned        count_seq;
     //     unsigned        global_sequence_from=0, global_sequence_to=-1u;
     //     bool            global_gaps_to_Ns=false;
     unsigned        all_lines_N;
     unsigned        curr_taxon;

     char read_mode;

     // The allowed symbols from the phylip manual:


     char ch;
     faststring line;
     std::vector<faststring> fsvec;
     faststring *p;

     // Clear the oject and initialize it with the supplied data type
     clear(datatype); // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX

     ch = infile.peek();

     // If the first char is a white space or a digit, we expect this
     // is the line with the number of taxa and seq. positions.
     if ( isspace(ch) || isdigit(ch) )
     {
       infile.getline(line);
       split(fsvec, line);

       if (fsvec.size()==2 && fsvec[0].isAnUnsigned() && fsvec[1].isAnUnsigned() )
       {
	 taxaNum = fsvec[0].ToUnsigned();
	 posNum  = fsvec[1].ToUnsigned();

/* 	 std::cout << "Number of taxa     " << taxaNum << std::endl; */
/* 	 std::cout << "Number of residues " << posNum  << std::endl; */
       }
     }

     //     std::vector<faststring*> sqs;
     std::vector<faststring*> all_lines;

/*      if (taxaNum == 0 && posNum == 0) */
/*      { */
/*        faststring *pp = new faststring; */
/*        *pp = line; */
/*        all_lines.push_back(pp); */
/*      } */

     // We read the complete file into a vector of strings.
     p = new faststring();
     infile.getline(*p);
     p->removeSpacesFront();
     p->removeSpacesBack();

     while (!infile.fail() )
     {
       all_lines.push_back(p);
       p = new faststring();
       infile.getline(*p);
       p->removeSpacesFront();
       p->removeSpacesBack();
     }

     // The problem is: The phylip format is very difficult
     // to parse. Main problem is to distinguish between
     // sequential and interleaved format.
     // The interleaved format is distinguished as follows:
     // After each block there must be a blank line.
     // So if we find blank lines which are not the last line
     // the format is interleaved. Otherwise it is sequential.

     // Another problem is relaxed versus strict format.
     // In the related format the sequence name and the sequence are separated by white spaces.
     // In the strict format the sequence name is less or sequal to 10 character.
     // Example: Name123456A MILLV
     // Can be interpreted as sequence AMILLV or as sequence MILLV with the two different names.
     // In the strit format the following would be legal:
     // NAME123456MILLVWSAL LADDPKHIMV
     // The sequence name would be NAME123456 the rest is the sequence.
     // Note that the format allows spaces in the sequence.
     // The two versions cannot be distinguished unambiguously!!
     // Pragmatic approach: In the presence of one white space island, we assume the relaxed format.
     // Otherwise the strict format.

     // Since empty lines are allowed at the end of the file in all formats we remove them.
     // They cannot be indicative of the interleaved format.

     unsigned i;
     all_lines_N = all_lines.size();

     //     std::cerr << all_lines[all_lines_N-1]->length() << std::endl;
          
     while (all_lines[all_lines_N-1]->length() == 0)
     {
       --all_lines_N;
       delete all_lines[all_lines_N];
       std::vector<faststring*>::iterator it = all_lines.end();
       --it;
       all_lines.erase(it);
     }

     unsigned number_blank_lines=0;

     for (i=0; i<all_lines_N; ++i)
     {
       if (all_lines[i]->length() == 0)
       {
	 ++number_blank_lines;
       }
     }
     
     if (all_lines_N == 0 || number_blank_lines == all_lines_N) // Nothing to do???
       return 0;
  
     if (number_blank_lines == 0) // sequential format
     {
       read_mode = 's';
       // Example:
       //a          ATGATTCAAC CTCAGACCCT TTTAAATGTA GCAGATAACA GTGGAGCTCG 
       //           AAAATTGATG 
       //b          ATGATTCAAC CTCAGACCCA TTTAAATGTA GCGGATAACA GCGGGGCTCG 
       //           AGAATTGATG 

       // This can only be interpreted, if the number
       // of taxa and positions has been specified in
       // the first line.
       // What makes it more complicated
       // is the fact, that the phylip format tolerates
       // numbers and other stuff at the beginning of pure sequence lines.

       faststring  seq_taxon_name;
       faststring  seq_in_one_line;
       //       unsigned    current_length = 0;
       unsigned    curr_line_num=0;

       faststring  *currLine;


       if (taxaNum == 0 || posNum == 0)
       {
	 // delete all faststrings in all_lines
	 for (unsigned ind=0; ind < all_lines.size(); ++ind)
	   delete all_lines[ind];
	 return -3;
       }


       // For all taxa
       //
       for (curr_taxon = 0; curr_taxon < taxaNum; ++curr_taxon)
       {
	 if (curr_line_num >= all_lines_N)
	 {
	   std::cerr << "Parse error while reading the input file:\nUnexpected end of file in line: " << curr_line_num <<  std::endl;
	   std::cerr << "Trying to read " << taxaNum << " taxa but found only " << curr_taxon << "!" <<  std::endl;

	   exit(-45);
	 }

	 currLine       = all_lines[curr_line_num];
	 const char* pos = currLine->c_str();
	 //	 currLineLen     = currLine->length();
	 seq_taxon_name.clear();
	 seq_in_one_line.clear();
	 
	 // Copy the sequence name
	 if (format_flag == 0)
	 {
	   for (i=0; i<10 && *pos != '\0'; ++i, ++pos)
	   {
	     seq_taxon_name.push_back(*pos);
	   }
	   if (i<10) // Fill name with spaces.
	     seq_taxon_name.append(10-i, ' ');
	 }
	 else
	 {
	   for (;*pos != '\0' && *pos != ' '; ++pos)
	   {
	     seq_taxon_name.push_back(*pos);
	   }
	 }

	 // Skip all spaces:
	 while (*pos == ' ' && *pos != '\0')
	   ++pos;

	 // Read data after sequence name:
	 if (*pos)
	 {
	   // Copy the sequence data:
	   {
	     // Find first character that shall be copied.
	     // pos should point to the first char of the sequence.
	     char c;

	     // Skip all non-allowed symbols - since we skiped spaces, there should be none.
	     c = *pos;
	     while (!is_phylip_aa_dna_symbol(c) && !is_allowed_non_ABC(c) )
	     {
	       ++pos;
	       c = *pos;
	     }
	   }
	   // Now pos points to the first char that should be copied:
	   while (*pos)
	   {
	     if (*pos != ' ')
	     {
	       seq_in_one_line.push_back(*pos);
	     }
	     ++pos;
	   }
	 } // END if (*pos) AND read data after sequence name.
 
	 // Read lines until we have the required number of residues:
	 // TODO: More error checking could be done here:
	 // We could check that only allowed symbols are added.
	 // By this we also might be able to detect an error earlier
	 while (seq_in_one_line.length() < posNum)
	 {
	   ++curr_line_num;
	   if (curr_line_num >= all_lines_N)
	   {
	     std::cerr << "Parse error while reading the input file. The file has been interpreted as being in sequential phylip format." << std::endl;
	     std::cerr << "Several problems can trigger this error: (i) Wrong number of residues in sequences compared to the number specified " << std::endl;
	     std::cerr << "in the file header. (ii) Wrong number of taxa (sequences) compared to the number specified in the file header." << std::endl;
	     std::cerr << "(iii) Missing blank line between blocks of data if this file should be in interleaved format." << std::endl;
	     exit(-44);
	   }
	   currLine       = all_lines[curr_line_num];
	   //	   const char* pos = currLine->c_str();
	   {
	     // Find first character that shall be copied.
	     char c;
	     pos = currLine->c_str();
	     c = *pos;
	     while (!is_phylip_aa_dna_symbol(c) && !is_allowed_non_ABC(c) )
	     {
	       ++pos;
	       c = *pos;
	     }
	   }
	   // Now pos points to the first char that should be copied:
	   // Here we are in the sequential section, reading all residues until
	   // we have filled our sequence.
	   while (*pos)
	   {
	     if (*pos != ' ')
	     {
	       seq_in_one_line.push_back(*pos);
	     }
	     ++pos;
	   }
	 }
	 // If we are here, the sequence should have been copied to
	 // seq_in_one_line

	 seq = new CSequence_Mol (datatype);
	 seq->set_taxon_and_sequence(seq_taxon_name, seq_in_one_line, datatype);
	 add_seq(seq);

	 // The next line should start with a new sequence:
	 ++curr_line_num;

       } // End for loop - for all taxa

       if (curr_line_num != all_lines_N)
       {
	 std::cerr << "Parse error while reading the input file: This can have several reasons:\n"
                      "(i) The number of taxa found in the first block does not match the number specified in the file header.\n"
	              "(ii) Since only one data block was found, this file is interpreted as sequential format. If this should be a file\n"
                      "in interleaved format, blank line must be added between blocks of data.\n"
                      "(iii) It could also be that the number of residues in the sequences is larger than specified in the header, so that addtional\n"
	              "lines of data are interpreted as additional taxa in the sequential format." << std::endl;
	 exit(-67);
       }

     }    // End sequential format
     else // Interleaved format
     {
       read_mode = 'i';
       // In the first block, the first 10 characters are
       // the sequence name. Then the first base starts
       // the sequence.
       // After the first block, the sequence beginns/continues at
       // the first base symbol in each line. Thus, numbers and spaces are ignored.
       // Sequence names are not allowed in interleaved blocks.

       // With format_flag == 1 we may have a different number of sequence name characters.
       // But then, no spaces are allowed in the sequnce name. 

       unsigned     curr_line_num=0;
       faststring   *currLine;

       //       unsigned current_length = 0;
       unsigned curr_taxon = 0;

       faststring *pseq_taxon_name   = NULL;
       faststring *pseq_in_one_line  = NULL;

       std::vector<faststring*> taxonnames;
       std::vector<faststring*> sequences;

       currLine        = all_lines[curr_line_num];
       const char* pos;

       // read first block (interleaved) - until we find a blank line:
       while ( currLine->length() != 0 ) // We read all lines of the first block. Counter: curr_line_num 
       {
	 pos = currLine->c_str();

	 // We will store all pointer to sequences in the vectors defined above,
	 // so we do not delete sequences here. We keep on creating new ones.
	 pseq_taxon_name  = new faststring;
	 pseq_in_one_line = new faststring;

	 // Copy the sequence name
	 if (format_flag == 0)
	 {
	   for (i=0; i<10 && *pos != '\0'; ++i, ++pos)
	   {
	     pseq_taxon_name->push_back(*pos);
	   }
	   if (i<10) // Fill name with spaces.
	     pseq_taxon_name->append(10-i, ' ');
	 }
	 else
	 {
	   for (;*pos != '\0' && *pos != ' '; ++pos)
	   {
	     pseq_taxon_name->push_back(*pos);
	   }
	 }

	 // Skip all spaces:
	 while (*pos == ' ' && *pos != '\0')
	   ++pos;

	 // Read data after sequence name:
	 if (*pos) // we neeed to read the rest of the line
	 {
	   // Copy the sequence data:
	   {
	     char c;
	     
	     c = *pos;
	     while (!is_phylip_aa_dna_symbol(c) && !is_allowed_non_ABC(c) )
	     {
	       ++pos;
	       c = *pos;
	     }
	   }
	   // Now pos points to the first char that should be copied:
	   while (*pos)
	   {
	     if (*pos != ' ')
	       pseq_in_one_line->push_back(*pos);
	     ++pos;
	   }
	 } // if (*pos) AND read data after sequence name.
	 
	 // Another taxon has been read: Let us store it:

	 taxonnames.push_back(pseq_taxon_name);
	 sequences.push_back(pseq_in_one_line);	 

	 //	 std::cerr << "Read partial sequence" << std::endl;
	 //	 std::cerr << pseq_taxon_name->c_str() << std::endl;
	 //	 std::cerr << pseq_in_one_line->c_str() << std::endl;


	 // Read next line with next taxon name
	 ++curr_taxon;
	 ++curr_line_num;
	 if (curr_line_num >= all_lines_N)
	   break;
	 // If we break, currLine will be set to its new value further down
	 currLine       = all_lines[curr_line_num];
       } // End while loop over all non blank lines in first block (interleaved)

       // All remaining blocks still need to be read.

       if (taxaNum == 0 && posNum == 0)
	 taxaNum = curr_taxon;

       // Check that the number of taxa that might had been specified agrees with the
       // number of taxa we found now.
       if (taxaNum != 0 && taxaNum != curr_taxon) // do we have a problem?
       {
	 std::cerr << "Parse error while reading the input file:\nThe number of taxa in the first block of the input file ("
		   << curr_taxon  << ") does not match the number of taxa specified in the header of the phylip file (" << taxaNum << ") " 
		   << std::endl;
	 exit(-55);
       }

       // This line should be blank: Skip blank lines between the blocks
       while (curr_line_num < all_lines_N && all_lines[curr_line_num]->length()==0)
       {
	 ++curr_line_num;
       }

       // Now we need to read all remaining blocks:
       while (true)  // (curr_line_num < all_lines_N)
       {
	 // Read the next block:
	 for (curr_taxon=0; curr_taxon < taxaNum; ++curr_taxon)
	 {
	   // The current line should be the first line of this block
	   if (curr_line_num >= all_lines_N)
	   {
	     // We are within a block and ran out of lines:
	     std::cerr << "Parse error while reading the input file:\nUnexpected end of file. More data has been expected." << std::endl;
	     exit(-46);
	   }
	   
	   currLine       = all_lines[curr_line_num];

	   // Find first character that shall be copied.
	   char c;
	   pos = currLine->c_str();
	   c = *pos;
	   while (!is_phylip_aa_dna_symbol(c) && !is_allowed_non_ABC(c) )
	   {
	     ++pos;
	     c = *pos;
	   }

	   // Copy the sequence on this line to our data.
	   while (*pos)
	   {
	     if (*pos != ' ')
	       sequences[curr_taxon]->push_back(*pos);
	     ++pos;
	   }
	   ++curr_line_num;
	 } // End of for loop which reads the next block
	 
	 // The next line(s) should be blank or we should be at the end of the list of lines
	 if (curr_line_num < all_lines_N && all_lines[curr_line_num]->length()!=0)
	 {
	   std::cerr << "Parse error while reading the input file:\nThe phylip file was interpreted to be in the phylip format."
                        "Either the number of lines of data in one of the data blocks is not identical to the number of taxa\n"
	                "or a blank line might be missing between different data blocks in the interleaved format." << std::endl;
	   std::cerr << "This problem was detected in line " << curr_line_num+1 << " of the file, but it could have occured earlier." << std::endl;
	   exit(-77);
	 }

	 while (curr_line_num < all_lines_N && all_lines[curr_line_num]->length()==0)
	 {
	   ++curr_line_num;
	 }
	 if (curr_line_num >= all_lines_N) // we parsed all lines of the data set.
	 {
	   break;
	 }
	 // This should be the first line of the next block
	 currLine = all_lines[curr_line_num];
	 // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
       } // End while loop read all blocks

       if (taxaNum > 0 && posNum == 0)
	 posNum = sequences[0]->length();

       // Now we should have read the sequences
       // and we can delete the data we reserved.
       for (curr_taxon=0; curr_taxon < taxaNum; ++curr_taxon)
       {
	 seq = new CSequence_Mol (datatype);
	 seq->set_taxon_and_sequence(*taxonnames[curr_taxon], *sequences[curr_taxon], datatype);
	 add_seq(seq);
	 delete taxonnames[curr_taxon];
	 delete sequences[curr_taxon];
       }
     } // End interleaved format
     

     // Do some final error checking

     if (taxaNum == 0)
     {
       // delete all faststrings in all_lines
       for (unsigned ind=0; ind < all_lines.size(); ++ind)
	 delete all_lines[ind];
       return 0;
     }

     unsigned len1 = seqData[0]->length();

     if (len1 != posNum)
     {
       std::cerr << "Parse error while reading the input file:\n"
	         << "The file was interpreted in the " << (read_mode == 's'? "sequential":"interleaved") << " format.\n" 
	            "The number of positions specified in the file header is not identical to the number of positons\n"
	            "in the first sequence.\n"
	            "This can have multiple reasons:\n"
	            "In the interleaved format, make sure, taxon names only occur in the first block of sequences.\n"
	            "In the interleaved format a blank line must be found between all blocks of data." << std::endl;
       std::cerr << "The number of positions in the first sequence is: " << len1
		 << ", while the number specified in the header is " << posNum << "."<<  std::endl;

       //       std::cerr << seqData[0]->getSeqStr() << std::endl;

       exit(-66);
     }

     unsigned len2;
     // Check that the number of positions is correct and equal for all sequences: 
     for (curr_taxon=1; curr_taxon < taxaNum; ++curr_taxon)
     {
       seq = seqData[curr_taxon]; 
       len2 = seq->length();

       if (len1 != len2)
       {
	 std::cerr << "Parse error while reading the input file:\nThe number of residues in sequence " << curr_taxon+1 << " differs from the number of residues in sequence 1. The sequence lengths for sequence 1 and sequence " << curr_taxon+1 
		   <<  " are respectively " << len1 << " and " << len2 << "." << std::endl;
       exit(-67);
       }
     }

     // Check that the data type is equal for all sequences and that 
     CSequence_Mol::DataTypesEnum dt_tmp;
     
     datatype = seqData[0]->get_datatype();
     for (curr_taxon=1; curr_taxon < taxaNum; ++curr_taxon)
     {
       dt_tmp = seqData[curr_taxon]->get_datatype();
	
       if (datatype != dt_tmp)
       {
	 std::cerr << "Warning: Not all sequences seem to have the same data type.\n";
	 datatype = CSequence_Mol::mixed;
       }
     }
    
     char ambig0, ambig_tmp;
     ambig0 = seqData[0]->get_ambiguity_character();
     for (curr_taxon=1; curr_taxon < taxaNum; ++curr_taxon)
     {
       ambig_tmp = seqData[curr_taxon]->get_ambiguity_character();
	
       if (ambig0 != ambig_tmp)
       {
	 std::cerr << "Warning: Not all sequences seem to have the same ambiguity character.\n";
       }
     }
  


 
/*      // PrintMessage_cerr(seq.getName());PrintMessage_cerr("-\n");PrintMessage_cerr(seq.getFullName());PrintMessage_cerr("-\n");PrintMessage_cerr(seq.getDescription());PrintMessage_cerr("-\n"); */
/*      if ( infile.fail() && !infile.eof() ) */
/*      { */
/*        PrintMessage_cerr("\n\n"); */
/*        faststring errormsg; */
/*        errormsg   =  "An error occurred while reading the input file. It might not be a valid fasta file.\n"; */
/*        errormsg  +=  "File position: line "; */
/*        errormsg  +=  faststring(infile.line()).c_str(); */
/*        ErrorMessage(errormsg.c_str()); */
/*        PrintMessage_cerr("\n"); */
/*        flush_cerr(); */
/*        return -25; */
/*      } */
     
/*      if ( count_seq < global_sequence_from || count_seq > global_sequence_to ) */
/*      { */
/*        if (false) */
/*        { */
/* 	 PrintMessage_cerr("Skipping sequence "); */
/* 	 PrintMessage_cerr(faststring(count_seq).c_str()); */
/* 	 PrintMessage_cerr(": "); */
/* 	 PrintMessage_cerr(seq->getName()); */
/* 	 PrintMessage_cerr("\n"); */
/* 	 flush_cerr(); */
/*        } */
/*        continue; */
/*      } */
     
/*      PrintMessage_cerr("Processing sequence "); */
/*      PrintMessage_cerr(faststring(count_seq).c_str()); */
/*      PrintMessage_cerr(": "); */
/*      PrintMessage_cerr(seq->getName()); */
/*      PrintMessage_cerr("\n"); */
/*      flush_cerr(); */
     
/*      seqData.push_back(seq); */
/*    } // End while */
   
     unsigned j;
     
     for (j=0; j<all_lines.size(); ++j)
       delete all_lines[j];
     
     return 0;
   }  // End read_from_Phylip_File


   // Normal export to different formats:
   void ExportSequences(std::ostream &os, char format, unsigned interleaved_len)
   {
     if (taxaNum == 0)
       return;

     // vector of sequences with gaps and stars.
     // Why do we need a vector here?? 9.9.2912
     std::vector<faststring> vofs;

     unsigned i;

     for (i=0; i< taxaNum; ++i)
     {
       vofs.push_back(faststring('\0', (size_t)10000));
       seqData[i]->getSequence_fill_in_gaps_and_stars(vofs[i]);
     }

     if (format=='r') // relaxed phylip
     {
       std::vector<faststring> phynames;
       get_sequence_names(phynames);

       unsigned pos=0;
       unsigned i;

       for (i=0; i < taxaNum; ++i)
       {
	 phynames[i].replace_char(' ', '_');
       }
       os << "\t" << taxaNum << "\t" << posNum << std::endl;

       // write first block with sequence names:
       if (pos < posNum)
       {
	 for (i=0; i<taxaNum; ++i)
	 {
           os << phynames[i] << " ";
           os << vofs[i].substr(pos, interleaved_len) << std::endl;
	 }
	 os << std::endl; // Separate blocks
	 pos += interleaved_len;
       }

       while (pos < posNum)
       {
	 for (i=0; i<taxaNum; ++i)
	 {
           os << vofs[i].substr(pos, interleaved_len) << std::endl;
	 }
	 os << std::endl; // Separate blocks
	 pos += interleaved_len;
       }
     }
     else if (format=='p') // Phylip
     {
       std::vector<faststring> uphynames;
       unique_phylip_names(uphynames);

       unsigned pos=0;
       unsigned i;

       for (i=0; i<taxaNum; ++i)
       {
	 uphynames[i].fill_if_shorter(10, ' ');
       }

       os << "\t" << taxaNum << "\t" << posNum << std::endl;

       // write first block with sequence names:
       if (pos < posNum)
       {
	 for (i=0; i<taxaNum; ++i)
	 {
           os << uphynames[i] << " ";
           os << vofs[i].substr(pos, interleaved_len) << std::endl;
	 }
	 os << std::endl; // Separate blocks
	 pos += interleaved_len;
       }

       while (pos < posNum)
       {
	 for (i=0; i<taxaNum; ++i)
	 {
           os << vofs[i].substr(pos, interleaved_len) << std::endl;
	 }
	 os << std::endl; // Separate blocks
	 pos += interleaved_len;
       }
     }
     else if (format=='f')
     {
       unsigned pos=0;
       unsigned i;

       for (i=0; i<taxaNum; ++i)
       {
         os << ">" << seqData[i]->getFullName() << std::endl;
	 pos = 0;
	 while (pos < posNum)
	 {
           os << vofs[i].substr(pos, interleaved_len) << std::endl;
	   pos += interleaved_len;
	 }
       }
     }
     else if (format=='n')
     {
/*        Begin data; */
/*        Dimensions ntax=4 nchar=15; */
/*        Format datatype=dna symbols="ACTG" missing=? gap=-; */
/*        Matrix */
/* 	 Species1   atgctagctagctcg */
/* 	 Species2   atgcta??tag-tag */
/* 	 Species3   atgttagctag-tgg */
/* 	 Species4   atgttagctag-tag            */
/* 	 ; */
/*        End; */

       os << "Begin data;" << std::endl;
       os << "Dimensions ntax="<< taxaNum << " nchar="<< posNum << ";" << std::endl;
       os << "Format datatype="<< seqData[0]->type_as_string()
	 //  << " symbols=\" "<< <<  "\"" 
	 << " missing=" << ambig_char << " gap=-;" << std::endl;
       os << "matrix" << std::endl;
       
       unsigned pos=0;
       unsigned i;

       while (pos < posNum)
       {
	 for (i=0; i<taxaNum; ++i)
	 {
           os << "'" << seqData[i]->getName() << "'" << " ";
           os << vofs[i].substr(pos, interleaved_len) << std::endl;
	 }
	 pos += interleaved_len;
       }
       os << ";" << std::endl << "end;" << std::endl;

     }
     else
     {
       std::cerr << "Internal error in ExportSequences(std::ostream, char, unsigned):\n"
	            "Unkown value for format parameter. No output has been generated."
		 << std::endl;
     }
   }

   // Export do different formats with site-filter
   // Commented out since not finished and not tested.
/*    void ExportSequences(std::ostream &os, char format, unsigned interleaved_len, */
/* 			std::vector<unsigned> site_filter) */
/*    { */
/*      // TODO: Currently onyl supports fasta without interleaved_len */

/*      if (taxaNum == 0) */
/*        return; */

/*      // vector of sequences with gaps and stars. */
/*      // Why do we need a vector here?? 9.9.2912 */
/*      std::vector<faststring> vofs; */

/*      unsigned i; */

/*      for (i=0; i< taxaNum; ++i) */
/*      { */
/*        vofs.push_back(faststring('\0', (size_t)10000)); */
/*        seqData[i]->getSequence_fill_in_gaps_and_stars(vofs[i]); */
/*      } */

/*      if (format=='p') // Phylip */
/*      { */
/* /\*        std::vector<faststring> uphynames; *\/ */
/* /\*        unique_phylip_names(uphynames); *\/ */

/* /\*        unsigned pos=0; *\/ */
/* /\*        unsigned i; *\/ */

/* /\*        os << "\t" << taxaNum << "\t" << posNum << std::endl; *\/ */

/* /\*        while (pos < posNum) *\/ */
/* /\*        { *\/ */
/* /\* 	 for (i=0; i<taxaNum; ++i) *\/ */
/* /\* 	 { *\/ */
/* /\*            os << uphynames[i] << " "; *\/ */
/* /\*            os << vofs[i].substr(pos, interleaved_len) << std::endl; *\/ */
/* /\* 	 } *\/ */
/* /\* 	 os << std::endl; *\/ */
/* /\* 	 pos += interleaved_len; *\/ */
/* /\*        } *\/ */
/*      } */
/*      else if (format=='f') */
/*      { */
/*        unsigned pos=0; */
/*        unsigned i; */

/*        for (i=0; i<taxaNum; ++i) */
/*        { */
/*          os << ">" << seqData[i]->getFullName() << std::endl; */
/* 	 pos = 0; */
/* 	 faststring &seq_ref = vofs[i]; */
/* 	 while (pos < posNum) */
/* 	 { */
/* 	   if (site_filter[pos] != 0) */
/* 	     os << seq_ref[pos] << std::endl; */
/* 	   ++pos; */
/* 	 } */
/*        } */
/*      } */
/*      else if (format=='n') */
/*      { */
/* /\*        Begin data; *\/ */
/* /\*        Dimensions ntax=4 nchar=15; *\/ */
/* /\*        Format datatype=dna symbols="ACTG" missing=? gap=-; *\/ */
/* /\*        Matrix *\/ */
/* /\* 	 Species1   atgctagctagctcg *\/ */
/* /\* 	 Species2   atgcta??tag-tag *\/ */
/* /\* 	 Species3   atgttagctag-tgg *\/ */
/* /\* 	 Species4   atgttagctag-tag            *\/ */
/* /\* 	 ; *\/ */
/* /\*        End; *\/ */

/* /\*        os << "Begin data;" << std::endl; *\/ */
/* /\*        os << "Dimensions ntax="<< taxaNum << " nchar="<< posNum << ";" << std::endl; *\/ */
/* /\*        os << "Format datatype="<< seqData[0]->type_as_string() *\/ */
/* /\* 	 //  << " symbols=\" "<< <<  "\""  *\/ */
/* /\* 	 << " missing=" << ambig_char << " gap=-;" << std::endl; *\/ */
/* /\*        os << "matrix" << std::endl; *\/ */
       
/* /\*        unsigned pos=0; *\/ */
/* /\*        unsigned i; *\/ */

/* /\*        while (pos < posNum) *\/ */
/* /\*        { *\/ */
/* /\* 	 for (i=0; i<taxaNum; ++i) *\/ */
/* /\* 	 { *\/ */
/* /\*            os << "'" << seqData[i]->getName() << "'" << " "; *\/ */
/* /\*            os << vofs[i].substr(pos, interleaved_len) << std::endl; *\/ */
/* /\* 	 } *\/ */
/* /\* 	 pos += interleaved_len; *\/ */
/* /\*        } *\/ */
/* /\*        os << ";" << std::endl << "end;" << std::endl; *\/ */

/*      } */

/*    } */


   void ExportSequences_no_fill_in_ext_phylip_range(std::ostream &os, unsigned begin_exp, unsigned end_exp)
   {
     if (taxaNum == 0)
       return;

     if (end_exp > posNum || end_exp <= begin_exp)
       return;

     unsigned i;
     
     os << "\t" << taxaNum << "\t" << end_exp-begin_exp << std::endl;
     
     for (i=0; i<taxaNum; ++i)
     {
       os.write(seqData[i]->getName(), seqData[i]->getName_faststring().length());
       os.put(' ');
       os.write(seqData[i]->getSeqBegin()+begin_exp, end_exp-begin_exp);
       os.put('\n');
     }
   }


   void ExportSequences_no_fill_in_ext_phylip_range(std::ostream &os, unsigned *coords, unsigned num_coord_pairs)
   {
     if (taxaNum == 0)
       return;

     //     if (end_exp > posNum || end_exp <= begin_exp)
     //       return;

     unsigned i,j, N=0;
     unsigned begin_exp;
     unsigned end_exp;     


     for (j=0; j<num_coord_pairs; ++j)
     {
	 begin_exp = coords[2*j];
	 end_exp   = coords[2*j+1];
	 if (end_exp > posNum)
	   end_exp = posNum;
	 if (begin_exp < end_exp)
	   N += end_exp-begin_exp;
     }

     os << "\t" << taxaNum << "\t" << N << std::endl;
     
     for (i=0; i<taxaNum; ++i)
     {
       os.write(seqData[i]->getName(), seqData[i]->getName_faststring().length());
       os.put(' ');
       for (j=0; j<num_coord_pairs; ++j)
       {
	 begin_exp = coords[2*j];
	 end_exp   = coords[2*j+1];
	 if (end_exp > posNum)
	   end_exp = posNum;
	 if (begin_exp < end_exp)
	   os.write(seqData[i]->getSeqBegin()+begin_exp, end_exp-begin_exp);
       }
       os.put('\n');
     }
   }



   // Not yet well tested!!!!!
   void Export_Single_Sequence(std::ostream &os, const faststring& fullname,
			       char format, unsigned interleaved_len)
   {
     if (taxaNum == 0)
       return;

     // vector of sequences with gaps and stars.
     // Why do we need a vector here?? 9.9.2912
     std::vector<faststring> vofs;

     unsigned i;

     for (i=0; i< taxaNum; ++i)
     {
       vofs.push_back(faststring());
       if ( seqData[i]->getFullName() == fullname)
	 seqData[i]->getSequence_fill_in_gaps_and_stars(vofs[i]);
     }

     if (format=='p') // Phylip
     {
       std::vector<faststring> uphynames;
       unique_phylip_names(uphynames);

       unsigned pos=0;
       unsigned i;

       for (i=0; i<taxaNum; ++i)
	 uphynames[i].fill_if_shorter(10, ' ');

       os << "\t" << taxaNum << "\t" << posNum << std::endl;

       while (pos < posNum)
       {
	 for (i=0; i<taxaNum; ++i)
	 {
	   if (!vofs[i].empty())
	   {
	     os << uphynames[i] << " ";
	     os << vofs[i].substr(pos, interleaved_len) << std::endl;
	   }
	 }
	 os << std::endl;
	 pos += interleaved_len;
       }
     }
     else if (format=='f')
     {
       unsigned pos=0;
       unsigned i;

       for (i=0; i<taxaNum; ++i)
       {
	 if (!vofs[i].empty())
	 {
	   os << ">" << seqData[i]->getFullName() << std::endl;
	   pos = 0;
	   while (pos < posNum)
	   {
	     os << vofs[i].substr(pos, interleaved_len) << std::endl;
	     pos += interleaved_len;
	   }
	 }
       }
     }
     else if (format=='n')
     {
/*        Begin data; */
/*        Dimensions ntax=4 nchar=15; */
/*        Format datatype=dna symbols="ACTG" missing=? gap=-; */
/*        Matrix */
/* 	 Species1   atgctagctagctcg */
/* 	 Species2   atgcta??tag-tag */
/* 	 Species3   atgttagctag-tgg */
/* 	 Species4   atgttagctag-tag            */
/* 	 ; */
/*        End; */

       os << "Begin data;" << std::endl;
       os << "Dimensions ntax="<< 1 << " nchar="<< posNum << ";" << std::endl;
       os << "Format datatype="<< seqData[0]->type_as_string()
	 //  << " symbols=\" "<< <<  "\"" 
	 << " missing=" << ambig_char << " gap=-;" << std::endl;
       os << "matrix" << std::endl;
       
       unsigned pos=0;
       unsigned i;

       while (pos < posNum)
       {
	 for (i=0; i<taxaNum; ++i)
	 {
	   if (!vofs[i].empty())
	   {
	     os << seqData[i]->getName() << " ";
	     os << vofs[i].substr(pos, interleaved_len) << std::endl;
	   }
	 }
	 pos += interleaved_len;
       }
       os << ";" << std::endl << "end;" << std::endl;

     }


   }





   //   void            GetBaseFrequencies(unsigned, std::vector<unsigned>);
   
   //   void            initSplits();
   //   void            computeSplits();
   //   CSplits         *splits;
   //   CBasefreq       *basefreq;

   double compute_site_pattern_stat(bool filter_gap_abmig_symbols=true)
   {
     unsigned i, j;
     
     // Here we will store pointers to all sequences for fast access:
     const char **tax = new const char* [taxaNum];

     if (taxaNum == 0 || posNum == 0)
       return 0;

     // Store pointers to all sequences in the tax array - this will be faster
     for (j=0; j < taxaNum; ++j)
     {
       tax[j] = seqData[j]->getSeqStr();
     }

     faststring tmp;
     std::map<faststring, unsigned> m;

     bool contains_gap_ambig;
     char c;

     // For all positions in alignment
     for (i=0; i<posNum; ++i)
     {
       tmp.clear();
       contains_gap_ambig = false;

       // Old and slow version:
       // Copy all symbols to site pattern vector
/*        for (j=0; j < taxaNum; ++j) */
/*        { */
/* 	 c = tax[j][i]; */
/*  	 if (filter_gap_abmig_symbols) */
/*  	 { */
/* 	   // Check whether this is a position we do not want to consider. */
/*  	   if (datatype == CSequence_Mol::dna && (is_DNA_iupac_ambig(c) || c == '-') ) */
/*  	   { */
/* 	     contains_gap_ambig = true; */
/*  	   } */
/* 	   else */
/* 	   { */
/* 	     if (datatype == CSequence_Mol::protein && (is_aa_ambig(c) || c == '-') ) */
/* 	     { */
/* 	       contains_gap_ambig = true; */
/* 	     } */
/* 	   } */
/*  	 } // END if (filter_gap_abmig_symbols) */
/* 	 tmp.push_back(c); */
/*        } // END: For all taxa */

       // Faster version:
       if (filter_gap_abmig_symbols && datatype == CSequence_Mol::dna)
       {
	 for (j=0; j < taxaNum; ++j) // Copy all symbols to site pattern vector
	 {
	   c = tax[j][i];
	   if (is_DNA_iupac_ambig(c) || c == '-')
 	   {
	     contains_gap_ambig = true;
 	   }
	   tmp.push_back(c);
	 }
       }
       else if (filter_gap_abmig_symbols && datatype == CSequence_Mol::protein)
       {
	 for (j=0; j < taxaNum; ++j) // Copy all symbols to site pattern vector
	 {
	    c = tax[j][i];
	   if (is_aa_ambig(c) || c == '-')
 	   {
	     contains_gap_ambig = true;
 	   }
	   tmp.push_back(c);
	 }
       }
       else
       {
	 for (j=0; j < taxaNum; ++j) // Copy all symbols to site pattern vector
	 {
	    c = tax[j][i];
	    tmp.push_back(c);
	 } 
       }
       if (!contains_gap_ambig)  
	 add_or_count(m, tmp);
     } // END: for all positions in alignment

     std::map<faststring, unsigned>::iterator it, it_end;
     it     = m.begin();
     it_end = m.end();

     double res=0;
     double x;

     // Formula from "Statistical Tests of Models of DNA Substitution", Journal of Molecular Evolution, Springer-Verlag NewYork Inc. 1993
     // Nick Goldman
     // Formula: \sum   N_i*ln(N_i) - N*ln(N), where N_i are the pattern frequencies and N ist the total number of sites.

     while (it != it_end)
     {
       x = it->second;
       if (x > 0)
	 res += x*log(x);
       //       std::cout << x << " " << res << std::endl;
       ++it;
     }
     delete [] tax;

     x = posNum;
     res -= x*log(x);
     //     std::cout << x << " " << res << std::endl;

     return res;
   } // END compute_site_pattern_stat(bool filter_gap_abmig_symbols=true)


   double compute_site_pattern_stat(std::map<faststring, unsigned> &m, bool filter_gap_abmig_symbols=true)
   {
     unsigned i, j;
     m.clear();
     
     // Here we will store pointers to all sequences for fast access:
     const char **tax = new const char* [taxaNum];

     if (taxaNum == 0 || posNum == 0)
       return 0;

     // Store pointers to all sequences in the tax array - this will be faster
     for (j=0; j < taxaNum; ++j)
     {
       tax[j] = seqData[j]->getSeqStr();
     }

     faststring tmp;

     bool contains_gap_ambig;
     char c;

     // For all positions in alignment
     for (i=0; i<posNum; ++i)
     {
       tmp.clear();
       contains_gap_ambig = false;

       // Old and slow version:
       // Copy all symbols to site pattern vector
/*        for (j=0; j < taxaNum; ++j) */
/*        { */
/* 	 c = tax[j][i]; */
/*  	 if (filter_gap_abmig_symbols) */
/*  	 { */
/* 	   // Check whether this is a position we do not want to consider. */
/*  	   if (datatype == CSequence_Mol::dna && (is_DNA_iupac_ambig(c) || c == '-') ) */
/*  	   { */
/* 	     contains_gap_ambig = true; */
/*  	   } */
/* 	   else */
/* 	   { */
/* 	     if (datatype == CSequence_Mol::protein && (is_aa_ambig(c) || c == '-') ) */
/* 	     { */
/* 	       contains_gap_ambig = true; */
/* 	     } */
/* 	   } */
/*  	 } // END if (filter_gap_abmig_symbols) */
/* 	 tmp.push_back(c); */
/*        } // END: For all taxa */

       // Faster version:
       if (filter_gap_abmig_symbols && datatype == CSequence_Mol::dna)
       {
	 for (j=0; j < taxaNum; ++j) // Copy all symbols to site pattern vector
	 {
	   c = tax[j][i];
	   if (is_DNA_iupac_ambig(c) || c == '-')
 	   {
	     contains_gap_ambig = true;
 	   }
	   tmp.push_back(c);
	 }
       }
       else if (filter_gap_abmig_symbols && datatype == CSequence_Mol::protein)
       {
	 for (j=0; j < taxaNum; ++j) // Copy all symbols to site pattern vector
	 {
	    c = tax[j][i];
	   if (is_aa_ambig(c) || c == '-')
 	   {
	     contains_gap_ambig = true;
 	   }
	   tmp.push_back(c);
	 }
       }
       else
       {
	 for (j=0; j < taxaNum; ++j) // Copy all symbols to site pattern vector
	 {
	    c = tax[j][i];
	    tmp.push_back(c);
	 } 
       }
       if (!contains_gap_ambig)  
	 add_or_count(m, tmp);
     } // END: for all positions in alignment

     std::map<faststring, unsigned>::iterator it, it_end;
     it     = m.begin();
     it_end = m.end();

     double res=0;
     double x;

     while (it != it_end)
     {
       x = it->second;
       if (x > 0)
	 res += x*log(x);
       //       std::cout << x << " " << res << std::endl;
       ++it;
     }
     delete [] tax;

     x = posNum;
     res -= x*log(x);
     //     std::cout << x << " " << res << std::endl;

     return res;
   }



   void get_partial_pattern(unsigned pos, unsigned n, faststring &pat)
   {
     unsigned i;

     if (n > taxaNum)
       n = taxaNum;

     pat.clear();
     for(i=0; i<n; ++i)
       pat.push_back((seqData[i]->getSeqStr())[pos]);
   }



   // Does not yet provide the possibility to remove amig positions.
   void get_vector_of_site_pattern_frequencies(fastvector<unsigned> &spf)
   {
     unsigned i, j;
     
     const char **tax = new const char* [taxaNum];

     for (j=0; j < taxaNum; ++j)
     {
       tax[j] = seqData[j]->getSeqStr();
     }

     faststring tmp;
     std::map<faststring, unsigned> m;

     for (i=0; i<posNum; ++i)
     {
       tmp.clear();

       for (j=0; j < taxaNum; ++j)
       {
	 tmp.push_back(tax[j][i]);
       }
       add_or_count(m, tmp);
     }
     std::map<faststring, unsigned>::iterator it, it_end;
     it     = m.begin();
     it_end = m.end();

     while (it != it_end)
     {
       spf.push_back(it->second);
       ++it;
     }
     delete [] tax;
   }

   void clear(CSequence_Mol::DataTypesEnum     datatype_param = CSequence_Mol::unknown,
	      char ambig_char_param = '?')
   {
     taxaNum    = 0;
     posNum     = 0;
     ambig_char = ambig_char_param;
     datatype   = datatype_param;

     originalPosNumbers.clear();
     originalPosNumbers_supplied = false;

     int i, n=seqData.size();

     // Recently resolved memory leak:
     for (i=0; i<n; ++i)
         delete seqData[i];

     seqData.clear();
     sn_map.clear();
   }

   void remove_gap_only_sequences()
   {
     int i, n=seqData.size();

     // We do this in reversed order, since this avoids adapting i and n after
     // an element has been removed.    
     for (i=n-1; i>=0; --i)
     {
       if (seqData[i]->has_only_gaps() )
       {
	 remove_del(i);
       }
     }
   }

   void remove_ambig_only_sequences()
   {
     int i, n = seqData.size();
    
     // We do this in reversed order, since this avoids adapting i and n after
     // an element has been removed.
     for (i=n-1; i>=0; --i)
     {
       if (seqData[i]->has_only_ambigs() )
       {
	 remove_del(i);
       }
     }
   }

   void remove_gap_or_ambig_only_sequences()
   {
     int i, n = seqData.size();

     // We do this in reversed order, since this avoids adapting i and n after
     // an element has been removed.
     for (i=n-1; i>=0; --i)
     {
       //              std::cout << "Checking sequence: " << i+1 << std::endl;
       if (seqData[i]->has_only_ambigs_or_gap() )
       {
	 remove_del(i);
       }
     }
   }

   void steal_sequence(CSequences2 &seqs, unsigned i)
   {
     if (i < seqs.taxaNum)
     {
       add_seq(seqs.seqData[i]);
       seqs.remove(i);
     }
   }

   void remove_del(unsigned seq_num)
   {
     CSequence_Mol *p = seqData[seq_num];
     faststring sname = p->getName();

     --taxaNum;
     seqData.erase(seqData.begin()+seq_num);
     std::map<faststring, CSequence_Mol*>::iterator find_it = sn_map.find(sname);
     if (find_it != sn_map.end() )
       sn_map.erase(find_it);
     delete p;
   }

   CSequence_Mol* remove(unsigned seq_num)
   {
     CSequence_Mol *p = seqData[seq_num];
     faststring sname = p->getName();

     --taxaNum;
     seqData.erase(seqData.begin()+seq_num);
     std::map<faststring, CSequence_Mol*>::iterator find_it = sn_map.find(sname);
     if (find_it != sn_map.end() )
       sn_map.erase(find_it);
     return p;
   }

   void add_seq_to_alignment(CSequence_Mol::DataTypesEnum     datatype_param,
                             const faststring                  &fullname_param,
                             const faststring                  &seq_data_param,
			     char                              ambig_param='\0' // default auto
			    )
   {
     if (datatype != datatype_param)
     {
	 std::cerr << "Warning: Not all sequences seem to have the same data type.\n";
	 datatype = CSequence_Mol::mixed;
     }

     //     std::cout << "datatype_param: " << datatype_param  << std::endl;

     // Handle general ambig TODO XXXXXXXXXXXXXXXXXXXXXXXXXXXXX
     //     if (general_ambig == CSequence_Mol::unknown)

     CSequence_Mol   *seq;
     seq = new CSequence_Mol (datatype_param, ambig_param);
     //     std::cout << "Datatype after constructor: " << seq->get_datatype() << std::endl;
     seq->set_taxon_and_sequence(fullname_param, seq_data_param, datatype_param);
     //     std::cout << "Datatype after constructor+set_taxon_and_sequence: " << seq->get_datatype() << std::endl;




     ++taxaNum;
     if (posNum == 0)
       posNum = seq->length();
     else if (posNum != seq->length())
     {
       std::cerr << "Error: Sequence added has different length than existing sequences." << std::endl;
       delete seq;
       return;
     }

     add_seq(seq);

     // Further autodetect stuff???????


   }


   void alignment_without_gap_N_positions(CSequences2 &seqs)
   {
     seqs.clear(datatype, seqs.ambig_char);

     std::vector<bool> gap_N_positions(posNum, false);

     size_t i, j;
     const char *the_seq;
     char c;

     if (datatype == CSequence_Mol::dna)
     {
       for (i=0; i< taxaNum; ++i)
       {
	 the_seq = seqData[i]->getSeqStr();
	 for(j=0; j<posNum; ++j)
	 {
	   c = the_seq[j];
	   if (is_DNA_iupac_ambig(c) || c == '-')
	   {
	     gap_N_positions[j] = true;
	   }
	 }
       }
     }

     if (datatype == CSequence_Mol::protein)
     {
       for (i=0; i< taxaNum; ++i)
       {
	 the_seq = seqData[i]->getSeqStr();
	 for(j=0; j<posNum; ++j)
	 {
	   c = the_seq[j];
	   if (is_aa_ambig(c) || c == '-')
	   {
	     gap_N_positions[j] = true;
	   }
	 }
       }
     }

     faststring new_seq;

     for (i=0; i< taxaNum; ++i)
     {
       new_seq.clear();
       the_seq = seqData[i]->getSeqStr();

       for(j=0; j<posNum; ++j)
       {
	 if (!gap_N_positions[j])
	   new_seq.push_back(the_seq[j]);
       }
       seqs.add_seq_to_alignment(datatype, seqData[i]->getFullName(), new_seq, seqData[i]->get_ambiguity_character() );
     }
   }


   void alignment_without_gap_only_positions(CSequences2 &seqs)
   {
     //     std::cout << "Entering alignment_without_gap_only_positions" << std::endl;

     seqs.clear(datatype, ambig_char);

     // Count the number of gaps at all positions.
     // Gap only positions are those that have taxaNum gaps.
     std::vector<unsigned> gap_positions(posNum, 0);

     size_t i, j;
     const char *the_seq;
     char c;

     for (i=0; i< taxaNum; ++i)
     {
       the_seq = seqData[i]->getSeqStr();
       for(j=0; j<posNum; ++j)
       {
	 c = the_seq[j];
	 if (c == '-')
	 {
	   ++gap_positions[j];
	 }
       }
     }

     faststring new_seq_str;

     for (i=0; i< taxaNum; ++i)
     {
       new_seq_str.clear();
       the_seq = seqData[i]->getSeqStr();

       for(j=0; j<posNum; ++j)
       {
	 if (gap_positions[j] != taxaNum) // only gaps
	   new_seq_str.push_back(the_seq[j]);
       }
       seqs.add_seq_to_alignment(datatype,seqData[i]->getFullName(),new_seq_str,
				 seqData[i]->get_ambiguity_character());
     }
/*      std::cout << "DEBUGOUT-alignment_without_gap_only_positions" << std::endl; */
/*      seqs.print_DEBUG(std::cout,3); */
   }


   void replace_sequence_interval(const CSequences2 &s, unsigned pos1=0, unsigned pos2=-1u) 
   {
     determine_map_of_sequence_names();

     faststring name;
     unsigned i, N = s.seqData.size();
     std::map<faststring, CSequence_Mol*>::iterator find_it;
 
    if (taxaNum != s.taxaNum)
     {
 	 std::cerr << "Critical error in replace_sequence_interval: Number of taxa are not equal"
		   << std::endl;
	 exit(-1);       
     }

     for (i=0; i<N; ++i)
     {
       name = s.seqData[i]->getName();
       find_it = sn_map.find(name);
       if (find_it == sn_map.end())
       {
	 std::cerr << "Critical error in replace_sequence_interval: Sequence "
		   << name << " does not exist." << std::endl;
	 exit(-1); 
       }
       find_it->second->replace_part_of_sequence(*(s.seqData[i]), pos1, pos2);
     }
     if (seqData.size() > 0)
       posNum = seqData[0]->length();
   }

   void append_to_sequences(const CSequences2 &s) 
   {
     determine_map_of_sequence_names();

     faststring name;
     unsigned i,N=s.seqData.size();
     std::map<faststring, CSequence_Mol*>::iterator find_it;

     if (taxaNum != s.taxaNum)
     {
 	 std::cerr << "Critical error in append_sequences: Number of taxa are not equal"
		   << std::endl;
	 exit(-1);       
     }

     for (i=0; i<N; ++i)
     {
       name = s.seqData[i]->getName();
       find_it = sn_map.find(name);
       if (find_it == sn_map.end())
       {
	 std::cerr << "Critical error in append_sequences: Sequence "
		   << name << " does not exist." << std::endl;
	 exit(-1); 
       }
       find_it->second->append_sequence(*(s.seqData[i]));
     }
     if (seqData.size() > 0)
       posNum = seqData[0]->length();
   }

   void append_to_sequences(const faststring &site_pattern)
   {
     unsigned i;

     if ( taxaNum != site_pattern.size() )
     {
       std::cerr << "Error when calling the append_to_sequences method:\n"
	            "the number of sites in the pattern that shall be\n"
                    "appended the sequences ofject differs from the\n"
	            "number of sequences. This pattern will not be appended." << std::endl; 
       exit(-22);
     }
     for (i=0; i<taxaNum; ++i)
     {
       seqData[i]->append_residue_unchecked(site_pattern[i]);
     }
     ++posNum;
   }


   bool is_insertion_column(unsigned col_pos)
   {
     unsigned i_tax;
     char c;

     for (i_tax=0; i_tax<taxaNum; ++i_tax)
     {
       c = seqData[i_tax]->get_pos(col_pos);
       if (c == '.' || (c >='a' && c <= 'z') )
	 return true;
     }
     return false;
   }


   bool is_mainly_gap_pos(unsigned col_pos, double min_prop)
   {
     unsigned i_tax;
     char c;
     unsigned count_gap= 0;

     for (i_tax=0; i_tax<taxaNum; ++i_tax)
     {
       c = seqData[i_tax]->get_pos(col_pos);
       if (c == '-' || c == '.')
	 ++count_gap;
     }
     if ((float)count_gap/taxaNum >= min_prop)
       return true;
     else
       return false;
   }


   // Assumes that insertion positions have lower case characters
   // and they could have dots '.' instead of gaps.
   // The letter is not necessary, since insertion positions should
   // at least have one lower case character.
   // So dots could have been replaced by gaps already.
   unsigned get_next_insertion_column(unsigned col_pos)
   {
     for (; col_pos < posNum; ++col_pos)
       if (is_insertion_column(col_pos))
	 return col_pos;
     return col_pos;
   }

   unsigned get_next_insertion_or_mainly_gap_column(unsigned col_pos, float min_gap_prop)
   {
     for (; col_pos < posNum; ++col_pos)
       if (is_insertion_column(col_pos) || is_mainly_gap_pos(col_pos, min_gap_prop))
	 return col_pos;
     return col_pos;
   }

   // Assumes that non insertion positions have upper case characters.
   // If dots have already been converted to gaps, they are not
   // indicative of non insertion positions!!
   unsigned get_next_non_insertion_column(unsigned col_pos)
   {
     for (; col_pos < posNum; ++col_pos)
       if (!is_insertion_column(col_pos))
	 return col_pos;
	 //	 if (c == '.' || c >='a' && c <= 'z') // Not a non-insertion column
     return col_pos;
   }

   unsigned get_next_non_insertion_and_non_mainly_gap_column(unsigned col_pos, float  min_gap_prop)
   {
     for (; col_pos < posNum; ++col_pos)
       if (!(is_insertion_column(col_pos) || is_mainly_gap_pos(col_pos, min_gap_prop)) )
	 return col_pos;
	 //	 if (c == '.' || c >='a' && c <= 'z') // Not a non-insertion column
     return col_pos;
   }

   // RETURN VALUE NEEDS TO BE SET. DO NOT USE IN THIS FORM.
/*    void get_bp_entropies_nucleotide(double *ent) */
/*    { */
/*      const char **tax = new const char* [taxaNum]; */
/*      unsigned  j; */

/*      if (taxaNum == 0 || posNum == 0) */
/*        return; */

/*      // Store pointers to all sequences in the tax array - this will be faster */
/*      for (j=0; j < taxaNum; ++j) */
/*      { */
/*        tax[j] = seqData[j]->getSeqStr(); */
/*      } */
/*    } */

   // Computes relative amount of residues that are non gap and non ambig in this range:
   double get_non_ambig_non_gap_coverage_in_range(unsigned pos, unsigned pos_end, int flag=0)
   {
     const char **tax = new const char* [taxaNum];
     unsigned i, j;
     unsigned N = 0;


     if (taxaNum == 0 || posNum == 0)
       return -1; // Indicates an error

     // Store pointers to all sequences in the tax array - this will be faster
     if (flag == 0 || flag == 1) // 1kite all or general mode
     {
       N = taxaNum;
       for (j=0; j < taxaNum; ++j)
       {
	 tax[j] = seqData[j]->getSeqStr();
       }
     }
     else if (flag == 13) // Special 1kite flag: Filter hamstered reference taxa:
     {
#ifdef DEBUG
       std::cerr << "Called: CSequences2.h:get_non_ambig_non_gap_coverage_in_range with flag==1" << std::endl;
#endif

       unsigned jj=0;

       for (j=0; j < taxaNum; ++j)
       {
	 // In case we compute coverage values, we only have hamstered reference taxa left in 1kite:
	 if (is_1kite_hamstered_reference_seq_id(seqData[j]->getFullName() ) )
	 {
	   tax[jj] = seqData[j]->getSeqStr();
	   ++jj;
	 }
       }
       N = jj;
     }
     else
     {
       return -1; // Indicates an error
     }


     unsigned count_all=0;
     unsigned count_non_ambig_non_gap=0;
     double   res;
     char     c;

     if (datatype ==  CSequence_Mol::protein)
     {
       for (i =pos; i< pos_end && i < posNum; ++i)
       {
	 for (j=0; j< N; ++j)
	 {
	   c = tax[j][i];
	   ++count_all;
	   if ( !(is_aa_ambig(c) || c == '-') )
	     ++count_non_ambig_non_gap;
	 }
       }
     }
     else if (datatype ==  CSequence_Mol::dna)
     {
       for (i =pos; i< pos_end && i < posNum; ++i)
       {
	 for (j=0; j< taxaNum; ++j)
	 {
	   c = tax[j][i];
	   ++count_all;
	   if ( !(is_DNA_iupac_ambig(c) || c == '-' ) )
	     ++count_non_ambig_non_gap;
	 }
       }
     }
     res = count_non_ambig_non_gap;
     res /= count_all;

     //          std::cerr << count_non_ambig_non_gap << std::endl;
     //          std::cerr << count_all               << std::endl;

     return res;
   }



   void get_1kite_original_reference_taxa_vec(std::vector<faststring> &vec) const
   {
#ifdef DEBUG
       std::cerr << "Called: CSequences2.h:get_1kite_original_reference_taxa_vec" << std::endl;
#endif
     unsigned i;

     for (i=0; i < taxaNum; ++i)
     {
       if ( is_1kite_hamstered_reference_seq_id( seqData[i]->getFullName() ) )
       {
	 vec.push_back(seqData[i]->getFullName() );
       }
     }
   }

   // Determines all 1kite original reference taxa, sorts them and moves them to the top of the alignment.
   int move_1kite_original_reference_taxa_sorted_to_top()
   {
#ifdef DEBUG
     std::cerr << "Called: CSequences2.h:move_1kite_original_reference_taxa_sorted_to_top" << std::endl;
#endif

     std::vector<faststring> vec;
     get_1kite_original_reference_taxa_vec(vec);
     std::sort(vec.begin(), vec.end());

     int i;

     for (i=vec.size()-1; i>= 0; --i)
     {
       reorder_move_seq_to_top(vec[i]);
     }
     return vec.size();
   }


/*    int move_1kite_core_of_other_seqs_to_top(const CSequences2 &other_seqs) */
/*    { */
/* #ifdef DEBUG */
/*      std::cerr << "Called: CSequences2.h:move_1kite_core_of_other_seqs_to_top" << std::endl; */
/* #endif */

/*      std::vector<faststring> vec; */

/*      other_seqs.get_1kite_cores_vec(vec); */
     
/*      int i; */

/*      for (i=vec.size()-1; i>= 0; --i) */
/*      { */
/*        reorder_move_seq_to_top(vec[i]); */
/*      } */
/*      return vec.size(); */
/*    } */



   // Currently only implemented for datatype DNA
   //  Allowed flag values:
   // 0:  normal mode
   // 1:  Gap regions always get a distance of -3
   // 2:  Compute gap distance, i.e. count gap versus non-gap nucleotides.
   // The result is in the range from 0..1.
   // Ns are treated as mismatch. If this is not desired, use the other overwrite of this function below.
   double sequence_distance(unsigned s1, unsigned s2, unsigned begin_range, unsigned end_range, unsigned flag=0)
   {
     char   c1, c2;

     if (s1 >= taxaNum || s2 >= taxaNum)
       return -1;

     CSequence_Mol* seq1 = seqData[s1];
     CSequence_Mol* seq2 = seqData[s2];

     // Do some range checks:
     if (end_range < begin_range) // This is considered a major error and the caller should correct this problem.
     {
       std::cerr << "Error in sequence_distance: end_range < begin_range" << std::endl;
       exit(-33);
     }

     // Out of range checks are corrected silently.
     // It is not clear this is always the best method to deal with this.
     if ((end_range > seq1->length()))
     {
       end_range = seq1->length();
       if (begin_range > seq1->length())
	 begin_range = seq1->length();
     }
     if ((end_range > seq2->length()))
     {
       end_range = seq2->length();
       if (begin_range > seq2->length())
	 begin_range = seq2->length();
     }
     
     unsigned i = begin_range;

     unsigned len  = 0;
     unsigned diff = 0;

     if (flag < 2)
     {
       if (seq1->get_datatype() != CSequence_Mol::dna)
	 return -2;

       while (i != end_range)
       {
	 // TODO: DO WE NEED TOUPPER??
	 c1 = toupper(seq1->get_pos(i));
	 c2 = toupper(seq2->get_pos(i));

	 if (c1 != '-' && c2 != '-')
	 {
	   ++len;
	   if (c1 != c2 || c1 == 'N') // If they differ or if both are Ns we count this as a difference.
	   {
	     ++diff;
	   }
	 }
	 else
	 {
	   if (flag == 1)
	     return -3.0;
	 }
	 ++i;
       }
       if (len == 0)
       {
	 diff = 2;
	 len  = 1;
       }
       return (double)diff / (double)len;
     }
     else if (flag == 2)
     {
       //       len = end_range - begin_range;

       while (i != end_range)
       {
	 // TODO: DO WE NEED TOUPPER??
	 c1 = toupper(seq1->get_pos(i));
	 c2 = toupper(seq2->get_pos(i));

	 // We only count differences in the gap structur: Non-gap versus gap.
	 if ( (c1 == '-' && c2 != '-') || (c1 != '-' && c2 == '-') )
	 {
	   ++diff;
	 }
	 ++i;
       }
       return diff;
     }
     else
     {
       std::cerr << "Flag value undefined in CSequences2:sequence_distance(..);" << std::endl;
       exit(-127);
     }
   }




   // Currently only implemented for datatype DNA
   // Note: N positions are not counted if the other sequence contains a gap.
   // In this case only the gap and not the N is counted, since the sequence positions
   // are not comparable.
   // The function does to compute the distance but the length, differences, gaps
   // and Ns positions.
   // Note: length_comparable_nucleotieds + gaps_extensions + Ns must add up to
   //       the length of the range.
   void sequence_distance(unsigned s1, unsigned s2, unsigned begin_range, unsigned end_range,
			  unsigned &differences,  unsigned &length_comparable_nucleotieds,
			  unsigned &gap_openings, unsigned &gaps_extensions, unsigned &Ns)
   {
     unsigned i = begin_range;
     char     c1, c2;

     if (s1 >= taxaNum || s2 >= taxaNum)
     {
       differences = -1;
       return;
     }

     CSequence_Mol* seq1 = seqData[s1];
     CSequence_Mol* seq2 = seqData[s2];

     if (end_range < begin_range)
     {
       std::cerr << "Error in sequence_distance: end_range < begin_range" << std::endl;
       exit(-33);
     }

     if ((end_range - begin_range > seq1->length()) || (end_range - begin_range > seq2->length()) )
     {
       end_range = begin_range + seq1->length();
       if (end_range - begin_range > seq2->length())
	 end_range = begin_range + seq2->length();
     }

     unsigned local_len           = 0;
     unsigned local_numgapsopen   = 0;
     unsigned local_numgapext     = 0;
     unsigned local_differences   = 0;
     unsigned local_numNs         = 0;
     bool     in_gap        = false;

     if (seq1->get_datatype() != CSequence_Mol::dna)
     {
       differences = -2;
       return;
     }


     while (i != end_range)
     {
       // TODO: DO WE NEED TOUPPER??
       c1 = toupper(seq1->get_pos(i));
       c2 = toupper(seq2->get_pos(i));

       if (c1 == '-' || c2 == '-')
       {
	 if (!in_gap)
	 {
	   ++local_numgapsopen;
	   in_gap = true;
	 }
	 ++local_numgapext;
       }
       else
       {
	 in_gap = false;
	 if (c1 == 'N' || c2 == 'N')
	   ++local_numNs;
	 else
	 {
	   ++local_len;
	   if (c1 != c2)
	     ++local_differences;
	 }
       }
       ++i;
     }

     differences                   = local_differences;
     length_comparable_nucleotieds = local_len;
     gap_openings                  = local_numgapsopen;
     gaps_extensions               = local_numgapext;
     Ns                            = local_numNs;
   }


   // Currently only implemented for datatype DNA
   //  Allowed flag values:
   // 0: normal mode.
   // 1: Gap regions always get a distance of -3.
   void get_sequences_distances(unsigned begin_range, unsigned end_range, void (*call_back_distance)(short, short, double), unsigned flag=0)
   {
     // Compute all pairwise distances:
     short N = taxaNum;
     short i, j;
     double dist;
     
     for (i=0; i<N; ++i)
     {
       for (j=i+1; j<N; ++j)
       {
	 dist = sequence_distance(i, j, begin_range, end_range, flag);
	 call_back_distance(i,j,dist);
       }
     }
   }

   // Same as above, but for some distances, those that have a false value in vector
   // we assign an invalidly large value to the distance.
   // Here a value of 100 is assigned to those distances
   void get_sequences_distances(unsigned begin_range, unsigned end_range, void (*call_back_distance)(short, short, double), const std::vector<bool> &sequences_valid, unsigned flag=0)
   {
     // Compute all pairwise distances:
     unsigned N = taxaNum;
     unsigned i, j;
     double dist;
     
     if (sequences_valid.size() != N)
     {
       std::cerr << "Critical internal error: In function get_sequences_distances the size of the vector sequences_valid must be equal to the number of sequences in the sequences object." << std::endl;
       exit(-3);
     }

     for (i=0; i<N; ++i)
     {
       for (j=i+1; j<N; ++j)
       {
	 if (!sequences_valid[i] || !sequences_valid[j])
	 {
	   dist = 100;
	 }
	 else
	 {
	   dist = sequence_distance(i, j, begin_range, end_range, flag);
	 }
	 call_back_distance(i,j,dist);
       }
     }
   }


   CSequence_Mol* get_sequence_with_sequence_name_match(faststring &partial_name)
   {
     unsigned i;

     for (i=0; i < taxaNum; ++i)
     {
       if ( seqData[i]->getFullName_faststring().find(partial_name) != faststring::npos )
       {
	 return seqData[i];
       }
     }
     return NULL;
   }
   
   void get_alignment_as_copy(std::vector<faststring> &v)
   {
     v.clear();
     v.reserve(taxaNum+1);

     unsigned i;
     for (i=0; i<taxaNum; ++i)
     {
       v.push_back(seqData[i]->getSeq_faststring());
       v[i].c_str();  // Make buffer 0-teminated.
     }
   }

};



// inline CSequences::CSequences():posNum(0),datatype(0),splits(NULL){}


//inline string    CSequences::GetTaxonLabel(unsigned i)
//                             {return taxa.GetTaxonLabel(i);}
//inline void      CSequences::SetTaxonLabel(unsigned i, string s)
//                             {taxa.SetTaxonLabel(i,s);}






/* inline void PrintSpecial_TaxaSeq(ostream& os, */
/* 			  const std::vector<unsigned>& taxa_vec, */
/* 			  const std::vector<unsigned>& pos_vec, */
/* 			  const CTaxa       *ptaxa, */
/* 			  const CSequences2 *pseq, */
/* 			  const char*      symbolsList, */
/* 			  bool       printTaxonLabels, */
/* 			  unsigned   taxonLabelPrintWidth, */
/* 			  unsigned   maxTaxonLabelLength) */
/* { */
/*   std::vector<unsigned>::const_iterator taxa_vec_it; */
/*   std::vector<unsigned>::const_iterator taxa_vec_end = taxa_vec.end(); */

/*   std::vector<unsigned>::const_iterator pos_vec_it; */
/*   std::vector<unsigned>::const_iterator pos_vec_end = pos_vec.end(); */

/*   unsigned theTaxonNumber; */

/*   std::ios_base::fmtflags orig_ios_flag = os.flags(); */
/*   os.flags(orig_ios_flag | ios::left); */

/*   for (taxa_vec_it=taxa_vec.begin(); taxa_vec_it!=taxa_vec_end; ++taxa_vec_it) */
/*   { */
/*     theTaxonNumber = *taxa_vec_it; */

/*     if (printTaxonLabels) */
/*     { */
/*       os.width(maxTaxonLabelLength); */
/*       faststring tl = ptaxa->GetTaxonLabel(theTaxonNumber); */
/*       tl.shorten(maxTaxonLabelLength); */
/*       os << tl; */
/*       int morespaces = (int)taxonLabelPrintWidth-(int)maxTaxonLabelLength; */

/*       if ( morespaces > 0 ) */
/* 	{ */
/* 	  os.width(morespaces); */
/* 	  os << " "; */
/* 	} */
/*     } */

/*     for (pos_vec_it = pos_vec.begin(); pos_vec_it != pos_vec_end; ++pos_vec_it) */
/*     { */
/*       os << symbolsList[pseq->GetChar(theTaxonNumber,*pos_vec_it)]; */
/*     } */
/*     os << std::endl; */
/*   } */
/*   os.flags(orig_ios_flag); */
/* } */

/* inline  */
/* void PrintSpecial_PosNumbers(ostream& os, */
/* 			     std::vector<unsigned> pos_vec, // we get a real copy */
/* 			     const CSequences2 *pseq, */
/* 			     unsigned         taxonLabelPrintWidth, */
/* 			     bool             UseNexusComments) */
/* { */
/*   unsigned i; */
/*   unsigned numPos           = pos_vec.size(); */
/*   unsigned printNoZerosTill = numPos; */
/*   bool     printNoZero_set  = false; */

/*   // We determine the original position numbers: */
/*   for (i=0; i<numPos; ++i) */
/*     pos_vec[i] = pseq->GetOriginalPosNumber(pos_vec[i]); */

/*   unsigned maxDigits = ((unsigned)log10((double)pos_vec[numPos-1]))+1; */
/*   unsigned factor    = 1; */
/*   unsigned digit; */

/*   for (i=1; i<maxDigits; ++i) */
/*     factor *= 10; */

/*   for (; factor; factor /= 10, printNoZero_set = false) */
/*   { */
/*     if (UseNexusComments) */
/*     { */
/*       os.width(taxonLabelPrintWidth); */
/*       os << "[ "; */
/*     } */
/*     else */
/*     { */
/*       os.width(taxonLabelPrintWidth); */
/*       os << ' '; */
/*     } */
/*     for (i=0; i < numPos; ++i) */
/*     { */
/*       digit = pos_vec[i]/factor; */

/*       if ( i >= printNoZerosTill || digit > 0) */
/*       { */
/* 	if (!printNoZero_set) */
/* 	{ */
/* 	  printNoZerosTill = i; */
/* 	  printNoZero_set = true; */
/* 	} */
/* 	os << (char)(digit+'0'); */
/* 	pos_vec[i] -= factor*digit; */
/*       } */
/*       else */
/*       { */
/* 	os << ' '; */
/*       } */
/*     } */
/*     if (UseNexusComments) */
/*     { */
/*       os << "]" << std::endl; */
/*     } */
/*     else */
/*     { */
/*       os << std::endl; */
/*     } */
/*   } */
/* } */


/* inline void      CSequences::GetBaseFrequencies(unsigned ti, */
/* 						std::vector<unsigned> bf) */
/* { */
/*   unsigned pos; */
/*   for (pos = 0; pos < posNum; ++pos) */
/*     ++bf[seqData[ti][pos]] */
/* } */


#endif
