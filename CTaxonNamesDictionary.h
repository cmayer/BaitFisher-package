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
#ifndef  CTAXONNAMESDICTIONARY_H
#define  CTAXONNAMESDICTIONARY_H

#include "faststring2.h"
#include "CDistance_matrix.h"
#include "CSequences2.h"
#include "CFile/CFile2_1.h"
#include "CSeqNameList.h"
#include <set>
#include <list>
#include <map>
#include <vector>

class CTaxonNamesDictionary
{
  CSeqNameList *snl;
  CSeqNameList all;

  // Each index storen in the taxon_names_to_index map is equal to the index of the same taxon name in the taxon_names_vec.
  // This allows us to translate indices and names back and forth.

  std::vector<faststring>        taxon_names_vec;
  std::map<faststring, unsigned> taxon_names_to_index; 

 public: 

 CTaxonNamesDictionary(faststring directoryname, std::list<faststring> &list_of_fasta_file_names, unsigned short sequence_name_taxon_field_number, char field_delimiter):all(true)
 {
    // Collect all taxon names to build a taxon name dictionary
    // The collection of all sequence names is stored in the all object. 
    faststring full_transcript_file_name;

    std::list<faststring>::iterator file_names_it, file_names_it_end, file_names_it_begin;

    file_names_it_begin = file_names_it = list_of_fasta_file_names.begin();
    file_names_it_end = list_of_fasta_file_names.end();

    // For each fasta file
    while (file_names_it != file_names_it_end)
    {
      full_transcript_file_name = directoryname + "/" + *file_names_it;

      //    std::cout << "Reading sequence names in " << full_transcript_file_name << std::endl;
      snl = new CSeqNameList(full_transcript_file_name.c_str(), true);

      if (file_names_it == file_names_it_begin)
      {
	all.set_to(*snl, "all-names-union");
      }
      else
      {
	all.add_List_non_redundant(*snl);
      }
      ++file_names_it;
    }
    //  all.print(cout,0);

    // TAXON NAMES DICTIONARY: consists of a vector, and map
    // A set is used to store all names temporarily.

    std::set<faststring>           taxon_names_set;
    bool                           error_get_set_of_name_field;

    // Take the names of the sequences, divide it according to the delimiter '|' and 
    // extract field number sequence_name_taxon_field_number. Store all names in the set.
    error_get_set_of_name_field = all.get_set_of_name_field(sequence_name_taxon_field_number+1, field_delimiter, taxon_names_set);

    if (error_get_set_of_name_field)
    {
      std::cerr << "WARNING: At least one error occured when trying to read taxon names from the fasta sequence headers."
		<< " Usually this error occurs, when headers are not in the proper format, when the parameter file contains a wrong"
		<< " wrong sequence name field delimiter or if the taxon name field number is wrong." << std::endl;
      std::cerr << "This error will influence: TODO." << std::endl;
    }

    // Now we have a set that contains all taxon names found in all fasta files.

    std::set<faststring>::iterator it, it_end;
    it      = taxon_names_set.begin();
    it_end  = taxon_names_set.end();

    unsigned i;
    i=0;
    while (it != it_end)
    {
      // Since the taxon_names_vec was empty originally, the names have an index equal to i.
      taxon_names_vec.push_back(*it); // push all names into the vector.            unique-index -> name
      taxon_names_to_index[*it] = i;  // for each name store its index in the map.  name -> unique-index
      ++it;
      ++i;
    }
    // Now we can translate unique indices and taxon names back and forth.
  }

  void print(std::ostream &outfile)
  {
    outfile << "## Index <-> taxon name:" << std::endl;

    unsigned i, N = taxon_names_vec.size();
    for (i=0; i<N; ++i)
    {
      outfile << i << " "  << taxon_names_vec[i] << std::endl;
    }
  }


  // Main dictionary retrieval function.
  unsigned dictionary_get_index_of_taxonname(faststring taxonname)
  {
/*     map<faststring, unsigned>::iterator it, it_end; */
/*     it     = taxon_names_to_index.begin(); */
/*     it_end = taxon_names_to_index.end(); */
/*     cout << "taxon_names_to_index map" << endl; */
/*     while (it != it_end) */
/*     { */
/*       cout << it->first << "=>" << it->second << endl; */
/*       ++it; */
/*     } */


    std::map<faststring, unsigned>::iterator it_find = taxon_names_to_index.find(taxonname);

    if (it_find != taxon_names_to_index.end() )
      return it_find->second;
    else
      return -1u;
  }

  const std::vector<faststring> & get_taxon_names_vec()
  {
    return taxon_names_vec;
  }

  void print_taxon_names_to_index(std::ostream &os)
  {
    std::map<faststring, unsigned>::iterator it, it_end;
    it     = taxon_names_to_index.begin();
    it_end = taxon_names_to_index.end();

    os << "Content of taxon names to global index dictionary: " << std::endl;
    while (it != it_end)
    {
      os.width(30);
      os << it->first << " " << it->second << std::endl;
      ++it;
    }
  }

};


#endif
