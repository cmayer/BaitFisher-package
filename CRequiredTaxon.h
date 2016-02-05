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
#ifndef CREQUIREDTAXA_H
#define CREQUIREDTAXA_H

#include "faststring2.h"
#include <iostream>
#include <fstream>
#include <set>
#include <vector>
#include "CTaxonNamesDictionary.h"
#include "print_container.h"

class CRequiredTaxa
{
  std::vector<std::set<faststring> > required_groups;
  std::vector<std::set<unsigned> >   required_groups_u;
  CTaxonNamesDictionary &dictionary;
  bool empty;

  bool sets_intersect(const set<unsigned> &s1, const set<unsigned> &s2) const
  {
    //   cout << "Set1" << endl;
    //   print_set(cout, s1);
    //   cout << endl;
    //   cout << "Set2" << endl;
    //   print_set(cout, s2);
    //   cout << endl;

    unsigned size1 = s1.size();
    unsigned size2 = s2.size();

    unsigned size3 = macromin(size1, size2);

    // We allocate memory for the result set that will be stored in a vector.
    // Using allocated memory is more efficient than using an insert-iterator
    // that has more overhead.
    // Important: Later we have to shorten the vector to its actual size. 
    vector<unsigned> res(size3);
    vector<unsigned>::iterator it;
    it=set_intersection (s1.begin(), s1.end(), s2.begin(), s2.end(), res.begin());

/*     cout << "SET1: "; */
/*     print_container(cout, s1, "(", ",", ")\n"); */
/*     cout << "SET2: "; */
/*     print_container(cout, s2, "(", ",", ")\n"); */
/*     cout << "SET-RES1: "; */
/*     print_container(cout, res, "(", ",", ")\n");     */
/*     cout << "SIZE of intersection (A): " << res.size() << endl; */

    // Shorten the vector to the actual size of the intersection.
    // Otherwise we have trailing 0s from the constructor of the vector
    // and the size of res would always be size3.
    res.resize(it-res.begin());
/*     cout << "SIZE of intersection (B): " << res.size() << endl; */
    
/*     cout << "SET-RES2: "; */
/*     print_container(cout, res, "(", ",", ")\n");     */

    return (res.size()>0);
  }
  
  // Do not allow gaps:
  void parse_required_groups(faststring line, int &error)
  {
    if (line.empty())
    {
      cout << "LOG: Required taxa string is empty. No required taxa specified." << endl;
      empty = true;
    }

    //    dictionary.print(cout);

    set<faststring> empty_set;
    set<unsigned> empty_set_u;

    unsigned i1, i2, N1, N2;

    vector<faststring> splitter1;
    vector<faststring> splitter2;
  
    split_respect(splitter1, line, ", ");

    N1 = splitter1.size();

    for (i1=0; i1<N1; ++i1)
    {
      required_groups.push_back(empty_set);
      required_groups_u.push_back(empty_set_u);

      //    cout << splitter1[i] << endl;
      split(splitter2, splitter1[i1], "(, )");
      N2 = splitter2.size();
      for (i2=0; i2<N2; ++i2)
      {
	splitter2[i2].removeSpacesBack();
	splitter2[i2].removeSpacesFront();
	required_groups[i1].insert(splitter2[i2]);
	unsigned index = dictionary.dictionary_get_index_of_taxonname(splitter2[i2]);
	if (index == -1u)
	{
	  std::cerr << "ERROR: Could not find taxon name " << splitter2[i2]
		    << " in the list of all taxon names found in all given sequence files." << std::endl;
	  std::cerr << "Please check that the taxon name is spelled properly or remove the taxon name from the list of required"
		    << " taxa to proceed. Exiting." << std::endl;
	  error = 3;
	}
	required_groups_u[i1].insert(index);
      }
    }
  }


 public:
  // Error codes:
  // 0: No error,
  // 1: File does not contain a single line as required. We might have more or less than one line in the file.
  // 3: Sequence names not found in dictionary 

  // The dummy parameter is only used to distinguish this and the following constructor, since both have the same parameter list
  // except of this dummy parameter.
  /* Depricates. Required taxa have to be specified in the parameter file.
 CRequiredTaxa(faststring infilename, CTaxonNamesDictionary &dict, int &error, int dummy):dictionary(dict),empty(false)
  {
    std::ifstream is(infilename.c_str());

    error = 0;
    std::vector<faststring> file_as_vec;

    appendFileToVector(is, file_as_vec);

    unsigned i;
    i=0;
    while (i< file_as_vec.size() )
    {
      file_as_vec[i].removeSpacesFront();
      file_as_vec[i].removeSpacesBack();
      if (file_as_vec[i].empty() )
      {
	file_as_vec.erase(file_as_vec.begin()+i);
      }
      else
      {
	++i;
      }
    }
    if (file_as_vec.size() != 1)
    {
      error = 1;
      return;
    }
    parse_required_groups(file_as_vec[0], error);
  }
  */

   // Error codes:
  // 0: No error,
  // 1: -,  unused
  // 3: Sequence names not found in dictionary 
 CRequiredTaxa(faststring required_string, CTaxonNamesDictionary &dict, int &error):dictionary(dict),empty(false)
  {
    error = 0;
    parse_required_groups(required_string, error);
  }
 
  bool is_empty() const
  {
    return empty;
  }

  bool check_precense_of_all_groups_of_taxa(std::vector<faststring> taxon_vec) const
  {
    if (empty)
      return true;

    std::set<unsigned> taxon_set;

    unsigned i, N, index;
    for (i=0, N=taxon_vec.size(); i < N; ++i)
    {
      index = dictionary.dictionary_get_index_of_taxonname(taxon_vec[i]);
      taxon_set.insert(index);
    }

    // v_set must have a non vanishing intersection with all groups of taxa that are required.
    
    unsigned iset, Nset;
    bool all_taxa_present = true;

    for (iset=0, Nset=required_groups_u.size();iset < Nset; ++iset)
    {
      if ( !sets_intersect(taxon_set, required_groups_u[iset] ) )
      {
	//		cout << " Sets do not intersect: " << endl;
	//		print_set(cout, taxon_set);                cout << endl;
	//		print_set(cout, required_groups_u[iset]);  cout << endl;
	all_taxa_present = false;
	break;
      }
      
    }
    return all_taxa_present;
  }

  bool check_precense_of_all_groups_of_taxa(std::set<unsigned> taxon_ids) const
  {
    if (empty)
      return true;

    // taxon_ids must have a non vanishing intersection with all groups of taxa that are required.
    
    unsigned iset, Nset;
    bool all_taxa_present = true;

    for (iset=0, Nset=required_groups_u.size();iset < Nset; ++iset)
    {
      if ( !sets_intersect(taxon_ids, required_groups_u[iset] ) )
      {
	//		cout << " Sets do not intersect: " << endl;
	//		print_set(cout, taxon_set);                cout << endl;
	//		print_set(cout, required_groups_u[iset]);  cout << endl;
	all_taxa_present = false;
	break;
      }
    }
    return all_taxa_present;
  }

  


  void print(std::ostream &os, int mode = 0)
  {
    std::vector<std::set<faststring> >::iterator it_vec_faststrings;
    std::vector<std::set<unsigned> >::iterator   it_vec_unsigned;

    std::set<faststring>::iterator it_set_faststrings;
    std::set<unsigned>::iterator   it_set_unsigned;

    it_vec_faststrings = required_groups.begin();

    if (empty)
    {
      os << "CRequiredTaxa: empty" << endl;
      return;
    }

    while (it_vec_faststrings != required_groups.end())
    {
      os << "(";
      it_set_faststrings = it_vec_faststrings->begin();
      while (it_set_faststrings != it_vec_faststrings->end() )
      {
	os << *it_set_faststrings;
	++it_set_faststrings;
	if (it_set_faststrings != it_vec_faststrings->end())
	  os << ",";
      }
      ++it_vec_faststrings;
      os << "),";
      if (mode == 1)
	os << endl;
    }

    it_vec_unsigned  = required_groups_u.begin();

    os << "*********" << endl;

    while (it_vec_unsigned != required_groups_u.end())
    {
      os << "(";
      it_set_unsigned = it_vec_unsigned->begin();
      while (it_set_unsigned != it_vec_unsigned->end() )
      {
	os << *it_set_unsigned;
	++it_set_unsigned;
	if (it_set_unsigned != it_vec_unsigned->end())
	  os << ",";
      }
      ++it_vec_unsigned;
      os << "),";
      if (mode == 1)
	os << endl;
    }

  }
  

};

#endif
