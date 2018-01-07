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
// Uses modified comparison functions to make sort stable.

#include <iostream>
#include <fstream>
#include "faststring2.h"
#include <vector>
#include <set>
#include "Ctriple.h"
#include "print_container.h"
#include "CBlastParser.h"

// Introducing this define has several advantages:
// - simple switch to stable_sort and back
// - simple find all calls to sort in this file.
#define MYSORT  sort
// Instead of using stable_sort we now use modified comparison functions. So this switch should not be necessary any more.
// #define MYSORT  stable_sort

class CBaitRegionSelectorSets
{
  std::set<faststring>                               set_of_alignment_names;
  std::set<std::pair<faststring, unsigned> >         set_of_alignment_names_feature_num;
  std::set<Ctriple<faststring, unsigned, unsigned> > set_of_alignment_name_feature_num_start__ie_locus;

 public:
  void add(const faststring &str, unsigned i, unsigned j)
  {
    set_of_alignment_names.insert(str);
    set_of_alignment_names_feature_num.insert(std::make_pair(str, i));
    set_of_alignment_name_feature_num_start__ie_locus.insert(Ctriple<faststring, unsigned, unsigned>(str,i, j));
  }


  bool exists(const faststring &str) const
  {
    if (set_of_alignment_names.find(str) != set_of_alignment_names.end() )
      return true;
    else
      return false;
  }

  bool exists(const faststring &str, unsigned i) const
  {
    if (set_of_alignment_names_feature_num.find(std::make_pair(str, i)) != set_of_alignment_names_feature_num.end() )
      return true;
    else
      return false;
  }

  bool exists(const faststring &str, unsigned i, unsigned j) const
  {
    if (set_of_alignment_name_feature_num_start__ie_locus.find(Ctriple<faststring, unsigned, unsigned>(str,i,j)) != set_of_alignment_name_feature_num_start__ie_locus.end() )
      return true;
    else
      return false;
  }

  void print_all(std::ostream &os) const
  {
    os << "Selector string (alignment) argument:" << std::endl;
    print_container(os, set_of_alignment_names, "", "\n" , "\n");
    os << "Selector string (alignment) argument, unsigned (feature number): " << std::endl;
    print_container(os, set_of_alignment_names_feature_num, "", "\n" , "\n");
    os << "Selector string (alignment) argument, unsigned (feature number), unsigned (start): " << std::endl;
    print_container(os, set_of_alignment_name_feature_num_start__ie_locus, "", "\n" , "\n");
  }

  void print_num_genes_num_features(std::ostream &os, const char *s, bool print_num_featues) const
  {
    os << s << "The bait regions belong to this number of alignment files: " << set_of_alignment_names.size() << std::endl;
    os << s << "Number of bait regions:                                    " << set_of_alignment_name_feature_num_start__ie_locus.size() << std::endl;
    if (print_num_featues)
      os << s << "Number of features                                       " << set_of_alignment_names_feature_num.size() << std::endl;
  }
};

class CBait
{
  faststring bait;
  double     CG;
  double     max_dist;
  unsigned   tiling_index;  // 1 based index in tiling array. 1: first tile; 2: second tile, etc.

 public:

  CBait(faststring &str)
  {
    faststring pre, post1, post2;

    // The error messages below could be improved. It would be a good idea to tell the user the line in the file which causes the error.
    // At the moment the user is being told the bait string or the number that causes the problem. While this should help in most cases to find the problem
    // it would be much more convenient for the user to have the line number in the file.

    str.divide_at(' ', pre, post2);
    if (!pre.isAnUnsigned())
    {
      std::cerr << "Error when parsing the bait file: Has it been edited manually?" << std::endl;
      std::cerr << "Before each bait string BaitFilter expects a number indicating the index of the tiling design column. This number could not be found for the following bait:"
		<< std::endl;
      std::cerr << str << std::endl;
      std::cerr << "Exiting." << std::endl;
      exit(-4);
    }
    tiling_index = pre.ToUnsigned();
    post2.divide_at(' ', pre, post1);
    bait = pre;
    post1.divide_at(' ', pre, post2);
    if (!pre.isADouble())
    {
      std::cerr << "Error when parsing the bait file: Has it been edited manually?" << std::endl;
      std::cerr << "After each bait string BaitFilter expects two double numbers containing the baits' CG content and its maximum distance to the MSA.\n"
		<< "Expecting a double number for the CG content but found: " << pre 
		<< std::endl;
      std::cerr << "Exiting." << std::endl;
      exit(-4);
    }
    if (!post2.isADouble())
    {
      std::cerr << "Error when parsing the bait file: Has it been edited manually?" << std::endl;
      std::cerr << "After each bait string BaitFilter expects two double numbers containing the baits' CG content and its maximum distance in the MSA.\n"
		<< "Expecting a double number for the distance value but found: " << post2 
		<< std::endl;
      std::cerr << "Exiting." << std::endl;
      exit(-4);
    }

    CG       = pre.ToDouble();
    max_dist = post2.ToDouble();
  }
  
  faststring& get_bait()
  {
    return bait;
  }

  unsigned get_tiling_index()
  {
    return tiling_index;
  }

  void print(std::ostream &os)
  {
    os << tiling_index  << " " << bait << " " << CG << " " << max_dist;
  }

  double get_max_distance()
  {
    return max_dist;
  }

  double get_CG_content()
  {
    return CG;
  }

}; // END CBait

// Forward declaration of friend functions of the CBaitLocus class. Avoids friend injection problems:
class CBaitLocus;
bool num_baits_lessThanCBaitLocus(CBaitLocus *a, CBaitLocus *b);
bool num_sequences_in_region_lessThanCBaitLocus(CBaitLocus *a, CBaitLocus *b);
bool gene_name__feature_num_lessThanCBaitLocus(CBaitLocus *a, CBaitLocus *b);

bool num_baits_greaterThanCBaitLocus(CBaitLocus *a, CBaitLocus *b);
bool num_sequences_in_region_greaterThanCBaitLocus(CBaitLocus *a, CBaitLocus *b);
bool gene_name__feature_num_greaterThanCBaitLocus(CBaitLocus *a, CBaitLocus *b);

bool gene_name__feature_num_lessThanOrEqualCBaitLocus(CBaitLocus *a, CBaitLocus *b);

// bool original_file_name__feature_num_start_lessThanCBaitLocus(CBaitLocus *a, CBaitLocus *b);

bool locus_start_lessThanCBaitLocus(CBaitLocus *a, CBaitLocus *b);


// Information about a bait region, i.e. a locus with a tiling design and its baits.
// Basically this is the information stored in one line of a BaitFisher output file.

class CBaitLocus
{
  faststring original_file_name;
  faststring gene_name;
  unsigned   feature_num;
  unsigned   num_features;
  unsigned   start;
  unsigned   num_sequences_in_region;
  unsigned   num_baits;
  unsigned   bait_length;
  unsigned   number_of_tiles;

  std::vector<CBait*> bait_vec;

 public:

  CBaitLocus(faststring &str)
  {
    std::vector<faststring> container;
    container.reserve(7);
    split(container, str, "\t|#");
    //split(container, str, "\t");

    if (container.size() != 8)
    {
      std::cerr << std::endl;
      std::cerr << "An error occurred while parsing the bait file. The line of input that causes the problem is:\n";
      std::cerr << str << std::endl;
      std::cerr << "The line does not contain 8 tab delimited fields as required." << std::endl;
      exit(-1);
    }

    container[0].removeSpacesFront();
    container[0].removeSpacesBack();
    container[1].removeSpacesFront();
    container[1].removeSpacesBack();
    container[2].removeSpacesFront();
    container[2].removeSpacesBack();
    container[3].removeSpacesFront();
    container[3].removeSpacesBack();
    container[4].removeSpacesFront();
    container[4].removeSpacesBack();
    container[5].removeSpacesFront();
    container[5].removeSpacesBack();
    container[6].removeSpacesFront();
    container[6].removeSpacesBack();
    container[7].removeSpacesFront();
    container[7].removeSpacesBack();

    original_file_name      = container[0];
    gene_name               = container[1];
    feature_num             = container[2].ToUnsigned();
    num_features            = container[3].ToUnsigned();
    start                   = container[4].ToUnsigned();
    num_sequences_in_region = container[5].ToUnsigned();
    num_baits               = container[6].ToUnsigned();

    std::list<faststring> b;
    split(b, container[7], ",");

    std::list<faststring>::iterator it, it_end;

    CBait *p_tmp;

    it     = b.begin();
    it_end = b.end();

    number_of_tiles = 0;
    // For all baits in this bait region:
    while(it != it_end)
    {
      p_tmp = new CBait(*it);
      bait_vec.push_back(p_tmp);
      if (p_tmp->get_tiling_index() > number_of_tiles)
	number_of_tiles = p_tmp->get_tiling_index();

/*       if (0) */
/*       { */
/* 	std::cerr << "Adding bait:" << std::endl; */
/* 	bait_vec.back()->print(std::cerr); */
/* 	std::cerr << std::endl; */
/*       } */

      ++it;
    }
    if (bait_vec.size() != 0)
    {
      bait_length = bait_vec[0]->get_bait().size();
    }
    else
    {
      bait_length = 0;
    }
  }

  ~CBaitLocus()
  {
    std::vector<CBait *>::iterator it, it_end;
    it     = bait_vec.begin();
    it_end = bait_vec.end();

    while (it != it_end)
    {
      delete *it;
      ++it;
    }
  }

  bool has_feature_information()
  {
    return (num_features == 0);
  }

  faststring &get_original_file_name()
  {
    return original_file_name;
  }

  const faststring &get_gene_name() const
  {
    return gene_name;
  }

  void get_vector_of_baits(unsigned tiling_index, std::vector<CBait*> &v)
  {
    std::vector<CBait*>::iterator it, it_end;
    it     = bait_vec.begin();
    it_end = bait_vec.end();

    v.clear();
    v.reserve(bait_vec.size());

    if (tiling_index == 0) // all
    {
      while (it != it_end)
      {
	v.push_back(*it);
	++it;
      }
    }
    else
    {
      while (it != it_end)
      {
	if ((**it).get_tiling_index() == tiling_index)
	{
	  v.push_back(*it);
	}
	++it;
      }
    }
  }

  unsigned get_feature_num() const
  {
    return feature_num;
  }

  unsigned get_start() const
  {
    return start;
  }
  
  unsigned get_number_of_tiles() const
  {
    return number_of_tiles;
  }

  unsigned get_bait_length() const
  {
    return bait_length;
  }

  faststring get_bait_region_name_string()
  {
    return original_file_name + "|" + gene_name + "|" + faststring(feature_num) + "|" + faststring(start);
  }

  // TODO: Change this to FILE * instead of using ostream.
  // Known formats are 'f', 's' for the 2 parameter version of this function.
  void print(std::ostream &os, char format)
  {
    std::vector<CBait *>::iterator i1, i2, i3;
    //    faststring name;
    unsigned count = 1;

    i1 = bait_vec.begin();
    i2 = bait_vec.end();

    if (format == 'f') // fasta
    {
      //      name = ">" + original_file_name + "|" + gene_name + "|" + faststring(feature_num) + "|" + faststring(start);

      while (i1 != i2)
      {
	os << ">" << get_bait_region_name_string()  << "|" << (**i1).get_tiling_index(); 
	os         << "|" << count    << std::endl;
	os << (**i1).get_bait()       << std::endl;
	++i1;
	++count;
      }
      
    } // END if (format == 'f') // fasta
    else if (format == 's') // standard - same as input
    {
      i3 = i2;
      --i3;
      os << original_file_name << "\t"
	 << gene_name << "\t"
	 << feature_num << "\t"
	 << num_features  << "\t"
	 << start  << "\t"
	 << num_sequences_in_region  << "\t"
	 << num_baits << "\t";

      // Works for 0, 1, or more elements in vector:
      if (i1 != i2) // Should always be the case, but be cautious
      {
	while (true)
	{
	  
	  (**i1).print(os);
	  if (i1 == i3)
	    break;
	  os << ",";
	  ++i1;
	  ++count;
	}
      }
      os << std::endl;
    }
  }

  // TODO: Change this to FILE * instead of using ostream.
  // Known formats are 'A' for the 3 parameter version of this function.
  void print(std::ostream &os, char format, const char *probeIDprefix)
  {
    std::vector<CBait *>::iterator i1, i2, i3;
    //    faststring name;

    i1 = bait_vec.begin();
    i2 = bait_vec.end();

    if (format == 'A') // Agilent-Germany format:
    {
      unsigned t_index, t_index_old = 0;
      unsigned count = 1;

      while (i1 != i2)
      {
	t_index = (**i1).get_tiling_index();
	
	if (t_index != t_index_old)
	{
	  t_index_old = t_index;
	  count = 1;
	}

	os << gene_name     << "_" << feature_num << '\t'
	   << probeIDprefix
	   << gene_name     << "_" 
	   << feature_num   << "_" 
	   << num_features  << "_"
	   << start         << "_"
	   << t_index       << "_"
	   << count                              << '\t' 
	   << (**i1).get_bait()
	   << "\t1"        << std::endl;
	++i1;
	++count;
      }
    }
  }



  void print_stats(std::ostream &os)
  {
    os << "Locus: " << gene_name << " " << feature_num << " " << start << " Num baits: " << bait_vec.size() << " bait-length :" << bait_length << std::endl;
  }

  unsigned get_num_baits()
  {
    return num_baits;
  }

  unsigned get_num_sequences()
  {
    return num_sequences_in_region;
  }

  // Can be used to compute mean values together with get_num_baits()
  double get_sum_max_distances()
  {
    int     i, N=num_baits;
    double  sum=0;

    for (i=0; i<N; ++i)
    {
      sum += bait_vec[i]->get_max_distance();
    }
    return sum;
  }

  double get_overall_max_dist()
  {
    int     i, N = num_baits;
    double  overall_max  = -1000;  // True values must be in the range 0..1
    double  dist;

    for (i=0; i<N; ++i)
    {
      dist = bait_vec[i]->get_max_distance();
      if (dist > overall_max)
	overall_max = dist;
    }
    return overall_max;
  }

  // Can be used to compute mean values together with get_num_baits()
  double get_sum_CG()
  {
    int     i, N=num_baits;
    double  sum=0;

    for (i=0; i<N; ++i)
    {
      sum += bait_vec[i]->get_CG_content();
    }
    return sum;
  }


  friend bool num_baits_lessThanCBaitLocus(CBaitLocus *a, CBaitLocus *b)
  {
    // First criterion:
    if (a->num_baits != b->num_baits)
      return a->num_baits < b->num_baits;

    // Second criterion: (opposite to first, therefore we use > and not <:
    if (a->num_sequences_in_region != b->num_sequences_in_region)
      return a->num_sequences_in_region > b->num_sequences_in_region;

    // Otherwise we sort by starting position. This makes the program output independent of
    // different versions of the sort algorithm. They are guaranteed to differ in different alignments files.
    return a->start < b->start;
  }

  // Not used:
  friend bool num_sequences_in_region_lessThanCBaitLocus(CBaitLocus *a, CBaitLocus *b)
  {
    // For loci based on an equal number of sequences, we prefer the one with less needed baits.
    // Note: This comparator sorts in the order greater is better, less is worse.
    // Better loci come later in a sorted list.

    // First criterion:
    if (a->num_sequences_in_region != b->num_sequences_in_region)
      return a->num_sequences_in_region < b->num_sequences_in_region;

    // Second criterion: (opposite to first, therefore we use > and not <:
    if (a->num_baits != b->num_baits)
      return a->num_baits > b->num_baits;
   
    // Otherwise we sort by starting position. This makes the program output independent of
    // different versions of the sort algorithm. They are guaranteed to differ in different alignments files.
    return a->start < b->start;
  }

  // First criterion:  gene_name
  // Second criterion: feature num
  // Third  criterion: locus start coordinate
  friend bool gene_name__feature_num_lessThanCBaitLocus(CBaitLocus *a, CBaitLocus *b)
  {
    // First criterion
    if (a->gene_name != b->gene_name)
      return a->gene_name < b->gene_name;

    // Second criterion
    if (a->feature_num != b->feature_num)
      return a->feature_num < b->feature_num;
    
    // Third criterion
    return a->start < b->start;
  }

/*   // not used: */
/*   friend bool num_baits_greaterThanCBaitLocus(CBaitLocus *a, CBaitLocus *b) */
/*   { */
/*    // First criterion: */
/*     if (a->num_baits != b->num_baits) */
/*       return a->num_baits > b->num_baits; */

/*     // Second criterion: (for the sequence coverage, more is better , therefore we use > and not <: */
/*     if (a->num_sequences_in_region != b->num_sequences_in_region) */
/*       return a->num_sequences_in_region < b->num_sequences_in_region; */

/*     // Otherwise we sort by starting position. This makes the program output independent of */
/*     // different versions of the sort algorithm. They are guaranteed to differ in different alignments files. */
/*     return a->start < b->start; */
/*   } */

  friend bool num_sequences_in_region_greaterThanCBaitLocus(CBaitLocus *a, CBaitLocus *b)
  {
    // First criterion: more sequence (coverage)
    if (a->num_sequences_in_region != b->num_sequences_in_region)
      return a->num_sequences_in_region > b->num_sequences_in_region;

    // Second criterion: less baits, opposite comparison
    if (a->num_baits != b->num_baits)
      return a->num_baits < b->num_baits;

   // Otherwise we sort by starting position. This makes the program output independent of
    // different versions of the sort algorithm. They are guaranteed to differ in different alignments files.
    return a->start < b->start;
  }

  friend bool gene_name__feature_num_greaterThanCBaitLocus(CBaitLocus *a, CBaitLocus *b)
  {
    if (a->gene_name != b->gene_name)
      return a->gene_name > b->gene_name;

    if (a->feature_num != b->feature_num)
      return a->feature_num > b->feature_num;
    
    return a->start > b->start;
  }


  friend bool locus_start_lessThanCBaitLocus(CBaitLocus *a, CBaitLocus *b)
  {
    // First criterion:
    if (a->gene_name != b->gene_name)
      return a->gene_name < b->gene_name;

    // Second criterion:
    return a->start < b->start;
  }

/*   friend bool original_file_name__feature_num_start_lessThanCBaitLocus(CBaitLocus *a, CBaitLocus *b) */
/*   { */
/*     // First criterion: */
/*     if (a->original_file_name != b->original_file_name) */
/*       return a->original_file_name < b->original_file_name; */

/*     // Second criterion: */
/*     if (a->feature_num != b->feature_num) */
/*       return a->feature_num < b->feature_num; */
    
/*     // Third criterion: */
/*     return a->start < b->start; */
/*   } */

}; // END CBaitLocus



// Contains usually multiple loci of bait regions. 
// Usually, loci come from multiple alignments, genes and features.

class CBaitLoci_collection
{
  std::vector<CBaitLocus*>                                       baitLoci_vec;
  std::map<Ctriple<faststring, unsigned, unsigned>, CBaitLocus*> baitLoci_map;
  bool                                                           delete_in_destructor;
  bool                                                           has_feature_information_bool;

  //  unsigned                                                       bait_length; // Will be set to 0 in constructor and will be adjusted to its true value in the add method.
  //                                                                              // The fact that we have one value for a collection of baits implies, that we assume that all baits in all loci have the same length.

  // Internal add 
  void add_locus(CBaitLocus *p)
  {
    baitLoci_vec.push_back(p);
    baitLoci_map.insert(std::pair<Ctriple<faststring, unsigned, unsigned>, CBaitLocus*>(Ctriple<faststring, unsigned, unsigned>(p->get_original_file_name(), p->get_feature_num(), p->get_start()), p) );

    if (p->has_feature_information())
      has_feature_information_bool = true;
    //    unsigned new_bait_length = p->get_bait_length();

/*     if (bait_length == 0) */
/*       bait_length = new_bait_length; */
/*     else */
/*     { */
/*       if (bait_length != new_bait_length) */
/*       { */
/* 	std::cerr << "ERROR: It is assumed that all baits in all bait loci have identical lenghts." << std::endl;  */
/* 	exit(-4); */
/*       } */
/*     } */
  }


 public:
 CBaitLoci_collection(std::ifstream &is):delete_in_destructor(true), has_feature_information_bool(false) //, bait_length(0)
  {
    faststring line;

    getline(is, line);
    while (!is.eof())
    {
      line.removeSpacesBack();
      line.removeSpacesFront();
      if (line[0] != '#' && !line.empty() )
	add_locus(new CBaitLocus(line) );
      getline(is, line);
    }
    if (!is_partially_ordered() )
    {
      std::cerr << "Note: Bait loci in input file are not ordered or genes appear twice." << std::endl;
      std::cerr << "Note: Lets try to order them." << std::endl;
      sort_by_gene_name__feature_num(baitLoci_vec.begin(), baitLoci_vec.end());
    }
  }


 CBaitLoci_collection(CBaitLoci_collection &source, char genes_or_features, char num_baits_or_sequences):delete_in_destructor(false), has_feature_information_bool(false) //, bait_length(0)
  {
/*     { */
/*       std::cout << "++++++++++++++" << std::endl; */
/*       std::cout << "Full set: DEBUG OUT" << std::endl; */
/*       source.print_all(std::cout, 's'); */
/*       std::cout << "++++++++++++++" << std::endl; */
/*     } */


    if (!source.is_partially_ordered() )
    {
      //      std::cout << "Source loci have to be sorted partially:" << std::endl;
      source.sort_within_genes(gene_name__feature_num_lessThanCBaitLocus);
      if (!source.is_partially_ordered() )
	if (global_verbosity >= 1)
	{
	  std::cerr << "WARNING: Sorting NOT successful." << std::endl;
	}
/*       else */
/* 	std::cout << "Sorting was successful." << std::endl; */
    }

/*     { */
/*       std::cout << "**************" << std::endl; */
/*       std::cout << "Full set: DEBUG OUT" << std::endl; */
/*       source.print_all(std::cout, 's'); */
/*       std::cout << "**************" << std::endl; */
/*     } */

    // Sort such the the best loci is first in each "group" of loci.
    if (genes_or_features == 'g')
    {
      //      std::cout << "Sorting within genes" << std::endl;
      if (num_baits_or_sequences == 'b')
	source.sort_within_genes(num_baits_lessThanCBaitLocus);
      else if (num_baits_or_sequences == 's')
	source.sort_within_genes(num_sequences_in_region_greaterThanCBaitLocus);
    }
    else if (genes_or_features == 'f')
    {
      //      std::cout << "Sorting within features" << std::endl;
      if (num_baits_or_sequences == 'b')
	source.sort_within_featues(num_baits_lessThanCBaitLocus);
      else if (num_baits_or_sequences == 's')
	source.sort_within_featues(num_sequences_in_region_greaterThanCBaitLocus);
    }

/*     { */
/*       std::cout << "##############" << std::endl; */
/*       std::cout << "Full set: DEBUG OUT" << std::endl; */
/*       source.print_all(std::cout, 's'); */
/*       std::cout << "##############" << std::endl; */
/*     } */

    // Add best to this:
    std::vector<CBaitLocus*>::iterator it, it_end, it_tmp;
    it     = source.baitLoci_vec.begin();
    it_end = source.baitLoci_vec.end();

    if (genes_or_features == 'g')
    {
      while (it != it_end)
      {
	add_locus(*it);
	//	baitLoci_vec.push_back(*it);
	it = source.find_end_this_gene(it);
      }
    }
    else if (genes_or_features == 'f')
    {
      while (it != it_end)
      {
	//	std::cout << "Next it:" << std::endl;
	//	(**it).print(std::cout, 's');
	//	std::cout << "==" << std::endl;
	add_locus(*it);
	//	baitLoci_vec.push_back(*it);
	it = source.find_end_this_feature(it);
      }
    }
    // We restore the original sorting of the source CBaitLoci_collection
    sort_within_genes(gene_name__feature_num_lessThanCBaitLocus);
  } // END CBaitLoci_collection(CBaitLoci_collection &source, char genes_or_features, char num_baits_or_sequences)


 CBaitLoci_collection(const CBaitLoci_collection &source, const CBaitRegionSelectorSets &lss, char mode):delete_in_destructor(false), has_feature_information_bool(false) // , bait_length(0)
 {
   unsigned i, N;
   N=source.baitLoci_vec.size();

   if (mode == 'b')      // Remove all loci from ALIGNMENT if a double hit occurs in this set
   {
     for (i=0; i<N; ++i)
     {
       CBaitLocus *p = source.baitLoci_vec[i];

       if (!lss.exists(p->get_original_file_name()) )
       {
	 add_locus(p);
	 //	 baitLoci_vec.push_back(p);
       }
     }
   }
   else if (mode == 'B') // Remove all loci from FEATURE if a double hit occurs in this set
   {
     for (i=0; i<N; ++i)
     {
       CBaitLocus *p = source.baitLoci_vec[i];

       if (!lss.exists(p->get_original_file_name(), p->get_feature_num()) )
       {
	 add_locus(p);
       }
     }
   }
   if (mode == 'x')      // Remove Locus if a double hit occurs in it
   {
     for (i=0; i<N; ++i)
     {
       CBaitLocus *p = source.baitLoci_vec[i];

       if (!lss.exists(p->get_original_file_name(), p->get_feature_num(), p->get_start()) )
       {
	 add_locus(p);
       }
     }
   }
   else if (mode == 'C') // In this mode we simply copy the loci to this collector.
   {
     for (i=0; i<N; ++i)
     {
       CBaitLocus *p = source.baitLoci_vec[i];
       add_locus(p);
     }
   }
   else
   {
     std::cerr << "Error: Unknown mode in constructor: CBaitLoci_collection(const CBaitLoci_collection &source, const CBaitRegionSelectorSets &lss, char mode): " << std::endl;
     std::cerr << "Bad mode character: " << mode << std::endl;
     exit(-1);
   }
 }


 CBaitLoci_collection(const CBaitLoci_collection &source, CBlast_parser &bp, double minimum_coverage):delete_in_destructor(false), has_feature_information_bool(false) //, bait_length(source.bait_length)
  {
    unsigned    tile_index, number_of_tiles;
    faststring  recognition_string;
    double      best_coverage_this_tiling_stack;
    bool        bait_region_has_tiling_stack_with_insufficient_coverage;
    int         max_hit_length, bait_length_this_bait_region;
    unsigned    count, N;

    std::vector<CBaitLocus*>::const_iterator it, it_end;

    it     = source.baitLoci_vec.begin();
    it_end = source.baitLoci_vec.end();

    N      = source.baitLoci_vec.size();
    count  = 0;

    // for each bait region
    while (it != it_end)
    {
      ++count;
      number_of_tiles                 = (**it).get_number_of_tiles();
      bait_length_this_bait_region    = (**it).get_bait_length();;
/*             std::cout << "Number of tiles for: " << (**it).get_bait_region_name_string() */
/* 		<< " " << number_of_tiles << std::endl; */

      if (bait_length_this_bait_region == 0)
      {
	// We have a problem. Something seems to have gone wrong. This should not happen.
	// We catch this potential error, since it leads to division by 0 and a crash would be expected.
	std::cerr << "The bait length has been determined to be 0 for the bait region: " << std::endl;  
	(**it).print(std::cerr, 's');
	std::cerr << "This bug should be reported. It should not occur. It might be caused by malformatted input files,"
	             " but the error should be caught earlier than here." << std::endl;
	exit(-5);
      }

      if (global_verbosity >= 3)
      {
	std::cout << "Handling region: " << count << " outof " << N << std::endl;
      }

      // Within each bait region determine the maximum hit coverage of all baits
      // For each tiling stack of baits:

      // Each tiling stack has to be checked individually.
      // A single tiling stack with insufficient query coverage should result in the exclusion of this bait region.
      bait_region_has_tiling_stack_with_insufficient_coverage = false;
      for (tile_index=1; tile_index <= number_of_tiles; ++tile_index)
      {
	// Do something with (**it):
	// Determine string to identify blast hits of baits of this bait region and tiling_index against the genome. 
	recognition_string = (**it).get_bait_region_name_string() + "|" + tile_index + "|";  

	//	std::cout << "T1: " << time(NULL) << std::endl;
	max_hit_length = bp.get_max_hit_query_span__matching_beginning2(recognition_string, global_verbosity);
	//	std::cout << "T2: " << time(NULL) << std::endl;

	best_coverage_this_tiling_stack = (double) max_hit_length/bait_length_this_bait_region;

	if (global_verbosity >= 50)
	{
	  std::cout << "Tile: " << tile_index << " " << recognition_string << " " << max_hit_length << "/" << bait_length_this_bait_region << "=" 
                    << best_coverage_this_tiling_stack << std::endl;
	}

	if (best_coverage_this_tiling_stack < 0)
	{
	  if (global_verbosity >= 50)
	  {
	    std::cout << "No hit for recognition string: " << recognition_string << " --> "
		      << 0.00 << std::endl;
	  }
	}
	else
	{
	  if (global_verbosity >= 50)
	  {
	    std::cout << "Determined best coverage of this tiling design column: " << recognition_string << " --> "
		      << best_coverage_this_tiling_stack << std::endl;
	  }
	}

	if (best_coverage_this_tiling_stack < minimum_coverage)
	  bait_region_has_tiling_stack_with_insufficient_coverage = true;
      }

      if (!bait_region_has_tiling_stack_with_insufficient_coverage)
      {
	add_locus(*it);
      }
      ++it;
    }
  }

  ~CBaitLoci_collection()
  {
    if (delete_in_destructor)
    {
      std::vector<CBaitLocus*>::iterator it, it_end;
      it     = baitLoci_vec.begin();
      it_end = baitLoci_vec.end();

      while (it != it_end)
      {
	delete *it;
	++it;
      }
    }
  }

  bool has_feature_information() const
  {
    return has_feature_information_bool;
  }

  void get_map_of_baits_to_locus_info(faststring original_file_name_match, unsigned tiling_index, std::map<faststring, CBaitLocus *> &bait_map)
  {
    //    sort_by_original_file_name__feature_num_start();

    std::vector<CBaitLocus*>::iterator it_Locus, it_Locus_end;
    it_Locus     = baitLoci_vec.begin();
    it_Locus_end = baitLoci_vec.end();

    std::vector<CBait *> tmp_bait_vec;

    unsigned i, N;

    std::map<faststring, CBaitLocus *>::iterator find_it;

    while (it_Locus != it_Locus_end)
    {
      if ( (**it_Locus).get_original_file_name().find(original_file_name_match) != faststring::npos)
      {
	(**it_Locus).get_vector_of_baits(tiling_index, tmp_bait_vec);

	N=tmp_bait_vec.size();
	for (i=0; i<N; ++i)
	{
	  faststring &thebait = tmp_bait_vec[i]->get_bait();
	  find_it = bait_map.find(thebait);
	  if (find_it == bait_map.end())  // Insert, if not present yet.
	    bait_map[thebait] = *it_Locus;
	}
      }
      ++it_Locus;
    }
  }

/*   void clean_NULL_elements_form_baitLoci_vec() */
/*   { */
/*     std::vector<CBaitLocus*>::iterator it_end, it1, it2; */
/*     unsigned count = 0; */
    
/*     it1 = it2 = baitLoci_vec.begin(); */
/*     it_end    = baitLoci_vec.end(); */

/*     // Advance it2 until we have the first NULL element: */
/*     while (*it2 != NULL && it2 != it_end) */
/*     { */
/*       ++it1; */
/*       ++it2; */
/*       ++count; */
/*     } */

/*     // if it2 is not at the end, it will point to a NULL element */
/*     while (it2 != it_end) */
/*     { */
/*       // Find next non-NULL element: */
/*       while (*it2 == NULL && it2 != it_end) */
/*       { */
/* 	++it2; */
/*       } */
/*       if (it2 == it_end) */
/* 	break; */
/*       //Copy: */
/*       *it1 = *it2; */
/*       ++count; */
      
/*       ++it1; */
/*       ++it2; */
/*     } */

/*     // Resize vector: */
/*     baitLoci_vec.resize(count, 0); */

/*   } */


/*   void sort_by_num_baits(std::vector<CBaitLocus*>::iterator it, std::vector<CBaitLocus*>::iterator it_end) */
/*   { */
/*     MYSORT (it, it_end, num_baits_lessThanCBaitLocus); */
/*   } */

/*   void sort_by_num_sequences(std::vector<CBaitLocus*>::iterator it, std::vector<CBaitLocus*>::iterator it_end) */
/*   { */
/*     MYSORT (it, it_end, num_sequences_in_region_greaterThanCBaitLocus); */
/*   } */

  void sort_by_gene_name__feature_num(std::vector<CBaitLocus*>::iterator it, std::vector<CBaitLocus*>::iterator it_end)
  {
    MYSORT (it, it_end, gene_name__feature_num_lessThanCBaitLocus);
  }

  unsigned get_number_of_loci() const
  {
    return baitLoci_vec.size();
  }


  void print_stats(std::ostream &os)
  {
    unsigned numBaits = 0;
    
    os << "Number of loci: " << baitLoci_vec.size() << std::endl;

    std::vector<CBaitLocus*>::iterator it, it_end;
    it     = baitLoci_vec.begin();
    it_end = baitLoci_vec.end();

    while (it != it_end)
    {
      (**it).print_stats(os);
      numBaits += (**it).get_num_baits();
      ++it;
    }
    os << "Total number of baits: " << numBaits << std::endl;
    //    os << "Bait length:           " << bait_length  << std::endl;
  }


  //////////////////
  // We now impelmented this in a constructor to be consistent with other filtering steps.
  //////////////////
/*   void filter_hit_coverage(double minimum_coverage, CBlast_parser &bp) */
/*   { */
/*     unsigned    tile_index, number_of_tiles; */
/*     faststring  recognition_string; */
/*     double      best_coverage_this_tiling_stack; */
/*     bool        bait_region_has_tiling_stack_with_insufficient_coverage; */

/*     std::vector<CBaitLocus*>::iterator it, it_end; */

/*     it     = baitLoci_vec.begin(); */
/*     it_end = baitLoci_vec.end(); */

/*     // for each bait region */
/*     while (it != it_end) */
/*     { */
/*       number_of_tiles = (**it).get_number_of_tiles(); */
/* /\*             std::cout << "Number of tiles for: " << (**it).get_bait_region_name_string() *\/ */
/* /\* 		<< " " << number_of_tiles << std::endl; *\/ */

/*       // Within each bait region determine the maximum hit coverage of all baits */
/*       // For each tiling stack of baits: */

/*       // Each tiling stack has to be checked individually. */
/*       // A single tiling stack with insufficient query coverage should result in the exclusion of this bait region. */
/*       bait_region_has_tiling_stack_with_insufficient_coverage = false; */
/*       for (tile_index=1; tile_index <= number_of_tiles; ++tile_index) */
/*       { */
/* 	// Do something with (**it): */
/* 	// Determine string to identify blast hits of baits of this bait region and tiling_index against the genome.  */
/* 	recognition_string = (**it).get_bait_region_name_string() + "|" + tile_index + "|";   */
	
/* 	best_coverage_this_tiling_stack = (double) bp.get_max_hit_query_span__matching_beginning1(recognition_string)/bait_length; */

/* 	if (best_coverage_this_tiling_stack < 0) */
/* 	{ */
/* 	  std::cout << "No hit for recognition string: " << recognition_string << " --> " */
/* 		    << 0.00 << std::endl; */
/* 	} */
/* 	else */
/* 	{ */
/* 	  std::cout << "Determined best coverage of this tiling stack: " << recognition_string << " --> " */
/* 		    << best_coverage_this_tiling_stack << std::endl; */
/* 	} */

/* 	if (best_coverage_this_tiling_stack < minimum_coverage) */
/* 	  bait_region_has_tiling_stack_with_insufficient_coverage = true; */
/*       } */

/*       if (bait_region_has_tiling_stack_with_insufficient_coverage) */
/*       { */
/* 	if (delete_in_destructor) */
/* 	  delete *it; */
/* 	*it = NULL; */
/*       } */
/*     } */

/*     clean_NULL_elements_form_baitLoci_vec(); */
/*   } */


  void print_all(std::ostream &os, char format)
  {
    std::vector<CBaitLocus*>::iterator it, it_end;
    it     = baitLoci_vec.begin();
    it_end = baitLoci_vec.end();

    // TODO: In case we have no alignment cutting we should use an alternative heading in the output file.
    //       This is not urgent, just more elegant. It requires some work to get the information to this point.
    if (format == 's')
    {
      //  if (alignment_cutting)
      os << "# Original-file-name\tGene-name\tFeature-number\tNumber-of-features-of-specified-type-for-this-gene-in-gff-file\tStart-coordinate-of-bait-region-in-alignment\tNumber-of-sequences-in-region-the-baits-are-based-on\tNumber-of-baits-constructed-for-this-bait-region\tList of baits - for each bait the tiling design column is given, then the bait string, its CG content and its maximum distance to the sequences in the alignment." << std::endl;
      //  else
      //    os << "# Corresponding sequence name\tGene-name\tn.a.\tn.a.\tStart-coordinate-of-bait-region-in-alignment\t Number-of-sequences-in-region-the-baits-are-based-on\tNumber-of-baits-constructed-for-this-bait-region\tList of baits - for each bait the tiling design column is given, then the bait string, its CG content and its maximum distance to the sequences in the alignment." << std::endl;
    }

    if ( !(format == 's' || format == 'f' ) )
    {
      std::cerr << "ERROR: Unknown internal format " << format << " found in call to CBaitLoci_collection::print_all(2 parameters)" << std::endl;
      exit(-24);
    }

    while (it != it_end)
    {
      (**it).print(os, format);
      ++it;
    }
  }


  void print_all(std::ostream &os, char format, const char *probeIDprefix)
  {
    std::vector<CBaitLocus*>::iterator it, it_end;
    it     = baitLoci_vec.begin();
    it_end = baitLoci_vec.end();

    if (format == 'A')
    {
      os << "TargetID\tProbeID\tSequence\tReplication" << std::endl;
    }
    else
    {
      std::cerr << "ERROR: Unknown internal format " << format << " found in call to CBaitLoci_collection::print_all(3 parameters)" << std::endl;
      exit(-25);
    }

    while (it != it_end)
    {
      (**it).print(os, format, probeIDprefix);
      ++it;
    }
  }


  void fill_locus_selector(CBaitRegionSelectorSets &ls) const
  {
    std::vector<CBaitLocus*>::const_iterator it, it_end;
    it     = baitLoci_vec.begin();
    it_end = baitLoci_vec.end();

    while (it != it_end)
    {
      ls.add( (**it).get_gene_name(),
	      (**it).get_feature_num(),
	      (**it).get_start() );
      ++it;
    }
  }


  void print_with_start_stepwidth(std::ostream &os, char format, unsigned start, unsigned step)
  {
    // TODO: In case we have no alignment cutting we should use an alternative heading in the output file.
    //       This is not urgent, just more elegant. It requires some work to get the information to this point.
    if (format == 's')
    {
      //  if (alignment_cutting)
      os << "# Original-file-name\tGene-name\tFeature-number\tNumber-of-features-of-specified-type-for-this-gene-in-gff-file\tStart-coordinate-of-bait-region-in-alignment\tNumber-of-sequences-in-region-the-baits-are-based-on\tNumber-of-baits-constructed-for-this-bait-region\tList of baits - for each bait the tiling design column is given, then the bait string, its CG content and its maximum distance to the sequences in the alignment." << std::endl;
      //  else
      //    os << "# Corresponding sequence name\tGene-name\tn.a.\tn.a.\tStart-coordinate-of-bait-region-in-alignment\t Number-of-sequences-in-region-the-baits-are-based-on\tNumber-of-baits-constructed-for-this-bait-region\tList of baits - for each bait the tiling design column is given, then the bait string, its CG content and its maximum distance to the sequences in the alignment." << std::endl;
    }

    unsigned i,N = baitLoci_vec.size();

    i = start-1;  // start is 1 based, i is 0 based. Therefore we have to subtract 1.
    while (i < N)
    {
      baitLoci_vec[i]->print(os, format);
      i += step;
    }
  }

  void print_with_start_stepwidth_multi_gene_feature_aware(std::ostream &os, char format, unsigned start_count, unsigned offset_count,
							   std::vector<CBaitLocus*>::iterator it_start_this_alignment)
  {
    if (format == 's')
    {
      //  if (alignment_cutting)
      os << "# Original-file-name\tGene-name\tFeature-number\tNumber-of-features-of-specified-type-for-this-gene-in-gff-file\tStart-coordinate-of-bait-region-in-alignment\tNumber-of-sequences-in-region-the-baits-are-based-on\tNumber-of-baits-constructed-for-this-bait-region\tList of baits - for each bait the tiling design column is given, then the bait string, its CG content and its maximum distance to the sequences in the alignment." << std::endl;
      //  else
      //    os << "# Corresponding sequence name\tGene-name\tn.a.\tn.a.\tStart-coordinate-of-bait-region-in-alignment\t Number-of-sequences-in-region-the-baits-are-based-on\tNumber-of-baits-constructed-for-this-bait-region\tList of baits - for each bait the tiling design column is given, then the bait string, its CG content and its maximum distance to the sequences in the alignment." << std::endl;
    }

    std::vector<CBaitLocus*>::iterator it_tmp, it;

    it = it_start_this_alignment;

    if ( (**it).has_feature_information() )
      it_tmp = find_end_this_feature(it);
    else
      it_tmp = find_end_this_gene(it);

    // First we determine the starting point for the iterator it "at start_count"
    unsigned i=0, startdiff, laststart, start;

    laststart = (**it).get_start();
    while (i < start_count) // We won't walk past the end of this feature  or alignment
    {
      ++it;
      if (it == it_tmp)
	break;
      // If we had to skip possible starting points, we try to correct for this.
      start = (**it).get_start();
      startdiff = start - laststart; // Difference of indices, should be at least 1.
      i += startdiff; 
      laststart = start;
    }

    if (it != it_tmp) // we should never have it == it_tmp, but if we have, we will crash, so we better check.
    {
      // Loop over all loci in steps of offset_count
      while (true) // We won't walk past the end of this feature or alignment. We break inside the loop if it == it_tmp.
      {
	// This locus will be printed.
	(*it)->print(os, format);

	i = 0; // Count the number of loci we skip in next iteration.	
	// Move to next locus we want to print:
	// We move by one locus position. (The first move.)  
	// We can handle the first move outside the loop below since offset_count must be > 0.
	++it;
	if (it == it_tmp)
	  break;

	// If we had to skip possible starting points, we try to correct for this.
	start = (**it).get_start();
	startdiff = start - laststart;
	i += startdiff; // Diff is 1 or larger.
	laststart = start;
	  
	// Skip offset_count loci.
	while (it < it_tmp && i < offset_count)
	{
	  // Next iteration.
	  ++it;
	  if (it == it_tmp)
	    break;
	  start = (**it).get_start();
	  startdiff = start - laststart;
	  i += startdiff;
	  laststart = start;
	}
	  
	// If it == it_tmp, we have an iterator not pointing into this alignment/feature any more.
	// So we have to break immediately.
	if (it == it_tmp)
	  break;
	  
	// In the ideal case i == offset_count. Then we are at the next locus we wanted to count.
	// If i < offset_count we have not reached the next good locus in this alignment/gene/feature
	//   before we moved past the end. In this case we have it == it_tmp and we leave the
	//   loop with the next while check.
	// If i > offset_count, then we had to skip more loci than we wanted to skip, since not all
	// loci are valid starting points. In this case we add a penalty.
       
      } // END while true
    } // END if (it != it_tmp)
  } // END print_with_start_stepwidth_multi_gene_feature_aware



  std::vector<CBaitLocus*>::iterator find_end_this_feature(std::vector<CBaitLocus*>::iterator it)
  {
    unsigned   the_feature_num = (**it).get_feature_num();
    faststring the_gene_name    = (**it).get_gene_name();

    while (it != baitLoci_vec.end() )
    {
      if (the_feature_num != (**it).get_feature_num() || the_gene_name != (**it).get_gene_name())
      {
	break;
      }
      ++it;
    }
    return it;
  }

  std::vector<CBaitLocus*>::iterator find_end_this_gene(std::vector<CBaitLocus*>::iterator it)
  {
    faststring the_gene_name    = (**it).get_gene_name();

    while (it != baitLoci_vec.end() )
    {
      if (the_gene_name != (**it).get_gene_name())
      {
	break;
      }
      ++it;
    }
    return it;
  }

  void sort_within_featues(bool (*f)(CBaitLocus*,CBaitLocus*))
  {
    std::vector<CBaitLocus*>::iterator it, it_end, it_tmp;
    it     = baitLoci_vec.begin();
    it_end = baitLoci_vec.end();

    while (it != it_end)
    {
      // Find end of this gene:
      it_tmp = find_end_this_feature(it);
      //      std::cout << "Sort range: " << distance(baitLoci_vec.begin(), it) << " " << distance(baitLoci_vec.begin(), it_tmp) << std::endl;
      MYSORT(it,it_tmp, f);
      it = it_tmp;
    }
  }

  void sort_within_genes(bool (*f)(CBaitLocus*,CBaitLocus*))
  {
    std::vector<CBaitLocus*>::iterator it, it_end, it_tmp;
    it     = baitLoci_vec.begin();
    it_end = baitLoci_vec.end();

    while (it != it_end)
    {
      // Find end of this gene:
      it_tmp = find_end_this_gene(it);
      MYSORT(it,it_tmp, f);
      it = it_tmp;
    }
  }

/*   void sort_by_original_file_name__feature_num_start() */
/*   { */
/*     MYSORT(baitLoci_vec.begin(), baitLoci_vec.end(), original_file_name__feature_num_start_lessThanCBaitLocus); */
/*   } */


  bool is_ordered_within_genes()
  {
    std::vector<CBaitLocus*>::iterator it, it_end, it_tmp;
    it     = baitLoci_vec.begin();
    it_end = baitLoci_vec.end();

    if (it == it_end)
      return true;

    it_tmp = it;
    ++it_tmp;

    while (it_tmp != it_end)
    {
      if (!gene_name__feature_num_lessThanCBaitLocus(*it, *it_tmp) )
	return false;
      it = it_tmp;
      ++it_tmp;
    }
    return true;
  }


  // Determines whether feature numbers reappear in successive loci
  // if a different feature number was found in between.
  // This ordering is assumed when filtering for the best loci in a feature.
  bool is_partially_ordered()
  {
    std::set<faststring> set_genes_features;

    std::vector<CBaitLocus*>::iterator it, it_end;
    it     = baitLoci_vec.begin();
    it_end = baitLoci_vec.end();

    unsigned   the_feature_num;
    faststring the_gene_name;

    if (it == it_end)
      return true;

    the_feature_num  = (**it).get_feature_num();
    the_gene_name    = (**it).get_gene_name();

    faststring gene_feature_name;
    gene_feature_name = the_gene_name + "_" + the_feature_num;
    set_genes_features.insert(gene_feature_name);

    while (it != it_end)
    {
      if (the_feature_num != (**it).get_feature_num() || the_gene_name != (**it).get_gene_name() )
      {
	// new gene or exon:
	the_feature_num   = (**it).get_feature_num();
	the_gene_name     = (**it).get_gene_name();
	gene_feature_name = the_gene_name + "_" + the_feature_num;

	if ( set_genes_features.find(gene_feature_name) != set_genes_features.end() )
	{
	  return false;  // We have seen this gene/exon combination before. So we return false;
	}
	set_genes_features.insert(gene_feature_name);
      }
      ++it;
    }
    return true;
  }

  void get_seqs_baits_count(unsigned &NumSeqs, unsigned &NumBaits) const
  {
    std::vector<CBaitLocus*>::const_iterator it, it_end;
    it     = baitLoci_vec.begin();
    it_end = baitLoci_vec.end();
    
    unsigned local_NumSeqs  = 0;
    unsigned local_NumBaits = 0;

    while (it != it_end)
    {
      local_NumSeqs  += (**it).get_num_sequences();
      local_NumBaits += (**it).get_num_baits();
      ++it;
    }

    NumSeqs  = local_NumSeqs;
    NumBaits = local_NumBaits; 
  }

  // Used to evaluate a starting index in thinning mode.
  // Thinning can be done for multiple alignment files even though we originally had one file in mind.
  // At the moment we do not start at position 1 when coming to the next alignment.

  void count_seqs_baits_for_start_and_offset(unsigned start_count, unsigned offset_count,
					     unsigned &sum_seq,    unsigned &sum_baits)
  {
    unsigned i,N;
    N = baitLoci_vec.size();

    sum_seq = sum_baits = 0;

    for (i=start_count; i<N; i += offset_count)
    {
      sum_baits += baitLoci_vec[i]->get_num_baits();
      sum_seq   += baitLoci_vec[i]->get_num_sequences();
    }
  }

  // Used to evaluate a starting index in thinning mode.
  // Thinning can now be done for multiple alignment files even though we it was originally designed to work only for single files.
  // This function is called for the loci of each alignment individually.
  // This functions obtains the first locus it should consider and walks through all loci until the loci belong to the next alignment.
  // It passes back the start locus for the next alignment.
  //// At the moment we do not start at position 1 when coming to the next alignment.

  void count_seqs_baits_for_start_and_offset_multi_gene_feature_aware(unsigned alignment_counter,
								      std::vector<CBaitLocus*>::iterator &it_start_this_alignment,
								      std::vector<CBaitLocus*>::iterator &it_start_next_alignment,
								      bool     &finished,
								      unsigned start_count, unsigned offset_count,
								      unsigned &sum_seq,    unsigned &sum_baits,
								      unsigned &penalty_diff)
  {
    unsigned i;
    unsigned lastLocusStart, thisLocusStart, locusStartDiff;
    int penalty_diff_local = 0;

    sum_seq = sum_baits = 0;

    if (global_verbosity >= 50)
    {
      std::cerr << "Working on alignment number: " << alignment_counter << " start " << start_count << std::endl;
    }

    std::vector<CBaitLocus*>::iterator it, it_end_al;

    if (alignment_counter == 1)  // Since the calling function does not have access to the baitLoci_vec, the start is initialized in this function, if alignment_counter == 1.
    {
      it = it_start_this_alignment = baitLoci_vec.begin();
    }
    else
    {
      it = it_start_this_alignment;
    }

    if ( (**it).has_feature_information() )
      it_end_al = find_end_this_feature(it);
    else
      it_end_al = find_end_this_gene(it);

    
    // First we determine the start for iterator at start_count
    i=0;
    lastLocusStart = (**it).get_start();  // This information is used to keep track of already deleted starting positions of loci.
                                          // This variable is initialised with the locus staring coordinate of the first bait region starting position.

    while (i < start_count) // We won't walk past the end of this feature  or alignment
    {
      ++it;
      if (it == it_end_al)
	break;
      // If we had to skip possible starting points, we try to correct for this.
      thisLocusStart = (**it).get_start();
      locusStartDiff = thisLocusStart - lastLocusStart; // Difference of indices, should be at least 1. Is 0 for first bait region.
      i += locusStartDiff; 
      lastLocusStart = thisLocusStart;
    }

    //      std::cout << "Found start (index:start_count, index-found:i): " << start_count << " " << i << std::endl;

    // Differences are not penalised if they occur while searching for the starting position.
    // In the ideal case, i == start_count. i < start_count should not be possible.
    // However, if we moved passed the ideal start_count, we penalise this.
    if (i > start_count)
    {
      penalty_diff_local = i-start_count;
    }

    // If we moved past the end of the alignment, we found no valid starting point in this alignment.
    // This will be penalised. This penalty might superseed a previous penalty.
    if (it == it_end_al)
    {
      // The amount of the penalty could be a matter of debate. A small/moderate penalty is the
      // start_count itself. It penalises larger start counts more.
      penalty_diff_local = start_count;
    }
    else
    {
      // Loop over all loci in steps of offset_count
      while (true) // We won't walk past the end of this feature or alignment. We break inside the loop if it == it_end_al.
      {
	// This locus is scored.
	//	  std::cout << "Score locus (seq-coord:start,counter:i): " <<  (**it).get_start() << " " << i << std::endl;
	sum_baits += (**it).get_num_baits();
	sum_seq   += (**it).get_num_sequences();

	i = 0; // Count the number of loci we skip in next iteration.
	
	// Move to next loci we want to add:
	// We move by one locus position. (The first move.)  
	// We can handle the first move outside the loop below since offset_count must be > 0.
	++it;
	if (it == it_end_al)
	  break;

	// If we had to skip possible starting points, we try to correct for this.
	thisLocusStart = (**it).get_start();
	locusStartDiff = thisLocusStart - lastLocusStart;
	i += locusStartDiff; // Diff is 1 or larger.
	lastLocusStart = thisLocusStart;
	  
	// Skip offset_count loci.
	while (it < it_end_al && i < offset_count)
	{
	  // Next iteration.
	  ++it;
	  if (it == it_end_al)
	    break;
	  thisLocusStart = (**it).get_start();
	  locusStartDiff = thisLocusStart - lastLocusStart;
	  i += locusStartDiff;
	  lastLocusStart = thisLocusStart;
	}
	  
	// If it == it_end_al, we have an iterator not pointing this alignment/feature any more.
	// So we have break immediately.
	if (it == it_end_al)
	  break;
	  
	// In the ideal case i == offset_count. Then we are at the next locus we wanted to count.
	// If i < offset_count we have not reached the next good locus in this alignment/gene/feature
	//   before we moved past the end. In this case we have it == it_end_al and we leave the
	//   loop with the next while check.
	// If i > offset_count, then we had to skip more loci than we wanted to skip, since not all
	// loci are valid starting points. In this case we add a penalty.
	  
	if (i > offset_count)
	{
	  penalty_diff_local += (i-offset_count);
	}
      } // END while true
    } // END ELSE if (it == it_end_al)

    penalty_diff = penalty_diff_local;
    if (it == baitLoci_vec.end() )
      finished = true;
    else
      finished = false;

    it_start_next_alignment = it_end_al;

  } //END FUNCTION count_seqs_baits_for_start_and_offset_multi_gene_feature_aware


  void get_original_file_name_gene_name(unsigned alignment_counter,
			      std::vector<CBaitLocus*>::iterator &it_start_this_alignment,
			      faststring &P_original_file_name,
			      faststring &P_gene_name
			      )
  {
    if (alignment_counter == 1)
    {
      it_start_this_alignment = baitLoci_vec.begin();
    }

    P_original_file_name  = (**it_start_this_alignment).get_original_file_name();
    P_gene_name           = (**it_start_this_alignment).get_gene_name();
  }




  double get_max_dist_mean() const
  {
    double   sum=0;
    double   num=0;
    int     i, N = baitLoci_vec.size();
    
    for (i=0; i<N; ++i)
    {
      sum += baitLoci_vec[i]->get_sum_max_distances();
      num += baitLoci_vec[i]->get_num_baits();
    }

    //    std::cerr << "sum max distances: " << sum << std::endl;
    //    std::cerr << "num baits:         " << num << std::endl;

    return sum/num;
  }
  
  double get_overall_max_dist() const
  {
    double overall_max_dist = -1000;  // True values must be in the range 0..1.
    double om_dist;
    int    i, N             = baitLoci_vec.size();
    
    for (i=0; i<N; ++i)
    {
      om_dist = baitLoci_vec[i]->get_overall_max_dist();
      if (om_dist > overall_max_dist)
      {
	overall_max_dist = om_dist;
      }
    }
    return overall_max_dist;
  }


  double get_CG_mean() const
  {
    double   sum=0;
    double   num=0;
    int     i, N = baitLoci_vec.size();
    
    for (i=0; i<N; ++i)
    {
      sum += baitLoci_vec[i]->get_sum_CG();
      num += baitLoci_vec[i]->get_num_baits();
    }
    return sum/num;
  }

}; // END CBaitLoci_collection
