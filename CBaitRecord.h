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
#include "faststring2.h"
#include <vector>
#include <set>
#include "Ctriple.h"
#include "print_container.h"
#include "CBlastParser.h"



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
    os << "Selector string (alignment) arument, unsigned (feature number): " << std::endl;
    print_container(os, set_of_alignment_names_feature_num, "", "\n" , "\n");
    os << "Selector string (alignment) arument, unsigned (feature number), unsigned (start): " << std::endl;
    print_container(os, set_of_alignment_name_feature_num_start__ie_locus, "", "\n" , "\n");
  }

  void print_num_genes_num_features(std::ostream &os, const char *s) const
  {
    os << s << "Bait regions have been determined for this number of different:" << std::endl;
    os << s << "Alignment files                                             " << set_of_alignment_names.size() << std::endl;
    os << s << "Feature (is identical to number of files if not applicable: " << set_of_alignment_names_feature_num.size() << std::endl;
    os << s << "Bait regions:                                               " << set_of_alignment_name_feature_num_start__ie_locus.size() << std::endl;
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

    str.divide_at(' ', pre, post2);
    if (!pre.isAnUnsigned())
    {
      std::cerr << "Baits are preceeded by a number indicating the index of the tiling design column. This number could not be found for the following bait:"
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
      std::cerr << "Baits are proceeded by two double numbers containing the CG content and the maximum distance in its cluster.\n"
		<< "Expecting a double number for the CG content but found: " << pre 
		<< std::endl;
      std::cerr << "Exiting." << std::endl;
      exit(-4);
    }
    if (!post2.isADouble())
    {
      std::cerr << "Baits are proceeded by two double numbers containing the CG content and the maximum distance in its cluster.\n"
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

bool seq_name__feature_num_start_lessThanCBaitLocus(CBaitLocus *a, CBaitLocus *b);

bool locus_start_lessThanCBaitLocus(CBaitLocus *a, CBaitLocus *b);



class CBaitLocus
{
  faststring seq_name;
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

    seq_name                = container[0];
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

  faststring &get_seq_name()
  {
    return seq_name;
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
    return seq_name + "|" + gene_name + "|" + faststring(feature_num) + "|" + faststring(start);
  }

  // TODO: Change this to FILE * instead of using ostream.
  // Known formats are 'f', 's' for the 2 paramter version of this function.
  void print(std::ostream &os, char format)
  {
    std::vector<CBait *>::iterator i1, i2, i3;
    //    faststring name;
    unsigned count = 1;

    i1 = bait_vec.begin();
    i2 = bait_vec.end();

    if (format == 'f') // fasta
    {
      //      name = ">" + seq_name + "|" + gene_name + "|" + faststring(feature_num) + "|" + faststring(start);

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
      os << seq_name << "\t"
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
  // Known formats are 'A' for the 3 paramter version of this function.
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
    // For an equal number of baits, prefer a locus based on more sequences.
    if (a->num_baits == b->num_baits)
      return a->num_sequences_in_region > b->num_sequences_in_region;

    return a->num_baits < b->num_baits;
  }

  friend bool num_sequences_in_region_lessThanCBaitLocus(CBaitLocus *a, CBaitLocus *b)
  {
    // For loci based on an equal number of sequences, prefer the one with less needed baits.
    // Note: This comparator sorts in reversed order.
    if (a->num_sequences_in_region == b->num_sequences_in_region)
      return a->num_baits > b->num_baits;
      
    return a->num_sequences_in_region < b->num_sequences_in_region;
  }

  // First criterion:  gene_name
  // Second criterion: feature num
  // Third  criterion: locus start coordinate
  friend bool gene_name__feature_num_lessThanCBaitLocus(CBaitLocus *a, CBaitLocus *b)
  {
    if (a->gene_name != b->gene_name)
      return a->gene_name < b->gene_name;

    if (a->feature_num != b->feature_num)
      return a->feature_num < b->feature_num;
    
    return a->start < b->start;
  }

  friend bool num_baits_greaterThanCBaitLocus(CBaitLocus *a, CBaitLocus *b)
  {
    if (a->num_baits == b->num_baits)
      return  a->num_sequences_in_region < b->num_sequences_in_region;

    return a->num_baits > b->num_baits;
  }

  friend bool num_sequences_in_region_greaterThanCBaitLocus(CBaitLocus *a, CBaitLocus *b)
  {
    if (a->num_sequences_in_region == b->num_sequences_in_region)
      return a->num_baits < b->num_baits;

    return a->num_sequences_in_region > b->num_sequences_in_region;
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
    if (a->gene_name == b->gene_name)
      return a->start < b->start;
    return a->gene_name < b->gene_name;
  }

  friend bool seq_name__feature_num_start_lessThanCBaitLocus(CBaitLocus *a, CBaitLocus *b)
  {
    if (a->seq_name != b->seq_name)
      return a->seq_name < b->seq_name;

    if (a->feature_num != b->feature_num)
      return a->feature_num < b->feature_num;
    
    return a->start < b->start;
  }

}; // END CBaitLocus




class CBaitLoci_collection
{
  std::vector<CBaitLocus*>                                       baitLoci_vec;
  std::map<Ctriple<faststring, unsigned, unsigned>, CBaitLocus*> baitLoci_map;
  bool                                                           delete_in_destructor;
  //  unsigned                                                       bait_length; // Will be set to 0 in constructor and will be adjusted to its true value in the add method.
  //                                                                              // The fact that we have one value for a collection of baits implies, that we assume that all baits in all loci have the same length.

  // Internal add 
  void add_locus(CBaitLocus *p)
  {
    baitLoci_vec.push_back(p);
    baitLoci_map.insert(std::pair<Ctriple<faststring, unsigned, unsigned>, CBaitLocus*>(Ctriple<faststring, unsigned, unsigned>(p->get_seq_name(), p->get_feature_num(), p->get_start()), p) );

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
 CBaitLoci_collection(std::ifstream &is):delete_in_destructor(true) //, bait_length(0)
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


 CBaitLoci_collection(CBaitLoci_collection &source, char genes_or_features, char num_baits_or_sequences):delete_in_destructor(false) //, bait_length(0)
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
	  std::cout << "WARNING: Sorting NOT successful." << std::endl;
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
    sort_within_genes(gene_name__feature_num_lessThanCBaitLocus);
  } // END CBaitLoci_collection(CBaitLoci_collection &source, char genes_or_features, char num_baits_or_sequences)


 CBaitLoci_collection(const CBaitLoci_collection &source, const CBaitRegionSelectorSets &lss, char mode):delete_in_destructor(false) // , bait_length(0)
 {
   unsigned i, N;
   N=source.baitLoci_vec.size();

   if (mode == 'b')      // Remove all loci from ALIGNMENT if a double hit occurs in this set
   {
     for (i=0; i<N; ++i)
     {
       CBaitLocus *p = source.baitLoci_vec[i];

       if (!lss.exists(p->get_seq_name()) )
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

       if (!lss.exists(p->get_seq_name(), p->get_feature_num()) )
       {
	 add_locus(p);
	 //	 baitLoci_vec.push_back(p);
       }
     }
   }
   if (mode == 'x')      // Remove Locus is if a double hit occurs in it
   {
     for (i=0; i<N; ++i)
     {
       CBaitLocus *p = source.baitLoci_vec[i];

       if (!lss.exists(p->get_seq_name(), p->get_feature_num(), p->get_start()) )
       {
	 add_locus(p);
	 //	 baitLoci_vec.push_back(p);
       }
     }
   }
 }


 CBaitLoci_collection(const CBaitLoci_collection &source, CBlast_parser &bp, double minimum_coverage):delete_in_destructor(false) //, bait_length(source.bait_length)
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

    N = source.baitLoci_vec.size();
    count = 0;

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
	std::cerr << "The bait length has been determined to be 0 for the bait region: " << std::endl;  
	(**it).print(std::cerr, 's');
	std::cerr << "This bug should be reported. It should not occur. It might be caused by mal formated input files,"
	             " but the error should be caught earlier than here." << std::endl;
	exit(-5);
      }

      if (global_verbosity > 2)
      {
	std::cout << "Region: " << count << " of " << N << std::endl;
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

	if (global_verbosity > 50)
	{
	  std::cout << "Tile: " << tile_index << " " << recognition_string << " " << max_hit_length << "/" << bait_length_this_bait_region << "=" << best_coverage_this_tiling_stack << std::endl;
	}

	if (best_coverage_this_tiling_stack < 0)
	{
	  if (global_verbosity > 50)
	  {
	    std::cout << "No hit for recognition string: " << recognition_string << " --> "
		      << 0.00 << std::endl;
	  }
	}
	else
	{
	  if (global_verbosity > 50)
	  {
	    std::cout << "Determined best coverage of this tiling stack: " << recognition_string << " --> "
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


  void get_map_of_baits_to_locus_info(faststring seq_name_match, unsigned tiling_index, std::map<faststring, CBaitLocus *> &bait_map)
  {
    //    sort_by_seq_name__feature_num_start();

    std::vector<CBaitLocus*>::iterator it_Locus, it_Locus_end;
    it_Locus     = baitLoci_vec.begin();
    it_Locus_end = baitLoci_vec.end();

    std::vector<CBait *> tmp_bait_vec;

    unsigned i, N;

    std::map<faststring, CBaitLocus *>::iterator find_it;

    while (it_Locus != it_Locus_end)
    {
      if ( (**it_Locus).get_seq_name().find(seq_name_match) != faststring::npos)
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

/*     // Advance it2 until we have the first NUll element: */
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


  void sort_by_num_baits(std::vector<CBaitLocus*>::iterator it, std::vector<CBaitLocus*>::iterator it_end)
  {
    sort (it, it_end, num_baits_lessThanCBaitLocus);
  }

  void sort_by_num_sequences(std::vector<CBaitLocus*>::iterator it, std::vector<CBaitLocus*>::iterator it_end)
  {
    sort (it, it_end, num_sequences_in_region_greaterThanCBaitLocus);
  }

  void sort_by_gene_name__feature_num(std::vector<CBaitLocus*>::iterator it, std::vector<CBaitLocus*>::iterator it_end)
  {
    sort (it, it_end, gene_name__feature_num_lessThanCBaitLocus);
  }


  unsigned get_number_of_loci()
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
  // We nowimpelmented this in a constructor to be consistent with other filtering steps.
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

    if (format == 's')
    {
      os <<  "# Corresponding sequence name\tGene-name\tfeature-number\tnumber-of-features-of-specified-type-for-this-gene-in-gff-file\tstart-coordinate-of-bait-region-in-alignment\tnumber-of-sequences-in-region-the-baits-are-based-on\tnumber-of-baits-constructed-for-this-bait-region\tList of baits - each bait is proceed by the tiling stack number and followed by CG content and the maximum distance to the sequences in the alignment." << std::endl;
    }

    if ( !(format == 's' || format == 'f' ) )
    {
      std::cerr << "ERROR: Unkown internal format " << format << " found in call to CBaitLoci_collection::print_all(2 paramters)" << std::endl;
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
      std::cerr << "ERROR: Unkown internal format " << format << " found in call to CBaitLoci_collection::print_all(3 parameters)" << std::endl;
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
    if (format == 's')
    {
      os <<  "# Corresponding sequence name\tGene-name\tfeature-number\tnumber-of-features-of-specified-type-for-this-gene-in-gff-file\tstart-coordinate-of-bait-region-in-alignment\tnumber-of-sequences-in-region-the-baits-are-based-on\tnumber-of-baits-constructed-for-this-bait-region\tList of baits - each bait is proceed by the tiling stack number and followed by CG content and the maximum distance to the sequences in the alignment." << std::endl;
    }

    unsigned i,N = baitLoci_vec.size();

    i = start-1;  // start is 1 based, i is 0 based. Therefore we have to subtract 1.
    while (i < N)
    {
      baitLoci_vec[i]->print(os, format);
      i += step;
    }
  }

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
      sort(it,it_tmp, f);
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
      sort(it,it_tmp, f);
      it = it_tmp;
    }
  }

  void sort_by_seq_name__feature_num_start()
  {
    sort(baitLoci_vec.begin(), baitLoci_vec.end(), seq_name__feature_num_start_lessThanCBaitLocus);
  }


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
