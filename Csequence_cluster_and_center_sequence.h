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
#ifndef CSEQUENCE_CLUSTER_AND_CENTER_SEQUENCE
#define CSEQUENCE_CLUSTER_AND_CENTER_SEQUENCE

#include "DEBUG_STUFF.h"
#include <set>
#include <list>
//#include "CTaxonNamesDictionary.h" // Has been used in earlier version, in which clusters where reported with unique sequence number. This intermediate step was not convenient for the user and is not used any more.
#include <iostream>
#include <fstream>
#include "faststring2.h"
#include "basic-DNA-RNA-AA-routines.h"
#include <climits>
#include "statistic_functions.h"
#include "print_container.h"
//#include "csv-toolbox.h"
//#include "mydir.h"
//#include "sequence-cluster-result-parser.h"
#include "CSequences2.h"
#include "CSequence_Mol2_1.h"
#include "CDistance_matrix.h"
#include "CSeqNameList.h"
#include <vector>
#include <ctime>


// Prototypes of global functions implemented in csequence_cluster_and_center_sequence.cpp
void recode_range(char * b, char *e);
void back_recode_range(char * b, char *e);
void recode_MSA(char **sequences, int num_taxa, int num_pos);
void find_center_sequence(char **msa, char *consensus,
			    char *variant, char *isAmbig,
			    int len_variant,
			    int pos_start, int num_taxa,
			    int pos_start_msa, int pos_end_msa,
			    int *distances, faststring &best_variant,
			    int *p_best_maximum_distance);
void find_center_sequence_for_msa_exhaustive(char **msa, faststring &center_seq,
					       int num_taxa,
					       int pos_start_msa, int pos_end_msa,
					       double &mdist);

void find_center_sequence_for_msa_heuristic(char **msa, faststring &center_seq,
					      int num_taxa,
					      int pos_start_msa, int pos_end_msa,
					      double &mdist);


class Csequence_loci_cluster_collection;

class Csequence_cluster_and_center_sequence
{
  unsigned    pos_start; // 0 based
  unsigned    pos_end;   // 0 based, after last element of range

  double      max_dist_in_cluster;  // Maximum distance within this cluster
  double      max_dist_bait_to_msa;
                                    // This is not the threshold but the true maximum distance.

  std::set<short>     taxon_set;  // This is the set of sequence numbers in the cluster.
                                  // Sequence numbers are the numbers that can be used to access the sequences in the CSequences2 object
                                  // The sequences can be accessed from this class by following the *loci variable.
  faststring center_sequence;
  char       **msa;               // MSA of this cluster. In the assign function, the filter the sequences that belong to the cluster.

  Csequence_loci_cluster_collection *loci;  // Via the loci we have access to the CSequnces2 object and the CTaxonNamesDictionary object.
  double     CG_content;


 public:
 Csequence_cluster_and_center_sequence():pos_start(-1u-1), pos_end(-1u), max_dist_in_cluster(-1), msa(NULL), loci(NULL), CG_content(0)
  {}

  // Constructor: Never used, never tested.
/*  Csequence_cluster_and_center_sequence(unsigned          P_pos_start, */
/* 					 unsigned           P_pos_end, */
/* 					 double             P_max_dist_in_cluster, */
/* 					 std::set<short> P_taxon_set, */
/* 					 Csequence_loci_cluster_collection *P_loci): */
/*   pos_start(P_pos_start), pos_end(P_pos_end), max_dist_in_cluster(P_max_dist_in_cluster), max_dist_bait_to_msa(-1), */
/*     taxon_set(P_taxon_set), msa(NULL), loci(P_loci), CG_content(0) */
/*   {} */

  ~Csequence_cluster_and_center_sequence()
  {
    if (msa)
      delete [] msa;
    else
    {
      // We can get here if we destroy a default constructed object
      // that has not been assigned any data. That can happen.
    }
  }

  // The only function that is allowed to assign data to objects of this class.
  // The msa that is passed to this member function is not the msa that is being stored
  // in the object, since not all sequences belong to the cluster we deal with.
  void assign(unsigned           P_pos_start,
	      unsigned           P_pos_end,
	      double             P_max_dist_in_cluster,
	      std::set<short>    &P_taxon_set, // Contains the indices of the sequences in P_msa that are copied to the msa array.
	      Csequence_loci_cluster_collection *P_loci,
	      char **P_msa)
  {
    pos_start = P_pos_start;
    pos_end   = P_pos_end;
    max_dist_in_cluster  = P_max_dist_in_cluster;
    // We could swap the content of the two sets. This might be faster.
    taxon_set = P_taxon_set;
    loci      = P_loci;

    std::set<short>::iterator it, it_end;

    it     = taxon_set.begin();
    it_end = taxon_set.end();
    short i;

    if (msa)
    {
      std::cerr << "WARNING: An existing msa is overwritten. This behaviour is not expected. " << std::endl;
      delete msa;
    }

    msa = new char * [taxon_set.size()];

    i=0;
    while (it != it_end)
    {
      msa[i] = P_msa[*it];
      ++i;
      ++it;
    }
  }


  unsigned get_number_of_elements_in_cluster()
  {
    return taxon_set.size();
  }

  bool center_computed()
  {
    return !center_sequence.empty();
  }

  const faststring& get_center_sequence()
  {
    return center_sequence;
  }

  double get_CG()
  {
    return CG_content;
  }

  double get_max_dist_in_cluster()
  {
    return max_dist_in_cluster;
  }

  double get_max_dist_bait_to_msa()
  {
    return max_dist_bait_to_msa;
  }

  unsigned window_size()
  {
    return pos_end - pos_start;
  }

  // Implemented in Csequence_cluster_and_center_sequence.cpp
  void compute_center(char);

  void print(std::ostream &os)
  {
    os << "Cluster: " << pos_start+1 << " " << pos_end << " "
       << max_dist_in_cluster << " ";
    print_container(os, taxon_set, "{", ",", "}\n");

    if (msa != NULL)
    {
      print_msa_backrecoded(os);
    }

    if (center_sequence.empty())
    {
      os << " ## Center sequenuence has not been computed yet." << std::endl;
    }
    else
    {
      os << center_sequence << std::endl;
    }
  }

  void print_msa_backrecoded(std::ostream &os, int flag = 0)
  {
    std::set<short>::iterator it;

    it     = taxon_set.begin();

    unsigned i;
    os << "Cluster MSA: " << pos_start+1 << " " << pos_end << std::endl;
    for (i=0; i<taxon_set.size(); ++i, ++it)
    {
      faststring tmp( msa[i]+pos_start, msa[i]+pos_end);
      back_recode_range(tmp.begin(), tmp.end());
      os.width(4);
      os << *it << " " << tmp << std::endl;
    }
    if (flag == 1)
    {
      if (center_computed() )
      {
	os << "Cent " << center_sequence << std::endl;
	os << "Max dist: " << max_dist_bait_to_msa << std::endl;
      }
      else
	os << "Center sequence not computed." << std::endl; 
    }
  }

};


// Stores and handles all clusters for a specific bait window: 

class Csequence_loci_cluster_collection
{
  CSequences2                                       *seqs;
  unsigned                                          pos_start;      // 0 based
  unsigned                                          pos_end;        // 0 based, after last element of range
  double                                            dist_threshold; //
  //  CTaxonNamesDictionary                             &taxon_names_dict;
  std::vector<bool>                                 valid_sequences; // Does the sequence with given index have Ns or gaps in this range?
                                                                     // Indices are the indices in the CSequences2 object.

  char**                                            msa;             // Pointer to a recoded copy of the data
                                                                     // pos_start and pos_end are the coordinates
                                                                     // in these sequences.

  std::list<Csequence_cluster_and_center_sequence> list_of_clusters;
  unsigned                                           number_of_valid_sequences;



 public:
  CSequences2 *get_pointer_to_sequences()
  {
    return seqs;
  }

  //  CTaxonNamesDictionary& get_taxonNamesDictionary()
  //  {
  //    return taxon_names_dict;
  //  }

  // Constructor:
 Csequence_loci_cluster_collection(unsigned P_pos_start, unsigned P_pos_end, CSequences2 *P_seqs,
				   //				   CTaxonNamesDictionary &P_taxon_names_dict,
				   double P_dist_threshold,
				   char** P_msa):
  seqs(P_seqs), pos_start(P_pos_start), pos_end(P_pos_end), dist_threshold(P_dist_threshold),
    //taxon_names_dict(P_taxon_names_dict),
    valid_sequences(P_seqs->GetTaxaNum()), msa(P_msa)
  {
    //     std::cerr << "Call to constructor: Csequence_loci_cluster_collection "
    // 	      << pos_start+1 << " " << pos_end << " " << dist_threshold << std::endl;

    unsigned          i, N=seqs->GetTaxaNum();
    CSequence_Mol*    seq;

    number_of_valid_sequences = 0;

// IMPORTANT: If this code should be uncommented, code further down needs to be uncommented as well.
/*     if (global_VERBOSITY > 200) */
/*     { */
/*       std::cerr << "#1 " << "Constructor: Csequence_loci_cluster_collection, range: " << pos_start+1 << " " << pos_end << std::endl; */
/*       std::cerr << "#1 " << "Valid: "; */
/*     } */ 

    for (i=0; i<N; ++i)
    {
      seq = seqs->get_seq_by_index(i);
      
      // Has no gaps, Ns in range?
      // In principle we have already conducted this check in check_full_coverage_and_required_taxa_in_range.
      // This information is not shared a global scope, so it is not available any more.
      valid_sequences[i] = !seq->range_contains_gaps_or_Ns(pos_start, pos_end);
      if (valid_sequences[i])
      {
	//	std::cerr << i << ", ";
	++ number_of_valid_sequences;
      }
    }
    //    std::cerr << std::endl;

    // DEBUG code: Check output of back translation 
    if (global_VERBOSITY > 100)
    {
      std::cout << "#1 MSA:" << std::endl;
      print_msa_backrecoded(std::cout);
    }

    // Valid distances among sequences are in the range 0..1.
    // With a larger value we indicate a problem, e.g. Ns or gaps in the sequences.
    // In order not to cluster sequences that contain Ns or gaps, we 
    // set the dist_threshold to a maximum value of 10.
    // Distances with a value >10 will not be clustered.
    // A good value for distances of sequences that contain gaps or Ns is 100.
    if (dist_threshold > 1)
      dist_threshold = 10;
  }

  void print_msa_backrecoded(std::ostream &os)
  {
    unsigned          i, N=seqs->GetTaxaNum();

    os << "Exon-locus MSA:" << pos_start+1 << " " << pos_end << std::endl;
    for (i=0; i < N; ++i)
    {
      faststring tmp(	msa[i]+pos_start, msa[i]+pos_end);
      back_recode_range(tmp.begin(), tmp.end());
      os.width(4);
      if (valid_sequences[i])
      {
	//	std::cerr.setf(ios::right);
	os << i << "+ " << tmp << std::endl;
      }
      else
      {
	os << i << "- " << tmp << std::endl;
      }
    }
  }


  void print_msa_backrecoded_all(std::ostream &os, int flag = 0)
  {
    print_msa_backrecoded(os);

    std::list<Csequence_cluster_and_center_sequence>::iterator it, it_end;    
    it       = list_of_clusters.begin();
    it_end   = list_of_clusters.end();
    while (it != it_end)
    {
      it->print_msa_backrecoded(os, flag);
      ++it;
    }
  }


  ~Csequence_loci_cluster_collection()
  {
  }

  // Implemented in Csequence_cluster_and_center_sequence.cpp
  void cluster_locus(faststring &, const faststring &);

  // Implemented in Csequence_cluster_and_center_sequence.cpp

  void compute_center_for_all_clusters(char);

  void print(std::ostream &os)
  {
    std::list<Csequence_cluster_and_center_sequence>::iterator it, it_end;

    os << "Collection of clusters and their center sequences: "
       << pos_start+1 << " " << pos_end
       << " dist_threshold " << dist_threshold
       << std::endl;

    it     = list_of_clusters.begin();
    it_end = list_of_clusters.end();

    while (it != it_end)
    {
      it->print(os);
      ++it;
    }
  }

  void get_numbers_of_valid_sequences_number_of_clusters(unsigned short &valid_sequences, unsigned short &number_of_clusters)
  {
    valid_sequences = number_of_valid_sequences;
    number_of_clusters = list_of_clusters.size();
  }

  unsigned append_center_sequences(std::vector<faststring> &v, std::vector<double> &vCG, std::vector<double> &vMaxDist)
  {
    unsigned count = 0;
    std::list<Csequence_cluster_and_center_sequence>::iterator it, it_end;

    it     = list_of_clusters.begin();
    it_end = list_of_clusters.end();

    while (it != it_end)
    {
      v.push_back(  it->get_center_sequence() );
      vCG.push_back(it->get_CG() );

      //      std::cerr << "Getting max_dist_bait_to_msa: " << it->get_max_dist_bait_to_msa() << std::endl;

      vMaxDist.push_back(it->get_max_dist_bait_to_msa() );
      ++it;
      ++count;
    }
    return count;
  }


};



class CBait_region_information
{
  faststring     gene_name;
  short          exon_number;
  short          number_of_exons;
  unsigned       start_pos;
  //  unsigned short bait_length;
  //  unsigned short offset;
  unsigned short number_of_baits_in_region;

  std::vector<Csequence_loci_cluster_collection *> loci_this_bait_region;

  // The total_number_of_valid_bait_windows_over_all_sequences_this_bait_region is effectively,
  // the size of the evidence, or the number of sequences the baits are based on, after adding all contributions in a bait region.
  unsigned short total_number_of_valid_bait_windows_over_all_sequences_this_bait_region;
  unsigned short total_number_of_clusters_this_region;


 public:
 CBait_region_information(faststring &P_gene_name, short P_exon_number, short P_number_of_exons, unsigned P_start_pos, //    unsigned short P_bait_length, unsigned short P_offset, 
			  unsigned short P_number_of_baits_in_region):
  gene_name(P_gene_name), exon_number(P_exon_number), number_of_exons(P_number_of_exons),
    start_pos(P_start_pos), //bait_length(P_bait_length), offset(P_offset),
    number_of_baits_in_region(P_number_of_baits_in_region),
    total_number_of_valid_bait_windows_over_all_sequences_this_bait_region(0),
    total_number_of_clusters_this_region(0)
  {}

  void add_info__cluster_locus(Csequence_loci_cluster_collection *p,
			       unsigned short add_to_number_of_sequences,
			       unsigned short add_to_number_of_baits)
  {
    loci_this_bait_region.push_back(p);
    total_number_of_valid_bait_windows_over_all_sequences_this_bait_region += add_to_number_of_sequences;
    total_number_of_clusters_this_region                                      += add_to_number_of_baits;
  }


  void get_vector_of_bait_sequences(std::vector<faststring> &v, std::vector<unsigned> &bait_counts, std::vector<double> &vCG, std::vector<double> &vMaxDist)
  {
    v.clear();
    vCG.clear();
    vMaxDist.clear();

    if (total_number_of_clusters_this_region != 0)
    {
      v.reserve(total_number_of_clusters_this_region);
      vCG.reserve(total_number_of_clusters_this_region);
      vMaxDist.reserve(total_number_of_clusters_this_region);
    }
    else
    {
      v.reserve(10*number_of_baits_in_region);
      vCG.reserve(10*number_of_baits_in_region);
      vMaxDist.reserve(10*number_of_baits_in_region);
    }

    int i, N=loci_this_bait_region.size();
    unsigned count;
    for (i=0; i<N; ++i)
    {
      count = loci_this_bait_region[i]->append_center_sequences(v, vCG, vMaxDist);
      bait_counts.push_back(count);
    }
  }


  void print_baits_result(std::ostream &os, const faststring &filename)
  {
    std::vector<faststring> v;
    std::vector<double>     vCG;
    std::vector<double>     vMaxDist;
    std::vector<unsigned>   bait_counts;

    get_vector_of_bait_sequences(v, bait_counts, vCG, vMaxDist);

    os.setf(std::ios::fixed);
    os.precision(3);

    os << filename        << "\t"
       << gene_name       << "\t"
       << exon_number     << "\t"
       << number_of_exons << "\t"
       << start_pos+1     << "\t"
       << total_number_of_valid_bait_windows_over_all_sequences_this_bait_region << "\t"
       << total_number_of_clusters_this_region
       << "\t";

    //    print_container(os, bait_counts, "(",",",")\t");

    unsigned i, N;
    N = v.size();

    if (v.size() != vCG.size())
    {
      std::cerr << "Internal error: v and vCG have different sizes. Sizes are " << v.size() << " and " << vCG.size() << std::endl;
      exit(-1);
    }

    if (N == 0)
    {
      std::cerr << "ERROR: NUMBER OF BAITS IS 0 for this locus. NEEDS TO BE VERIFIED." << std::endl;
      exit(-11);
    }

    int       index_bait_counts=0;
    unsigned  bait_sum=0;

    bait_sum += bait_counts[index_bait_counts];
    ++index_bait_counts;

    for (i=0; i<N-1; ++i)
    {
      if (i == (bait_sum-1) )
      {
	os << index_bait_counts << " " << v[i] << " " << vCG[i] << " " << vMaxDist[i] << ",";

	 bait_sum += bait_counts[index_bait_counts];
	 ++index_bait_counts;
      }
      else
      {
	os << index_bait_counts << " " << v[i] << " " << vCG[i] << " " << vMaxDist[i] << ",";
      }
    }
    os  << index_bait_counts << " " << v[i] << " " << vCG[i] << " "  << vMaxDist[i] << std::endl << std::flush;
  } // END void print_baits_result(std::ostream &os, const faststring &filename)
}; // END class CBait_region_information


#endif
