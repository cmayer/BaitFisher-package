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
#include "Csequence_cluster_and_center_sequence.h"
#include "DEBUG_STUFF.h"

using namespace std;

CDistance_Cluster *p_global_distcluster;
bool               global_bad_dist_this_range;


void call_back_distance(short i1, short i2, double dist)
{
  if (dist < 0)
    global_bad_dist_this_range = true;
  else
    p_global_distcluster->add(i1, i2, dist);
}



void Csequence_loci_cluster_collection::cluster_locus(faststring & region_name, const faststring &filename)
{
  // TODO:
  // We have to check how many (valid) sequences we have.
  // If this number is just 0 or 1, the clustering will fail and we need a workaround.

  // The number of valid sequences has been determined in the constructor of this object.
  if (number_of_valid_sequences < 2)
  {
    if (number_of_valid_sequences == 0)
    {
      if (global_VERBOSITY >= VERBOSITY_MOREWARNINGS)
      {
	cout << "WARNING: Found 0 valid sequences in bait window in Csequence_loci_cluster_collection::cluster_locus()." << endl;
	cout << "         Alignment file:          " << filename    << endl;
	cout << "         Region/gene name:        " << region_name << endl;
	cout << "         Within gene coordinates: " << pos_start+1 << "-" << pos_end << endl;
      }
    }
    if (number_of_valid_sequences == 1)
    {
      // TODO: Report this only if 
      if (global_VERBOSITY  >= VERBOSITY_MOREWARNINGS)
      {
	cout << "Note: Found 1 valid sequence in bait window in Csequence_loci_cluster_collection::cluster_locus(). Clustering will be skipped. Bait will be constructed from single sequence." << endl;
	cout << "      Alignment file:          " << filename    << endl;
	cout << "      Region/gene name:        " << region_name << endl;
	cout << "      Within gene coordinates: " << pos_start+1 << "-" << pos_end << endl;
      }

      // TODO: Uncomment and adapt:
      Csequence_cluster_and_center_sequence empty_cluster_and_center_sequence;

      // Add an empty object
      list_of_clusters.push_back(empty_cluster_and_center_sequence);
      
      if (global_VERBOSITY > 2000)
       {
	 cout << "Object that will be assigned: *******" << endl;
	 list_of_clusters.back().print(cout);
	 cout << "*************************************" << endl;
       }
       // Fill the object
      set<short> seq_set;
      unsigned   i, N = valid_sequences.size();
      for (i=0; i < N; ++i)
      {
	if (valid_sequences[i])
	  seq_set.insert(i);
      }
      list_of_clusters.back().assign(pos_start, pos_end, 0, seq_set, this, msa);
      if (global_VERBOSITY > 2000)
      {
	cout << "Object after having been assigned:***" << endl;
	list_of_clusters.back().print(cout);
	cout << "*************************************" << endl;
      }
    } // END if (number_of_valid_sequences == 1)
    // Prints data for all clusters in this locus. Center sequences are not computed yet.
    if (global_VERBOSITY > 2000)
    {
      print(cout);
    }
    return;
  } // END if (number_of_valid_sequences < 2)

  // Create cluster object:
  CDistance_Cluster cluster_object;

  // We initialise the global pointer to the cluster_object.
  // This pointer is used by the global call_back_distance function,
  // which we pass to the CSequences2 class where it is used to
  // fill the cluster_object with the distances of the sequences.
  p_global_distcluster = &cluster_object;

  // One problem exists with the sequences that contain gaps or Ns in the window we investigate.
  // The clustering algorithm assumes that all sequence numbers are present during clustering,
  // so we cannot skip sequences. So. we have to add distances for all sequences in the sequences objects.
  // We also do not want to copy the sequences to a new CSequences2 object. Then we have to keep track
  // of the new sequence number, and we have global unique numbers, sequence numbers in the exon,
  // sequence numbers in the window.
  // The best compromise is to add very large distances to the clustering algorithm and to remove
  // the "non-present" sequences after the clustering has been done.
  // Due to the large distances we assign to all sequence pairs in which one sequence is not valid, non of
  // the invalid sequences will be clustered with any other sequence.

  // Valid distances among sequences are in the range 0..1.
  // With a larger value we indicate a problem, e.g. Ns or gaps in the sequences.
  // In order not to cluster sequences that contain Ns or gaps, we 
  // set the dist_threshold to a maximum value of 10.
  // Distances with a value >10 will not be clustered.
  // A good value for distances of sequences that contain gaps or Ns is 100.

  // Add distances:
  seqs->get_sequences_distances(pos_start, pos_end,
				call_back_distance,
				valid_sequences, 0); 
    
  if (global_VERBOSITY > 3000)
  {
    cout << "Sorted distances in this cluster object" << endl;
    cluster_object.print_sorted(cout);
    cout << endl;
  }

  if (global_VERBOSITY >= 900)
  {
    // Print filename, region name, within gene/alignment coordinate, maximum sequence distance
    // for this starting position/locus.
    cout << "Maximum distance:\t" << filename << "\t" << region_name << "\t" << pos_start+1 << "\t"
	 << cluster_object.max_distance() << endl;
  }

  // Do clustering:
  cluster_object.run_clustering(dist_threshold);
    
  if (global_VERBOSITY > 900)
  {
    faststring cluster_string;
    cluster_object.get_cluster_string2(cluster_string);
    cout << "Clustering string: " << cluster_string << endl;
  }

  // Retrieve results: Store all valid clusters in the list of clusters "list_of_clusters" of this locus.
  set<short> seq_set;
  //  set<short> seq_set_uniq_numbers;
  set<short>::iterator it_set, it_set_end;
  double max_dist_in_clu;
  short i;
  bool  valid_cluster;

  Csequence_cluster_and_center_sequence empty_cluster_and_center_sequence;
  
  i=0;
  // For all clusters at this locus:
  while ( cluster_object.get_cluster_by_index(i, seq_set, max_dist_in_clu) )
  {
    valid_cluster = true;
    // Determine whether this is a good set of sequences in this range.
    // Bad sequences should never cluster with others with our settings. (distance=100)
    // So bad sequences should be clusters with just one sequence.
    if (seq_set.size() == 1)
    {
      it_set     = seq_set.begin();
      unsigned seq_index = *it_set;

      valid_cluster = valid_sequences[seq_index];
    }
    // else - the cluster must be valid in our setup since clusters with more than 1 element cannot contain bad sequences.


    // Do something with the cluster:
    // Debug code - more important it shows how to extract
    // the resulting clustering data
    if (global_VERBOSITY > 3000)
    {
      it_set     = seq_set.begin();
      it_set_end = seq_set.end();

      cout << "Found cluster: (";
      while (it_set != it_set_end)
      {
	cout << *it_set;
	++it_set;
	if (it_set != it_set_end)
	  cout << ",";
      }
      cout << "):" << max_dist_in_clu;
      if (valid_cluster)
	cout << endl;
      else
	cout << "(invalid)" << endl;
    }

    // Copy the clustering data to the Csequence_loci_cluster_collection data

    if (valid_cluster)
    {
      // Too bad: emplace_back is a C++11 feature
      //      list_of_clusters.emplace_back(pos_start, pos_end, max_dist_in_clu, seq_set, this);

      // Add an empty object
      list_of_clusters.push_back(empty_cluster_and_center_sequence);

      if (global_VERBOSITY > 2000)
      {
            cout << "Object that will be assigned: *******" << endl;
            list_of_clusters.back().print(cout);
            cout << "*************************************" << endl;
      }
      // Fill the object
      list_of_clusters.back().assign(pos_start, pos_end, max_dist_in_clu, seq_set, this, msa);
      if (global_VERBOSITY > 2000)
      {
	cout << "Object after having been assigned:***" << endl;
	list_of_clusters.back().print(cout);
	cout << "*************************************" << endl;
      }
    }
    ++i;

  } // END while ( cluster_object.get_cluster_by_index(i, seq_set, max_dist_in_clu) )

  // Prints data for all clusters in this locus. Center sequences are not computed yet.
  if (global_VERBOSITY > 1000)
  {
    print(cout);
  }
} // END Csequence_loci_cluster_collection::cluster_locus()

///////////////////////////////////////////////////////
///////////////// CLUSTERING //////////////////////////
///////////////////////////////////////////////////////

// Helper functions:

// Indices:
// First index:  Row, i.e. sequence
// Second index: Col, i.e. sequence position

void recode_range(char * b, char *e)
{
  while (b != e)
  {
    *b = recode_DNA(*b);
    ++b;
  }
}

void back_recode_range(char * b, char *e)
{
  while (b != e)
  {
    *b = backrecode_DNA(*b);
    ++b;
  }
}

// inline bool is_valid_recode_DNA_range(const char *p1, const char *p2)
// {
//   while (p1 != p2)
//   {
//     if (!is_valid_recode_DNA(*p1))
//       return false;
//     ++p1;
//   }
//   return true;
// }


void recode_MSA(char **sequences, int num_taxa, int num_pos)
{
  int  row;

  for (row=0; row < num_taxa; ++row)
  {
    recode_range(sequences[row], sequences[row]+num_pos);
  }
}

// Requires recoded sequences as input.
// The consensus will also be recoded.
////////////////// UPDATED/////////////////
void find_consensus(char **msa, int num_taxa, unsigned pos_start_msa,
		    unsigned pos_end_msa, char *consensus, char *isAmbig)
{
  int  row, col, col_in_msa;
  char c;

  int len = pos_end_msa - pos_start_msa;

  for (col=0, col_in_msa = pos_start_msa; col < len; ++col, ++col_in_msa)
  {
    c = msa[0][col_in_msa];  // Take symbol from row=0
    for (row=1; row < num_taxa; ++row)
    {
      c = c | msa[row][col_in_msa];  // Add all symbols of all other rows.
    }
    consensus[col] = c;
    if (recode_is_DNA_ambig(c) )
    {
      isAmbig[col] = 1; // Keep track of amibig characters: set isAmbig to true
    }
    else
    {
      isAmbig[col] = 0; // Keep track of amibig characters: set isAmbig to false
    }
  }
}

// ONLY DEBUGGING OUTPUT:
// Should work unchanged. CHECK!!!!
void print_all_variants(char *consensus, char *variant, char *isAmbig, int pos, int len)
{
  while (pos < len && !isAmbig[pos])
    ++pos;

  if (pos == len)
  {
    faststring p = variant;
    back_recode_range(p.begin(), p.end());
    cout << p << endl;
  }
  else
  {
    if (consensus[pos]&1)
    {
      variant[pos] = 1;
      print_all_variants(consensus, variant, isAmbig, pos+1, len);
    }
    if (consensus[pos]&2)
    {
      variant[pos] = 2;
      print_all_variants(consensus, variant, isAmbig, pos+1, len);
    }
    if (consensus[pos]&4)
    {
      variant[pos] = 4;
      print_all_variants(consensus, variant, isAmbig, pos+1, len);
    }
    if (consensus[pos]&8)
    {
      variant[pos] = 8;
      print_all_variants(consensus, variant, isAmbig, pos+1, len);
    }
  }
} // END find_consensus(..)




// Important: msa, consensus, and variant are expected to be passed as recoded msa.
//
// pos1 ... pos2 is the range to work on. pos2 is the index after the end of the range.
//
////////////////// UPDATED/////////////////
void find_center_sequence_exhaustive(char **msa, char *consensus,
				     char *variant, char *isAmbig,
				     int len_variant,
				     int pos_start, int num_taxa, 
				     int pos_start_msa, int pos_end_msa,
				     int *distances, faststring &best_variant,
				     int *p_best_maximum_distance)
{
  int i;
  int len_remaining = pos_end_msa - pos_start_msa;
  int pos_end       = pos_start   + len_remaining;

  // Skip all bases for which isAmbig == 0
  while (pos_start < pos_end && !isAmbig[pos_start])
  {
    ++pos_start;
    ++pos_start_msa;
  }

  if (pos_start == pos_end)
  {
    // We have a new variant:
    // Check wether this is better than the best we found so far.

    // Determine maximum distance:
    int m_dist = maximum(distances, distances+num_taxa);

    if (global_VERBOSITY > 800)
    {
      cout << "Current sequence:\nDistance: " << m_dist << endl;
      faststring p = variant;
      back_recode_range(p.begin(), p.end());
      cout << p << endl;

      print_container(cout, distances, distances+num_taxa, "distances: ", ":", "\n\n");
    }

    if (m_dist < *p_best_maximum_distance) // We found a better sequence:
    {

      // Be careful in the following assignment: "=" does not work here, since variant is not 0 terminated!!!
      best_variant.assign(variant, 0, len_variant);
      *p_best_maximum_distance = m_dist;

      // Some debug output:
      if (global_VERBOSITY > 500)
      {
	cout << "*** New best variant: ***" << endl;
	faststring tmp = best_variant;
	back_recode_range(tmp.begin(), tmp.end() );
	cout << tmp << endl;
	cout << "Distance: " << *p_best_maximum_distance << endl << endl;
      }
    }
  }
  else
  {
    int *tmp_dist = new int [num_taxa];

    // If we get here, we have an non unique "consensus at pos_start.
    // So let us try all variants which are found in the msa.
    if (consensus[pos_start]&1)
    {
      variant[pos_start] = 1;

      // Determine new distances:
      memcpy(tmp_dist, distances, sizeof(int)*num_taxa);
      // Adjust distances for this new base:
      for (i=0; i<num_taxa; ++i)
      {
	if (variant[pos_start] != msa[i][pos_start_msa])
	  ++tmp_dist[i];
      }

      find_center_sequence_exhaustive(msa, consensus, variant, isAmbig, len_variant,
				      pos_start+1, num_taxa, pos_start_msa+1, pos_end_msa,
				      tmp_dist, best_variant, p_best_maximum_distance);
    }
    if (consensus[pos_start]&2)
    {
      variant[pos_start] = 2;

      // Determine new distances:
      memcpy(tmp_dist, distances, sizeof(int)*num_taxa);
      // Adjust distances for this new base:
      for (i=0; i<num_taxa; ++i)
      {
	if (variant[pos_start] != msa[i][pos_start_msa])
	  ++tmp_dist[i];
      }

      find_center_sequence_exhaustive(msa, consensus, variant, isAmbig, len_variant,
				      pos_start+1, num_taxa, pos_start_msa+1, pos_end_msa,
				      tmp_dist, best_variant, p_best_maximum_distance);
    }
    if (consensus[pos_start]&4)
    {
      variant[pos_start] = 4;

      // Determine new distances:
      memcpy(tmp_dist, distances, sizeof(int)*num_taxa);
      // Adjust distances for this new base:
      for (i=0; i<num_taxa; ++i)
      {
	if (variant[pos_start] != msa[i][pos_start_msa])
	  ++tmp_dist[i];
      }

      find_center_sequence_exhaustive(msa, consensus, variant, isAmbig, len_variant,
				      pos_start+1, num_taxa, pos_start_msa+1, pos_end_msa,
				      tmp_dist, best_variant, p_best_maximum_distance);
    }
    if (consensus[pos_start]&8)
    {
      variant[pos_start] = 8;

      // Determine new distances:
      memcpy(tmp_dist, distances, sizeof(int)*num_taxa);
      // Adjust distances for this new base:
      for (i=0; i<num_taxa; ++i)
      {
	if (variant[pos_start] != msa[i][pos_start_msa])
	  ++tmp_dist[i];
      }

      find_center_sequence_exhaustive(msa, consensus, variant, isAmbig, len_variant,
				      pos_start+1, num_taxa, pos_start_msa+1, pos_end_msa,
				      tmp_dist, best_variant, p_best_maximum_distance);
    }
    delete [] tmp_dist;
  } // END else of if (pos1 == pos2)
} // END void find_center_sequence_exhaustive

// distances must be initialised by caller.
void find_center_sequence_heuristic(char **msa,     char *consensus,
				    char *center,  char *isAmbig,
				    unsigned len_center,
				    unsigned pos_start, unsigned num_taxa, 
				    unsigned pos_start_msa, unsigned pos_end_msa,
				    unsigned *distances,
				    unsigned *p_best_maximum_distance)
{
  static char recode_A = recode_DNA('A');
  static char recode_C = recode_DNA('C');
  static char recode_G = recode_DNA('G');
  static char recode_T = recode_DNA('T');

  unsigned len_remaining = pos_end_msa - pos_start_msa;
  unsigned pos_end       = pos_start   + len_remaining;

  // pos_start walks through isAmbig and consensus.
  // pos_start_msa walks through the msa. 

  // The following static variables are initialised only when calling this function for the first time.
  static unsigned max_num_taxa = 20;
  static unsigned *dist_A = new unsigned [max_num_taxa];
  static unsigned *dist_C = new unsigned [max_num_taxa];
  static unsigned *dist_G = new unsigned [max_num_taxa];
  static unsigned *dist_T = new unsigned [max_num_taxa];

  if (num_taxa > max_num_taxa)
  {
    max_num_taxa = num_taxa; // max_num_taxa takes on the new maximum value. It will be used as size when reallocating new arrays.

    if (dist_A != NULL)
    {
      delete [] dist_A;
      delete [] dist_C;
      delete [] dist_G;
      delete [] dist_T;

      dist_A = new unsigned [max_num_taxa];
      dist_C = new unsigned [max_num_taxa];
      dist_G = new unsigned [max_num_taxa];
      dist_T = new unsigned [max_num_taxa];
    }
  }

  unsigned distances_sum_squares;
  unsigned dist_A_sum_squares;
  unsigned dist_C_sum_squares;
  unsigned dist_G_sum_squares;
  unsigned dist_T_sum_squares;

//   unsigned dist_max_A;
//   unsigned dist_max_C;
//   unsigned dist_max_G;
//   unsigned dist_max_T;
//   unsigned dist_max;

  unsigned j;

  // Search for first non-invariant site:
  while (pos_start < pos_end && !isAmbig[pos_start])
  {
    ++pos_start;
    ++pos_start_msa;
  }

  if (pos_start == pos_end) // This is the case if all sites are invariant.
  {
    return;
  }

  // pos_start is first non-invariant site:
  {
    unsigned A=0, C=0, G=0, T=0;
    char c;

    for (j=0; j<num_taxa; ++j)
    {
      if (msa[j][pos_start_msa] == recode_A )
	++A;
      else if (msa[j][pos_start_msa] == recode_C )
	++C;
      else if (msa[j][pos_start_msa] == recode_G )
	++G;
      else if (msa[j][pos_start_msa] == recode_T )
	++T;
    }

    if (A > C)
    {
      if (G > T)
      {
	if (A > G)
	  c = recode_A;
	else
	  c = recode_G;
      }
      else // (G <= T)
      {
	if (A > T)
	  c = recode_A;
	else
	  c = recode_T;
      }
    }
    else // (A <= C)
    {
      if (G > T)
      {
	if (C > G)
	  c = recode_C;
	else
	  c = recode_G;
      }
      else // (G <= T)
      {
	if (C > T)
	  c = recode_C;
	else
	  c = recode_T;
      }
    }
    center[pos_start] = c; // First position has been determined.
    distances_sum_squares = 0;
    //    dist_max = 0;

    for (j=0; j<num_taxa; ++j)
    {
      if (center[pos_start] != msa[j][pos_start_msa])
	++distances[j];
      distances_sum_squares += distances[j]*distances[j];
//       if (dist_max < distances[j])
// 	dist_max = distances[j];
    }
  } // pos_start is first non-invariant site:

  ++pos_start;
  ++pos_start_msa;

  while (pos_start < pos_end)
  {
    if (isAmbig[pos_start]) // He now search for best nucleotide in center.
    {
      // Try 'A'
      if (consensus[pos_start]&recode_A)
      {
	memcpy(dist_A, distances, sizeof(unsigned)*num_taxa);
	center[pos_start] = recode_A;
	dist_A_sum_squares = 0;
	for (j=0; j<num_taxa; ++j)
	{
	  if (center[pos_start] != msa[j][pos_start_msa])
	    ++dist_A[j];
	  dist_A_sum_squares += dist_A[j]*dist_A[j];
	}
      }
      else
	dist_A_sum_squares = -1u;
      // END // Try 'A'
      // Try 'C'
      if (consensus[pos_start]&recode_C)
      {
	memcpy(dist_C, distances, sizeof(unsigned)*num_taxa);
	center[pos_start] = recode_C;
	dist_C_sum_squares = 0;
	for (j=0; j<num_taxa; ++j)
	{
	  if (center[pos_start] != msa[j][pos_start_msa])
	    ++dist_C[j];
	  dist_C_sum_squares += dist_C[j]*dist_C[j];
	}
      }
      else
	dist_C_sum_squares = -1u;
 // END // Try 'C'
      // Try 'G'
      if (consensus[pos_start]&recode_G)
      {
	memcpy(dist_G, distances, sizeof(unsigned)*num_taxa);
	center[pos_start] = recode_G;
	dist_G_sum_squares = 0;
	for (j=0; j<num_taxa; ++j)
	{
	  if (center[pos_start] != msa[j][pos_start_msa])
	    ++dist_G[j];
	  dist_G_sum_squares += dist_G[j]*dist_G[j];
	}
      }
      else
	dist_G_sum_squares = -1u;
      // END // Try 'G'
      // Try 'T'
      if (consensus[pos_start]&recode_T)
      {
	memcpy(dist_T, distances, sizeof(unsigned)*num_taxa);
	center[pos_start] = recode_T;
	dist_T_sum_squares = 0;
	for (j=0; j<num_taxa; ++j)
	{
	  if (center[pos_start] != msa[j][pos_start_msa])
	    ++dist_T[j];
	  dist_T_sum_squares += dist_T[j]*dist_T[j];
	}
      }
      else
	dist_T_sum_squares = -1u;
      // END // Try 'T'

      // Which "Try" is the best?
      {
	// Smallest dist_x_sum_squares is the best:
	if (dist_A_sum_squares  < dist_C_sum_squares) // A is better than C
	{
	  if (dist_G_sum_squares  < dist_T_sum_squares) // G is better than T
	  {
	    if (dist_A_sum_squares  < dist_G_sum_squares) // A is better than G
	    {
	      // A is the best: 
	      center[pos_start] = recode_A;
	      memcpy(distances, dist_A, sizeof(unsigned)*num_taxa);
	    }
	    else
	    {
	      // G is the best:
	      center[pos_start] = recode_G;
	      memcpy(distances, dist_G, sizeof(unsigned)*num_taxa);
	    }
	  }
	  else  // T is better than or equal to G
	  {
	    if (dist_A_sum_squares  < dist_T_sum_squares) // A is better than T
	    {
	      // A is the best
	      center[pos_start] = recode_A;
	      memcpy(distances, dist_A, sizeof(unsigned)*num_taxa);
	    }
	    else
	    {
	      // T is the best
	      center[pos_start] = recode_T;
	      memcpy(distances, dist_T, sizeof(unsigned)*num_taxa);
	    }
	  }
	}
	else // C is better than or equal to A
	{
	  if (dist_G_sum_squares  < dist_T_sum_squares) // G is better than T
	  {
	    if (dist_C_sum_squares  < dist_G_sum_squares) // C is better than G
	    {
	      // C is the best: 
	      center[pos_start] = recode_C;
	      memcpy(distances, dist_C, sizeof(unsigned)*num_taxa);
	    }
	    else
	    {
	      // G is the best:
	      center[pos_start] = recode_G;
	      memcpy(distances, dist_G, sizeof(unsigned)*num_taxa);
	    }
	  }
	  else  // T is better than or equal to G
	  {
	    if (dist_C_sum_squares  < dist_T_sum_squares) // C is better than T
	    {
	      // C is the best
	      center[pos_start] = recode_C;
	      memcpy(distances, dist_C, sizeof(unsigned)*num_taxa);
	    }
	    else
	    {
	      // T is the best
	      center[pos_start] = recode_T;
	      memcpy(distances, dist_T, sizeof(unsigned)*num_taxa);
	    }
	  }
	} // END C is better than or equal to A
      } // END Which "Try" is the best?
    } // END ELSE (!isAmbig[pos_start])
    ++pos_start;
    ++pos_start_msa;
  } // END while (pos_start < pos_end)
  // We have the greedy center sequence.
 *p_best_maximum_distance = maximum(distances, distances+num_taxa);
} // END find_center_sequence_heuristic




// TODO ///////////////////

// Requires:
// Multiple sequence alignment that is already recoded.
// References to (empty) faststrings for the center.
// Number of taxa
// Start and end indices of range for which the center sequence shall be computed.
// The center sequence is returned in recoded form.
void find_center_sequence_for_msa_exhaustive(char **msa, faststring &center_seq,
					     int num_taxa,
					     int pos_start_msa, int pos_end_msa,
					     double &mdist)
{
  faststring cons;
  faststring isAmbig;
  faststring variant;
  faststring best_variant;
  int        best_maximum_distance = INT_MAX;
  // The equality len_variant = pos_end - pos_start is only valid before ever incrementing
  // pos_start, which is always true in this function.
  int        len_variant = pos_end_msa - pos_start_msa;
  int        pos_start   = 0;
  //  int        pos_end     = len_variant;


  cons.assign(len_variant,    ' ');
  isAmbig.assign(len_variant, ' ');

  // Here we can pass .data() since we do not assign it
  find_consensus(msa, num_taxa, pos_start_msa, pos_end_msa, cons.begin(), isAmbig.begin() );
  variant = cons;

  int *dist_array = new int [num_taxa];
  memset(dist_array, 0, num_taxa*sizeof(int));
//   for (i=0; i< num_taxa; ++i)
//   {
//     dist_array[i] = 0;
//   }

   //  print_all_variants(cons.data(), variant.data(), isAmbig.data(), 0, len);
  
//   if (global_VERBOSITY > 5)
//   {
//     cout << "Start time: " << time(NULL) << endl;
//   }

  find_center_sequence_exhaustive(msa, cons.data(), variant.data(), isAmbig.data(),
				  len_variant,
				  pos_start, num_taxa,
				  pos_start_msa, pos_end_msa,
				  dist_array, best_variant,
				  &best_maximum_distance);

//   if (global_VERBOSITY > 5)
//   {
//     cout << "End time: " << time(NULL) << endl;
//   }

  center_seq = best_variant;
  // TODO remove if is clear this cannot happen
  if (len_variant == 0)
  {
    cerr << "ERROR: len_variant == 0 in find_center_sequence_for_msa_exhaustive. Please report this bug." << endl;
    exit(-54);
  }
  mdist = (double)best_maximum_distance/len_variant;

  delete [] dist_array;
  // 
  //  back_recode_range(center_seq.begin(), center_seq.end() );
} // END find_center_sequence_for_msa_exhaustive

// Requires:
// Multiple sequence alignment that is already recoded.
// References to (empty) faststrings for the center.
// Number of taxa
// Start and end indices of range for which the center sequence shall be computed.
// The center sequence is returned in recoded form.
void find_center_sequence_for_msa_heuristic(char **msa, faststring &center_seq,
					    int num_taxa,
					    int pos_start_msa, int pos_end_msa,
					    double &mdist)
{
  faststring cons;
  faststring isAmbig;
  faststring variant;
  unsigned   best_maximum_distance = 0;
  // The equality len_variant = pos_end - pos_start is only valid before ever incrementing
  // pos_start, which is always true in this function.
  unsigned   len_variant = pos_end_msa - pos_start_msa;
  unsigned   pos_start   = 0;
  //  unsigned   pos_end     = len_variant;


  cons.assign(len_variant,    ' ');
  isAmbig.assign(len_variant, ' ');

  // Here we can pass .data() since we do not assign it
  find_consensus(msa, num_taxa, pos_start_msa, pos_end_msa, cons.begin(), isAmbig.begin() );
  variant = cons;

  unsigned *dist_array = new unsigned [num_taxa];
  memset(dist_array, 0, num_taxa*sizeof(int));
//   for (i=0; i< num_taxa; ++i)
//   {
//     dist_array[i] = 0;
//   }

   //  print_all_variants(cons.data(), variant.data(), isAmbig.data(), 0, len);

//   if (global_VERBOSITY > 5)
//   {
//     cout << "Start time: " << time(NULL) << endl;
//   }

  find_center_sequence_heuristic(msa, cons.data(), variant.data(), isAmbig.data(),
				 len_variant,
				 pos_start, num_taxa,
				 pos_start_msa, pos_end_msa,
				 dist_array, &best_maximum_distance);

//   if (global_VERBOSITY > 5)
//   {
//     cout << "End time: " << time(NULL) << endl;
//   }

//   for (i=0; i< num_taxa; ++i)
//   {
//     if (dist_array[i] > best_maximum_distance)
//       best_maximum_distance = dist_array[i];
//   }
//  best_maximum_distance = maximum(distances, distances+num_taxa);

  center_seq = variant;

  // TODO remove is clear this cannot happen
  if (len_variant == 0)
  {
    cerr << "ERROR: len_variant == 0 in find_center_sequence_for_msa_heuristic. Please report this bug." << endl;
    exit(-55);
  }

  mdist = (double)best_maximum_distance/len_variant;

  if (global_VERBOSITY > 100)
  {
    cout << "Reporting mdist for a center sequence with starting position: " << pos_start << " " << mdist << endl;
  }
  // 
  //  back_recode_range(center_seq.begin(), center_seq.end() );
  delete [] dist_array;
} // END find_center_sequence_for_msa_heuristic


void Csequence_cluster_and_center_sequence::compute_center(char center_computation_mode)
{
  // The msa is prepared for the center algorithm:
  // We still need: the "consensus" and the "ambig" strings.

  int num_sequences = taxon_set.size();

//   if (global_VERBOSITY > 5)
//   {
//     cout << "Computing center sequence: " << endl << flush;
//   }

  if (num_sequences == 1)
  {
    center_sequence.assign(msa[0]+pos_start, msa[0]+pos_end);
    max_dist_bait_to_msa = 0;
  }
  else // We have to compute a center sequence:
  {
    if (center_computation_mode == 0) // exhaustive
    {
      find_center_sequence_for_msa_exhaustive(msa, center_sequence, num_sequences, pos_start, pos_end, max_dist_bait_to_msa);
    }
    else if (center_computation_mode == 1) // heuristic
    {
      find_center_sequence_for_msa_heuristic(msa, center_sequence, num_sequences, pos_start, pos_end, max_dist_bait_to_msa);
    }
    else // Comparison mode
    {
      std::set<short>::iterator it;
      it     = taxon_set.begin();

      {
	unsigned i;
	cout << "The MSA for the result below:" << endl;
	cout << "Time: " << time(NULL) << endl;
	
	for (i=0; i<taxon_set.size(); ++i, ++it)
	{
	  faststring tmp( msa[i]+pos_start, msa[i]+pos_end);
	  back_recode_range(tmp.begin(), tmp.end());
	  cout.width(4);
	  cout << *it << ":" << tmp << std::endl;
	}
      }

      unsigned t1, t2, t3;
      faststring c1, c2;
      double md1, md2;
      t1 = time(NULL);
      find_center_sequence_for_msa_exhaustive(msa, c1, num_sequences, pos_start, pos_end, md1);
      t2 = time(NULL);
      find_center_sequence_for_msa_heuristic(msa, c2, num_sequences, pos_start, pos_end, md2);
      t3 = time(NULL);
    
      // This assignment must come before back recoding c1.
      // Otherwise we back recode twice, which is not correct.
      center_sequence = c1;

      back_recode_range(c1.begin(), c1.end());
      back_recode_range(c2.begin(), c2.end());

      cout << "Cent1:     " << c1    << endl;
      cout << "max_dist1: " << md1   << endl;
      cout << "Tdiff1:    " << t2-t1 << endl;
      
      cout << "Cent2:     " << c2    << endl;
      cout << "max_dist2: " << md2   << endl;
      cout << "Tdiff2:    " << t3-t2 << endl;


      max_dist_bait_to_msa = md1;

      if (md1 != md2)
      {
	cout << "Deviating distances found: " << md2-md1 << " " << (md2-md1)*120 << endl;
      }
      else 
      {
	cout << "Same distance." << endl;
      }
    } // END  else // Comparison mode
  } // END  else // We have to compute a center sequence:

  //  cout << "LOG: pos: " << pos_start << " mis " << max_dist_bait_to_msa*120 << endl;

  // We have to back recode the center sequence:
  back_recode_range(center_sequence.begin(), center_sequence.end() );
  
  if (!center_sequence.empty())
  {
    unsigned CG, AT;
    
    CG_AT_content_in_region(center_sequence.begin_str(), center_sequence.end_str(), AT, CG);
    double dCG=CG, dAT=AT;
    CG_content = (dCG / (dCG+dAT));
  }
}


void Csequence_loci_cluster_collection::compute_center_for_all_clusters(char center_computation_mode)
{
  std::list<Csequence_cluster_and_center_sequence>::iterator it_list     =  list_of_clusters.begin();
  std::list<Csequence_cluster_and_center_sequence>::iterator it_list_end =  list_of_clusters.end();

  //  unsigned t1, t2;

  while (it_list != it_list_end)
  {
//     t1 = time(NULL);
//     cout << "T1: " << t1 << endl;
    it_list->compute_center(center_computation_mode);
//     t2 = time(NULL);
//     cout << "T2: " << t2 << endl;
//     cout << "Tdiff: " << t2-t1 << endl;
    

    if (!it_list->center_computed() )
    {
      cerr << "ERROR: Internal error: Center sequence not computed even though it was requested! Number of sequences in cluster: " << it_list->get_number_of_elements_in_cluster() << endl;
    }

    ++it_list;
  }
  
}
