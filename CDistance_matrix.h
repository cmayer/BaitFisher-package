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
#ifndef CDISTANCE_MATRIX_H
#define CDISTANCE_MATRIX_H

#include <iostream>
#include <map>
#include "faststring2.h"
#include <list>
#include <iterator>
#include <set>


typedef unsigned long index_distance_map;

//************************************************
// Usage-Restrictions:
//
// Minimum number of entries:
// CDistanceCollection: Number of objects >= 0. Note: 0 and 1 Objects cannot be distinguished.
// CDistance_Cluster:   Number of objects > 1.  2 Objects: remains to be tested.

//************************************************


//************************************************
// Part I: Distance indices
// In order to index pairs of sequences we do the following:
// If we have two sequence numbers we assume they are smaller than short integers (16 bits).
// We can combine the two short numbers to a 32 bit index of both. This index is unique if
// we put the smaller sequence index into the higher bits and the larger index into the lower bits.
// Distance indices are of type unsigned and must be 32 bits in size in order to work.
//
// The following routines handle the index conversions:
// makeIndex creates an unsigned combined index from 2 short numbers.
// first_index and second_index extract the first and the second index of a pair.
//
// Note: Initial indices of sequnces: These should be the indices of the sequnces in the file.
//       They have to be unique in the range 0...taxaNum-1.
//       During the clustering we create new indices of cluster nodes.
//       In contrast to the distance indices, cluster indices are of type short.
//
// How the clustering is done:
// cluster_node is the main data structure in the clustering method.
// The hierachy of clustering steps are stored in the cluster_hierarchy,
// a vector of cluster nodes.
// The cluster_hierarchy has a very special form:
// The index of each node coincides with the this_nodes_index.
// This makes it efficient to obtain the cluster hierachy without
// having to search for the correct index.
// 
//
//==============
// Algorithm:
//==============
// We start by adding all pairwise distances to the dist_map.
// The indices of the pairwise distances are added to the current_cluster_indices
// as well as the all_cluster_indices.
// 
// Initially the pairwise distances are added as terminal nodes to the cluster_hierarchy.
// This is done by adding the indices of the pairwise distances to the cluster_hierarchy.
// Internally, they have children set to -1 and distance 0, indicating that they
// do not represent clustering steps.
// Note: The cluster_hierarchy is only used to keep track of the clustering steps.
//       The next step in the algorithm is determined from the sorted_distances,
//       a list of iterators to the current_cluster_indices. The list of 
//       sorted_distances keeps a list of iterators in accending order of the 
//       distances associated with the index. 
// Later we add cluster_nodes that have other cluster_nodes as children.
// Distances that have been clustered are marked for deletion and are erased permanently later. (For details, read the code.)
//************************************************

inline unsigned makeIndex(short i1, short i2)
{
  unsigned res;
  if (i1 > i2)
  {
    res = i2;
    res <<= 16;
    res += i1;
    return res;
  }
  else
  {
    res = i1;
    res <<= 16;
    res += i2;
    return res;
  }
}

inline short first_index(unsigned i)
{
  i >>= 16;
  return (short) i;
}

inline short second_index(unsigned i)
{
  return (short) i;
}



// CDistanceCollection collects distance of indexed objects, i.e. objects are referred to by their index.
// Distances are stored in a map, which is useful (i) if new objects are introduced regularly, or (ii) if not
// all distances are specified. 

class CDistanceCollection
{
 public:
  typedef std::map<index_distance_map, double> CDistanceCollection_map; 
  typedef std::map<index_distance_map, double>::iterator CDistanceCollection_map_iterator; 

 protected:
  CDistanceCollection_map dist_map;

 public:

  void clear()
  {
    dist_map.clear();
  }

  void add(short i1, short i2, double d)
  {
    dist_map[makeIndex(i1, i2)] = d;
  }

  void print(std::ostream &os)
  {
    CDistanceCollection_map_iterator it     = dist_map.begin();
    CDistanceCollection_map_iterator it_end = dist_map.end();

    unsigned u;
    double   d;

    while (it != it_end)
    {
      u = it->first;
      d = it->second;

      os << "d("<< first_index(u) << "," << second_index(u) << ")="<<d << std::endl; 
      ++it;
    }

  }

  double max_distance()
  {
    CDistanceCollection_map_iterator it     = dist_map.begin();
    CDistanceCollection_map_iterator it_end = dist_map.end();
    double d, d_max;

    if (it == it_end)
      return 0;

    d_max = it->second;

    while (it != it_end)
    {
      d = it->second;
      d_max = (d_max > d ? d_max:d);
      ++it;
    }
    return d_max;
  }


};

struct less_than_distance_iterator
{
  bool operator()(CDistanceCollection::CDistanceCollection_map_iterator it1, CDistanceCollection::CDistanceCollection_map_iterator it2)
  {
    return it1->second < it2->second; 
  }
};


struct cluster_node
{
  short this_nodes_index; // The index of this node in the cluster_hierarchy vector, but also the index of the node that has the following two children.
  short child1;
  short child2;
  double dist;

  bool unused;

  // constructor
  cluster_node(short index_this, short ch1, short ch2, double d):
  this_nodes_index(index_this),child1(ch1), child2(ch2), dist(d), unused(false)
  {}

  // constructor - terminal node
  // This node only represents itself and it has no children. So these indices are -1.
  cluster_node(short index):
  this_nodes_index(index),child1(-1), child2(-1), dist(0), unused(false)
  {}

  //  cluster_node(): unused(true)
  //  {}

  bool is_terminal()
  {
    return (child1 == -1);
  }

  bool is_unused()
  {
    return unused;
  }


};



class CDistance_Cluster : public CDistanceCollection
{
  std::list<CDistanceCollection_map_iterator> sorted_distances;

  // Sets of indices that are known:
  std::set<short> current_cluster_indices;  // Indices of objects added to clustering. (Not indices of distances.)
                                            // Keeps track of remaining indices that need to be clustered.
                                            // Indices that have been clustered are removed.
  std::set<short> all_cluster_indices;      // Keeps track of all cluster indices that have ever been created.

  std::vector< cluster_node > cluster_hierarchy; 
  unsigned                verbosity;


 public:
 
  CDistance_Cluster(unsigned p_verbosity=0):verbosity(p_verbosity)
  {}

  void clear(unsigned newverbosity=-1u)
  {
    sorted_distances.clear();
    current_cluster_indices.clear();
    all_cluster_indices.clear();
    cluster_hierarchy.clear();

    CDistanceCollection::clear();
 
    if (newverbosity != -1u)
      verbosity = newverbosity;
  }

  // At the end of the clustering this is the number of clusters.
  unsigned get_num_groups_in_current_cluster()
  {
    return current_cluster_indices.size();
  }


  void sort()
  {
    less_than_distance_iterator lti;
    sorted_distances.sort(lti);
  }

  // Add terminal node to cluster hierarchy
  void addto_cluster_hierarchy(short index)
  {
    cluster_hierarchy.push_back(cluster_node(index));

    // The index of this new node should be equal to index. 
    short i = cluster_hierarchy.size();

    if ((i-1) != index)
    {
      std::cerr << "Error: Indices must be passed consecutively and 0 based to the cluster hierarchy" << std::endl;
    }
  }

  // Add terminal node to cluster hierarchy
  void addto_cluster_hierarchy(short index, short child1, short child2, double d)
  {
    cluster_hierarchy.push_back(cluster_node(index, child1, child2, d));

    // The index of this new node should be equal to index. 
    short i = cluster_hierarchy.size();

    if ((i-1) != index)
    {
      std::cerr << "Error: Indices must be passed consecutively and 0 based to the cluster hierarchy" << std::endl;
    }
  }

  void add_singleton(short i)
  {
    current_cluster_indices.insert(i);
    all_cluster_indices.insert(i); 
  }


  void add(short i1, short i2, double d)
  {
    if (i1 == i2)
    {
      std::cerr << "Distances of indices to itself cannot be added to the distance matrix. They will always be assumed to be 0" << std::endl;
      exit(-2);
    }

    std::pair<CDistanceCollection_map_iterator,bool> ret;

    ret = dist_map.insert(std::make_pair(makeIndex(i1, i2), d));
    if (ret.second)
      sorted_distances.push_back(ret.first);

    current_cluster_indices.insert(i1);
    current_cluster_indices.insert(i2);

    all_cluster_indices.insert(i1);
    all_cluster_indices.insert(i2); 
  }

  // Do not enter the indices to the sets.
  // This can be done later.
  void add_partial(short i1, short i2, double d)
  {
    if (i1 == i2)
    {
      std::cerr << "Distances of indices to itself cannot be added to the distance matrix. They will always be assumed to be 0" << std::endl;
      exit(-2);
    }

    std::pair<CDistanceCollection_map_iterator,bool> ret;

    ret = dist_map.insert(std::make_pair(makeIndex(i1, i2), d));
    if (ret.second)
      sorted_distances.push_back(ret.first);
  }


  void print_sorted(std::ostream &os)
  {
    sort();

    std::list<CDistanceCollection_map_iterator>::iterator it, it_end;

    it     = sorted_distances.begin();
    it_end = sorted_distances.end();

    unsigned u;
    double   d;

    while (it != it_end)
    {
      u = (*it)->first;
      d = (*it)->second;

      os << "d("<< first_index(u) << "," << second_index(u) << ")="<<d << std::endl; 
      ++it;
    }
  }

  void print_current_cluster_indices(std::ostream &os)
  {
    std::set<short>::iterator it_cluster_indices     = current_cluster_indices.begin();
    std::set<short>::iterator it_cluster_indices_end = current_cluster_indices.end();

    while (it_cluster_indices != it_cluster_indices_end)
    {
      os << *it_cluster_indices << ",";
      ++it_cluster_indices;
    }
    os << std::endl;

    
  }
  
  
  void run_clustering(double distance_limit)
  {
    std::set<short>::iterator it_cluster_indices     = current_cluster_indices.begin();
    std::set<short>::iterator it_cluster_indices_end = current_cluster_indices.end();
    
    // Add all indices as terminal nodes in cluster_hierarchy.
    // The cluster_hierarchy is only used to keep track of the clustering steps.
    while (it_cluster_indices != it_cluster_indices_end)
    {
      addto_cluster_hierarchy(*it_cluster_indices);
      ++it_cluster_indices;
    }

    if(verbosity > 10)
    {
      faststring cluster_string;
      std::cerr << "Initial cluster string:" << std::endl;
      get_cluster_string2(cluster_string);
      std::cerr << cluster_string << std::endl;
    }

    // The next step in the algorithm is determined from the "sorted_distances",
    // a list of iterators to the current_cluster_indices. The list of 
    // sorted_distances keeps a list of iterators in accending order of the 
    // distances associated with the index.

    std::list<CDistanceCollection_map_iterator>::iterator it, it_end;

    // sorted_distances is a list of iterators of type CDistanceCollection_map_iterator.
    // We sort this iterator list such that the iterators to the dist_map (index->distance)
    // is sorted according to the distance. Smaller distances come first.

    sort();
    it     = sorted_distances.begin();
    it_end = sorted_distances.end();

    unsigned u;
    short    index1;
    short    index2;
    short    tmp_index;

    short    new_cluster_index;

    // We continue to cluster pairs of indices while we still find distances which are less than the threshold
    while (sorted_distances.size() > 0 && (*it)->second <= distance_limit)
    {
      // The pair with the smalles distance is the first one in the sorted_distances list.
      // We cluster the pair of indices it is pointing to:
      u = (*it)->first;
      index1 = first_index(u);
      index2 = second_index(u);
     
      // Obsolete:      cluster_hierarchy.push_back(std::make_pair(index1, index2));
      // Determine a new unique index for the next cluster node.
      new_cluster_index = all_cluster_indices.size();

      if (verbosity > 10)
	std::cerr << "Outer loop -- Clustering: " << index1 << ", " << index2 << " to " << new_cluster_index << std::endl; 

      // Eliminate all distances of index pairs for which one index is equal to index1 or index2.

      it_cluster_indices     = current_cluster_indices.begin();
      it_cluster_indices_end = current_cluster_indices.end();

      // Save within distance of new cluster:
      double within_distance = dist_map[makeIndex(index1, index2)];

      // Add new cluster to hierarchy.
      addto_cluster_hierarchy(new_cluster_index, index1, index2, within_distance);
      
      // Eliminate distance pair  index1, index2:
      dist_map[makeIndex(index1, index2)] = -1; // Mark as unused

      if (verbosity > 10)
	std::cerr << "Mark unused: " << index1 << " , " << index2 << std::endl;  

      // Move through all indices that are not in the new cluster.
      while (it_cluster_indices != it_cluster_indices_end)
      {
	tmp_index = *it_cluster_indices;
	
	if (verbosity > 10)
	  std::cerr << "**Loop index: " << tmp_index << std::endl;

	if (tmp_index != index1 && tmp_index != index2)
	{
	  if (verbosity > 10)
	    std::cerr << "**Treat index: " << tmp_index << std::endl;
	
	  double dist1 = dist_map[makeIndex(index1, tmp_index)];
	  double dist2 = dist_map[makeIndex(index2, tmp_index)];

	  double max_dist  = dist1;
	  //	  short  max_index = index1;

	  if (dist2 > max_dist)
	  {
	    //	    max_index = index2;
	    max_dist = dist2;
	  }

	  // Add distance between new cluster and tmp_index
	  add_partial(tmp_index, new_cluster_index, max_dist);
	  /////// 	  dist_map[makeIndex(tmp_index, new_cluster_index)] = max_dist;
	  // The distance between tmp_index (!=index1, != index2) and index1 and index2 are not needed any more.
	  dist_map[makeIndex(tmp_index, index1)] = -1; // Mark as unused
	  dist_map[makeIndex(tmp_index, index2)] = -1; // Mark as unused

	  if (verbosity > 10)
	  {
	    std::cerr << "New distance: " << tmp_index << " , " << new_cluster_index << " " << max_dist << std::endl;  
	    std::cerr << "Mark unused:  " << tmp_index << " , " << index1 << std::endl;  
	    std::cerr << "Mark unused:  " << tmp_index << " , " << index2 << std::endl;  
	  }
	}
	++it_cluster_indices;
      } // END  while (it_cluster_indices != it_cluster_indices_end)

      // Remove index1 and index2 from cluster indices. Add new cluster index.

      if (verbosity > 10)
      {
	std::cerr << "Insert new cluster index: " << new_cluster_index << std::endl; 
	std::cerr << "Erase index: " << index1 << std::endl; 
	std::cerr << "Erase index: " << index2 << std::endl; 
      }

      // Now we add the new cluster index to the sets.
      all_cluster_indices.insert(new_cluster_index);
      current_cluster_indices.insert(new_cluster_index);

      // We remove the indices we clustered from the current_cluster_indices
      current_cluster_indices.erase(index1);
      current_cluster_indices.erase(index2);

      // Adjust sorted_distances:
      sort();

      it     = sorted_distances.begin();
      it_end = sorted_distances.end();
     
      // TODO: can be done a bit faster:
      // Delete items marked for deletion.
      // Delete them from sorted_distances as well as from the distance_map.
      // Remember: Distances are sorted from small to large, so the distances
      //           marked for deletion which have distance -1 come first.
      //           If we passed those with distance -1 we are done. (see loop condition)

      while (it != it_end && (*it)->second == -1)
      {
	unsigned u = (*it)->first;

	if (verbosity > 10)
	  std::cerr << "Erase permanently: " << first_index(u) << " " << second_index(u) << std::endl; 

	dist_map.erase(*it);            // it points to the iterator in the dist_map. That is the object we want to remove. 
	sorted_distances.erase(it);

	// Next iteration step:
	// The next distance with distance -1, if present, should be at the beginning again.
	it = sorted_distances.begin();
	// All iterators that are not removed should keep their validity.
	//	it_end = sorted_distances.end();
      }

      // This completes this clustering step:
      if(verbosity > 10)
      {
	faststring cluster_string;
	std::cerr << "Cluster string after this step: " << new_cluster_index << std::endl;
	get_cluster_string2(cluster_string);
	std::cerr << cluster_string << std::endl;
      }

      // See above:
      // Iteration statement: 	it     = sorted_distances.begin();

    } // END while (sorted_distances.size() > 1 && (*it)->second < distance_limit)
    
  } // run_clustering(double distance_limit)


  void append_to_cluster_string(short node, faststring &str, const std::vector<faststring> &names)
  {
    short child1 = cluster_hierarchy[node].child1;
    short child2 = cluster_hierarchy[node].child2;
    double within_distance = cluster_hierarchy[node].dist;
    bool isterminal = cluster_hierarchy[node].is_terminal();
    bool use_names = false;

    if (names.size() > 0)
      use_names = true;

    if (cluster_hierarchy[node].this_nodes_index != node)
    {
      std::cerr << "Error in cluster hierarchy" << std::endl;
      exit(-1);
    }

    if (isterminal)
    {
      str.append("(");  // Bracket terminal nodes. This makes parsing easier.
      if (use_names)
      {
	str.append( names[node] );
      }
      else
      {
	str.append(faststring(node));
      }
      str.append(")");
    }
    else
    {
      str.append("(");
      append_to_cluster_string(child1, str, names);
      str.append(",");
      append_to_cluster_string(child2, str, names);
      str.append("):");
      str.append(faststring(within_distance));
    }
  }


  void get_cluster_string(faststring &str, const std::vector<faststring> &names)
  {
    std::set<short>::iterator it_cluster_indices     = current_cluster_indices.begin();
    std::set<short>::iterator it_cluster_indices_end = current_cluster_indices.end();
    short    tmp_index;
    short    count = 0;

    unsigned remaining_cluster_indices = current_cluster_indices.size();

    if (remaining_cluster_indices > 1)
      str.append("{");
    else
      str.append("(");
    
    while (it_cluster_indices != it_cluster_indices_end)
    {
      if (count > 0)
	str.append(",");
	tmp_index = *it_cluster_indices;
	append_to_cluster_string(tmp_index, str, names);
	++it_cluster_indices;
	++count;
    } // END  while (it_cluster_indices != it_cluster_indices_end)

    if (remaining_cluster_indices > 1)
      str.append("}:");
    else
      str.append("):");

    str.append(faststring(max_distance()));
  }


  void append_to_cluster_string2(short node, faststring &str)
  {
    short child1 = cluster_hierarchy[node].child1;
    short child2 = cluster_hierarchy[node].child2;
    double within_distance = cluster_hierarchy[node].dist;
    bool isterminal = cluster_hierarchy[node].is_terminal();

    if (cluster_hierarchy[node].this_nodes_index != node)
    {
      std::cerr << "Error in cluster hierarchy" << std::endl;
      exit(-1);
    }

    if (isterminal)
    {
      str.append("(");  // Bracket terminal nodes. This makes parsing easier.
      str.append(faststring(node));
      str.append(")");
    }
    else
    {
      str.append("(");
      append_to_cluster_string2(child1, str);
      str.append(",");
      append_to_cluster_string2(child2, str);
      str.append("):");
      str.append(faststring(within_distance));
    }
  }


  // As above, but uses indices as names.
  void get_cluster_string2(faststring &str)
  {
    std::set<short>::iterator it_cluster_indices     = current_cluster_indices.begin();
    std::set<short>::iterator it_cluster_indices_end = current_cluster_indices.end();
    short    tmp_index;
    short    count = 0;

    unsigned remaining_cluster_indices = current_cluster_indices.size();

    if (remaining_cluster_indices > 1)
      str.append("{");
    else
      str.append("(");
    
    // Each (remaining) index in current_cluster_indices
    // specifies its own cluster which we now print
    // in hierarchical form.
    // The current_cluster_indices are the starting point.
    // Its sub-clusters will be obtained from the cluster_hierarchy
    // vector.

    while (it_cluster_indices != it_cluster_indices_end)
    {
      if (count > 0)
	str.append(",");
	tmp_index = *it_cluster_indices;
	append_to_cluster_string2(tmp_index, str);
	++it_cluster_indices;
	++count;
    } // END  while (it_cluster_indices != it_cluster_indices_end)

    if (remaining_cluster_indices > 1)
      str.append("}:");
    else
      str.append("):");

    str.append(faststring(max_distance()));
  }



  void print_debug(std::ostream &os)
  {
    os << "***Debug-stats************************************"                << std::endl;
    os << "sorted_distances.size()        " << sorted_distances.size()        << std::endl;
    os << "current_cluster_indices.size() " << current_cluster_indices.size() << std::endl;
    os << "all_cluster_indices.size()     " << all_cluster_indices.size()     << std::endl;
    os << "cluster_hierarchy.size()       " << cluster_hierarchy.size()       << std::endl;
    os << "dist_map.size()                " << dist_map.size()                << std::endl;


    os << "Sorted distances:" << std::endl;
    print_sorted(os);

    os << "Remaining cluster indices:" << std::endl;
    print_current_cluster_indices(os);
    os << "**************************************************"                << std::endl;    

  }

  // How to write an adapter to get the clustering result:
  // The vector current_cluster_indices contains all "umbrella" indices.
  // There is one cluster for each entry in this vector.
  // Some entries might not be clustered. Remember: We stopped clustering
  // if the remaining distances are above the threshold.
  // Starting from the indices in current_cluster_indices, we obtain
  // the subclusters from following the cluster_hierarchy.

  // Gets one cluster by index. Returns false if the index is out of range.
  // Returns true if in range.

  void add_to_set(short node, std::set<short> &cluster, double &max_distance)
  {
    short child1 = cluster_hierarchy[node].child1;
    short child2 = cluster_hierarchy[node].child2;
    double within_distance = cluster_hierarchy[node].dist;
    bool isterminal = cluster_hierarchy[node].is_terminal();

    if (cluster_hierarchy[node].this_nodes_index != node)
    {
      std::cerr << "Error in cluster hierarchy" << std::endl;
      exit(-1);
    }

    if (isterminal)
    {
      cluster.insert(node);
    }
    else
    {
      add_to_set(child1, cluster, max_distance);
      add_to_set(child2, cluster, max_distance);
      if (within_distance > max_distance)
	max_distance = within_distance;

      //      str.append(faststring(within_distance));
    }
  }

  // The index must be 0 based.
  bool get_cluster_by_index(short i, std::set<short> &cluster, double &max_dist)
  {
    std::set<short>::iterator it_cluster_indices     = current_cluster_indices.begin();
    std::set<short>::iterator it_cluster_indices_end = current_cluster_indices.end();
   
    advance(it_cluster_indices, i);

    if (it_cluster_indices == it_cluster_indices_end)
      return false;

    cluster.clear();
    max_dist = 0;
    add_to_set(*it_cluster_indices, cluster, max_dist);

    return true;
  }




};

#endif
