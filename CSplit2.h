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
#ifndef CSplitH
#define CSplitH

#include "fast-dynamic-bitset/fast-dynamic-bitset.h"
#include <iostream>
#include <vector>
#include <fstream>
//#include "number2cstring.h"

class CSplit
{
 private:
  fast_dynamic_bitset split;

 public:
  // Bitset with a specfied number of bits, here taxa
  CSplit(unsigned long size):split(size) {
    split.clear();
  }

 CSplit(const fast_dynamic_bitset& b):split(b){}

 CSplit(const CSplit& b):split(b.split){}
    //  CSplit(unsigned long, const vector<unsigned long>&);
    //  CSplit(unsigned long, unsigned long*, unsigned long);

  friend bool operator ==(const CSplit&, const CSplit&);
  friend bool operator <(const CSplit&, const CSplit&);
  friend bool lessThan_mod_flip_friend(const CSplit&, const CSplit&);

  fast_dynamic_bitset& get_split()
  {
    return split;
  }

  void reset()
  {
    split.reset();
  }

  void set(unsigned i)
  {
    split.set(i);
  }

  void reset(unsigned i)
  {
    split.reset(i);
  }

  bool test(unsigned i) const
  {
    return split.test(i);
  }

  void set(std::vector<int> vec)
  {
    int i=0, n=vec.size();

    for (; i<n; ++i)
    {
      split.set(vec[i]);
    }
  }

  void flip()
  {
    split.flip();
  }

  void flip(unsigned i)
  {
    split.flip(i);
  }


  unsigned count_taxa_in_ingroup()
  {
    return split.count();
  }

  unsigned size()
  {
    return split.size();
  }

  bool compatible_with(const CSplit& s_b) const
  {
    CSplit s_a_flipped(split);
    CSplit s_b_flipped(s_b.split);
    s_a_flipped.split.flip();
    s_b_flipped.split.flip();

    if( ((split & s_b.split).none())             || ((split & s_b_flipped.split).none())            || 
	((s_a_flipped.split & s_b.split).none()) || ((s_a_flipped.split & s_b_flipped.split).none()) )
      return true;
    else
      return false;
  }

  // Returns true if all 1s in (*this) are also 1s in s_b
  bool is_this_a_subset_of_the_parameter(const CSplit& s_b) const
  {
    CSplit s_b_flipped(s_b.split);
     s_b_flipped.split.flip();

     // Some debug code:
     //     std::cout << "*this & parameter_flipped: " << (split & s_b_flipped.split) << std::endl;

     return  (split & s_b_flipped.split).none();
  }

  // Returns true if all 1s in s_b are also 1s in (*this)
  bool is_parameter_a_subset_of_this(const CSplit& s_b) const
  {
    CSplit s_a_flipped(split);
     s_a_flipped.split.flip();

     // Some debug code:
     //     std::cout << "*this flipped & parameter: " << (s_a_flipped.split & s_b.split) << std::endl;

     return  (s_a_flipped.split & s_b.split).none();
  }

  void print(std::ostream& os) const
  {
    os << split;
  }

  const std::string as_string_lowBitsFirst() const
  {
    return split.as_string_lowBitsFirst();
  }

  const std::string as_string_highBitsFirst() const
  {
    return split.as_string_highBitsFirst();
  }

};


/* // Neue Reihenfolge der Argumente */
/* CSplit::CSplit(unsigned long len, const vector<unsigned long>& vec):split(len) { */
/*   vector<unsigned long>::const_iterator it; */
/*   it = vec.begin(); */
/*   while( it != vec.end() ) { */
/*     split.set(*it);       //indexerror */
/*     ++it; */
/*   } */
/* } */

/* // Neue Reihenfolge der Argumente */
/* CSplit::CSplit(unsigned long len, unsigned long *u, unsigned long array_len):split(len) { */
/*   split.clear(); */
/*   int i = 0; */
/*   while(i < array_len) { */
/*     split.set(u[i]);      //indexerror */
/*     i++; */
/*     *u++; */
/*   } */
/* } */





inline bool operator ==(const CSplit& split_a, const CSplit& split_b) {
  if ( split_a.split == split_b.split ) {
    return true;
  } else {
    CSplit split_b_flipped(split_b);
    split_b_flipped.split.flip();
    if( split_a.split == split_b_flipped.split ) {
      return true;
    } else {
      return false;
    }
  }
}

inline bool operator <(const CSplit& split_a, const CSplit& split_b)
{
  if ( split_a.split < split_b.split ) {
    return true;
  }
  else {
    return false;
  }
}

/* inline bool lessThan_mod_flip_friend(const CSplit& split_a, const CSplit& split_b) */
/* { */
/*   CSplit s_a_flipped(split_a); */
/*   CSplit s_b_flipped(split_b); */

/*   s_a_flipped.split.flip(); */
/*   s_b_flipped.split.flip(); */

/*   //  return (split_a < split_b) || (s_a_flipped < s_b_flipped); */
/*   return (s_a_flipped < s_b_flipped); */
/* } */

/* struct lessThan_mod_flip */
/* { */
/*   bool operator()(const CSplit& split_a, const CSplit& split_b) const */
/*   { */
/*     return lessThan_mod_flip_friend(split_a, split_b); */
/*   } */
/* }; */

//
// Split with branch length:
//
class CSplit_l : public CSplit
{
  double b_length;

 public:

 CSplit_l(unsigned long u):CSplit(u), b_length(-1)
 {}

 CSplit_l(const fast_dynamic_bitset& b, double bl = -1):CSplit(b), b_length(bl) {}

 CSplit_l(const CSplit& b,  double bl = -1):CSplit(b), b_length(bl){}

 CSplit_l(const CSplit_l& b):CSplit(b), b_length(b.b_length){}

  

  void set_b_length(double b_l)
  {
    b_length = b_l;
  }

  double get_b_length() const
  {
    return b_length;
  }


};

#endif
