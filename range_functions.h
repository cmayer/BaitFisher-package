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
#ifndef RANGE_FUNCTIONS_H
#define RANGE_FUNCTIONS_H

#include <iostream>
#include <list>
#include <cstdlib>
#include "faststring2.h"
#include <cerrno>
#include <iterator>
#include <vector>

#define MACRO_MAX(x,y) ((x)<(y) ? (y) : (x))
#define MACRO_MIN(x,y) ((x)<(y) ? (x) : (y))

// This is a range class for unsigned numbers not for double numbers, which would
// be interesting in its own right.
// In principle one could make a template class out of this.


//
// General conventions:
//  - ranges are undirected
//  - start <= end
//  - end is position after last position in range
//  - ranges can be 0 or 1 based. Be careful!!!!
//  - the preferred way to represent a vanishing range is (0,0).
// TODO: See comments in functions for which this convention does not apply.
// TODO: Several funcions work independent of a convention here.
//




class  CRangeList;

class Crange
{

  friend class CRangeList;

  unsigned b;
  unsigned e;

 public:

 Crange(unsigned pb, unsigned pe):b(pb), e(pe) 
  {
    order();
  }

  void order() // According to the conventions above: b <= e
  {
    if (b > e)
    {
      unsigned tmp = b;
      b   = e;
      e   = tmp;
    }
  }


  friend unsigned overlap(Crange r1, Crange r2)
  {
    unsigned overlap_start = MACRO_MAX(r1.b, r2.b);
    unsigned overlap_end   = MACRO_MIN(r1.e, r2.e);
    
    if (overlap_end > overlap_start)
      return overlap_end - overlap_start;
    else
      return 0;
  }
  
  friend Crange overlap_range(Crange r1, Crange r2)
  {
    unsigned overlap_start = MACRO_MAX(r1.b, r2.b);
    unsigned overlap_end   = MACRO_MIN(r1.e, r2.e);
    
    if (overlap_end > overlap_start)
    {
      Crange res(overlap_start, overlap_end);
      return res;
    }
    else
      return Crange(0,0);
  }
  
  friend void overlap_range(Crange r1, Crange r2, Crange &res1, Crange &res_over, Crange &res2)
  {
    unsigned overlap_start = MACRO_MAX(r1.b, r2.b);
    unsigned overlap_end   = MACRO_MIN(r1.e, r2.e);
    unsigned res1_start    = MACRO_MIN(r1.b, r2.b);
    unsigned res2_end      = MACRO_MAX(r1.e, r2.e);

    if (overlap_end > overlap_start)
    {
      res_over.assign(overlap_start, overlap_end);
      res1.assign(res1_start, overlap_start);
      res2.assign(overlap_end, res2_end);
    }
    else
    {
      res_over.assign(0, 0);
      if (r1 < r2)
      {
	res1 = r1;
	res2 = r2;
      }
      else
      {
	res1 = r2;
	res2 = r1;
      }
    }
  }

  friend double overlap_proportion(Crange r1, Crange r2)
  {
    unsigned overlap_start = MACRO_MAX(r1.b, r2.b);
    unsigned overlap_end   = MACRO_MIN(r1.e, r2.e);

    if (overlap_end > overlap_start)
    {
      double overl = overlap_end - overlap_start;
      double l     = MACRO_MIN(r1.len(), r2.len());

      return (overl/l);
    }
    else
      return 0.0;
  }


  friend bool adjacent_or_overlap(Crange r1, Crange r2)
  {
    unsigned overlap_start = MACRO_MAX(r1.b, r2.b);
    unsigned overlap_end   = MACRO_MIN(r1.e, r2.e);
    
    if (overlap_end >= overlap_start)
      return true;    // == means they are adjacent 
    else
      return false;
  }

  Crange()
  {
    b = e = 0;
  }

  friend void swap(Crange &r1, Crange &r2)
  {
    Crange tmp = r2;
    r2 = r1;
    r1 = tmp;
  }

  friend void swap(Crange *r1, Crange *r2)
  {
    Crange tmp = *r2;
    *r2 = *r1;
    *r1 = tmp;
  }

  bool empty() const
  {
    return (b == e);
  }

  void limit_to_range(Crange limit_r)
  {
    // Nothing to be done in this case:
    if (is_r2_Subrange_of_r1(limit_r, *this) )
      return;

    b = MACRO_MAX(limit_r.b, b);
    e = MACRO_MIN(limit_r.e, e);
    if (b >= e)
      b = e = 0;
  }

  friend bool operator==(Crange r1, Crange r2)
  {
    return (r1.b == r2.b && r1.e == r2.e);
  }

  void assign(unsigned pb, unsigned pe)
  {
    b = pb;
    e = pe;
    order();
  }

  // true if r1 starts before r2
  // true if r1 starts with r2 and if r1 is shorter
  // r1    *******
  // r2    ***********   r1 < r2 is true
  friend bool operator<(Crange r1, Crange r2)
  {
    if (r1.b != r2.b)
      return r1.b < r2.b;
    return r1.e < r2.e;
  }

  friend bool equalRange(Crange r1, Crange r2)
  {
    return (r1 == r2);
  }


  friend bool is_r2_Subrange_of_r1(Crange r1, Crange r2)
  {
    return (r2.b >= r1.b && r2.e <= r1.e);
  }

  friend Crange range_intersection(Crange r1, Crange r2)
  {
    Crange res( MACRO_MAX(r1.b,r2.b), MACRO_MIN(r1.e,r2.e));
    if (res.b >= res.e)
      res.e = res.b = 0;

    return res;
  }

  friend bool before(Crange r1, Crange r2)
  {
    return r1.e <= r2.b;
  }

  friend bool after(Crange r1, Crange r2)
  {
    return r1.b >= r2.e;
  }

  // formerly called range_union, which was a misnomer.
  //  friend  ......

  friend Crange range_span(Crange r1, Crange r2)
  {
    Crange res( MACRO_MIN(r1.b,r2.b), MACRO_MAX(r1.e,r2.e));
    return res;
  }


  void shift (unsigned s)
  {
    b +=s;
    e +=s;
  }

  unsigned begin() const
  {
    return b;
  }

  unsigned end() const
  {
    return e;
  }

  unsigned inclusive_end() const
  {
    return e-1;
  }

  unsigned len() const
  {
    return e-b;
  }

  bool point_is_in_range(unsigned x) const
  {
    if (x>=b && x < e)
      return true;
    else
      return false;
  }

  bool point_is_behind_range(unsigned x) const
  {
    if (x >= e)
      return true;
    else
      return false;
  }

  bool point_is_before_range(unsigned x) const
  {
    if (x < b)
      return true;
    else
      return false;
  }

  void print_with_brackets(std::ostream &os) const
  {
    os << "(" << b << "," << e << ")"; 
  }

  void print_with_brackets_inclusive_end(std::ostream &os) const
  {
    os << "(" << b << "," << (e-1) << ")"; 
  }

  void print_minus_deliminated_inclusive_end(std::ostream &os) const
  {
    os << b << "-" << (e-1); 
  }

  void transform_start_end(int mult_start, int add_start, int mult_end, int add_end)
  {
    b = b*mult_start + add_start;
    e = e*mult_end   + add_end;
  }

  void transform_coord(std::vector<unsigned> &coord_trans)
  {
    // We subtract 1 from e before we look up this index in the field, since the end index of a range is 1 behind its last element.
    if ( e-1 >= coord_trans.size() )
    {
      std::cerr << "Index out of range: " << e << std::endl;
      exit(-33);
    }

    b = coord_trans[b];
    e = coord_trans[e-1]+1;   // We transform the coord of the last element and add 1.
  }


}; // class Crange




/* class CRangeList_array */
/* { */
/*   Crange      limit; */
/*   faststring  rl_array; */

/*   CRangeList_array(Crang l):limit(l), faststring('\0', l.end()-l.begin()+1) */
/*   {} */

/*   Crange   get_limit() */
/*   { */
/*     return limit; */
/*   } */

/*   unsigned get_limit_begin() */
/*   { */
/*     return limit.begin(); */
/*   } */

/*   unsigned get_limit_end() */
/*   { */
/*     return limit.end(); */
/*   } */
  


/* }; */




// NOTE: Code referring to array is either experimental or unfinished.
//
//

class CRangeList
{
  std::list<Crange> rl;
  bool              sorted;
  unsigned          lower_bound;
  unsigned          upper_bound; // Accoring to the convetion this position comes after the last position in this list of ranges.
  unsigned          count; // Introduced since the empty() function could be
                           // via size()==0, whereas size could be O(n)
                           // which might be very inefficient.
  // rl_array is not fully implemented: In the future it could be used to store "in" and "out of" the range list
  // by setting chars to 1 or 0, which would make many functions more efficient.
  faststring        *rl_array;


 public:
  typedef std::list<Crange>::iterator CRangeList_iterator;
  typedef std::list<Crange>::const_iterator CRangeList_iterator_const;

  // Default constructor:
  // Caution: lower_bound and upper_bound are set to 0, since they are needed by
  //          some functions. This is better than arbitrary values.
  //          even though one could say the values are not correct.
  //          But arbitrary values are not more correct than 0. 
  //          Take care that all member functions set these variables
  //          if the first element of a list is added.

 CRangeList():sorted(false), lower_bound(0), upper_bound(0), count(0), rl_array(NULL)
  {}

  // Copy constructor
 CRangeList(const CRangeList &a):rl(a.rl), lower_bound(a.lower_bound), upper_bound(a.upper_bound), count(a.count), rl_array(NULL)
  {
    if (a.rl_array != NULL)
    {
      rl_array = new faststring(*a.rl_array);
    }
  }

  void  clear_rl_array()
  {
    if (rl_array)
    {
      delete rl_array;
      rl_array = NULL;
    }
  }

  ~CRangeList()
  {
    clear_rl_array();
  }


 public:
  void add(unsigned pb, unsigned pe)
  {
    Crange new_r (pb, pe);

    if (new_r.empty())
      return;

    clear_rl_array();

    // If the new element is > than the last element in the list,
    // sorting is not affected, otherwise "it is false".
    // Note: We cannot set sorted to true, if the rest of the list is not sorted.

    if (count == 0)
    {
      sorted = true;
      lower_bound = new_r.b;
      upper_bound = new_r.e;
      rl.push_back(new_r);
    }
    else
    {
      if ( !(rl.back() < new_r) )
	sorted = false;
      // else - do not change sorted status

      if (lower_bound > new_r.b)
	lower_bound = new_r.b;
      if (upper_bound < new_r.e)
	upper_bound = new_r.e;

      rl.push_back(new_r);
    }
    ++count;
  }

  CRangeList_iterator_const begin() const
  {
    return rl.begin();
  } 

  CRangeList_iterator_const begin()
  {
    return rl.begin();
  } 

  CRangeList_iterator_const end() const 
  {
    return rl.end();
  } 

  CRangeList_iterator_const end()
  {
    return rl.end();
  } 

  void create_reset_rl_array()
  {
    if (rl_array)
    {
      rl_array->assign(upper_bound-lower_bound, '\0');
    }
    else
    {
      rl_array = new faststring('\0', upper_bound-lower_bound);      
    }
  }


  void fill_rl_array()
  {
    if (rl_array == NULL)
      create_reset_rl_array();

    CRangeList_iterator_const it     = rl.begin();
    CRangeList_iterator_const it_end = rl.end();

    unsigned pos1, pos2;

    // foreach range:
    while (it != it_end)
    {
      pos1 = it->begin();
      pos2 = it->end();

      while (pos1 < pos2)
      {
	rl_array[pos1] = 1;
	++pos1;
      }
      ++it;
    }
  }

  void add(Crange new_r)
  {
    // If the new element is > than the last element in the list,
    // sorting is not affected, otherwise "it is false".
    // Note: We cannot set sorted to true, if the rest of the list is not sorted.

    clear_rl_array();

    if (new_r.empty())
      return;

    if (count == 0)
    {
      sorted = true;
      lower_bound = new_r.b;
      upper_bound = new_r.e;
      rl.push_back(new_r);
    }
    else
    {
      if ( !(rl.back() < new_r) )
	sorted = false;
      // else - do not change sorted status

      if (lower_bound > new_r.b)
	lower_bound = new_r.b;
      if (upper_bound < new_r.e)
	upper_bound = new_r.e;

      rl.push_back(new_r);
    }
    ++count;
  }


  void add(const CRangeList &a)
  {
    CRangeList_iterator_const it     = a.rl.begin();
    CRangeList_iterator_const it_end = a.rl.end();

    while (it != it_end)
    {
      if ( !it->empty() )
	add(*it);
      ++it;
    }
  }

  bool add(faststring str)
  {
    size_t pos1=0, pos2, len=str.size();
    unsigned c1,c2;

    while (pos1 < len && str[pos1] != ';')
    {
      str.skip_symbol(pos1,' ');

      errno = 0;
      c1 = str.ToUnsigned(pos1, pos2);
      //      std::cerr << "c1 " << c1 << std::endl;
      if (errno != 0)
      {
	return false;
      }
      str.skip_symbol(pos2,' ');

      if (str[pos2] == ',' || str[pos2] == ';') // One number range.
      {
	add(c1,c1+1);
      }
      else if (str[pos2] == '-') // Two number range
      {
	pos1 = pos2+1;  // Be careful '-' is read as minus
	str.skip_symbol(pos1,' ');
	//	std::cerr << "pos1 " << pos1 << std::endl;
	errno = 0;
	c2 = str.ToUnsigned(pos1, pos2);
	//	std::cerr << "c2 " << c2 << std::endl;
	if (errno != 0)
	{
	  return false;
	}
	//	std::cerr << "pos2 " << pos2 << std::endl;
	str.skip_symbol(pos2,' ');
	if (str[pos2] > '\0' && str[pos2] <= '\7')
	{
	  unsigned tmp = c1;
	  unsigned delta = str[pos2]-'\0';
	  if (errno != 0)
	  {
	    return false;
	  }
	  ++pos2;
	  str.skip_symbol(pos2,' ');
	  
	  while (tmp < c2)
	  {
	    add(tmp, tmp+1);
	    tmp += delta;
	  }
	}
	else
	{
	  add(c1, c2+1);
	}
      }
      else
      {
	return false;
      }
      pos1 = pos2;
      
      if (str[pos1] == ',') // we are not done yet - the list continues
      {
	++pos1;
	str.skip_symbol(pos1,' ');
	if (str[pos1] == ';')

	  return false;
      }
    }
    return true;
  }


  unsigned size() const
  {
    return count;
  }

  // Assignment operator: Use standard assigment operator. No pointer content needs to be copied
  //                      This is still true since rl_array is not ready to use yet.
  //                      Currently there is no need/no plan to implement it.
         
  void clear()
  {
    rl.clear(); // Will call the destructor for all objects.
    sorted = false;
    clear_rl_array();
    // I am not sure this is the best way to assign values to the bounds.
    // Bounds are not well defined for an empty list.
    lower_bound = upper_bound = 0;
    count = 0;
  }


  void sort()
  {
    if (!sorted)
    {
      rl.sort();
      sorted = true;
    }
  }


  // We expect the use to check for count == 0
  Crange bounds()
  {
    // There is one problem here: The last entry need to have the "highest coordinate"
    // 100-150, 200-400, 250-300

    Crange res(lower_bound, upper_bound);
    return res;
  }

  // We expect the use to check for count == 0
  unsigned get_lower_bound()
  {
    return lower_bound;
  }

  // We expect the use to check for count == 0
  unsigned get_upper_bound()
  {
    return upper_bound;
  }


  // Find an island of overlaping ranges, starting at it_beg.
  // The initial value of it_end (the second parameter) is overwritten,
  // and is therefore irrelevant.
  // When leaving this function, it_end points to the range behind the island.
  // This is also true, if the island is the last island in the list.
  // In this case, it_end == to list_of_sat_ptr.end()

  unsigned findOverlappingIsland(CRangeList_iterator_const  it_beg,
                                 CRangeList_iterator_const &it_end) const
  {
    unsigned count=1;

    it_end = it_beg;
    ++it_end;

    while (it_end != rl.end())
    {
      if (overlap(*it_beg, *it_end))   // If they overlap, we move it_end,                                                                                
      {
        ++count;
        ++it_end;
      }
      else // if they do not overlap any more, we need not be out of the island.                                                                            
      {    // Let us move it_beg, which should now be behind or equal to it_end.                                                                            
        ++it_beg;
        if (it_beg == it_end) // this is the case of the last two sats do not overlap                                                                       
          break;              // In this case it_end points one behind the last island.                                                                     
      }
    }
    return count;
  }

  unsigned findAdjacentIsland(CRangeList_iterator_const  it_beg,
			      CRangeList_iterator_const &it_end) const
  {
    unsigned count=1; // The current range is the Nr 1.

    it_end = it_beg;
    ++it_end;

    while (it_end != rl.end())
    {
      if (adjacent_or_overlap(*it_beg, *it_end))   // If they overlap, we move it_end,                                                                                
      {
        ++count;
        ++it_end;
      }
      else // if they do not overlap any more, we need not be out of the island.                                                                            
      {    // Let us move it_beg, which should now be behind or equal to it_end.                                                                            
        ++it_beg;
        if (it_beg == it_end) // this is the case of the last two sats do not overlap                                                                       
          break;              // In this case it_end points one behind the last island.                                                                     
      }
    }
    return count;
  }


  void findRangeWithMaxExtToRight(CRangeList_iterator_const  it_beg,
				  CRangeList_iterator_const  it_end,
				  CRangeList_iterator_const &it_extends_most) const
  {
    it_extends_most = it_beg;
    if (it_beg == it_end)
      return;

    ++it_beg;
    while (it_beg != it_end)
    {
      if ((*it_beg).e > (*it_extends_most).e )
        it_extends_most = it_beg;
      ++it_beg;
    }
  }


  Crange range_of_island(CRangeList_iterator_const    it_beg,
			 CRangeList_iterator_const   &it_end_of_island)
  {
    Crange res;
    CRangeList_iterator_const  it_extends_most;

    if (!sorted)
    {
      std::cerr << "Error: range_of_island requires the list of ranges to be sorted." << std::endl;
      exit(0);
    }

    res.b = it_beg->b;

    findOverlappingIsland(it_beg, it_end_of_island);
    findRangeWithMaxExtToRight(it_beg, it_end_of_island, it_extends_most);
    
    res.e = it_extends_most->e;

    return res;
  }

  Crange range_of_adjacentisland(CRangeList_iterator_const   it_beg,
				 CRangeList_iterator_const   &it_end_of_island)
  {
    Crange                     res;
    CRangeList_iterator_const  it_extends_most;

    if (!sorted)
    {
      std::cerr << "Internal error: range_of_adjacentisland requires the list of ranges to be sorted." << std::endl;
      exit(0);
    }

    res.b = it_beg->b;

    findAdjacentIsland(it_beg, it_end_of_island);
    findRangeWithMaxExtToRight(it_beg, it_end_of_island, it_extends_most);
    
    res.e = it_extends_most->e;

    return res;
  }


  // coverage_len allows to specify a range - double check that feature.
  double coverage_rel(Crange bds)
  {
    CRangeList rl_tmp(*this);

    rl_tmp.limit_to_range(bds);

    double bds_len = bds.len();
    
    return rl_tmp.coverage_len(bds)/bds_len;
  }

  unsigned coverage_len2(Crange bds)
  {
    CRangeList rl2;
    get_consolidated_range_list(rl2);

    rl2.limit_to_range(bds);
    return rl2.coverage_raw();
  }

  unsigned coverage_raw()
  {
    unsigned count=0;

    CRangeList_iterator_const it, it_end;

    it     = rl.begin();
    it_end = rl.end();

    while (it != it_end)
    {
      count += it->len();
      ++it;
    }
    return count;
  }

  // Side effect: *this will be sorted
  unsigned coverage_len(Crange bds) // bounds in which to compute the coverage
  {
    unsigned cover = 0;

    CRangeList_iterator_const it, it_end, end_island;
    Crange r1, r2;

    //    std::cerr << "Computing coverage for this: " << std::endl;
    //    print_bracketed_list(std::cerr);
    //    std::cerr << std::endl;

    sort();

    it     = rl.begin();
    it_end = rl.end();

    while (it != it_end && it->e < bds.b)
      ++it;

    if (it == it_end)
      return 0;

    // it points to first range after start of bds
    r1 = range_of_island(it, end_island);
    cover += overlap(r1, bds);
    it = end_island;
    //    std::cerr << "First overlap: " << cover << std::endl;

    while (it != it_end && it->b < bds.e)
    {
      r1 = range_of_island(it, end_island);
      it = end_island;
      if (r1.e > bds.e)
      {
	cover += overlap(r1, bds);
	break;
      }
      else
	cover += r1.len();
      //      std::cerr << "next value: " << cover << std::endl;
    }

    return cover;
  }

  // Side effect: *this will be sorted
  unsigned coverage_len() //  without bounds
  {
    unsigned cover = 0;

    CRangeList_iterator_const it, it_end, end_island;
    Crange r1;

    //    std::cerr << "Computing coverage for this: " << std::endl;
    //    print_bracketed_list(std::cerr);
    //    std::cerr << std::endl;

    sort();

    it     = rl.begin();
    it_end = rl.end();

    while (it != it_end)
    {
      r1 = range_of_island(it, end_island);
      it = end_island;

      cover += r1.len();
      //      std::cerr << "next value: " << cover << std::endl;
    }

    return cover;
  }

  unsigned coverage_len_loop() //  without bounds
  {
    unsigned cover = 0;
    Crange r1, b;

    //    std::cerr << "Computing coverage for this: " << std::endl;
    //    print_bracketed_list(std::cerr);
    //    std::cerr << std::endl;

    if (count == 0)
      return 0;

    // Required by point_is_in_RangeList
    sort();

    CRangeList_iterator_const it1 = rl.begin();

    b = bounds(); //Upper bound of this range is element after range.
    unsigned i, b_end;

    b_end = b.end();
    for (i=b.begin(); i<b_end; ++i)
    {
      if (point_is_in_RangeList(i, it1) )
      {
	++cover;
      }
    }
    return cover;
  }



  void test1()
  {
    CRangeList_iterator_const  it1 = rl.begin();
    CRangeList_iterator_const  it2 = rl.end();
    CRangeList_iterator_const  it3, it4;
 
    while (it1 != it2)
    {
      Crange r = range_of_island(it1, it3);

      it1 = it3;

      std::cout << "Next range: " << std::endl;
      r.print_with_brackets(std::cout);
      std::cout << std::endl;
    }
  }

  void print_bracketed_list(std::ostream &os) const
  {
    CRangeList_iterator_const  it_beg = rl.begin();
    CRangeList_iterator_const  it_end = rl.end();

    if (it_beg == it_end)
      return;

    while (1) //(it_beg != it_end)
    {
      it_beg->print_with_brackets(os);
      ++it_beg;
      if (it_beg == it_end)
	return;
      os << ",";
    }
  }

  void print_bracketed_list_inclusive_ends(std::ostream &os) const
  {
    CRangeList_iterator_const  it_beg = rl.begin();
    CRangeList_iterator_const  it_end = rl.end();

    if (it_beg == it_end)
      return;

    while (1) //(it_beg != it_end)
    {
      it_beg->print_with_brackets_inclusive_end(os);
      ++it_beg;
      if (it_beg == it_end)
	return;
      os << ",";
    }
  }


  void print_deliminated_inclusive_end(std::ostream &os, const char *delim) const
  {
    CRangeList_iterator_const  it_beg = rl.begin();
    CRangeList_iterator_const  it_end = rl.end();

    if (it_beg == it_end)
      return;

    while (1) //(it_beg != it_end)
    {
      it_beg->print_minus_deliminated_inclusive_end(os);
      ++it_beg;
      if (it_beg == it_end)
	return;
      os << delim;
    }
  }

  void print_deliminated_inclusive_end(std::ostream &os, char delim= ',') const
  {
    CRangeList_iterator_const  it_beg = rl.begin();
    CRangeList_iterator_const  it_end = rl.end();

    if (it_beg == it_end)
      return;

    while (1) //(it_beg != it_end)
    {
      it_beg->print_minus_deliminated_inclusive_end(os);
      ++it_beg;
      if (it_beg == it_end)
	return;
      os << delim;
    }
  }


  // For this CRangeList compute a consolidated range list:
  // Adjacent or overlapping ranges are merged in the consolidated list
  void get_consolidated_range_list(CRangeList &crl)
  {
    sort();
    crl.clear();

    CRangeList_iterator_const  it1, it2;

    it1 = it2 = rl.begin();

    while (it2 != rl.end() )
    {
      crl.add(range_of_adjacentisland(it1, it2) );
      //      crl.sort();
      it1 = it2;
    }
  }


  void remove_empty()  // Now we try to avoid to add empty ranges in the first place
  {
    CRangeList_iterator  it1;

    it1 = rl.begin();
    
    while (it1 != rl.end() )
    {	
      if ( it1->empty() )
      {
	it1 = rl.erase(it1);
	--count;
      }
      else
	++it1;
    }
  }


  // Side effect: Removes empty ranges form list:
  void limit_to_range(Crange limit_r)
  {
    CRangeList_iterator  it1;

    clear_rl_array();

    it1 = rl.begin();
    
    if (sorted)
    {
      //      std::cout << "sorted mode" << std::endl;
      while (it1 != rl.end() && before(*it1, limit_r) )
      {	
	// remove ranges before limiting range:
	++it1;
      }
      if (it1 != rl.begin())
	rl.erase(rl.begin(), it1);
      // count will be reset at the end!!

      while (it1 != rl.end() && !after(*it1, limit_r) )
      {
	it1->limit_to_range(limit_r);
	
	if ( it1->empty() )
	{
	  it1 = rl.erase(it1);
	  // count will be reset at the end!!
	}
	else
	  ++it1;
      }

      if (it1 != rl.end() )
      {	
	// remove ranges after limiting range:
	rl.erase(it1, rl.end());
	// count will be reset at the end!!
      }
    }
    else
    {
      while (it1 != rl.end() )
      {
	it1->limit_to_range(limit_r);
	
	if ( it1->empty() )
	{
	  it1 = rl.erase(it1);
	  // count will be reset at the end!!
	}
	else
	  ++it1;
      }
    }
    // This could be made a bit more efficient.
    // This solution however is most secure.
    count = rl.size();
  }


  friend CRangeList range_union(Crange a, Crange b)
  {
    CRangeList tmp, res;

    tmp.add(a);
    tmp.add(b);

    tmp.get_consolidated_range_list(res);
    
    return res;
  }

  friend CRangeList cut_out_range(Crange a, Crange b)
  {
    Crange r1, r2, r3;
    CRangeList res;

    overlap_range(a,b,r1, r2, r3);

    res.add(r1);
    res.add(r3);
    return res;
  }


  void set_to_range_list_union(const CRangeList &a, const CRangeList &b)
  {
    CRangeList tmp(a);
    tmp.add(b);

    clear();
    tmp.get_consolidated_range_list(*this);
  }

  void set_to_range_list_intersection(CRangeList &a, CRangeList &b)
  {
    clear();

    CRangeList ca, cb;

    a.get_consolidated_range_list(ca);
    b.get_consolidated_range_list(cb);

    CRangeList_iterator_const it_a, it_b;

    it_a = ca.rl.begin();
    it_b = cb.rl.begin();

    while (it_a != ca.rl.end() )
    {
      // While the range it_b lies before it_a:
      while ( it_b != cb.rl.end() && before(*it_b, *it_a) )
      {
	++it_b;
      }
      // While the range it_a overlaps with it_b:
      while ( it_b != cb.rl.end() && overlap(*it_a, *it_b) )
      {
	add(range_intersection(*it_a, *it_b) );
	++it_b;
      }
      if (it_b != cb.rl.begin())
	--it_b;
      ++it_a;
    }
  }

  // Side effect: list of ranges will be sorted
  bool point_is_in_RangeList(unsigned x)
  {
    sort();
    CRangeList_iterator_const it1, it2;

    it1 = rl.begin();
    it2 = rl.end();

    // Skip all ranges which have an end value that lies before x
    while (it1 != it2 && x >= it1->e)
    {
      ++it1;
    }

    // Test all ranges which have a beginning that lies before x
    while (it1 != it2 && x >= it1->b)
    {
      if (it1->point_is_in_range(x))
	return true;
      ++it1;
    }
    // All remaining ranges have a beginning that lies behind x:
    return false;
  }

  // Requirement: list of ranges must be sorted
  bool point_is_in_RangeList(unsigned x, CRangeList_iterator_const &it1) const
  {
    if (!sorted)
    {
      std::cerr << "When calling point_is_in_RangeList, the range list must be sorted." << std::endl;
      exit(-11);
    }

    CRangeList_iterator_const it2;
    it2 = rl.end();

    // Skip all ranges which have an end value that lies before x
    while (it1 != it2 && x >= it1->e)
    {
      ++it1;
    }

    // Test all ranges which have a beginning that lies before x
    while (it1 != it2 && x >= it1->b)
    {
      if (it1->point_is_in_range(x))
	return true;
      ++it1;
    }
    if (it1 != rl.begin() )
      --it1; // We might have moved too far.

    // All remaining ranges have a beginning that lies behind x:
    return false;
  }


  void export_as_gff(std::ostream &os, const char *seq, const char* source, const char *type, char direction)
  {
    CRangeList_iterator_const it     = rl.begin();
    CRangeList_iterator_const it_end = rl.end();

    while (it != it_end)
    {
      os << seq    << "\t"
	 << source << "\t"
	 << type   << "\t"
	 << it->b  << "\t"
	 << it->e  << "\t"
	 << ".\t"
	 << direction << "\t"
	 << ".\t"
	 << "ID=" << seq << "_" << it->b << "_" << it->e << std::endl;
      ++it;
    }
  }

  // Side effect: CRangeList rl is being sorted and set to a consolidated list.
  void set_to_complement_of(CRangeList &rl, Crange limit)
  {
    CRangeList crl;
    Crange r1, r2, r3;

    // Clear this CRangeList. Otherwise we add to a probably non empty list.
    clear();

    rl.get_consolidated_range_list(crl);

    crl.sort(); // Should already be sorted. Just to be sure.
    crl.limit_to_range(limit);

    CRangeList_iterator_const it     = crl.rl.begin();
    CRangeList_iterator_const it_end = crl.rl.end();

    r3 = limit;
    while (it != it_end)
    {
      // r2 is just a dummy for the overlapping part.
      overlap_range(r3, *it, r1, r2, r3);
      //      std::cerr << "Adding "; r1.print_with_brackets(std::cerr); std::cerr << std::endl;
      add(r1);
      ++it;
    }
    add(r3);
  }


  unsigned* get_range_list_field()
  {
    unsigned* p = new unsigned [2*count];

    CRangeList_iterator_const it     = rl.begin();
    CRangeList_iterator_const it_end = rl.end();
    unsigned i=0;
   
    sort();

    while (it != it_end)
    {
      p[i]   = it->begin();
      p[i+1] = it->end();
      ++i;
      ++it;
    }
    return p;
  }

  friend bool less_than_CRangeList(const CRangeList &a, const CRangeList &b)
  {
     CRangeList_iterator_const it_a     = a.rl.begin();
     CRangeList_iterator_const it_b     = b.rl.begin();
     CRangeList_iterator_const it_a_end = a.rl.end();
     CRangeList_iterator_const it_b_end = b.rl.end();


     while (it_a != it_a_end && it_b != it_b_end && *it_a == *it_b )
     {
       ++it_a;
       ++it_b;
     }
     // return true if a is shorter than b
     if (it_a != it_a_end || it_b != it_b_end)
       return a.size() < b.size();

     return (*it_a < *it_b);
  }

  friend bool less_than_CRangeList_pointer(const CRangeList *a, const CRangeList *b)
  {
    return less_than_CRangeList(*a, *b);
  }

  void  transform_start_end(int mult_start, int add_start, int mult_end, int add_end)
  {
    CRangeList_iterator it     = rl.begin();
    CRangeList_iterator it_end = rl.end();

    while (it != it_end)
    {
      it->transform_start_end(mult_start, add_start, mult_end, add_end);
      ++it;
    }
  }

  void transform_coord(std::vector<unsigned> &coord_trans)
  {
    CRangeList_iterator it     = rl.begin();
    CRangeList_iterator it_end = rl.end();

    while (it != it_end)
    {
      it->transform_coord(coord_trans);
      ++it;
    }
  }

}; // End class CRangeList



inline void coverage_values_loop(Crange limit,
				 CRangeList &rl1,
				 CRangeList &rl2,
				 CRangeList &rl3,
				 CRangeList &rl4,

				 unsigned &len_gff1,
				 unsigned &len_gff2,
				 unsigned &len_gff3,
				 unsigned &len_gff4,
				 unsigned &in1not3not4,
				 unsigned &in2not3not4,
				 unsigned &in1in2not3not4,
				 unsigned &in3in4not1not2)
{
  unsigned i;
  bool    in1, in2, in3, in4;

  len_gff1=0;
  len_gff2=0;
  len_gff3=0;
  len_gff4=0;
  in1not3not4=0;
  in2not3not4=0;
  in1in2not3not4=0;
  in3in4not1not2=0;

  unsigned en=limit.end();

  CRangeList::CRangeList_iterator_const it1, it2, it3, it4;

  it1 = rl1.begin();
  it2 = rl2.begin();
  it3 = rl3.begin();
  it4 = rl4.begin();
  
  rl1.sort();
  rl2.sort();
  rl3.sort();
  rl4.sort();

  for (i=limit.begin(); i < en; ++i)
  {
    //    std::cerr << "Checking " << i << std::endl;

    in1 = rl1.point_is_in_RangeList(i, it1);
    in2 = rl2.point_is_in_RangeList(i, it2);
    in3 = rl3.point_is_in_RangeList(i, it3);
    in4 = rl4.point_is_in_RangeList(i, it4);

    if (in1) ++len_gff1;
    if (in2) ++len_gff2;
    if (in3) ++len_gff3;
    if (in4) ++len_gff4;

    if (in1 && !in3 && !in4)         ++in1not3not4;
    if (in2 && !in3 && !in4)         ++in2not3not4;

    if ( in1 &&  in2 && !in3 && !in4)  ++in1in2not3not4;
    if (!in1 && !in2 &&  in3 &&  in4)  ++in3in4not1not2;
  }
}



// This function is basically used for testing purposes. Its logic is straight forward.

inline  void coverage_values_loop(Crange limit,
				  CRangeList rl1,  // no side effect
				  CRangeList rl2,  // no side effect
				  unsigned &limit_len,
				  unsigned &len_rl1,
				  unsigned &len_rl2,
				  unsigned &len_union,
				  unsigned &len_intersection)
{
  unsigned i;
  bool    in1, in2;

  len_rl1 = 0;
  len_rl2 = 0;
  len_intersection = 0;
  len_union        = 0;
  limit_len = limit.end() - limit.begin();

  unsigned en=limit.end();

  CRangeList::CRangeList_iterator_const it1, it2;

  it1 = rl1.begin();
  it2 = rl2.begin();

  rl1.sort();
  rl2.sort();

  for (i=limit.begin(); i < en; ++i)
  {
    //    std::cerr << "Checking " << i << std::endl;

    in1 = rl1.point_is_in_RangeList(i, it1);
    in2 = rl2.point_is_in_RangeList(i, it2);

    if (in1)
    {
      ++len_rl1;
      ++len_union;
      if (in2)
      {
	--len_union;
	++len_intersection;
      }
    }
    if (in2)
    {
      ++len_rl2;
      ++len_union;
    }
  }
}

inline  void coverage_values(Crange limit,
			     CRangeList rl1,  // no side effect
			     CRangeList rl2,  // no side effect
			     unsigned &limit_len,
			     unsigned &len_rl1,
			     unsigned &len_rl2,
			     unsigned &len_union,
			     unsigned &len_intersection)
{
  rl1.limit_to_range(limit);
  rl2.limit_to_range(limit);

  //  std::cout << "DEBUG: ";
  //  rl1.print_bracketed_list(std::cout);
  //  std::cout << std::endl;

  CRangeList un, inter;
  
  un.set_to_range_list_union(rl1,rl2);
  inter.set_to_range_list_intersection(rl1,rl2);
  limit_len = limit.len();

  len_rl1   = rl1.coverage_len();
  len_rl2   = rl2.coverage_len();
  len_union = un.coverage_len();
  len_intersection = inter.coverage_len();
}


inline unsigned overlap(unsigned a1, unsigned a2,
		 unsigned b1, unsigned b2)
{
  unsigned overlap_start = MACRO_MAX(a1, b1);
  unsigned overlap_end   = MACRO_MIN(a2, b2);
 
  if (overlap_end > overlap_start)
    return overlap_end - overlap_start;
  else
    return 0;
}


inline bool equalRange(unsigned a1, unsigned a2,
		unsigned b1, unsigned b2)
{
  return (a1 == b1 && a2 == b2);
}


inline bool is_r2_Subrange_of_r1(unsigned a1, unsigned a2, unsigned b1, unsigned b2)
{
  return (b1 >= a1 && b2 <= a2);
}


inline void range_intersection(unsigned a1, unsigned a2,
			unsigned b1, unsigned b2,
			unsigned &c1, unsigned &c2)
{
  c1 = MACRO_MAX(a1,b1);
  c2 = MACRO_MIN(a2,b2);
  if (c1 > c2)  // Empty intersection
    c1=c2=0;
}

// formerly called range_union. Only if the ranges overlap we have the true union.
inline void range_span(unsigned a1, unsigned a2,      // First range
		       unsigned b1, unsigned b2,      // Second range
		       unsigned &c1, unsigned &c2)    // Result range
{
  c1 = MACRO_MIN(a1,b1);
  c2 = MACRO_MAX(a2,b2);
}



#endif
