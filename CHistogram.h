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
#ifndef CHISTOGRAM_H
#define CHISTOGRAM_H

// This is the CHistogram class written by Christoph Mayer, Forschungsmuseum Alexander Koenig, Bonn, Germany.
// Discription:
// This class allowes different binning schemes:
// 1) User supplied number of bins together with lower and upper value.
// IMPORTANT: The lower value is inclusive, the upper value is exclusive.
//            Therefore, a value equal to the upper value cannot be added
//            to the histrogram data.
//            If this is a problem, increase the upper value, or take
//            a different constructor. See below.
// 2) Vector or range of data. From this the lower and upper bound
//    are determined automatically. Number fo bins or automatic binning
//    if the number of bins parameter has certain negative values.
// 3) Discrete data. In this case the data is not binned when it is added,
//    but stored in a map. The binning is performed upon request.




#include <iostream>
#include <vector>
#include <map>
#include <cmath>
#include <cctype>
#include <iterator>
#include <cstdlib>
#include <algorithm>
#include "statistic_functions.h"

#define EPSS  0.00000000001;








inline bool greater_than_pair(const std::pair<double, unsigned> &a, const std::pair<double, unsigned> &b )
{
  return a.second > b.second;
}




inline void add_or_count(std::map<double, unsigned> &m, double &x)
{
  std::map<double,unsigned>::iterator it;

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




class CHistogram

{
 private:
  double   a;   // Inclusive range minimum
  double   b;   // Exclusive range maximum
  unsigned bins;
  double   step;
  unsigned entries;

  std::vector<unsigned> data;

  bool     discrete_data;  // is determined in initialisation. Switched on with -10.
  std::map<double, unsigned> discrete_hist_data;

 private:
  bool add_intern_non_discrete(double x)
  {
    // Be careful: If x==b this is allowed, but it places
    // this hit into the 

    if (x < a || x >= b)
    {
      if (x == b)
	std::cerr << "Warning: CHistrogram::add_intern_non_discrete data.siz\n"
		  << "   Data points that are added must be exclusive the upper bound of the given range.\n\n";
      return false;
    }

    ++entries;

    unsigned bin = (x-a)/step;
    //    std::cerr << "HH " << bin << " " << data.size() << std::endl;

    ++data[bin];
    return true;
  }

  // Adds the value to the discete data map.
  // Number of bins as well as a and b are not changed.
  bool add_intern_discrete_no_update(double x)
  {
    ++entries;
    add_or_count(discrete_hist_data, x);
    return true;
  }

  void update_a_b_bins_discrete_data()
  {
    // Remember: The discrete data is stored in an automatically sorted map.
    // Therefore, the bounds are determined from the first and last element.

    // Determine a, b and bins:
    a = discrete_hist_data.begin()->first;
    b = discrete_hist_data.rbegin()->first;

    bins = discrete_hist_data.size();
  }


 public:
 void init_with_vector(const std::vector<double> &vals, int bins_param=0)
 {
   init_with_range(vals.begin(), vals.end(), vals.size(), bins_param);
 }


 template <typename T>
   void init_with_range(T         it_beg,         // pointer or iterator to first element
			T         it_end,         // pointer or iterator to element behind last element
			unsigned  num_range_elements, // number of elements in range
			int       bins_param=0)       // number of bins or if <=0 specifies the method
                                                      // to determine this number
 {
   double mi, ma, me, sd;

   discrete_data = false;

   if (num_range_elements == 0) // Fill in some dummy values so that everything is well defined. Pratically this does not make sense, but should prevent other classes using this data from crashing.
   {
     a = 0;
     b = 1;
     data.push_back(0);
     bins = 1;
     step = 1;
     
     return;
   }
   
   mi = ma = 0;  // This code silences the used uninitialized warning of some compilers. It is not necessary, since the two variables are initialized in
                 // range_min_max as the name suggests.

   range_min_max(it_beg, it_end, mi, ma);

   // Special case: Only one category and automatic number of bins, but not discrete data
   // In principle this is a special case of discrete data but not explicitly mentioned.
   if (bins_param != -10 && bins_param <= 0 && mi== ma)
     bins_param = 1;

   

   // Formulas for auto bins number: http://en.wikipedia.org/wiki/Histogram

   if (bins_param == 0) // Sturges formula -- can perform poorly if n < 30.
   {
     bins = ceil(log(num_range_elements)/log(2)+1);
   }
   else if (bins_param == -1) // Scotts formula
   {
     range_mean_sd(it_beg, it_end, me, sd);
     bins = 3.5*sd/pow(num_range_elements, 1./3.);
     if (bins < 1)
       bins = 1;
   }
   else if (bins_param == -2) // sqrt formula
   {
     bins = sqrt(num_range_elements);
   }
   else if (bins_param == -3) // Freedman-Diaconis' formula // Not implemented
   {
     bins = 1;
   }
   else if (bins_param == -4) // minimization of risk function // Not implemented
   {
     bins = 1;
   }
   else if (bins_param == -10) // Each value gets its own bin
   {
     discrete_data = true;
   }
   else if (bins_param > 0)
   {
     bins = bins_param;
   }
   else
   {
     std::cerr << "Error: the init_param parameter is negative, but its value does not specify a method. This needs to be corrected." << std::endl;
     exit(-133);
   }

   // We have minimum (mi), maximum (ma) and number of bins.
   // Now we determine the step size as well as lower and upper bound (a,b).

   if (!discrete_data) // Add non discrete values:
   {
     if (bins > 1)
     {
       step = (ma-mi)/(bins-1.0);
       a = mi - 0.5*step;
       b = ma + 0.5*step;
     }
     else
     {
       step = 1.0; // fantasy value
       a = mi - 0.5*step;
       b = ma + 0.5*step;
     }

/*      std::cerr << "--a    " << a << std::endl; */
/*      std::cerr << "--b    " << b << std::endl; */
/*      std::cerr << "--bins " << bins << std::endl; */
/*      std::cerr << "--step " << step << std::endl; */


     data = std::vector<unsigned>(bins, 0);

     while (it_beg != it_end)
     {
       if (!add_intern_non_discrete(*it_beg) )
       {
	 std::cerr << "Internal error: Value out of range when building the histgram.\n" << std::endl;
	 exit(-144);
       }
       ++it_beg;
     }
   }
   else // Add dicrete values
   {
     step = 0;
     // Add all values to map:
     while (it_beg != it_end)
     {
       add_intern_discrete_no_update(*it_beg);
       ++it_beg;
     }
     update_a_b_bins_discrete_data();
   }


 }

// Constructor with: lower and upper value as well as the number of
// equidistant bins.
// Here we do not allow to determine the number of bins automatically for the simple reason:
// The number of necessary bins can only be determined if the number of data points is known
// which is not assumed here.
CHistogram(double a_param, double b_param, int bins_param):
 a(a_param), b(b_param), bins(bins_param), entries(0), data(bins_param, 0)
 {
    if (bins_param < 1)
    {
      std::cerr << "Critical error in constructor of CHistogram. Parameter bins_param must be positive here.\n" << std::endl;
      exit(-143);
    }

    discrete_data = false;
    step = (b-a)/bins;
 
    // Adds no data.
 }
  
 // Constructor with: range or vector as well as the bins_param
 // if bins_param < 1: then the parameter specifies the method to determine the number
 //                    of bins automatically.
 //                    0: Sturges formula -- can perform poorly if n < 30.
 //                   -1: Scotts formula
 //                   -2: sqrt formula
 //                   -3: Freedman-Diaconis' formula // Not implemented
 //                   -4: minimization of risk function // Not implemented
 // if bins_param > 1: this specifies the actual number of bins

 CHistogram(const std::vector<double> &vals, int bins_param=0):entries(0)
 {
   init_with_range(vals.begin(), vals.end(), vals.size(), bins_param);
 }

 CHistogram(const std::vector<float> &vals, int bins_param=0):entries(0)
 {
   init_with_range(vals.begin(), vals.end(), vals.size(), bins_param);
 }

 CHistogram(const std::vector<unsigned> &vals, int bins_param=0):entries(0)
 {
   init_with_range(vals.begin(), vals.end(), vals.size(), bins_param);
 }

 template <typename T>
 CHistogram(T it_beg,
	    T it_end,
	    int bins_param=0):entries(0)
 {
   init_with_range(it_beg, it_end, std::distance(it_beg, it_end), bins_param);
 }

  double get_lower()
  {
    return a;
  }

  double get_upper()
  {
    return b;
  }

  unsigned get_bins()
  {
    return bins;
  }

  double get_binsize()
  {
    return step;
  }

  // This function should be reimplemented in a cleaner way!!!
  // The use of the data field to store disc
  const std::vector<unsigned> &get_histogram_data()
  {
    if (discrete_data) // data needs to be prepared:
    {
      data.clear();

      std::map<double, unsigned>::iterator it_beg, it_end;

      it_beg = discrete_hist_data.begin(); 
      it_end = discrete_hist_data.end();

      while (it_beg != it_end)
      {
	data.push_back(it_beg->second);
	++it_beg;
      }
    }
    return data; 
  }

  void get_bin_coords(std::vector<double> &v)
  {
/*     std::cerr << "a    " << a << std::endl; */
/*     std::cerr << "b    " << b << std::endl; */
/*     std::cerr << "bins " << bins << std::endl; */
/*     std::cerr << "step " << step << std::endl; */


    unsigned i;

    v.clear();

    if (!discrete_data)
    { 
      for (i=0; i<bins; ++i)
      {
	v.push_back(a+(i+0.5)*step);
      }
    }
    else
    {
      std::map<double, unsigned>::iterator it_beg, it_end;

      it_beg = discrete_hist_data.begin(); 
      it_end = discrete_hist_data.end();

      while (it_beg != it_end)
      {
	v.push_back(it_beg->first);
	++it_beg;
      }
    }
  }


  //XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  // old code - needs revision. See method to compute coordinates
/*   void print_bin_coord_hist_data(std::ostream &os, bool normal=true) */
/*   { */

/*     unsigned i; */
/*     double   tmp; */

/*     // Not very elegant to change this for os - should be changed */

/*     os.setf(std::ios::fixed); */
/*     os.precision(6); */

/*     if (!discrete_data) */
/*     { */
/*       for (i=0; i<bins; ++i) */
/*       { */
/* 	if (normal) */
/* 	  tmp = (double) data[i]/entries/step; */
/* 	else */
/* 	  tmp = (double) data[i]; */
/* 	os << (a+(i+0.5)*step) << "\t"; */
/* 	os << tmp << std::endl; */
/*       } */
/*     } */
/*     else */
/*     { */
/*       std::map<double, unsigned>::iterator it_beg, it_end; */

/*       it_beg = discrete_hist_data.begin();  */
/*       it_end = discrete_hist_data.end(); */

/*       while (it_beg != it_end) */
/*       { */
/* 	if (normal) */
/* 	  tmp = (double) it_beg->second/entries/step; */
/* 	else */
/* 	  tmp = (double) it_beg->second; */
/* 	os << it_beg->first  << "\t" << it_beg->second << std::endl; */

/* 	++it_beg; */
/*       } */
/*     } */
/*   } */


  // Returns false if out of range error, else true.
  bool add(double x)
  {
    if (discrete_data)
    {
      add_intern_discrete_no_update(x);
      update_a_b_bins_discrete_data();
      return true; // Can never fail
    }
    else
    {
      return add_intern_non_discrete(x);
    }
  }

  bool add(std::vector<double> v)
  {
    int i, n=v.size();

    if (!discrete_data) // Non-discrete data:
    {
      for (i=0; i<n; ++i)
	if (!add_intern_non_discrete(v[i]))
	  return false;
    }
    else // discrete data:
    {
      for (i=0; i<n; ++i)
      {
	add_intern_discrete_no_update(v[i]);
      }
      update_a_b_bins_discrete_data();
    }
    return true;
  }

  template<typename T> // T is some iterator or pointer type
  bool add(T b, T e) // defines range
  {
    if (!discrete_data) // Discrete data:
    { 
      while (b != e)
      {
	if (!add_intern_non_discrete(*b))
	  return false;
	++b;
      }
    }
    else
    {
      while (b != e)
      {
	add_intern_discrete_no_update(*b);
	++b;
      }
      update_a_b_bins_discrete_data();
    }
    return true;
  }

  bool is_discrete_data()
  {
    return discrete_data;
  }

  unsigned num_entries()
  {
    return entries;
  }

  //XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
  // old code - needs revision.
/*   bool get_dominant_category_and_num(double &cat, unsigned &num) */
/*   { */
/*     if (discrete_data) */
/*     { */
/*       std::map<double, unsigned>::iterator it_beg, it_end; */

/*       it_beg = discrete_hist_data.begin();  */
/*       it_end = discrete_hist_data.end(); */
      
/*       if (it_beg == it_end) */
/*       { */
/* 	return false;  */
/*       } */

/*       cat = it_beg->first; */
/*       num = it_beg->second; */

/*       while (it_beg != it_end) */
/*       { */
/* 	if (it_beg->second > num) */
/* 	{ */
/* 	  cat = it_beg->first; */
/* 	  num = it_beg->second; */
/* 	} */
/* 	++it_beg; */
/*       } */
/*     } */
/*     else // Non discrete data: */
/*     { */
/*       double   max_key; */
/*       unsigned max_value; */

/*       unsigned i, n; */

/*       n = data.size(); */
/*       i = 0; */
      
/*       if (n == 0) */
/*       { */
/* 	return false;  */
/*       } */

/*       max_key   = a + step/2; */
/*       max_value = data[0]; */

/*       while (i != n) */
/*       { */
/* 	if (data[i] > num) */
/* 	{ */
/* 	  max_key   = a + (i + 0.5)*step; */
/* 	  max_value = data[i]; */
/* 	} */
/* 	++i; */
/*       } */
/*     } */
/*     return true; */
/*   } */

/*   // Should be tested before it is used. */
/*   void dominant_category(double &value, double &prop) */
/*   { */
/*     double   cat; */
/*     unsigned num; */

/*     if (get_dominant_category_and_num(cat, num)) */
/*     { */
/*       value = cat; */
/*       prop  = (double)num/entries;  */
/*     } */
/*     else */
/*     { */
/*       value = 0; */
/*       prop  = -1; */
/*     } */
/*   } */



  void categories_sorted_by_dominance(std::vector<std::pair<double, unsigned> > &cat_val_pair_vec)
  {
    if (discrete_data)
    {
      std::map<double, unsigned>::iterator it_beg, it_end;

      it_beg = discrete_hist_data.begin(); 
      it_end = discrete_hist_data.end();
      
      if (it_beg == it_end)
      {
	return; 
      }

      while (it_beg != it_end)
      {
	double   cat;
	unsigned num;

	cat = it_beg->first;
	num = it_beg->second;

	cat_val_pair_vec.push_back(std::make_pair(cat, num));

	++it_beg;
      }
    }
    else // Non discrete data:
    {
      double   key;
      unsigned value;

      unsigned i, n;

      n = data.size();
      i = 0;
      
      if (n == 0)
      {
	return; 
      }

      while (i != n)
      {
	key = a + step/2;
	value = data[0];

	cat_val_pair_vec.push_back(std::make_pair(key, value));
	++i;
      }
    }

    sort(cat_val_pair_vec.begin(), cat_val_pair_vec.end(), greater_than_pair);
  }

};

#endif
