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
#ifndef STATISTIC_FUNCTIONS_H
#define STATISTIC_FUNCTIONS_H

#include <iostream>
#include <vector>
#include <map>
#include <cmath>
#include <cctype>
#include <iterator>
#include <cstdlib>
#include <algorithm>
#include <climits>

#define EPSS  0.00000000001;


#define macro_min(x,y) ((x)<(y) ? (x) : (y))
#define macro_max(x,y) ((x)>(y) ? (x) : (y))
#define macro_mean_double(x,y) (((x)+(y))/2.0)


template<typename T>
size_t vec_mean_sd(const std::vector<T> &vec,
		   double &mean, double &sd)
{
  size_t     i, n=vec.size();
  double  sum_x  = 0;
  double  sum_xx = 0;

  if (n==0)
    return 0;

  for (i=0; i<n; ++i)
  {
    sum_x  +=  (double)vec[i];
    sum_xx +=  (double)vec[i]*vec[i];
  }
  if (n!=0)
    mean = sum_x/(double)n;
  else
    mean = 0;
  if (n>1)
    sd   = sqrt( (sum_xx-sum_x*mean)/(double)(n-1.0) );
  else
    sd   = 0;
  if (sd < 0.000001)
    sd = 0;
  return n;
}

template<typename T>
size_t vec_mean_sd(const std::vector<T> &vec,
		   double &mean, double &sd,
		   T &sum, T &sum_of_squares)
{
  size_t     i, n=vec.size();
  T  sum_x  = 0;
  T  sum_xx = 0;
  T  elem;

  if (n==0)
    return 0;

  for (i=0; i<n; ++i)
  {
    elem = vec[i];
    sum_x  +=  elem;
    sum_xx +=  elem*elem;
  }
  mean = (double)sum_x/(double)n;
  if (n>1)
    sd   = sqrt( (sum_xx-sum_x*mean)/(double)(n-1.0) );
  else
    sd   = 0;
  if (sd < 0.000001)
    sd = 0;

  sum            = sum_x;
  sum_of_squares = sum_xx;

  return n;
}

template<typename T>
void median_range(std::vector<T> &vec, size_t f, size_t l, double &median, size_t &index, bool &is_datum)
{
  size_t s = f+l;
  size_t m = s/2;

  if (s%2 == 0)
  {
    median   = vec[m];
    index    = m;
    is_datum = true;
  }
  else
  {
    median   = (vec[m]+vec[m+1])/2;
    index    = m;
    is_datum = false;
  }
}

// Method 2 routines:

// Prerequisite: Type T must be a type for which it makes sense to multiply it with a double to obtain a double.
// Compute quartiles:
// First quartile, second quartile=median, and third quartile
template<typename T>
size_t vec_median_quartile_sort_method2(std::vector<T> &vec, double &Q1, double &Q2, double &Q3)
{
  size_t     n=vec.size();
  size_t     index, dummy;
  bool       is_datum;

  std::sort(vec.begin(), vec.end());

  median_range(vec, 0, n-1, Q2, index, is_datum);

  if (is_datum)
  {
    median_range(vec, 0,     index, Q1, dummy, is_datum);
    median_range(vec, index, n-1,   Q3, dummy, is_datum);
  }
  else
  {
    median_range(vec, 0,       index, Q1, dummy, is_datum);
    median_range(vec, index+1, n-1,   Q3, dummy, is_datum);
  }
  return n;
}

template<typename T>
size_t vec_median_quartile_method2(std::vector<T> vec, double &Q1, double &Q2, double &Q3)
{
  return vec_median_quartile_sort_method2(vec, Q1, Q2, Q3);
}

//
// Computes outlier bounds for a vector.
// Side effect: vec gets sorted. More efficient than copying the vector.
template<typename T>
size_t vec_median_quartile_outlier_bounds_sort_method2(std::vector<T> &vec, double &Q1, double &Q2, double &Q3, double &O_lower, double &O_upper)
{
  int    res = vec_median_quartile_sort_method2(vec, Q1, Q2, Q3);
  double IQR = Q3 - Q1;
    
  O_lower = Q1 - 1.5 * IQR;
  O_upper = Q3 + 1.5 * IQR;

  return res;
}

// Computes outlier bounds for a vector.
// Same vec_median_quartile_outlier_bounds_sort but without sorting the vector. Slighly slower, since the vector has to be copied.
template<typename T>
size_t vec_median_quartile_outlier_bounds_method2(std::vector<T> vec, double &Q1, double &Q2, double &Q3, double &O_lower, double &O_upper)
{
  int    res = vec_median_quartile_sort_method2(vec, Q1, Q2, Q3);
  double IQR = Q3 - Q1;
    
  O_lower = Q1 - 1.5 * IQR;
  O_upper = Q3 + 1.5 * IQR;

  return res;
}

template<typename T>
void vec_mark_outlier_mehod2(std::vector<T> &vec, std::vector<bool> &outlier)
{
  outlier.clear();

  double Q1, Q2, Q3, O_lower, O_upper;
  vec_median_quartile_outlier_bounds_method2(vec, Q1, Q2, Q3, O_lower, O_upper);

  unsigned i, N = vec.size();
  outlier.reserve(N);
  for (i=0; i<N; ++i)
  {
    if (vec[i] < O_lower || vec[i] > O_upper)
      outlier.push_back(true);
    else
      outlier.push_back(false);
  }
  std::cout << "Qi: "  << Q1 << " " << Q2 << " " << Q3 << std::endl;
  std::cout << "Bounds: "  << O_lower << " " << O_upper << std::endl;
}

// Method 3 routines:

// Prerequisite: Type T must be a type for which it makes sense to multiply it with a double to obtain a double.
// Compute quartiles:
// First quartile, second quartile=median, and third quartile
template<typename T>
size_t vec_median_quartile_sort_method3(std::vector<T> &vec, double &Q1, double &Q2, double &Q3)
{
  size_t     n=vec.size();
  size_t     m = n/4;
  size_t     r = n%4;

  std::sort(vec.begin(), vec.end());

  if (n==0)
    return 0;

  if (n==1)
  {
    Q1 = vec[0];
    Q2 = vec[0];
    Q3 = vec[0];
    return 1;
  }

  if (r == 0)
  {
    Q1 = 0.5*vec[  m-1] + 0.5*vec[  m];
    Q2 = 0.5*vec[2*m-1] + 0.5*vec[2*m];
    Q3 = 0.5*vec[3*m-1] + 0.5*vec[3*m];
  }
  else if (r == 1)
  {
    //    std::cerr << "m: " << m << std::endl;

    Q1 = 0.25*vec[ m-1] + 0.75*vec[  m];
    Q2 = vec[2*m];
    Q3 = 0.75*vec[ 3*m] + 0.25*vec[3*m+1];
  }
  else if (r == 2)
  {
    //    std::cerr << "m: " << m << std::endl;

    Q1 = vec[m];
    Q2 = 0.50*vec[2*m] + 0.50*vec[2*m+1];
    Q3 = vec[3*m+1];
  }
  else // if (r == 3)
  {
    //    std::cerr << "m: " << m << std::endl;

    Q1 = 0.75*vec[  m] + 0.25*vec[  m+1];
    Q2 = vec[2*m+1];
    Q3 = 0.25*vec[3*m+1] + 0.75*vec[3*m+2];
  }
  return n;
}

template<typename T>
size_t vec_median_quartile_method3(std::vector<T> vec, double &Q1, double &Q2, double &Q3)
{
  return vec_median_quartile_sort_method3(vec, Q1, Q2, Q3);
}

// Computes outlier bounds for a vector.
// Side effect: vec gets sorted. More efficient than copying the vector.
template<typename T>
size_t vec_median_quartile_outlier_bounds_sort_method3(std::vector<T> &vec, double &Q1, double &Q2, double &Q3, double &O_lower, double &O_upper)
{
  int    res = vec_median_quartile_sort_method3(vec, Q1, Q2, Q3);
  double IQR = Q3 - Q1;
    
  O_lower = Q1 - 1.5 * IQR;
  O_upper = Q3 + 1.5 * IQR;

  return res;
}

// Computes outlier bounds for a vector.
// Same vec_median_quartile_outlier_bounds_sort but without sorting the vector. Slighly slower, since the vector has to be copied.
template<typename T>
size_t vec_median_quartile_outlier_bounds_method3(std::vector<T> vec, double &Q1, double &Q2, double &Q3, double &O_lower, double &O_upper)
{
  int    res = vec_median_quartile_sort_method3(vec, Q1, Q2, Q3);
  double IQR = Q3 - Q1;
    
  O_lower = Q1 - 1.5 * IQR;
  O_upper = Q3 + 1.5 * IQR;

  return res;
}

template<typename T>
void vec_mark_outlier_mehod3(std::vector<T> &vec, std::vector<bool> &outlier)
{
  outlier.clear();

  double Q1, Q2, Q3, O_lower, O_upper;
  vec_median_quartile_outlier_bounds_method3(vec, Q1, Q2, Q3, O_lower, O_upper);

  unsigned i, N = vec.size();
  outlier.reserve(N);
  for (i=0; i<N; ++i)
  {
    if (vec[i] < O_lower || vec[i] > O_upper)
      outlier.push_back(true);
    else
      outlier.push_back(false);
  }
  std::cout << "Qi: "  << Q1 << " " << Q2 << " " << Q3 << std::endl;
  std::cout << "Bounds: "  << O_lower << " " << O_upper << std::endl;
}

template<typename T>
size_t vec_min_max(const std::vector<T> &vec, T &min, T &max)
{
  size_t  i, n=vec.size();

  if (n==0)
    return 0;

  max = min = vec[0];

  for (i=0; i<n; ++i)
  {
    if (vec[i] < min)
      min = vec[i];

    if (vec[i] > max)
      max = vec[i];
  }
  return n;
}


// Determines the mean and sd for a range.
template <typename T>
size_t range_mean_sd(T it_beg,
		     T it_end,
		     double &mean, double &sd)
{
  double  sum_x  = 0;
  double  sum_xx = 0;
  size_t n=0;
  
  if (it_beg == it_end)
  {
    return 0;
  }

  while (it_end != it_beg)
  {
    sum_x  +=  (double)*it_beg;
    sum_xx +=  (double)*it_beg* *it_beg;
    ++it_beg;
    ++n;
  }
  mean = sum_x/(double)n;
  if (n>1)
    sd   = sqrt( (sum_xx-sum_x*mean)/(double)(n-1.0) );
  else // if (n==1)
    sd   = 0;
  if (sd < 0.000001)
    sd = 0;
  return n;
}


template<typename T>
size_t range_min_max(T it_beg,
		     T it_end,
		     double &min, double &max)
{
  size_t  n=0;

  if (it_beg == it_end)
    return 0;

  max = min = *it_beg;

  while (it_end != it_beg)
  {
    if (*it_beg< min)
      min = *it_beg;

    if (*it_beg > max)
      max = *it_beg;
    ++n;
    ++it_beg;
  }
  return n;
}

template <typename T>
void vec_sum_sum_squares(std::vector<T> v, T& sum, T&sum_squares)
{
  sum=0;
  sum_squares=0;

  unsigned i, n=v.size();

  for (i=0; i<n; ++i)
  {
    sum         += v[i];
    sum_squares += v[i]*v[i];
  }
}

// More functions: Single result as return value, using pointers to specify ranges. 
template <typename T>
T minimum(T*b, T*e)
{
  T m;
  if (b>=e)
    return INT_MIN;
  m = *b;
  ++b;
  while (b!=e)
  {
    if (*b < m)
    {
      m = *b;
    }
    ++b;
  }
  return m;
}

template <typename T>
T maximum(T*b, T*e)
{
  T m;
  if (b>=e)
    return INT_MIN;
  m = *b;
  ++b;
  while (b!=e)
  {
    if (*b > m)
    {
      m = *b;
    }
    ++b;
  }
  return m;
}





#endif
