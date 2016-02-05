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
#ifndef CTRIPLE_H
#define CTRIPLE_H

#include <iostream>

//class Ctriple;
//bool Ctriple::operator<(Ctriple &a, Ctriple &b)

template <typename T1, typename T2, typename T3>
class Ctriple
{
 private:
  T1 t1;
  T2 t2;
  T3 t3;

 public:
 Ctriple(const T1 &pt1,const T2 &pt2,const T3 &pt3):t1(pt1),t2(pt2),t3(pt3)
  {}

  T1 first()  { return t1; }
  T2 second() { return t2; }
  T3 third()  { return t3; }

  friend bool operator<(const Ctriple<T1, T2, T3> &a, const Ctriple<T1, T2, T3> &b)
  {
    if (a.t1 != b.t1)
      return a.t1 < b.t1;
    if (a.t2 != b.t2)
      return a.t2 < b.t2;
    return a.t3 < b.t3;
  }

  void print(std::ostream &os) const
  {
    os << "(" << t1 << "," << t2 << "," << t3 << ")";
  }

};


template <typename T1, typename T2, typename T3>
inline std::ostream &operator<<(std::ostream &os, const Ctriple<T1, T2, T3> &t)
{
  t.print(os);
  return os;
}









#endif


