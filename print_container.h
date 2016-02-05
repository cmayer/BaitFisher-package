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
#ifndef PRINT_CONTAINER_H
#define PRINT_CONTAINER_H

#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <list>
#include <map>

template <typename T, typename S>
std::ostream& operator<<(std::ostream &os, std::pair<T,S> p)
{
  os << p.first << "=>" << p.second;
  return os;
}


template <typename T>
void print_container(std::ostream &os, T it, T it_end,
		     const char *pre, const char *delim, const char *post)
{
  if (it == it_end)
    return;

  os << pre;
  os << *it;
  ++it;
  while (it != it_end)
  {
    os << delim << *it;
    ++it;
  }
  os << post;
}

template<typename T>
void print_container(std::ostream &os, T c,
		     const char *pre, const char *delim, const char *post)
{
  print_container(os, c.begin(), c.end(), pre, delim, post);
}

#endif
