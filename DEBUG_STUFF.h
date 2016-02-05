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
#ifndef DEBUG_STUFF_H
#define DEBUG_STUFF_H

#define DEBUG
//#define VERBOSITY 11

// Variant 1:
extern unsigned global_VERBOSITY;

//Variant 2:
// Also:
// comment out: unsigned global_VERBOSITY = 0; in bait-fisher.cpp
//              remove & in front of verbosity in parse_parameter_file(..) in bait-fisher.cpp 



#define  VERBOSITY_ERROR_EXIT   0
#define  VERBOSITY_WARNINGS     1
#define  VERBOSITY_MOREWARNINGS 2
#define  VERBOSITY_PROGRESS     3
#define  VERBOSITY_MOREPROGRESS 4
#define  VERBOSITY_COORDINATES  8
#define  VERBOSITY_FILEINFO     10
#define  VERBOSITY_DEBUG        20

             
//#define global_VERBOSITY 0
#endif
