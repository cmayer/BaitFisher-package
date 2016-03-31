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
#ifndef GLOBAL_TYPES_AND_PARAMETERS_H
#define GLOBAL_TYPES_AND_PARAMETERS_H

#include <cstdio>
#include <iostream>
#include <vector>
#include <string>
#include <cstring>
#include "faststring2.h"
#include <climits>
#include <fstream>

#define DEBUG

#define PROGNAME "Bait-Filter"
#define VERSION  "1.0.5"

extern faststring                       global_bait_filename;     //
extern char                             global_mode;              //
extern faststring                       global_blast_exe;         // 
extern faststring                       global_output_filename;

extern double                           global_blast_min_hit_coverage_of_bait;
extern double                           global_blast_max_first_hit_evalue;
extern double                           global_blast_max_second_hit_evalue;
extern faststring                       global_blast_db;  //

extern faststring                       global_blast_extra_commandline;
extern faststring                       global_blast_evalue_commandline;

extern unsigned                         global_thinning_step_width;

extern bool                             global_use_GUI;
extern unsigned                         global_conversion_mode;

extern faststring                       global_ProbeID_prefix;

extern unsigned                         global_verbosity;

#define macromax(x,y) ((x)<(y) ? (y) : (x))
#define macromin(x,y) ((x)<(y) ? (x) : (y))

#ifdef  DEBUG
#define DEBUGOUT1(x)        std::cerr << x                << std::endl;
#define DEBUGOUT2(x,y)      std::cerr << x << y           << std::endl;
#define DEBUGOUT3(x,y,z)    std::cerr << x << y << z      << std::endl;
#define DEBUGOUT4(x,y,z,w)  std::cerr << x << y << z << w << std::endl;
#else
#define DEBUGOUT1(x)
#define DEBUGOUT2(x,y)
#define DEBUGOUT3(x,y,z)
#define DEBUGOUT4(x,y,z,w)
#endif


void good_bye_and_exit(FILE *of, int);
void init_param();
void read_and_init_parameters(int argc, char** argv);
void print_parameters(FILE*, const char *s);

#endif
