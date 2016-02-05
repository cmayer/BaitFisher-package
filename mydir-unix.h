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
#ifndef MYDIR_H
#define MYDIR_H

#include "faststring2.h"


bool file_exists(const char *);
bool file_or_dir_exists(const char *);
bool move_file_or_dir(const char*, const char *);

bool make_directory(const char *);
bool change_dir(const char *);

faststring next_unique_file_or_dir_name(faststring, unsigned *z=NULL);
void       next_file_name_in_series(faststring &, unsigned *z=NULL);

faststring get_working_directory();
faststring get_working_directory2(char *);

bool      get_file_names_in_dir(std::list<faststring>&, faststring &path_name);
bool      get_file_names_in_dir(std::list<faststring>&, faststring &path_name, faststring &extension_filter);
bool      get_file_names_in_dir(std::list<faststring>&, faststring &path_name, std::list<faststring> &extension_filter);

long long get_file_size(faststring path_name);
int       file_or_dir_type(const char *fname);

//faststring pwd();

#endif


