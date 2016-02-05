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
#include <iostream>
#include <fstream>
#include <unistd.h>    // Might be available for mingw. I have not verified this. (??)
#include <sys/stat.h>  // Should be available for mingw. I have not verified this. (??)
#include <sys/types.h>
#include "faststring2.h"
#include <stdio.h>
#include <dirent.h>


#include "mydir-unix.h"

using namespace std;

bool file_exists(const char *fname)
{
  std::ifstream is;
  is.open(fname);
  if (is.fail())
  {
    return false;
  }
  is.close();
  return true;
}



bool file_or_dir_exists(const char *fname)
{
  struct stat st;
  if(stat(fname,&st) == 0)
  {
    return true;
  }
  else
    return false;
}

// Determines the file type:
// 0: does not exist
// 1: regular file
// 2: directory
// 5: other

int file_or_dir_type(const char *fname)
{
  struct stat st;
  if(stat(fname,&st) != 0)
  {
    return 0;
  }

  // Link detection does not seem to work properly on my mac.
//   if (st.st_mode & S_IFLNK)
//   {
//     if (st.st_mode & S_IFREG)
//       return 3;
//     if (st.st_mode & S_IFDIR)
//       return 4;
//   }

  if (st.st_mode & S_IFREG)
    return 1;
  if (st.st_mode & S_IFDIR)
    return 2;

  return 5;
}

// bool directory_exists(const char *dname)
// {
//   struct stat st;
//   if(stat(dname,&st) == 0)
//   {
//     return true;
//   }
//   else
//     return false;
// }


bool change_dir(const char *dname)
{
  return (chdir(dname)==0);
}

// Remember: This does not work like the mv on linux. The new_name must be complete. It cannot end in an existing directory.
bool move_file_or_dir(const char* old_name, const char *new_name)
{
  int res = rename(old_name, new_name);

  return (res == 0);
}


bool make_directory(const char *dname)
{
  return (mkdir(dname, 0700) == 0); // Default: Only user is allowed to use this folder.
}

void next_file_name_in_series(faststring &str, unsigned *z)
{
  faststring fname_without_extension;
  faststring putative_version;

  str.backwards_search_divide_at_break_at('_','/',fname_without_extension, putative_version);

  unsigned   version = 0;
  //  cerr << "putative_version: " << putative_version << endl;

  if (!putative_version.empty() && putative_version.isAnUnsigned() )
  {
    version = putative_version.ToUnsigned(); 
    ++version;
  }
  str = fname_without_extension + "_" + faststring(version);
  //  cerr << "next str in series: " << str << endl;
  if (z != NULL)
  {
    *z = version;
  }
}


//   faststring path = str.get_path();
//   faststring name = str.get_FileOrDirName();
//   faststring major_name = name;
//   unsigned   version = 0;
//   //  cerr << path << "/" << name << endl;
//   //  cerr << major_name << endl;

//   unsigned pos_underscore = major_name.rfind('_');
//   if (pos_underscore != faststring::npos)
//   {
//     faststring putative_version(major_name.begin()+pos_underscore+1, major_name.end());
//     //    cerr << "putative_version: " << putative_version << endl;
//     if (putative_version.isAnUnsigned())
//     {
//       version = putative_version.ToUnsigned();
//       major_name.shorten(pos_underscore);
//       //      cerr << major_name << endl;
//     }
//     ++version;
//   }
//   str = path + "/" + major_name + "_" + faststring(version);
//   //  cerr << "next str in series: " << str << endl;
//   if (z != NULL)
//   {
//     *z = version;
//   }

faststring next_unique_file_or_dir_name(faststring str, unsigned *z)
{
  faststring fname_without_extension;
  faststring putative_version;

  str.backwards_search_divide_at_break_at('_','/',fname_without_extension, putative_version);

 
  unsigned   version = 0;
  //  cerr << "putative_version: " << putative_version << endl;

  if (!putative_version.empty() && putative_version.isAnUnsigned() )
  {
    version = putative_version.ToUnsigned(); 
    ++version;
  }
  str = fname_without_extension + "_" + faststring(version);
  //  cerr << "next str in series: " << str << endl;
  if (z != NULL)
  {
    *z = version;
  }

  while (file_or_dir_exists(str.c_str()))
  {
    ++version;
    str = fname_without_extension + "_" + faststring(version);
  }
  return str;
}

// faststring next_unique_file_or_dir_name(faststring str, unsigned *z)
// {
//   faststring path = str.get_path();
//   faststring name = str.get_FileOrDirName();
//   faststring major_name = name;
//   unsigned   version = 0;
//   //  cerr << path << "/" << name << endl;
//   //  cerr << major_name << endl;


//   unsigned pos_underscore = major_name.rfind('_');
//   if (pos_underscore != faststring::npos)
//   {
//     faststring putative_version(major_name.begin()+pos_underscore+1, major_name.end());
//     //   cerr << "putative_version: " << putative_version << endl;

//     if (putative_version.isAnUnsigned())
//     {
//       version = putative_version.ToUnsigned();
//       major_name.shorten(pos_underscore);
//       //      cerr << major_name << endl;
//     }
//   }
//   ++version;
//   str = path + "/" + major_name + "_" + faststring(version);
//   //  cerr << "str: " << str << endl;


//   while (file_or_dir_exists(str.c_str()))
//   {
//     ++version;
//     str = path + "/" + major_name + "_" + faststring(version);
//   }

//   if (z != NULL)
//   {
//     *z = version;
//   }
//   return str;  
// }

// This version does not seem to work:
// faststring get_working_directory()
// {
//   char buf[400];

//   getcwd(buf, 400);
//   return faststring(buf);
// }

// faststring get_working_directory2(char* argv0)
// {
// }


bool      get_file_names_in_dir(std::list<faststring>& container, faststring &path_name)
{
  faststring filepath;
  DIR *dp;
  struct dirent *dirp;
  struct stat filestat;

  path_name.removeSpacesBack(" /");

  dp = opendir( path_name.c_str() );
  if (dp == NULL)
  {
    return false;
  }

  while ((dirp = readdir( dp )))
  {
    filepath = path_name + "/" + dirp->d_name;

    // If the file is a directory (or is in some way invalid) we'll skip it 
    if (stat( filepath.c_str(), &filestat )) continue;
    if (S_ISDIR( filestat.st_mode ))         continue;

    container.push_back(dirp->d_name);
  }
  closedir( dp );

  return true;
}

bool     get_file_names_in_dir(std::list<faststring>& container, faststring &path_name, faststring &extension_filter)
{
  faststring filepath, extension;
  DIR *dp;
  struct dirent *dirp;
  struct stat filestat;

  path_name.removeSpacesBack(" /");

  dp = opendir( path_name.c_str() );
  if (dp == NULL)
  {
    return false;
  }

  while ((dirp = readdir( dp )))
  {
    filepath = path_name + "/" + dirp->d_name;

    // If the file is a directory (or is in some way invalid) we'll skip it 
    if (stat( filepath.c_str(), &filestat )) continue;
    if (S_ISDIR( filestat.st_mode ))         continue;

    extension = filepath.filename_extension();

    if (extension_filter == extension)
	container.push_back(dirp->d_name);
  }
  closedir( dp );

  return true;  
}

bool     get_file_names_in_dir(std::list<faststring>& container, faststring &path_name, std::list<faststring> &extension_filter)
{
  faststring filepath, extension;
  DIR *dp;
  struct dirent *dirp;
  struct stat filestat;

  std::list<faststring>::iterator it, it1 = extension_filter.begin();
  std::list<faststring>::iterator it2     = extension_filter.end();

  path_name.removeSpacesBack(" /");

  dp = opendir( path_name.c_str() );
  if (dp == NULL)
  {
    return false;
  }

  while ((dirp = readdir( dp )))
  {
    filepath = path_name + "/" + dirp->d_name;

    // If the file is a directory (or is in some way invalid) we'll skip it 
    if (stat( filepath.c_str(), &filestat )) continue;
    if (S_ISDIR( filestat.st_mode ))         continue;

    extension = filepath.filename_extension();

    it = it1;
    while (it != it2)
    {
      if (*it == extension)
	container.push_back(dirp->d_name);
      ++it;
    }
  }
  closedir( dp );

  return true;  
}

long long get_file_size(faststring path_name)
{
  //  cout << path_name.c_str() << endl;
  struct stat filestatus;
  if (stat( path_name.c_str(), &filestatus )==0)
  {
    return filestatus.st_size;
  }

  // This value is returned if no file size could be determined.
  // E.g. if the file does not exist.
  return -1;
}





