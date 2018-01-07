/*  BaitFisher (version 1.2.8) a program for designing DNA target enrichment baits
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
#ifndef CFILE_H
#define CFILE_H

#include <fstream>
#include <string>
#include "../faststring2.h"

// Should not exceed 1 000 000 without making sure that the stack is large enough.
#ifndef BUFFERSIZE
#define BUFFERSIZE 1000000
#endif
#define UNDOs      1


// Changes:
// 22.08.2009: Added getline function for faststring. Disadvantage: additional dependence on faststring.h
// 07.01.2011: Changed Version to CFile2_1.h
// 07.01.2011: Now uses the faststring2.h

// Idea: ios::   set Buffer size

// TODO: Old mac format not supported.
//       This requires to allow two successive calls to the internal ungetchar command.


class CFile : private std::ifstream
{
 private:
  // There is a problem here: The stack size is limited. One should consider to move this to the heap.
  // allocate it with malloc, new or declare it as static - the latter is not good in a class since a
  // static class member is the same for all instances of the class, which we do not want here.

  char      buffer[BUFFERSIZE];
  char*     buffer_end;
  char*     buffer_pos;
  
  char      __status;
  unsigned  __line;

  //  bool      __eof;
  //  bool      __openFailed;

  unsigned fill_buffer(unsigned overlap)
  {
    unsigned          good_overlap = buffer_end - buffer;
    std::streamsize   n;

    if (good_overlap < overlap)
    {
      overlap = good_overlap;
    }

    if (overlap > 0)
      std::memmove(buffer, buffer_end - overlap, overlap);

    std::ifstream::read(buffer + overlap, BUFFERSIZE - overlap);
    n = std::ifstream::gcount();

    if ( n == 0 )
    {
      __status |=  __eof_flag;     // Set eof flag
      __status &= ~__good_flag;    // Unset good flag
      __status |=  __fail_flag;    // Setting the fail flag is not always correct. Needs to be unset in routines that read more than one char. 
    }

    buffer_pos = buffer+overlap;
    buffer_end = buffer_pos + n;

    return overlap;
  }

  char getchar_intern()  // should only be done if we did not fail yet!!!! This would allow us to recover!!!
  {
    if ( buffer_pos == buffer_end )
    {
      fill_buffer(UNDOs);
      if (__status & __fail_flag)
	return '\0';
    }
    return *buffer_pos++;
  }

  void ungetchar_intern()  // should only be done if we did not fail yet!!!! This would allow us to recover!!!
  {
    if (buffer_pos != buffer)
      --buffer_pos;
    else
    {
      __status |=  __fail_flag;      // Set fail flag.
      __status &= ~__good_flag;     // Unset good flag.
    }
  }

 public:

  enum {__eof_flag = 1, __good_flag = 2,  __fail_flag = 4, __bad_flag = 8, __fail_reason1 = 16, __fail_reason2 = 32, __fail_reason3 = 64,
	__fail_reason4 = 128};

  void open(const char *name)
  {
    ffopen(name);
  }

  void open(std::string name)
  {
    ffopen(name);
  }

  void ffopen(const char *name)
  {
    __status = __good_flag;

    std::ifstream::open(name);

    if ( std::ifstream::fail() )
    {
      __status |=  __fail_flag;      // Set fail flag.
      __status &= ~__good_flag;     // Unset good flag.
    }

    __line        = 1;
    buffer_end    = buffer;
    buffer_pos    = buffer;
  }

  void ffopen(std::string name)
  {
    ffopen(name.c_str());
  }

  void ffclose()
  {
    std::ifstream::close();
  }

  void close()
  {
    ffclose();
  }

  bool exists()  // -- deprecated -- do not use this any more - check fail() instead !!!!!!!!
  {
    return !fail(); 
  }

  bool fail()
  {
    return (__status & __fail_flag);
  }

  bool good()
  {
    return (__status & __good_flag);
  }

  unsigned line()
  {
    return __line;
  }

  bool eof()
  {
    return (__status & __eof_flag);
  }

  char status()
  {
    return __status;
  }

  bool fail_reason1()
  {
    return (__status & __fail_reason1);
  }

  bool fail_reason2()
  {
    return (__status & __fail_reason2);
  }

  bool fail_reason3()
  {
    return (__status & __fail_reason3);
  }

  bool fail_reason4()
  {
    return (__status & __fail_reason4);
  }

  void rewind()
  {
    clear();
    std::ifstream::seekg (0, std::ios::beg);
    //    clear();
  }

  void clear(char s = __good_flag)
  {
    __status = s;
    std::ifstream::clear();
  }

  void ungetchar()
  {
    if (*(buffer_pos-1) == '\n' && !(__status & __fail_flag))
      --__line;
    ungetchar_intern();
  }

  void ignore(int delim = -1)
  {
    while (__status == __good_flag && getchar() != delim ){}
  }

  char peekchar()
  {
    char c = getchar();
    ungetchar();
    return c;
  }

  char peek()
  {
    return peekchar();
  }

  char getchar()
  {
    register char c;

    c = getchar_intern();

    if ( c < 14 && !(__status & __fail_flag) )
    {
      if (  c == '\r' )
      {
	// Overwrite the last reading position that contained the \r with \n
	*(buffer_pos-1) = '\n';           // Should always be valid! Works for 1 Undo
	c = getchar_intern();
	if ( c != '\n' && !(__status & __fail_flag) )
	{
	  // 	std::cerr << "Old mac file format currently not supported." << std::endl;
	  ungetchar();    /* old mac format, else dos format     */
	}
	c = '\n';
	++__line;
      }
      else if ( c == '\n')
	++__line;
    }
    return c;  
  }


  char getrawchar()
  {
    return getchar_intern();
  }

  void getline(faststring& str, char delim='\n')
  {
    char c = getchar();
    str.clear();
    while ( c != delim && !(__status & __fail_flag) )
    {
      str.push_back(c);
      c = getchar();
    }
    if ((__status & __fail_flag) && str.size() > 0)
      __status &= ~__fail_flag; // Unset fail flag by using & on the complement of the fail flag;
  }

  void getline(std::string& str, char delim='\n')
  {
    char c = getchar();
    str="";
    while ( c != delim && !(__status & __fail_flag) )
    {
      str.push_back(c);
      c = getchar();
    }
    if ((__status & __fail_flag) && str.size() > 0)
      __status &= ~__fail_flag; // Unset fail flag by using & on the complement of the fail flag;
  }

  void getline(char* cstr, unsigned n, char delim='\n')
  {
    char     c = getchar();
    unsigned i = 0;

    while ( !(__status & __fail_flag ) && i < n-1 &&  c != delim)
    {
      cstr[i] = c;
      ++i;
      c = getchar(); 
    }
    if ((__status & __fail_flag) && i > 0)
      __status &= ~__fail_flag; // Unset fail flag by using & on the complement of the fail flag;
    cstr[i] = '\0';
  }

  void readFileIntoString(faststring &str)
  {
    char c;
    
    str.clear();
    
    c = getchar();
    while (!(__status & __fail_flag))
    {
      str.push_back(c);
      c = getchar();
    }
  }
  

  char lastchar()
  {
    char c;

    if ( buffer_pos != buffer_end )
      c = *(buffer_pos-1);
    else
    {
      __status |= __fail_flag;      // Set fail flag.
      __status &= ~__good_flag;     // Unset good flag.
      return 0;
    }

    /* It could be a '\r' in mac format */
    if (c == '\r')
      c = '\n';

    return c;
  }

  char lastrawchar()
  {
    char c;

    if ( buffer_pos != buffer_end )
      c = *(buffer_pos-1);
    else
    {
      __status |= __fail_flag;      // Set fail flag.
      __status &= ~__good_flag;     // Unset good flag.
      return 0;
    }
    return c;
  }


};




#endif
