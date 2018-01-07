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
#ifndef FASTSTRING_H
#define FASTSTRING_H

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstring>
#include <cctype>
#include <climits>
#include <stdarg.h>
#include <vector>
#include <cmath>
#include <cerrno>
#include <list>
#include <cfloat>
#include <climits>

#define _MINBUFFER_CAPACITY 4

// Uncomment when using a 64bit windows machine. 
// #define _WIN64

// Includes the symbols which occur in the blosum matrices: 
const char toupper_lookup[] = {
  '\0',  '\1',  '\2',  '\3',  '\4',  '\5',  '\6',  '\7', '\10', '\11', //   0- 9
 '\12', '\13', '\14', '\15', '\16', '\17', '\20', '\21', '\22', '\23', //  10-19
 '\24', '\25', '\26', '\27', '\30', '\31', '\32', '\33', '\34', '\35', //  20-29
 '\36', '\37', '\40', '\41', '\42', '\43', '\44', '\45', '\46', '\47', //  30-39
 '\50', '\51', '\52', '\53', '\54', '\55', '\56', '\57', '\60', '\61', //  40-49
 '\62', '\63', '\64', '\65', '\66', '\67', '\70', '\71', '\72', '\73', //  50-59
 '\74', '\75', '\76', '\77','\100','\101','\102','\103','\104','\105', //  60-69
'\106','\107','\110','\111','\112','\113','\114','\115','\116','\117', //  70-79
'\120','\121','\122','\123','\124','\125','\126','\127','\130','\131', //  80-89
'\132','\133','\134','\135','\136','\137','\140','\101','\102','\103', //  90-99
'\104','\105','\106','\107','\110','\111','\112','\113','\114','\115', // 100-109
'\116','\117','\120','\121','\122','\123','\124','\125','\126','\127', // 110-119
'\130','\131','\132','\173','\174','\175','\176','\177',               // 120-127
// extended ascii code:
'\200','\201',
'\202','\203','\204','\205','\206','\207','\210','\211','\212','\213',
'\214','\215','\216','\217','\220','\221','\222','\223','\224','\225',
'\226','\227','\230','\231','\232','\233','\234','\235','\236','\237',
'\240','\241','\242','\243','\244','\245','\246','\247','\250','\251',
'\252','\253','\254','\255','\256','\257','\260','\261','\262','\263',
'\264','\265','\266','\267','\270','\271','\272','\273','\274','\275',
'\276','\277','\300','\301','\302','\303','\304','\305','\306','\307',
'\310','\311','\312','\313','\314','\315','\316','\317','\320','\321',
'\322','\323','\324','\325','\326','\327','\330','\331','\332','\333',
'\334','\335','\336','\337','\340','\341','\342','\343','\344','\345',
'\346','\347','\350','\351','\352','\353','\354','\355','\356','\357',
'\360','\361','\362','\363','\364','\365','\366','\367','\370','\371',
'\372','\373','\374','\375','\376','\377'
};

const char tolower_lookup[] = {
  '\0',  '\1',  '\2',  '\3',  '\4',  '\5',  '\6',  '\7', '\10', '\11', //   0- 9
 '\12', '\13', '\14', '\15', '\16', '\17', '\20', '\21', '\22', '\23', //  10-19
 '\24', '\25', '\26', '\27', '\30', '\31', '\32', '\33', '\34', '\35', //  20-29
 '\36', '\37', '\40', '\41', '\42', '\43', '\44', '\45', '\46', '\47', //  30-39
 '\50', '\51', '\52', '\53', '\54', '\55', '\56', '\57', '\60', '\61', //  40-49
 '\62', '\63', '\64', '\65', '\66', '\67', '\70', '\71', '\72', '\73', //  50-59
 '\74', '\75', '\76', '\77','\100','\141','\142','\143','\144','\145', //  60-69
'\146','\147','\150','\151','\152','\153','\154','\155','\156','\157', //  70-79
'\160','\161','\162','\163','\164','\165','\166','\167','\170','\171', //  80-89
'\172','\133','\134','\135','\136','\137','\140','\141','\142','\143', //  90-99
'\144','\145','\146','\147','\150','\151','\152','\153','\154','\155', // 100-109
'\156','\157','\160','\161','\162','\163','\164','\165','\166','\167', // 110-119
'\170','\171','\172','\173','\174','\175','\176','\177',               // 120-127
  // extended ascii code:
'\200','\201',
'\202','\203','\204','\205','\206','\207','\210','\211','\212','\213',
'\214','\215','\216','\217','\220','\221','\222','\223','\224','\225',
'\226','\227','\230','\231','\232','\233','\234','\235','\236','\237',
'\240','\241','\242','\243','\244','\245','\246','\247','\250','\251',
'\252','\253','\254','\255','\256','\257','\260','\261','\262','\263',
'\264','\265','\266','\267','\270','\271','\272','\273','\274','\275',
'\276','\277','\300','\301','\302','\303','\304','\305','\306','\307',
'\310','\311','\312','\313','\314','\315','\316','\317','\320','\321',
'\322','\323','\324','\325','\326','\327','\330','\331','\332','\333',
'\334','\335','\336','\337','\340','\341','\342','\343','\344','\345',
'\346','\347','\350','\351','\352','\353','\354','\355','\356','\357',
'\360','\361','\362','\363','\364','\365','\366','\367','\370','\371',
'\372','\373','\374','\375','\376','\377'
};



#define macromin(x,y) ((x)<(y) ? (x) : (y))
//#define macromax(x,y) ((x)>(y) ? (x) : (y))

class faststring;

template <typename Container>
int
split (Container &l, const faststring &s, char const * const ws = "\r\n\t\v\f ");

template <typename Container>
int
split_at (Container &l, const faststring &s, const faststring &delim);

template <typename Container>
int
split_respect (Container &l, const faststring &s, char const * const ws = "\r\n\t\v\f ");

template <typename Container>
int
split_strict (Container &l, const faststring &s, char const * const ws = " \t\n");



namespace {
  inline bool is_in_symbol_list__(char c, const char* wstr="\r\n\t\v\f ")
  {
    return (strchr(wstr,c) != NULL);
  }

  inline bool is_in_symbol_list__(char c, const char* wstr, size_t n)
  {
    return (memchr(wstr,c,n) != NULL);
  }

  const char * find_fundamental ( const char *b1, const char *e1,
				  const char *b2, const char *e2)
  {
    for (; b1 < e1 ; ++b1)
    {
      const char * tb1 = b1;
      const char * tb2 = b2;
      for (; tb2 < e2 && tb1 < e1; ++tb2, ++tb1)
      {
	if (*tb1 != *tb2)
	  break;
      }
      if (tb2 == e2) // we have a match
      {
	return b1;
      }
    }
    return NULL;
  }

  const char * rfind_fundamental ( const char *rb1, const char *re1,
				   const char *rb2, const char *re2)
  {
    for (; rb1 > re1 ; --rb1)
    {
      const char * trb1 = rb1;
      const char * trb2 = rb2;
      for (; trb2 > re2 && trb1 > re1; --trb2, --trb1)
      {
	if (*trb1 != *trb2)
	  break;
      }
      if (trb2 == re2) // we have a match
      {
	return trb1+1;
      }
    }
    return NULL;
  }


} // End namespace




class faststring
{
 public:
  // Shadow the "global" size_t typedef for this class.
  typedef unsigned size_t;
  static const size_t npos = UINT_MAX;

  // Possibly even more efficient:
  //    Use pointer _end, _capacity_end instead of _len, _capacity


 private:
  char       *_buf;
  size_t     _len;
  size_t     _capacity;   // Capacity of container, i.e number of bytes available
                           // for writing. If a capacity is requested by the user,
                           // an additional byte is reserved in order to be able to append
                           // the additional \0 which is appended by some function calls.
                           // Consequently, the capacity reported to user is _capacity-1.
 
protected:

 public:
  ~faststring ()
  {
    if (_buf)
      std::free(_buf);
  }

  faststring ():_buf(NULL),  _len(0), _capacity(0) // :_len(0), _capacity(a)
  {
    //    _buf = NULL(char *)malloc(_MINBUFFER_CAPACITY);
  }

  faststring (const char *str)
  {
    _len      = strlen(str);
    _capacity = _len+1;
    _buf      = (char *)malloc(_capacity); 
    if (_buf == 0)
    {
      std::cerr << "Critical error: malloc failed when requesting " << _capacity << " bytes in faststring constructor." << std::endl;
      exit(-1);
    }
    memcpy(_buf, str, _len);
    
  }

  faststring (const char *str_begin, const char *str_end)
  {
    if (str_end > str_begin)
    {
      _len      = str_end - str_begin;
      _capacity = _len + 1;
      _buf      = (char *)malloc(_capacity);
    if (_buf == 0)
    {
      std::cerr << "Critical error: malloc failed when requesting " << _capacity << " bytes in faststring constructor." << std::endl;
      exit(-1);
    } 
      memcpy(_buf, str_begin, _len);
    }
    else
    {
      _buf = NULL;
      _len = _capacity = 0;
    }
  }

  faststring (const faststring &a):_len(a._len),_capacity(a._len+1)
  {
    if (a._buf == NULL)
    {
      _buf = 0;
      _len = 0;
      _capacity = 0;
      return;
    }

    _buf = (char *)malloc(_capacity);
    if (_buf == 0)
    {
      std::cerr << "Critical error: malloc failed when requesting " << _capacity << " bytes in faststring constructor." << std::endl;
      exit(-1);
    }
    memcpy(_buf, a._buf, _len);
  }

 faststring(char c, size_t n=1):_len(n),_capacity(n+1)
  {
    _buf = (char *)malloc(_capacity);
    if (_buf == 0)
    {
      std::cerr << "Critical error: malloc failed when requesting " << _capacity << " bytes in faststring constructor." << std::endl;
      exit(-1);
    }

    memset(_buf, c, _len);
    /*
    size_t i;
    for (i=0; i<n; ++i)
    _buf[i] = c; */
  }

 faststring(int i)
  {
    _capacity = (size_t)(sizeof(i)*8*.31+3);   // bits * 1/log_2(10) + 3 for ciel and sign and \0
    _buf = (char *)malloc(_capacity);
    if (_buf == 0)
    {
      std::cerr << "Critical error: malloc failed when requesting " << _capacity << " bytes in faststring constructor." << std::endl;
      exit(-1);
    }
    snprintf(_buf, _capacity, "%d", i);
    _len = strlen(_buf);
    //    std::cout << "Capacity: " << _capacity << std::endl;
  }

 faststring(unsigned i)
  {
    _capacity = (size_t)(sizeof(i)*8*.31+3);   // bits * 1/log_2(10) + 3 for ciel and sign and \0
    _buf = (char *)malloc(_capacity);
    if (_buf == 0)
    {
      std::cerr << "Critical error: malloc failed when requesting " << _capacity << " bytes in faststring constructor." << std::endl;
      exit(-1);
    }
	snprintf(_buf, _capacity, "%d", i);
    _len = strlen(_buf);
    //    std::cout << "Capacity: " << _capacity << std::endl;
  }

 // Prints the unsigned i into the string. The number appears right-justified.
 // len is the minimum length of the field. If the number needs more digits it well get them.
 // If the number need less than len digits they will be left padded with the fill char.
 // Example: Calling constructor with 999, '#', 5 will result in the string "##999", for 999, '#', 2 we get "999". 
 faststring(unsigned i, char fill, int len)
  {
    int num_digits = (i>0) ? ((int) std::log10 ((double)i))+1 : 1; // We take care not to pass 0 to log10  

    _capacity = num_digits+1;
    if (len >= num_digits)
      _capacity = len+1;
    _buf = (char *)malloc(_capacity);
    if (_buf == 0)
    {
      std::cerr << "Critical error: malloc failed when requesting " << _capacity << " bytes in faststring constructor." << std::endl;
      exit(-1);
    }
    snprintf(_buf,_capacity, "%0*u", len, i); // left padded with zeros (0)
    // Replace left padded zeros with fill char

    size_t pos = 0;
    size_t e   = _capacity-2;
    while (pos < e && _buf[pos]=='0')
    {
      _buf[pos] = fill;
      ++pos;
    }
    ++e;
    while (pos < e && _buf[pos] )
      ++pos;
    _len = pos;
    //    std::cout << "Capacity: " << _capacity << std::endl;
  }

  faststring(long i)
  {
    _capacity = (size_t)(sizeof(i)*8*.31+3);   // bits * 1/log_2(10) + 3 for ciel and sign and \0
    _buf = (char *)malloc(_capacity);
    if (_buf == 0)
    {
      std::cerr << "Critical error: malloc failed when requesting " << _capacity << " bytes in faststring constructor." << std::endl;
      exit(-1);
    }
    snprintf(_buf, _capacity, "%ld", i);
    _len = strlen(_buf);
  }

  faststring(unsigned long i)
  {
    _capacity = (size_t)(sizeof(i)*8*.31+3);   // bits * 1/log_2(10) + 3 for ciel and sign and \0
    _buf = (char *)malloc(_capacity);
    if (_buf == 0)
    {
      std::cerr << "Critical error: malloc failed when requesting " << _capacity << " bytes in faststring constructor." << std::endl;
      exit(-1);
    }
    snprintf(_buf, _capacity,"%lu", i);
    _len = strlen(_buf);
  }

  faststring(float f, int pres=6)
  {
    _capacity = FLT_MAX_10_EXP + pres + 3; // +3 sign, decimal point and \0
    _buf = (char *)malloc(_capacity);
    if (_buf == 0)
    {
      std::cerr << "Critical error: malloc failed when requesting " << _capacity << " bytes in faststring constructor." << std::endl;
      exit(-1);
    }
    snprintf(_buf, _capacity,"%.*f", pres, f);
    _len = strlen(_buf);
  }

  faststring(double f, int pres=14)
  {
    _capacity = DBL_MAX_10_EXP + pres + 3; // +3 sign, decimal point and \0
    _buf = (char *)malloc(_capacity);
    if (_buf == 0)
    {
      std::cerr << "Critical error: malloc failed when requesting " << _capacity << " bytes in faststring constructor." << std::endl;
      exit(-1);
    }
    snprintf(_buf,_capacity, "%.*f", pres, f);
    _len = strlen(_buf);
  }

  // In the following function, the parameter bool b is not used directly. It has been introduced to distinguish this variant of the constructor
  // from the variants declared above.
  faststring(double f, int pres, bool b)
  {
    (void) b;

    _capacity = DBL_MAX_10_EXP + pres + 5; // +5 sign, decimal point, E, sign and \0
    _buf = (char *)malloc(_capacity);
    if (_buf == 0)
    {
      std::cerr << "Critical error: malloc failed when requesting " << _capacity << " bytes in faststring constructor." << std::endl;
      exit(-1);
    }
    snprintf(_buf,_capacity, "%.*g", pres, f);
    _len = strlen(_buf);

    //    std::cerr << "faststring len: " << _len << std::endl;
  }

  faststring &operator+=(const faststring &s)
  {
    this->append(s);
    return *this;
  }

  void push_back(char c)
  {
    if (_capacity == _len)
    {
      if (_buf == 0)
      {
	_capacity = _MINBUFFER_CAPACITY;
	_buf = (char *)malloc(_capacity);
	if (_buf == 0)
	{
	  std::cerr << "Critical error: malloc failed when requesting " << _capacity << " bytes in faststring push_back." << std::endl;
	  exit(-1);
	}
      }
      else
      {
	_capacity <<= 1;
	_buf = (char *)realloc(_buf, _capacity);
	if (_buf == 0)
	{
	  std::cerr << "Critical error: realloc failed when requesting " << _capacity << " bytes in faststring push_back." << std::endl;
	  exit(-1);
	}
      }
    }
    _buf[_len] = c;
    ++_len;
  }



  void pop_back()
  {
    if (_len >0)
      --_len;
  }

  char get_last()
  {
    if (_len >0)
      return _buf[_len-1];
    else
      return '\0';
  }

  char get_last_and_pop()
  {
    if (_len >0)
    {
      --_len;
      return _buf[_len];
    }
    else
      return '\0';
  }

  // Sometimes strings have to end with a certain character.
  // If this is present, nothing will be done.
  // If this is not present, the character will be appended.
  void ensure_string_ends_with(char c)
  {
    if (get_last() != c)
    {
      push_back(c);
    }
  }

/*   void ensure_string_does_not_end_with(char c) */
/*   { */
/*     while (_len > 0 && get_last() == c) */
/*     { */
/*       pop_back(); */
/*     } */
/*   } */

// chomp in perl: (Removing a trailing \0 is not allowed.
  void pop_trailing_characters(char c)
  {
    // get_last() and pop_back() both ensure that _len>0
    if (c == '\0') // Problem: get_last returns \0 if the len=0, so we would end up with an infinite loop if we remove \0 in an empty string.
      return;

    while (get_last() == c)
    {
      pop_back();
    }
  }

  const char *c_str()
  {
    if (_capacity == 0) // We reserve 1 byte for the empty cstring.
    {
      _capacity = 1;
      _len = 0;
      _buf      = (char *)malloc(_capacity);
      if (_buf == 0)
      {
	std::cerr << "Critical error: malloc failed when requesting " << _capacity << " bytes in faststring c_str()." << std::endl;
	exit(-1);
      }
      *_buf = '\0';
      return _buf;
    }

    if (_capacity == _len)
    {
      ++_capacity;
      _buf = (char *)realloc(_buf, _capacity);
    }

    //  if (_buf[_len] == '\0')
    //    return _buf;

    _buf[_len] = '\0';
    return _buf;
  }


  char *begin() const
  {
    return _buf;
  }

  char *begin_str() const
  {
    return _buf;
  }

  char *data() const
  {
    return _buf;
  }

  char *end() const
  {
    return _buf+_len;
  }

  char *end_str() const
  {
    return _buf+_len;
  }

  char *rbegin() const
  {
    return _buf+_len-1;
  }


  char *rend() const
  {
    return _buf-1;
  }

  bool empty() const
  {
    return (_len == 0);
  }

  void reverse()
  {
    char *i1, *i2;
    char c;

    i1 = begin();
    i2 = rbegin();

    while (i1<i2)
    {
      c = *i1;
      *i1 = *i2;
      *i2 = c;
      ++i1;
      --i2;
    }
  }


  //=================================================
  // assign family:
  //=================================================
  // string& assign ( const string& str );
  // string& assign ( const string& str, size_t pos, size_t n );
  // string& assign ( const char* s, size_t n );
  // string& assign ( const char* s );
  // string& assign ( size_t n, char c );
  // template <class InputIterator>
  //   string& assign ( InputIterator first, InputIterator last );

  faststring& assign(const char *str)
  {
    _len = strlen(str);
    reserve(_len);
    memcpy(_buf, str, _len);
    return *this;
  }

  // n is the maximum number of chars that will be assigned.
  faststring& assign(const char *str, size_t n)
  {
    _len = n;
    size_t ll = strlen(str);
    if (_len > ll)  // if n > ll, we have to reduce _len to ll
      _len = ll;
    reserve(_len);
    memcpy(_buf, str, _len);
    return *this;
  }

  faststring& assign(const char *str_begin, const char *str_end)
  {
    if (str_end > str_begin)
    {
      _len = str_end - str_begin;
      reserve(_len);
      memcpy(_buf, str_begin, _len);
    }
    else
    {
      clear();
    }
    return *this;
  }

  faststring& assign(const faststring &str)
  {
    _len = str._len;
    reserve(_len);
    memcpy(_buf, str._buf, _len);
    return *this;
  }

  faststring& assign(const faststring &str, size_t pos, size_t n=npos )
  {
    if (n == npos || n+pos > str._len)
      n = str._len-pos;
    _len = n;
    reserve(_len+1);
    memcpy(_buf, str._buf+pos, n);
    return *this;
  }

  // Could be made more efficient. TODO
  faststring& assign( size_t n, char c )
  {
    clear();
    reserve(n);
    size_t i;
    for (i=0; i<n; ++i)
      push_back(c);
    return *this;  
  }
  



   //=================================================
   // The family of append functions:
   //=================================================
  // basic_string& append(const basic_string& s)	basic_string	 Append s to *this.
  // basic_string& append(const basic_string& s, 
  //                      size_type pos, size_type n)
  // basic_string	 Append a substring of s to *this.
  // basic_string& append(const charT* s)	basic_string	 Append s to *this.
  // basic_string& append(const charT* s, size_type n)	basic_string	 Append the first n characters of s to *this.
  // basic_string& append(size_type n, charT c)	basic_string	 Append n copies of c to *this.
  // template <class InputIterator>
  // basic_string& append(InputIterator first, InputIterator last)
  // basic_string	 Append a range to *this.


  // step imp-fin
  faststring& append(const char *s)
  {
    int old_len = _len;
    int s_len   = strlen(s);
    int new_len = old_len + s_len;

    _len = new_len;

    reserve(new_len);
    memcpy(_buf+old_len, s, s_len);
    return *this;
  }

  // step imp-fin
  faststring& append(const faststring &s)
  {
    int old_len = _len;
    int s_len   = s._len;
    int new_len = old_len + s_len;

    _len = new_len;

    reserve(new_len);
    memcpy(_buf+old_len, s._buf, s_len);
    return *this;
  }

  // step imp-fin
  faststring& append(const faststring& s, size_t pos, size_t n=npos)
  {
    if (n==npos || pos+n > s.length() )
      n = s.length()-pos;

    int old_len = _len;
    int new_len = old_len + n;

    _len = new_len;

    reserve(new_len);
    memcpy(_buf+old_len, s._buf+pos, n);
    return *this;
  }

  // step imp-fin
  faststring& append(const char* s, size_t n)
  {
    int old_len = _len;
    int new_len = old_len + n;

    _len = new_len;

    reserve(new_len);
    memcpy(_buf+old_len, s, n);
    return *this;
  }

  // step imp-fin
  faststring& append(size_t n, char c)
  {
    size_t old_len = _len;
    size_t new_len = old_len + n;

    reserve(new_len);
    size_t i;
    for (i=0; i<n; ++i)
      push_back(c);
    return *this;  
  }

  // step imp-fin
  faststring& append(const char * first_it, const char * end_it)
  {
    size_t old_len = _len;
    size_t n       = end_it-first_it;
    size_t new_len = old_len + n;

    _len = new_len;

    reserve(new_len);
    memcpy(_buf+old_len, first_it, n);
    return *this;  
  }



   //---------------------------------------
   // Side effects:  Invalidates pointers.
   //---------------------------------------
  void reserve(size_t s)
  {
    if (_capacity >= s)
      return;

    ++s;   // Reserve one more char in case we have to append a \0 at end of string later.
           // Well this could be seen as a wast of memory.
    if (_buf == NULL)
    {
      _buf = (char *)malloc(s);
      if (_buf == 0)
      {
	std::cerr << "Critical error: Trying again to malloc " << s << " bytes in faststring reserve()." << std::endl;
	_buf = (char *)malloc(s);
      }
      if (_buf == 0)
      {
	std::cerr << "Critical error: malloc failed when requesting " << s << " bytes in faststring reserve()." << std::endl;
	exit(-1);
      }
    }
    else
    {
      _buf = (char *)realloc(_buf, s);
      if (_buf == 0)
      {
	std::cerr << "Critical error: realloc failed when requesting " << s << " bytes in faststring reserve()." << std::endl;
	exit(-1);
      }
    }
    _capacity = s;
  }



  // Returns the part of the string following the last occurrence of delim.
  // If the character does not occur, the full string is returned.
  // The second parameter can be used to specify a minimum number of characters that need to
  // be found before the extension. fpos_minimum == 1: The extension cannot start at index 0.
  faststring get_string_extension_or_all(char delim, size_t fpos_minimum = 0)
  {
    size_t fpos = this->rfind(delim);

    if (fpos < fpos_minimum)
      fpos = npos;

    if (fpos == npos)
      return *this;
    else
    {
      faststring res(this->begin()+fpos+1, this->end());
      return res;
    }    
  }

  // Returns the part of the string following the last occurrence of delim.
  // If the character does not occur, the empty string is returned.
  // The second parameter can be used to specify a minimum number of characters that need to
  // be found before the extension. fpos_minimum == 1: The extension cannot start at index 0.
  faststring get_string_extension(char delim, size_t fpos_minimum = 0)
  {
    size_t fpos = this->rfind(delim);

    if (fpos < fpos_minimum)
      fpos = npos;

    if (fpos == npos)
      return faststring();
    else
    {
      faststring res(this->begin()+fpos+1, this->end());
      return res;
    }
  }


  // Returns the part of the string before the last occurrence of delim.
  // If the character does not occur, the empty string is returned.
  // The second parameter can be used to specify a minimum number of characters that need to
  // be found before the extension. fpos_minimum == 1: The extension cannot start at index 0.
  faststring get_string_before_extension(char delim, size_t fpos_minimum = 0)
  {
    size_t fpos = this->rfind(delim);

    if (fpos < fpos_minimum)
      fpos = npos;

    if (fpos == npos)
      return faststring();
    else
    {
      faststring res(this->begin(), this->begin()+fpos);
      return res;
    }    
  }

  // Returns the part of the string before the last occurrence of delim.
  // If the character does not occur, the full string is returned.
  // The second parameter can be used to specify a minimum number of characters that need to
  // be found before the extension. fpos_minimum == 1: The extension cannot start at index 0.
  faststring get_string_before_extension_or_all(char delim, size_t fpos_minimum = 0)
  {
    size_t fpos = this->rfind(delim);

    if (fpos < fpos_minimum)
      fpos = npos;

    if (fpos == npos)
      return *this;
    else
    {
      faststring res(this->begin(), this->begin()+fpos);
      return res;
    }    
  }


  // Filename functions:
  // Depricated: should be removed soon. Wrong name.
  faststring filename_base()
  {
    faststring res(*this);
    res.shorten(res.rfind('.'));
    return res;
  }

  faststring filename_dirname_and_basename()
  {
    faststring res = filename_dirname();
    res += "/";
    res += filename_basename();
    return res;
  }

  // Acts like the linux command does, but accepts only one argument.
  faststring filename_basename()
  {
    if (size() == 0) // If string is empty before removing anything we return ".".
      return ".";

    faststring tmp = *this;
    tmp.removeSpacesBack(" /"); // Removes trailing spaces and "/" characters

    if (tmp.size() == 0) // Note: Multiple //// are collapsed to "/"
      return "/";

    faststring res = tmp.get_string_extension_or_all('/');
    return res;
  }

  // Acts like the linux command does, but accepts only one argument.
  faststring filename_dirname()
  {
    if (size() == 0) // If string is empty before removing anything we return ".".
      return ".";

    faststring tmp = *this;
    tmp.removeSpacesBack(" /"); // Removes trailing spaces and "/" characters

    if (tmp.size() == 0) // Note: Multiple //// are collapsed to "/"
      return "/";

    faststring before_delim = tmp.get_string_before_extension('/');

    if (before_delim.empty() )
    {
      if (tmp[0] == '/')  // Example: /usr -> dirname must be "/". If the only / is the root, this is the dirname
	return "/";

      return ".";
    }
    return before_delim;
  }


  // Deprecated. - Replaced by filename_dirname() which mimics the dirname command in linux much better.
  faststring filename_path()
  {
    faststring res(*this);
    size_t fpos = res.rfind('/');

    if (fpos == npos)
      return ".";
    else
    {
      res.shorten(fpos);
      return res;
    }
  }

  faststring filename_basename_without_extension()
  {
    faststring tmp  = filename_basename();
    return tmp.get_string_before_extension_or_all('.',1);
  }

  faststring filename_dirname_and_basename_without_extension()
  {
    faststring tmp  = filename_dirname();
    tmp += "/";
    tmp += filename_basename_without_extension();
    return tmp;
  }



  // Returns the last part of the filename, i.e. the part after the last relevant "/".
  // If this sting contains the basename this is the base name with extension.
  // If the string ends with a directory, this is the directory.
  // Some special cases are catched.
  // Deprecated and replaced by the filename_basename function.

/*   faststring filename_FileOrDirName() */
/*   { */
/*     unsigned fpos = this->rfind('/'); */

/*     if (fpos == npos) */
/*       return *this; */
/*     else */
/*     { */
/*       faststring res(this->begin()+fpos+1, this->end()); */
/*       return res; */
/*     } */
/*   } */

  faststring filename_extension()
  {
    // If we do not restrict the search to the basename, dots in the path can cause problems
    // if the file does not have a file extension.
    faststring tmp = filename_basename();
    return tmp.get_string_extension('.', 1);
  }

  // END filename functions

  // Searches the string backwards for delim. Divides the string at delim if delim is found in backward direction
  // before break_at. So if in forward direction delim occurs only before break_at or if delim does not occur, the 
  // string cannot be divided in the requested way. In this case pre gets the string und post is empty.

  bool backwards_search_divide_at_break_at(char delim, char break_at, faststring &pre, faststring &post) const
  {
    size_t len = length();
    size_t i = len-1;
    int    mode = 0;

    while (i > 0)
    {
      if (get_unckecked(i) == delim) // Yea, we found delim behind any occurrence of break_at
      {
	mode = 2;
	break;
      }
      if (get_unckecked(i) == break_at) // We could not find the delim before break_at, so we have to stop.
      {
	mode = 1;
	break;
      }

      --i;
    }

    if (mode == 2)
    {
      pre  = substr(0, i);
      post = substr(i+1, npos);
      return true;
    }
    else // if mode == 0 or 1
    {
      pre = *this;
      post = "";
      return false;
    }
  }

  bool divide_at(char delim, faststring &pre, faststring &post) const
  {
    size_t pos = find(delim,0);
    if (pos == npos)
    {
      post.clear();
      pre = *this;
      return false;
    }
    post = substr(pos+1, npos);
    pre  = substr(0, pos);
    return true;
  }

  bool divide_at_first_of(const char *delims, faststring &pre, faststring &post) const
  {
    size_t pos = find_first_of(delims,0);
    if (pos == npos)
      return false;
    post = substr(pos+1, npos);
    pre  = substr(0, pos);
    return true;
  }

  bool divide_at_last_of(const char *delims, faststring &pre, faststring &post) const
  {
    size_t pos = find_last_of(delims);
    if (pos == npos)
      return false;
    post = substr(pos+1, npos);
    pre  = substr(0, pos);
    return true;
  }

  bool divide_at(faststring delim, faststring &pre, faststring &post) const
  {
    size_t pos = find(delim,0);
    if (pos == npos)
      return false;
    post = substr(pos+delim.size(), npos);
    pre  = substr(0, pos);
    return true;
  }

  size_t size() const
  {
    return _len;
  }

  size_t length() const
  {
    return _len;
  }

  size_t capacity() const
  {
    return _capacity; // One char is reserved for the \0 at end of string
  }

  void set_unckecked(size_t pos, char c)
  {
    _buf[pos] = c;
  }

  bool check_pos(size_t pos)
  {
    return (pos < _len);
  }

  bool set(size_t pos, char c)
  {
    if ( check_pos(pos) )
    {
      _buf[pos] = c;
      return true;
    }
    else
    {
      return false;
    }
  }


  void get_unckecked(size_t pos, char &c) const
  {
    c = _buf[pos];
  }

  char get_unckecked(size_t pos) const
  {
    return _buf[pos];
  }

  char & at(size_t pos)
  {
    return _buf[pos];
  }

  bool get(size_t pos, char &c)
  {
    if ( check_pos(pos) )
    {
      c = _buf[pos];
      return true;
    }
    else
    {
      return false;
    }
  }


  faststring& ToLower()
  {
    tolower();
    return *this;
  }

  faststring& ToUpper()
  {
    toupper();
    return *this;
  }

  void toupper_slow()
  {
    char *b = begin();
    char *e = end();

    while (b != e)
    {
      *b = std::toupper(*b);
      ++b;
    }
  }

  void toupper()
  {
    unsigned char *b = (unsigned char *) begin();
    unsigned char *e = (unsigned char *) end();

//    char *b = (char *) begin();
//    char *e = (char *) end();

    while (b != e)
    {
      *b = toupper_lookup[*b];
      ++b;
    }
  }

  void erase_front(size_t n)
  {
    if (n>0)
    {
      if (n > _len)
      {
	_len = 0;
	return;
      }
      else
      {
	_len -=n;
      }
      memmove(_buf, _buf+n, _len);
    }
  }

  void erase_back(size_t n)
  {
    if (n > _len)
      _len = 0;
    else
      _len -= n;
  }

  void tolower_slow()
  {
    char *b = begin();
    char *e = end();

    while (b != e)
    {
      *b = std::tolower(*b);
      ++b;
    }
  }

  void tolower()
  {
    unsigned char *b = (unsigned char *) begin();
    unsigned char *e = (unsigned char *) end();

    while (b != e)
    {
      *b = tolower_lookup[*b];
      ++b;
    }
  }

  void shorten(size_t new_len)
  {
    _len = macromin(_len, new_len);
  }

  void shorten_to_first_occurrence_of(char c)
  {
    size_t pos = find(c,0);
    shorten(pos);
  }

  // See also resize for a different version of this function.
  // Appends multiple characters to fill the string if it is shorter than wanted. 
   void fill_if_shorter ( size_t n, char c )
   {
     if (n <= _len)
       return;
     append(n-_len,c);
   }

  faststring &operator=(const faststring &a)
  {
    _len = a._len;

    //    std::cerr << "_len: " << _len << std::endl;

    reserve(_len+1);
    memcpy(_buf, a._buf, _len);
    return *this;
  }

  friend bool operator==(const faststring &a, const faststring &b)
  {
    if (a._len != b._len)
      return false;
    return (strncmp(a._buf, b._buf, a._len) == 0);
  }

  friend bool operator!=(const faststring &a, const faststring &b)
  {
    if (a._len != b._len)
      return true;
    return (strncmp(a._buf, b._buf, a._len) != 0);
  }

  // Has been corrected in 03.2011
  friend bool operator<(const faststring &a, const faststring &b)
  {
    int l   =  macromin(a._len, b._len );
    //     int res =  (strncmp(a._buf, b._buf, l) < 0);
    short res =  strncmp(a._buf, b._buf, l);

    if (res == 0)
      return (a._len < b._len );
    //     return res; // Was wrong in earlier versions.
    return (res < 0);
  }

  // Caution: strcmp assumes \0 terminated strings.
  //          Since we do not have \0 terminated strings we
  //          use strncmp with n such that we do not pass the ends
  //          of the strings.
  friend int fstrcmp(const faststring &a, const faststring &b)
  {
    int   l   =  macromin(a._len, b._len );
    short res =  strncmp(a._buf, b._buf, l);

    if (res == 0)  // We only compared the first l chars.
    {
      if (a._len < b._len)
	return -1;
      if (a._len > b._len)
	return 1;
      return 0;
    }
    return res;
  }

  friend int fstrncmp(const faststring &a, const faststring &b, int n)
  {
    int   l   = macromin(a._len, b._len );
    int   m   = macromin(l, n);

    short res =  strncmp(a._buf, b._buf, m);

    // We do not have \0 terminated strings, so we have to restrict the initial
    // comparison to the length m.

    // If n is larger than the length of any of the two strings
    // we need to check the lengths.

    if (n <= l)   // If n is smaller or equal to the shorter string, we have a result.
      return res; // This can be <,>, ==
    // else // if n > l

    // In this case, n is larger than one of the two string lengths or larger than both.
    // Still in this case we could have res!=0 among the first m string positions.
    // If res!=0, this is the result.
  
    if (res != 0)
      return res;

    // else if res == 0, the two strings are equal among the first m positions
    // and thus they are equal over the whole length of the shorter string.

    // If n is larger than both strings, the longer string is greater.
    // If n is larger than only one string, the longer string is greater.
    // => in all cases the shorter string is less.
    //    if lengths are equal, we have to return 0

    if (a._len < b._len)
      return -1;
    if (a._len > b._len)
      return 1;
    return 0;
  }


  // Has been corrected in 03.2011
  friend bool operator>(const faststring &a, const faststring &b)
  {
    int l   =  macromin(a._len, b._len );
    //    int res =  (strncmp(a._buf, b._buf, l) > 0);
    short res =  strncmp(a._buf, b._buf, l);
    
    if (res == 0)
      return (a._len > b._len );
    return (res > 0);
  }

  void removeSpacesFront(const char* delims="\r\n\t\v\f ")
   {
     size_t i;

     for (i=0; i < _len && is_in_symbol_list__(_buf[i], delims); ++i)
     {}
     erase_front(i);
   }


   void removeSpacesBack(const char* delims="\r\n\t\v\f ")
   {
     size_t i;
     size_t n_renomve=0;

     if (_len > 0)
     {
       for (i=_len-1; is_in_symbol_list__(_buf[i], delims); --i)
       {
	 ++n_renomve;
	 if (i==0)
	   break;
       }
       erase_back(n_renomve);
     }
   }

   void removeSpaces(const char* delims="\r\n\t\v\f ")
   {
     char *i, *j, *e;
     e = end();
     i = begin();

     // Skip symbols we do not remove
     while (i != e && !is_in_symbol_list__(*i, delims))
       ++i;

     if (i == e) // Nothing was done.
       return;

     j = i;
     ++i;

     while (i != e)
     {
       if (is_in_symbol_list__(*i, delims)) // remove this
       {
	 ++i;
       }
       else
       {
	 *j = *i;
	 ++i;
	 ++j;
       }
     }
     _len -= i-j;
   }


   void remove_symbol(char c)
   {
     char *i, *j, *e;
     e = end();
     i = begin();

     // Skip symbols we do not remove
     while (i != e && *i != c)
       ++i;

     if (i == e) // Nothing was done.
       return;

     j = i;
     ++i;

     while (i != e)
     {
       if (*i == c) // remove this
       {
	 ++i;
       }
       else
       {
	 *j = *i;
	 ++i;
	 ++j;
       }
     }
     _len -= i-j;
   }

   /*
   size_t count_symbol(char sym) const
   {
     size_t res=0;
     const char *pos_it = begin();
     const char *end_it = end();

     while (pos_it != end_it)
     {
       if (*pos_it == sym)
	 ++res;

       ++pos_it;
     }
     return res;
   }
   */

   size_t countChar(char sym) const
   {
     size_t res=0;
     const char *pos_it = begin();
     const char *end_it = end();

     while (pos_it != end_it)
     {
       if (*pos_it == sym)
	 ++res;

       ++pos_it;
     }
     return res;
   }

   void skip_symbols(size_t &pos, const char* delims="\r\n\t\v\f ") const
   {
     const char *theend = end();
     const char *it=begin()+pos;

     while (it != theend && is_in_symbol_list__(*it, delims) )
       ++it;
     pos = it-begin();
   }


   void skip_symbol(size_t &pos, char c) const
   {
     const char *theend = end();
     const char *it=begin()+pos;

     while (it != theend && *it == c)
       ++it;
     pos = it-begin();
   }
  


   void skip_spaces(size_t &pos) const
   {
     skip_symbols(pos);
   }

   void replace_char(char origChar, char newChar)
   {
           char *pos_it = begin();
     const char *end_it = end();

     while (pos_it != end_it)
     {
       if (*pos_it == origChar)
	 *pos_it = newChar;
       ++pos_it;
     }
   }

   int replace_char_count(char origChar, char newChar)
   {
     unsigned    count  = 0;
           char *pos_it = begin();
     const char *end_it = end();

     while (pos_it != end_it)
     {
       if (*pos_it == origChar)
       {
	 ++count;
	 *pos_it = newChar;
       }
       ++pos_it;
     }
     return count;
   }


/*----------------------------------------------------------------------------------------------------------------------
|	Shortens stored string to n - 3 characters, setting the last three characters "...". If string is already <=
|	n characters in length, this function has no effect. This is useful when it is desirable to show some of the
|	contents of a string, even when the string will not fit in its entirety into the space available for displaying it.
|*/
   faststring& ShortenTo(size_t n)
   {
     if (n <= 3 || _len <= n )
       return *this;

     // _len > n && n>3 now:
     _len = n;
     _buf[n-1] = '.';
     _buf[n-2] = '.';
     _buf[n-3] = '.';
     return *this;
   }
   
   // TODO: Use pointers to make this more efficient.
   bool isADouble() const
   {
     size_t i=0;

     bool hasDigits    = false;
     bool hasExp       = false;
     bool hasExpDigits = false;

     while (i< _len && isspace(_buf[i]))
       ++i;

     // I am not sure this line makes it more efficient
     //     if (i >= _len)
     //       return false;

     if (i < _len && (_buf[i] == '-' || _buf[i] == '+') )
       ++i;

     while (i < _len && isdigit(_buf[i]) )
     {
       ++i;
       hasDigits = true;
     }

     if (i < _len && _buf[i] == '.')
       ++i;

     while (i < _len && isdigit(_buf[i]) )
     {
       ++i;
       hasDigits = true;
     }

     if ( i < _len && (_buf[i] == 'e' || _buf[i] == 'E') )
     {
       hasExp = true;
       ++i;
     }

     while (i < _len && isdigit(_buf[i]) )
     {
       ++i;
       hasExpDigits = true;
     }

     while (i < _len && isspace(_buf[i]) )
       ++i;

     // We passed all characters that could be interpreted as part of the double number.
     // If we have not reached the end of the string, the string contains characters that are not
     // part of a double number. If this is the case we return false.
     if (i < _len) 
       return false;

     if (hasDigits)
     {
       if (hasExp && hasExpDigits)
	 return true;
       if (!hasExp && !hasExpDigits)
	 return true;
       return false;
     }
     else
       return false;
   }

   // TODO: Use pointers to make this more efficient.
   bool isAnInt() const
   {
     size_t i=0;
     bool hasDigits = false;
     
     while (i < _len  && _buf[i] == ' ')
       ++i;

     if (i < _len && (_buf[i] == '-' || _buf[i] == '+') ) 
       ++i;
    
     while (i < _len && isdigit(_buf[i]))
     {
       ++i;
       hasDigits = true;
     }

     while (i < _len && _buf[i] == ' ')
       ++i;

     if (i < _len)
       return false;
     else
       if (hasDigits)
	 return true;
       else
	 return false;
   }

   bool isAnUnsigned() const
   {
     size_t i=0;
     bool hasDigits = false;
     
     while (i < _len && _buf[i] == ' ')
       ++i;

     if (i < _len && _buf[i] == '+')
       ++i;
    
     while (i < _len && isdigit(_buf[i]))
     {
       ++i;
       hasDigits = true;
     }

     while (i < _len && _buf[i] == ' ')
       ++i;

     if (i < _len)
       return false;
     else
       if (hasDigits)
	 return true;
       else
	 return false;
   }
   
   //---------------------------------------
   // Side effects:  Invalidates pointers.
   //---------------------------------------
   faststring &AddQuotes(char quote_char = '\'')
   {
     reserve(_len + 2);
     memmove(_buf+1, _buf, _len);
     _buf[0]      = quote_char;
     _buf[_len+1] = quote_char;
     _len += 2;

     return *this;
   }

   void consume(faststring &s)
   {
     _capacity = s._capacity;
     _len      = s._len;
     _buf      = s._buf;

     s._buf = NULL;
     s._capacity = s._len = 0;
   }

   void swap(faststring &s)
   {
      char       *tmp_buf      = _buf;
      size_t      tmp_len      = _len;
      size_t      tmp_capacity = _capacity;

      _buf      = s._buf;
      _len      = s._len;
      _capacity = s._capacity;

      s._buf      = tmp_buf;
      s._len      = tmp_len;
      s._capacity = tmp_capacity;
   }


   faststring &doubleChar(char c)
   {
     faststring tmp;
     
     size_t count_c = countChar(c);

     if (count_c > 0)
     {
       tmp.reserve(_len + count_c);
       
       size_t i;
       size_t j;
       for (i=0, j=0; i<_len; ++i, ++j)
       {
	 if (c == (tmp._buf[j] = _buf[i]) )
	 {
	   ++j;
	   tmp._buf[j] = _buf[i];
	 }
       }
       tmp._len = _len + count_c;

       consume(tmp);
     }
     return *this;
   }


   int PrintF(const char *format, ...)
   {
     char *tmp;

     // Create and initialise list of optional arguments: 
     va_list argList;

#if !defined (WIN32) && !defined (_WIN64)  
     va_start(argList, format);
     int N = vasprintf(&tmp, format, argList);
     va_end(argList);
#else // In case we are on a WIN32 or _WIN64 machine
     va_start(argList, format);
     tmp = new char [5000];
     int N = vsnprintf(tmp, 5000, format, argList);
     va_end(argList);
#endif
     _buf = tmp;
     _capacity = N+1;
     _len = N;
     return N;
   }


   // Convert to unsigned long
   unsigned ToUnsigned()
   {
     char *end;
     return strtoul(c_str(), &end, 0);
   }

   // Make sure you check errno yourself
   unsigned	ToUnsigned(size_t  beg_pos,
			   size_t& end_pos)
   {
     const char     *beg = c_str()+beg_pos;
           char     *end;
           unsigned result;

     //     errno = 0;
     result  = strtoul(beg, &end, 0);
     //     if (errno != 0)
     //     {
     //       std::cerr << "Parse error in ToUnsigned." << std::endl;
     //     }
     end_pos = end-begin_str();
     return result;
   }



   // Convert to unsigned long
   unsigned long ToUnsignedLong()
   {
     char *end;
     return strtoul(c_str(), &end, 0);
   }

   // Convert to int
   int ToInt()
   {
     char *end;
     return strtol(c_str(), &end, 0);
   }

   // Convert to long
   long ToLong()
   {
     char *end;
     return strtol(c_str(), &end, 0);
   }

   long	ToLong(size_t  beg_pos,
	       size_t& end_pos)
   {
     const char   *beg = c_str()+beg_pos;
           char   *end;
           long   result;

     result  = strtol(beg, &end, 0);
     end_pos = end-c_str();
     return result;
   }


   // Convert to double
   double ToDouble()
   {
     char *end;
     return strtod(c_str(), &end);
   }

   double ToDouble(size_t  beg_pos,
		   size_t& end_pos)
   {
     const char   *beg = c_str()+beg_pos;
           char   *end;
           double result;

     result  = strtod(beg, &end);
     end_pos = end-c_str();
     return result;
   }

   double ToDouble(size_t pos,
		   char **end)
   {
     double result;

     result  = strtod(_buf+pos, end);
     return result;
   }


   bool EqualsCaseSensitive(const faststring &s)
   {
     if (_len != s._len)
       return false;

    return (strncmp(_buf, s._buf, _len) == 0);
   }

   bool EqualsCaseInsensitive(const faststring &s) const
   {
     if (_len != s._len)
       return false;

     faststring a(*this);
     faststring b(s);

     a.tolower();
     b.tolower();

    return (strncmp(a._buf, b._buf, _len) == 0);
   }

   //=================================================
   // The family of compare functions:
   //=================================================
   // int compare ( const string& str ) const;
   // int compare ( const char* s ) const;
   // int compare ( size_t pos1, size_t n1, const string& str ) const;
   // int compare ( size_t pos1, size_t n1, const char* s) const;
   // int compare ( size_t pos1, size_t n1, const string& str, size_t pos2, size_t n2 ) const;
   // int compare ( size_t pos1, size_t n1, const char* s, size_t n2) const;
   // ------------------------------------------------- all implemented

   // Result: >: +1, <: -1, == 0
   int compare(const faststring &s)
   {
     size_t l = macromin(_len, s._len);
     int    res = strncmp(_buf, s._buf, l);

     if (res != 0)
       return res;

     if (_len == s._len)
       return 0;
     else
       return (int)(_len > s._len);
   }

   // Result: >: +1, <: -1, == 0
   int compare(const char *s)
   {
     size_t   s_len = strlen(s);
     size_t   l     = macromin(_len, s_len);
     int      res   = strncmp(_buf, s, l);

     if (res != 0)
       return res;

     if (_len == s_len)
       return 0;
     else
       return (int)(_len > s_len);
   }


   // Result: >: +1, <: -1, == 0
   int compare(const char *t1, const char*t2, const char *s1, const char *s2)
   {
     // Out of range check:
     if (t1 > t2 || s1 > s2)
       return -2;

     if (t1 < _buf || t2 > _buf+_len)
       return -3;

     return compare_unchecked(t1,t2, s1, s2);
   }


   int compare ( size_t pos1, size_t n1, const faststring& str ) const
   {
     if (pos1 > _len || pos1 + n1 > _len)
       return -3;

     return compare_unchecked(_buf+pos1, _buf+pos1+n1, str._buf, str._buf+str._len);
   }

   int compare ( size_t pos1, size_t n1, const char* s) const
   {
    if (pos1 > _len || pos1 + n1 > _len)
       return -3;

    const char *e = s + strlen(s);

    return compare_unchecked(_buf+pos1, _buf+pos1+n1, s, e);
   }

   int compare ( size_t pos1, size_t n1, const faststring& str, size_t pos2, size_t n2 ) const
   {
     if (pos1 > _len || pos1 + n1 > _len)
       return -3;
     if (pos2 > str._len || pos2 + n2 > str._len)
       return -4;

     return compare_unchecked(_buf+pos1, _buf+pos1+n1, str._buf+pos2, str._buf+pos2+n2);
   }

   int compare ( size_t pos1, size_t n1, const char* s, size_t n2) const
   {
     size_t s_len = strlen(s);

     if (pos1 > _len || pos1 + n1 > _len)
       return -3;
     if (n2 > s_len)
       return -4;

     return compare_unchecked(_buf+pos1, _buf+pos1+n1, s, s+n2);
   }


   // Result: >: +1, <: -1, == 0
   int compare_unchecked(const char *t1, const char*t2, const char *s1, const char *s2) const
   {
     size_t t_len = t2-t1; 
     size_t s_len = s2-s1; 

     size_t l = macromin(t_len, s_len);
     int    res = strncmp(t1, s1, l);

     if (res != 0)
       return res;

     if (t_len == s_len)
       return 0;
     else
       return (int)(t_len > s_len);
   }


   unsigned char * p_str(unsigned char *buffer) const // 
   {
     if (_len > 255)
     {
       *buffer = '\0';
       return 0;
     }

     memcpy(buffer+1, _buf, _len);
     buffer[0] = (unsigned char)_len;
     return buffer;
   }

   void AsBinary(unsigned long l)
   {
     short bits = sizeof(l)*8;
     reserve(bits);
     _len = bits;

     short i;
     for (i=(short)_len-1; i >= 0; --i)
     {
       // is l odd
       _buf[i] = (char)((char)(l%2)+'0');
       l/=2;
     }
   }

   void AsBinary(unsigned long l, unsigned digits)
   {
     reserve(digits);
     _len = digits;

     short i;
     for (i=(short)_len-1; i >= 0; --i)
     {
       // is l odd
       _buf[i] = (char)((char)(l%2)+'0');
       l/=2;
     }
   }


   void AsHex(unsigned long l)
   {
     short hexes = sizeof(l)*2;
     char remainder;

     reserve(hexes);
     _len = hexes;

     short i;
     for (i=(short)_len-1; i >= 0; --i)
     {
       remainder = l%16;

       if (remainder < 10)
	 _buf[i] = remainder+'0';
       else
       {
	 _buf[i] = remainder+55;
       }
       l /= 16;
     }
   }

   void AsHex(unsigned long l, unsigned digits)
   {
     char remainder;

     reserve(digits);
     _len = digits;

     short i;
     for (i=(short)_len-1; i >= 0; --i)
     {
       remainder = l%16;

       if (remainder < 10)
	 _buf[i] = remainder+'0';
       else
       {
	 _buf[i] = remainder+55;
       }
       l /= 16;
     }
   }

   faststring substr(size_t pos, size_t n=npos) const
   {
     const char * b = _buf+pos;
     const char * e = _buf+pos+n;

     if (b > _buf+_len)
       b = _buf+_len;

     if (n==npos || e > _buf+_len)
     {
       e = _buf+_len;
     }
     
     faststring res(b, e);
     return res;
   }
   
   //=================================================
   // The family of erase functions:
   //=================================================
   // string& erase ( size_t pos = 0, size_t n = npos );
   // iterator erase ( iterator position );
   // iterator erase ( iterator first, iterator last );
   // ------------------------------------------------- all implemented

   faststring &erase ( size_t pos, size_t n = npos )
   {
     char *p1 = _buf+pos;
     char *p2 = _buf+pos+n;

     if (p1 >  _buf+_len)
       p1 = _buf+_len;

     if (n==npos || p2 > _buf+_len)
       p2 = _buf+_len;
     
     return erase(p1, p2);
   }

   faststring &erase ( char *p1)
   {
     return erase(p1, p1+1);
   }

   faststring &erase ( char *p1, char *p2 )
   {
     char *tmp1 = p1;
     char *tmp2 = p2;
     char *e    = _buf + _len;

     if (p1 >= p2)
       return *this;

     if (p1 < _buf)
       p1 = _buf;

     if (p2 > e)
       p2 = e;

     // Copy chars to the left within the string.
     while (tmp2 < e)
     {
       *tmp1 = *tmp2;
       ++tmp1;
       ++tmp2;
     }

     // Reduce length
     _len -= (p2-p1);

     return *this;
   }

  void erase()
  {
    _len = 0;
  }

  void clear()
  {
    _len = 0;
  }

  void reset()
  {
    if (_buf)
    {
      free(_buf);
      _buf = NULL;
      _capacity = 0;
      _len = 0;
    }
  }

   //=================================================
   // The family of insert functions:
   //=================================================
   // - string& insert ( size_t pos1, const string& str );
   // - string& insert ( size_t pos1, const string& str, size_t pos2, size_t n );
   // - string& insert ( size_t pos1, const char* s, size_t n);
   // - string& insert ( size_t pos1, const char* s );
   // - string& insert ( size_t pos1, size_t n, char c );
   // - iterator insert ( iterator p, char c );
   // - void insert ( iterator p, size_t n, char c );
   //=================================================

   //---------------------------------------
   // Side effects:  Invalidates pointers.
   //---------------------------------------
   faststring &insert(size_t pos, const char *p1, const char *p2)
   {
     if (pos > _len || p1 >= p2)
       return *this;


     size_t N = p2-p1; // Number of characters to be inserted
     reserve(_len+N);
     
     //     char *b    = _buf;
     char *e    = _buf+_len; // This is still the old _len. We have reserved new space but we have not used it yet.
     char *tmp  = e-1;       // Pointer to the last char. We start at the last char which we move first.
     char *tmpN = tmp+N;

     // Now, after having reserved new memory, we can work with pointers in _buf
     // again.

     char *p_pos = _buf + pos;

     // xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

     // Create the space needed:
     while (tmp >= p_pos )
     {
       *tmpN = *tmp;
       --tmp;
       --tmpN;
     }
     //     *tmpN = *tmp;
     _len += N;

     replace_intern(p_pos, p1, p2);
     return *this;
   }

   faststring &insert(char *p_pos, const char *p1, const char *p2)
   {
     return insert(p_pos-_buf, p1, p2); 
   }


   faststring &insert(size_t pos1, const faststring &str)
   {
     return insert(pos1, str.begin(), str.end()); 
   }

   faststring& insert ( size_t pos1, const faststring& str, size_t pos2, size_t n )
   {
     char *b = str._buf + pos2;
     char *e = str._buf + pos2+n;

     if (b > str.end() || e > str.end() )
       return *this;

     return insert(pos1, b, e); 
   }


   faststring& insert ( size_t pos1, const char* s, size_t n)
   {
     size_t s_len = strlen(s);
     
     if (s_len < n)
       return *this;
     
     return insert(pos1, s, s+n);
   }

   faststring& insert( size_t pos1, const char* s)
   {
     size_t s_len = strlen(s);
     
     return insert(pos1, s, s+s_len);
   }

   faststring& insert( size_t pos1, size_t n, char c )
   {
     faststring str(c, n);

     return insert(pos1, str.begin(), str.end() );
   }

   faststring& insert( size_t pos1, char c )
   {
     return insert(pos1, &c, &c+1 );
   }

   faststring& insert( char *p, size_t n, char c )
   {
     faststring str(c, n);

     return insert(p-_buf, str.begin(), str.end() );
   }

   char* insert( char *p, char c )
   {
     insert(p-_buf, &c, &c+1 );
     return p;
   }
   

   //=================================================
   // The family of replace functions:
   //=================================================
   // - string& replace ( size_t pos1, size_t n1,   const string& str ); 
   // - string& replace ( iterator i1, iterator i2, const string& str ); 
   // - string& replace ( size_t pos1, size_t n1, const string& str, size_t pos2, size_t n2 ); 
   // - string& replace ( size_t pos1, size_t n1,   const char* s, size_t n2 ); 
   // - string& replace ( iterator i1, iterator i2, const char* s, size_t n2 ); 
   // - string& replace ( size_t pos1, size_t n1,   const char* s ); 
   // - string& replace ( iterator i1, iterator i2, const char* s ); 
   // - string& replace ( size_t pos1, size_t n1,   size_t n2, char c ); 
   // - string& replace ( iterator i1, iterator i2, size_t n2, char c ); 
   // - template<class InputIterator> string& replace ( iterator i1, iterator i2, InputIterator j1, InputIterator j2 );

   // replace_intern assumes that ranges have already been checked!!!
   // This simple version does no error checking and simply overwrites
   // part of the this object.
   void replace_intern(char *t1, const char *p1, const char *p2)
   {
     while (p1 < p2)
     {
       *t1 = *p1;
       ++t1;
       ++p1;
     }
   }


   //TODO: test
   // Overwrite range must "fit" into target string.
   // Otherwise, nothing is done.
   void overwrite(char *t1, char *p1, char *p2)
   {
     size_t N = p2-p1;

     if (t1 < _buf || t1+N > _buf+_len)
       return;

     replace_intern(t1, p1, p2);
   }

   // erases the range from pos1 to pos1+n1 and then inserts str
   // if pos1 is not a valid position in the string, nothing will be done.

   faststring& replace ( size_t pos1, size_t n1,   const faststring& str ) 
   {
     if (pos1 > _len)
       return *this;

     // Number of chars we want/can to "remove":
     size_t real_n1 = macromin(n1, _len-pos1);

     replace_unchecked (_buf+pos1, _buf+pos1 + real_n1, str.begin(), str.length() );

     return *this;
   }

   faststring& replace ( char *i1, char *i2, const faststring& str )
   {
     if (i1 < _buf || i2 > _buf+_len)
       return *this;

     replace_unchecked(i1, i2, str.begin(), str.length() );
     return *this;    
   }

   faststring& replace ( size_t pos1, size_t n1, const faststring& str, size_t pos2, size_t n2 )
   {
     if (pos1 >= _len)
       return *this;

     if (pos2 >= str._len)
       return *this;

     // Number of chars we can "remove":
     size_t real_n1 = macromin(n1, _len-pos1);

     // Number of chars we can "insert":
     size_t real_n2 = macromin(n2, str._len-pos2);

     replace_unchecked(_buf+pos1,_buf+pos1+real_n1, str.begin()+pos2, real_n2);

     return *this;
   }
 
   faststring& replace ( size_t pos1, size_t n1,   const char* s, size_t n2 )
   {
     if (pos1 >= _len)
       return *this;

     size_t s_len = strlen(s);

     //     if (n2 <= s_len )
     //      return *this;

     // Number of chars we can "remove":
     size_t real_n1 = macromin(n1, _len-pos1);

     // Number of chars we can "insert":
     size_t real_n2 = macromin(n2, s_len);

     replace_unchecked(_buf+pos1,_buf+pos1+real_n1, s, real_n2);

     return *this;   
   }
 
   faststring& replace ( char *i1, char *i2, const char* s, size_t n2 )
   {
     if (i1 < _buf     || i2 > _buf+_len)
       return *this;

     size_t s_len = strlen(s);

     if (n2 > s_len)
       n2 = s_len;

     replace_unchecked (i1, i2, s, n2);
     return *this;
   }

   faststring& replace(char *i1, char *i2, const char* s1, const char * s2)
   {
     if (i1 < _buf     || i2 > _buf+_len || s2 < s1)
       return *this;

     replace_unchecked (i1, i2, s1, s2-s1);
     return *this;    

   }

   void replace_unchecked ( char *i1, char *i2, const char* s, size_t n2 )
   {
     // Number of chars we want/can to "remove":
     size_t n_remove = i2-i1;
     
     // Number of chars we want to insert:
     size_t n_insert = n2;

     if (n_remove == n_insert) // simple overwrite:
     {
       replace_intern(i1, s, s+n_insert);
     }
     else if (n_remove > n_insert) // We have too much space and need to erase some characters.
     {
       erase(i1, i1+n_remove - n_insert);
       replace_intern(i1, s, s+n_insert);
     }
     else	   // We have too little space and need to insert some characters.
     {
       size_t n_newspace = n_insert-n_remove;
       
       insert(i1-_buf, s, s + n_newspace);
       replace_intern(i1 + n_newspace, s+n_newspace, s+n_insert);
     }
   }
 
   faststring& replace ( size_t pos1, size_t n1,   const char* s )
   {
     if (pos1 >= _len)
       return *this;

     size_t s_len = strlen(s);

     //     if (n2 <= s_len )
     //      return *this;

     // Number of chars we can "remove":
     size_t real_n1 = macromin(n1, _len-pos1);

     // Number of chars we can "insert":
     size_t real_n2 = s_len;

     replace_unchecked(_buf+pos1,_buf+pos1+real_n1, s, real_n2);

     return *this;    
   }
 
   faststring& replace ( char *i1, char *i2, const char* s )
   {
     if (i1 < _buf || i2 > _buf+_len)
       return *this;

     size_t s_len = strlen(s);

     replace_unchecked(i1, i2, s, s_len);

     return *this; 
   }
 
   faststring& replace ( size_t pos1, size_t n1, size_t n2, char c )
   {
     if (pos1 >= _len)
       return *this;

     size_t real_n1 = macromin(n1, _len-pos1);
   
     faststring tmp(c, n2);

     replace_unchecked(_buf+pos1, _buf+pos1+real_n1, tmp.begin(), n2);

     return *this; 
   }
 
   faststring& replace ( char *i1, char *i2, size_t n2, char c )
   {
     if (i1 < _buf || i2 > _buf+_len)
       return *this;
     
     faststring tmp(c, n2);

     replace_unchecked(i1, i2, tmp.begin(), n2);
     
     return *this; 
   } 

   //=================================================


   static size_t max_size() {  return npos; }

   // Shorten string to n characters. If its length is less than n, we append c
   // as often as necessary to get a length of n.
   void resize ( size_t n, char c )
   {
     if (n < _len)
       shorten(n);
     else if (n > _len)
       append(n-_len,c);
   }


   void resize ( size_t n )
   {
     resize(n, char());
   }

   char & operator[](size_t i)
   {
     return _buf[i];
   }

   char operator[](size_t i) const
   {
     return _buf[i];
   }

   size_t copy ( char* s, size_t n, size_t pos = 0) const
   {
     if (pos > _len)
       return -2;

     int l = n;
     if (n+pos > _len)
       l = _len-pos;

     memcpy(s, _buf+pos, l);
     
     return l;
   }

//=================================================
// The family of find functions:
//=================================================
// size_t find ( const string& str, size_t pos = 0 ) const;
// size_t find ( const char* s, size_t pos, size_t n ) const;
// size_t find ( const char* s, size_t pos = 0 ) const;
// size_t find ( char c, size_t pos = 0 ) const;

   size_t find ( const faststring& str, size_t pos = 0 ) const
   {
     const char *p = find_fundamental(begin()+pos, end(), str.begin(), str.end() );
     
     if (p != NULL)
       return p-begin();
     else
       return npos;
   }

   size_t find ( const char* s, size_t pos, size_t n ) const
   {
     const char *p = find_fundamental(begin()+pos, end(), s, s+n );
     
     if (p != NULL)
       return p-begin();
     else
       return npos;
   }

   size_t find ( const char* s, size_t pos=0) const
   {
     const char *p = find_fundamental(begin()+pos, end(), s, s+strlen(s) );
     
     if (p != NULL)
       return p-begin();
     else
       return npos;
   }

   // There is one good reason determining the position and not the pointer:
   // - The position is still valid after a realloc, the pointer is not.
   size_t find(char c, size_t pos=0) const
   {
     char *b = begin()+pos;
     char *e = end();
     
     while (b < e && *b != c)
       ++b;
     if (b == e)
       return npos;
     else
       return b-begin();
   }

/*   size_t find_forward(char c, size_t pos=0) const */
/*   { */
/*     char *b = begin()+pos; */
/*     char *e = end(); */

/*     while (b < e && *b != c) */
/*       ++b; */
/*     if (b >= e) */
/*       return -1u; */
/*     else */
/*       return b-begin(); */
/*   } */


//=================================================
// The family of rfind functions:
//=================================================
// size_t rfind ( const string& str, size_t pos = npos ) const;
// size_t rfind ( const char* s, size_t pos, size_t n ) const;
// size_t rfind ( const char* s, size_t pos = npos ) const;
// size_t rfind ( char c, size_t pos = npos ) const;

   size_t rfind ( const faststring& str, size_t pos = npos ) const
   {
     if (pos >= _len)
       pos = _len-1;

     const char *p = rfind_fundamental(begin()+pos,rend(),str.rbegin(),str.rend());
     
     if (p != NULL)
       return p-begin();
     else
       return npos;
   }

   size_t rfind ( const char* s, size_t pos, size_t n ) const
   {
     if (pos >= _len)
       pos = _len-1;

     const char *p = rfind_fundamental(begin()+pos, rend(), s+n-1, s-1 );

     if (p != NULL)
       return p-begin();
     else
       return npos;
   }

   size_t rfind ( const char* s, size_t pos=npos) const
   {
     if (pos >= _len)
       pos = _len-1;

     const char *p = rfind_fundamental(begin()+pos, rend(), s+strlen(s)-1, s-1 );
     
     if (p != NULL)
       return p-begin();
     else
       return npos;
   }

  size_t rfind(char c, size_t pos=npos) const
  {
    if  (pos >= _len)
      pos = _len-1;

    char *b = begin()+pos;
    char *e = begin()-1;

    while (b > e && *b != c)
      --b;
    if (b == e) // should be == if nothing was found
      return npos;
    else       // should be b > e if char was found
      return b-begin();
  }


//=================================================
// The family of find_first_of functions:
//=================================================
/* size_t find_first_of ( const string& str, size_t pos = 0 ) const; */
/* size_t find_first_of ( const char* s, size_t pos, size_t n ) const; */
/* size_t find_first_of ( const char* s, size_t pos = 0 ) const; */
/* size_t find_first_of ( char c, size_t pos = 0 ) const; */

/* Consider the strcspn function here */

  size_t find_first_of( const faststring& str, size_t pos = 0 ) const
  {
    const char *b = begin()+pos;
    const char *e = end();

    while (b < e)
    {
      if (is_in_symbol_list__(*b, str.begin(), str.length() ))
	return b-begin();
      ++b;
    }
    return npos;
  }

  size_t find_first_of( const char *str, size_t pos, size_t n ) const
  {
    const char *b = begin()+pos;
    const char *e = end();

    while (b < e)
    {
      if (is_in_symbol_list__(*b, str, n))
	return b-begin();
      ++b;
    }
    return npos;
  }

  size_t find_first_of( const char *str, size_t pos = 0 ) const
  {
    const char *b = begin()+pos;
    const char *e = end();

    while (b < e)
    {
      if (is_in_symbol_list__(*b, str) )
	return b-begin();
      ++b;
    }
    return npos;
  }

 size_t find_first_of ( char c, size_t pos = 0 ) const
 {
   return find(c, pos);
 }

//=================================================
// The family of find_last_of functions:
//=================================================
/* size_t find_last_of ( const string& str, size_t pos = npos ) const; */
/* size_t find_last_of ( const char* s, size_t pos, size_t n ) const; */
/* size_t find_last_of ( const char* s, size_t pos = npos ) const; */
/* size_t find_last_of ( char c, size_t pos = npos ) const; */

 size_t find_last_of( const faststring& str, size_t pos = npos ) const
 {
   if (pos > _len)
     pos = _len;

   const char *e = begin()-1;
   const char *b = begin()+pos-1;
   
   while (b > e)
   {
     if (is_in_symbol_list__(*b, str.begin(), str.length() ) )
       return b-begin();
     --b;
   }
   return npos;
 }


 size_t find_last_of( const char *str, size_t pos, size_t n ) const
 {
   if (pos > _len)
     pos = _len;

   const char *e = begin()-1;
   const char *b = begin()+pos-1;

   while (b > e)
   {
     if (is_in_symbol_list__(*b, str, n))
	return b-begin();
     --b;
   }
   return npos;
 }


 size_t find_last_of( const char *str, size_t pos = npos ) const
 {
   if (pos > _len)
     pos = _len;

   const char *e = begin()-1;
   const char *b = begin()+pos-1;
   
   while (b > e)
   {
     if (is_in_symbol_list__(*b, str))
       return b-begin();
     --b;
   }
   return npos;
 }
 
 size_t find_last_of ( char c, size_t pos = npos ) const
 {
   return rfind(c, pos);
 }

//=================================================
// The family of find_first_not_of functions:
//=================================================
/* size_t find_first_not_of ( const string& str, size_t pos = 0 ) const; */
/* size_t find_first_not_of ( const char* s, size_t pos, size_t n ) const; */
/* size_t find_first_not_of ( const char* s, size_t pos = 0 ) const; */
/* size_t find_first_not_of ( char c, size_t pos = 0 ) const; */

/* Also consider the strspn function for the implementation !!!!!!!!! */

 size_t find_first_not_of( const faststring& str, size_t pos = 0 ) const
 {
   const char *b = begin()+pos;
   const char *e = end();
   
   while (b < e)
   {
     if (!is_in_symbol_list__(*b, str.begin(), str.length() ))
       return b-begin();
      ++b;
   }
   return npos;
 }

  size_t find_first_not_of( const char *str, size_t pos, size_t n ) const
  {
    const char *b = begin()+pos;
    const char *e = end();

    while (b < e)
    {
      if (!is_in_symbol_list__(*b, str, n))
	return b-begin();
      ++b;
    }
    return npos;
  }

  size_t find_first_not_of( const char *str, size_t pos = 0 ) const
  {
    const char *b = begin()+pos;
    const char *e = end();

    while (b < e)
    {
      if (!is_in_symbol_list__(*b, str))
	return b-begin();
      ++b;
    }
    return npos;
  }

 size_t find_first_not_of ( char c, size_t pos = 0 ) const
 {
   const char *b = begin()+pos;
   const char *e = end();

   while (b < e)
   {
     if (*b != c)
       return b-begin();
      ++b;
   }
   return npos;
 }

//=================================================
// The family of find_last_not_of functions:
//=================================================
/* size_t find_last_not_of ( const string& str, size_t pos = npos ) const; */
/* size_t find_last_not_of ( const char* s, size_t pos, size_t n ) const; */
/* size_t find_last_not_of ( const char* s, size_t pos = npos ) const; */
/* size_t find_last_not_of ( char c, size_t pos = npos ) const; */

 size_t find_last_not_of( const faststring& str, size_t pos = npos ) const
 {
   if (pos > _len)
     pos = _len;

   const char *e = begin()-1;
   const char *b = begin()+pos-1;
   
   while (b > e)
   {
     if (!is_in_symbol_list__(*b, str.begin(), str.length() ))
       return b-begin();
     --b;
   }
   return npos;
 }


 size_t find_last_not_of( const char *str, size_t pos, size_t n ) const
 {
   if (pos > _len)
     pos = _len;

   const char *e = begin()-1;
   const char *b = begin()+pos-1;

   while (b > e)
   {
     if (!is_in_symbol_list__(*b, str, n))
	return b-begin();
     --b;
   }
   return npos;
 }


 size_t find_last_not_of( const char *str, size_t pos = npos ) const
 {
   if (pos > _len)
     pos = _len;

   const char *e = begin()-1;
   const char *b = begin()+pos-1;
   
   while (b > e)
   {
     if (!is_in_symbol_list__(*b, str))
       return b-begin();
     --b;
   }
   return npos;
 }
 
 size_t find_last_not_of ( char c, size_t pos = npos ) const
 {
   if (pos > _len)
     pos = _len;

   const char *e = begin()-1;
   const char *b = begin()+pos-1;

   while (b > e)
   {
     if (*b != c)
       return b-begin();
      --b;
   }
   return npos;
 }

 void collapse_mono_symbol_repeats(char c)
 {
   faststring tmp;

   const char *p1, *p2;
   p2 = p1 = begin();
   ++p2;

   while (p2 != end() )
   {
     if (*p1 != *p2 || *p1 != c)
       tmp.push_back(*p1);
     p1 = p2;
     ++p2;
   }
   // When leaving the loop there is one character remaining that we need to append to tmp:
   tmp.push_back(*p1);
   swap(tmp);
 }

 void unquote()
 {
   removeSpacesBack();
   removeSpacesFront();
   if (*begin() == '\"' && *(end()-1) == '\"' )
   {
     this->erase_front(1);
     this->erase_back(1);

/*      this->erase(begin()); */
/*      this->erase(end()-1); */
   }
   else if (*begin() == '\'' && *(end()-1) == '\'' )
   {
/*        this->erase(begin()); */
/*        this->erase(end()-1); */

     this->erase_front(1);
     this->erase_back(1);
   }
 }

 // Unquote string with first char s and last char e
 void unquote(char s, char e)
 {
   removeSpacesBack();
   removeSpacesFront();
   if (*begin() == s && *(end()-1) == e )
   {
     this->erase_front(1);
     this->erase_back(1);
   }
 }



 bool is_CAPITALIZED_abbreviation_of(const faststring &keyword) const
 {
   size_t i   = 0;
   size_t max = size();
   
   // max cannot be greater than keyword.size()
   if ( max > keyword.size() )
     return false;
   
   // We start comparison in upper case region. Lets see how far we can go.
   while (i < max && std::toupper( _buf[i] ) == keyword[i] )
   {
     ++i;
   }
   
   // First we check whether both strings have the same length and we matched all chars.
   // This means that we matched exactly the upper case region of keyword, which is the
   // complete keyword.
   if ( i == max && max == keyword.size() )
     return true;
   
   // If we get here, i must be smaller than keyword.size() since we can only have
   // * i == max and i == keyword.size(): Not possible after previous if 
   // * i <  max and i == keyword.size(): Not possible since max <= keyword.size()
   // * i == max and i < keyword.size(): 
   // * i <  max and i < keyword.size(): 
   
   // Thus, keyword[i] is valid and we use it to check whether we are in the upper case region
   // (were all !alpha() should also be considered as upper case chars in this context!!). 
   // We simply check whether we are in the !islower() region.
   
   if ( !islower( keyword[i]) )    // We are still in non-lower case region
     return false;
   
   // If we get here, we are in the in the lower-case region
   // There are two possibilities:
   // - i == max, 
   if (i==max)
     return true;
   // - our string is longer than the upper case region, so we have to continue
   //   our comparison.
   
   // Let us move through lower case region
   while (i < max && ((*this)[i] == keyword[i] || std::tolower(_buf[i]) == keyword[i]) )
     ++i;
   
   // Could we match all letters in str?
   if (i==max)
     return true;
   else
     return false;
 }

 // get the next token and remove it from the called object.
 faststring& get_next_token(faststring &str, const char * ws = "\r\n\t\v\f ")
 {
   const size_t  n = size();
         size_t  i = 0;

     str.clear();

     // eat leading whitespace/delimiter -- this should be more efficient than
     // calling removeSpacesFront since we only need to erase in the front of the string once. 
     while ( (i < n) && ( is_in_symbol_list__( (*this)[i], ws)) ) ++i;

     if (i == n)   // nothing left so re remove all spaces and the string will be empty
     {
       erase();
     }
     else
     {
       // find end of word
       size_t  j = i+1;
       while ( (j < n) && (!is_in_symbol_list__( (*this)[j], ws)) ) ++j;

       str = substr(i, j);   // the next token
       erase(0, j);                 // remove token from called object 
     }

     return *this;
 }

 void set_to_string_with_values(std::vector<double> v, faststring delim= " ", int precision=14)
 {
   int i, n=v.size();
   clear();
   
   if (n == 0)
     return;
   
   append(faststring(v[0], precision));
   removeSpacesBack();

   for (i=1; i<n; ++i)
   {
     append(delim);
     append(faststring(v[i], precision));
     removeSpacesBack();
   }
 }

 void set_to_string_with_values(std::vector<int> v, faststring delim= " ")
 {
   int i, n=v.size();
   clear();
   
   if (n == 0)
     return;
   
   append(faststring(v[0]));
   removeSpacesBack();

   for (i=1; i<n; ++i)
   {
     append(delim);
     append(faststring(v[i]));
     removeSpacesBack();
   }
 }


 // Could be made more efficient.
 // A delimiter could be given as parameter.
 int get_vector_of_values(std::vector<double> &v)
 {
   int i, n;
   int bad_values=0;
   v.clear();
   
   std::vector<faststring> vs;
   split(vs, *this);
   n = vs.size();

   for(i=0; i<n; ++i)
   {
     if (vs[i].isADouble() )
       v.push_back(vs[i].ToDouble());
     else
       ++bad_values;
   }
   return bad_values;
 }

 // Funktioniert noch nicht richtig:
/*  void search_replace_first(faststring s1, faststring s2, unsigned pos=0) */
/*  { */
/*    unsigned find_pos = find(s1, pos); */

/*    std::cerr << "find_pos: " << find_pos << std::endl; */
/*    std::cerr << "s1.length(): " <<  s1.length()  << std::endl; */


/*    if (find_pos != npos) */
/*      replace(find_pos, s1.length(), s2 ); */
/*  } */

 void search_replace_all(faststring s1, faststring s2)
 {
   std::vector<faststring> v;
   split_at(v, *this, s1);

   size_t i, N=v.size();
   faststring res;
   
   if (N==1)
     return;  // Just for efficiency. Otherwise we would copy the string v[0] back to this.

   for (i=0; i<(N-1); ++i)
   {
     res += v[i];
     res += s2;
   }
   res += v[i];

   *this = res;
 }



};

inline std::ostream &operator<<(std::ostream &out, faststring  &str)
{
  out << str.c_str();
  return out;
}

inline std::ostream &operator<<(std::ostream &out, const faststring  &str)
{
  faststring fs = str;
  out << fs.c_str();
  return out;
}

inline faststring operator+(faststring a, faststring b)
{
  a += b;
  return a;
}

inline std::istream &getline(std::istream &is, faststring &str)
{
  int c;

  str.clear();

  c = is.get();
  if (  c == '\r' )
  {
    c = is.get();
    if ( c != '\n' )
    {
      is.putback(c);
    }
    c = '\n';
  }
  if (c == '\n')
    return is;

  while (c != EOF)
  {
    str.push_back(c);
    c = is.get();
    if (  c == '\r' )
    {
      c = is.get();
      if ( c != '\n' )
      {
	is.putback(c);
      }
      c = '\n';
      return is;
    }
    if (c == '\n')
    {
      return is;
    }
  }
  return is;
}

/* Can be optimised ?? */
/* Function is declared inline since it is a global function in a header file. It is not declared */
/* static in order not to get compiler warnings due to unused functions. */
//
// A note on efficiency: This is not the most efficient way to read a file.
// Indeed the many push_back calls can lead to many reallocation events for the string buffer.
// TODO: Could be made more efficient by measuring the file and by calling reserve.
inline void readFileIntoString(std::ifstream &is, faststring &str)
{
  int c;

  str.clear();

  c = is.get();
  while (true)
  {
    if (  c == '\r' )
    {
      str.push_back('\n');
      c = is.get();
      if ( c == '\n')
      {
	c = is.get();
      }
    }
    else if (c == '\n')
    {
      str.push_back('\n');
      c = is.get();
    }
    else
    {
      str.push_back(c);
      c = is.get();
    }
    if (c == EOF)
      return;
  }
}

//***********************************************************************
// Reads the content of a file and APPENDs it to he vector.
// NOTE: A vector<faststring> can be extemely inefficient.
//       When expanding the size of the vector, the vector class might
//       destruct and reconstruct all items, when reallocating them.
//***********************************************************************
inline void appendFileToVector(std::ifstream &is, std::vector<faststring> &v)
{
  faststring line;

  getline(is, line);
  while (is)
  {
    v.push_back(line);
    getline(is, line);
  }
}

inline void appendFileToList(std::ifstream &is, std::list<faststring> &l)
{
  faststring line;

  getline(is, line);
  while (is)
  {
    l.push_back(line);
    getline(is, line);
  }
}




  
//*************************************************************
// Tokenizer
//*************************************************************

// Note on efficiency:
// vector<faststring> can be inefficient as a container in the split functions.
// Expanding the vector leads not only to a relocation of the memory but also to
// the construction of copies and destruction of all existing elements.
// For many splitting events, a list should be favoured.


template <typename Container>
int
split (Container &l, const faststring &s, char const * const ws)
{
  l.clear();
    const size_t  S = s.size();
          size_t  i = 0;

    while (i < S) {
        // eat leading whitespace
        while ((i < S) && (is_in_symbol_list__(s[i],ws)))  ++i;
        if (i == S)  return l.size();  // nothing left but WS

        // find end of word
        size_t  j = i+1;
        while ((j < S) && (!is_in_symbol_list__(s[j],ws)))  ++j;

        // add word
        l.push_back(s.substr(i,j-i));

        // set up for next loop
        i = j+1;
    }
    return l.size();
}


template <typename Container>
int
split_at (Container &l, const faststring &s, const faststring &delim)
{
  size_t pos1=0;
  size_t pos2=0;

  while (pos2 != faststring::npos)
  {
    pos2 = s.find(delim, pos2); // Start of next delimiter

    if (pos2 >= s.size() )
      break;
    l.push_back(s.substr(pos1, pos2-pos1) );
    pos1 = pos2+delim.size(); // Start of next token
    pos2 = pos1;              // Start for next search
  }
  l.push_back(s.substr(pos1, pos2) );
  return l.size();
}


//
// Tokenizer which does not split braces.
//

inline void
skip_until_ws_respect(const faststring &s, size_t &pos, char const * const ws = "\r\n\t\v\f ")
{
  char c;
  
  unsigned count_parent = 0;
  unsigned count_brace  = 0;
  unsigned count_curly  = 0;
  unsigned count_sq     = 0;
  unsigned count_dq     = 0;

  size_t size = s.size();
      
  while (pos < size)
  {
    c = s[pos];
    if (c == '(')   { ++count_parent; }
    else if (c == ')')   { --count_parent; }
    else if (c == '[')   { ++count_brace;  }
    else if (c == ']')   { --count_brace;  }
    else if (c == '{')   { ++count_curly;  }
    else if (c == '}')   { --count_curly;  }
    else if (c == '\'')  { count_sq = (count_sq+1)%2;   }
    else if (c == '\"')  { count_dq = (count_dq+1)%2;   }
    
    if (is_in_symbol_list__(c, ws) && count_parent==0 && count_brace==0 && count_curly==0 && count_sq==0 && count_dq==0)
      return;
    
    ++pos;
  }
}

template <typename Container>
int
split_respect (Container &l, const faststring &s, char const * const ws)

{
  l.clear();
    const size_t  S = s.size();
          size_t  i = 0;

    while (i < S) {
        // eat leading whitespace
        while ((i < S) && (is_in_symbol_list__(s[i],ws)))  ++i;
        if (i == S)  return l.size();  // nothing left but WS

        // find end of word
        size_t  j = i;  // We can't skip the first char which can be a brace.
	skip_until_ws_respect(s, j, ws);

        // add word
        l.push_back(s.substr(i,j-i));

        // set up for next loop
        i = j+1;
    }
    return l.size();
}

// Multiple successive delimiter result in multiple hits.
template <typename Container>
int
split_strict (Container &l, const faststring &s, char const * const ws)
{
  l.clear();
  const size_t  S = s.size();
        size_t  i = 0;

  while (i < S) {
    // find end of word
    size_t  j = i;
    while ((j < S) && (!is_in_symbol_list__(s[j],ws)))  ++j;

    // add word
    l.push_back(s.substr(i,j-i));

    // set up for next loop
    i = j+1;

    if (i==S)  // Only true of last string is empty since otherwise i==S+1
    {
      l.push_back("");
    }
  }
  return l.size();
}

template <typename Container>
int
split_strict_respect (Container &l, const faststring &s, char const * const ws)
{
  l.clear();
  const size_t  S = s.size();
        size_t  i = 0;

  while (i < S) {
    // find end of word
    size_t  j = i;  // We can't skip the first char which can be a brace.
    skip_until_ws_respect(s, j, ws);

    // add word
    l.push_back(s.substr(i,j-i));

    // set up for next loop
    i = j+1;

    if (i==S)  // Only true of last string is empty since otherwise i==S+1
    {
      l.push_back("");
    }
  }
  return l.size();
}



inline bool less_than_pointer_to_faststring(const faststring *a, const faststring *b)
{
  return *a < *b;
}

struct less_than_pointer_to_faststring_struct
{
  bool operator()(const faststring* s1, const faststring* s2) const
  {
    return *s1 < *s2;
  }
};











#endif

/* String operations: */
/* I do not know what to do with this function: get_allocator	 Get allocator (public member function) */




