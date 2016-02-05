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

#ifndef EASYSTRING_H
#define EASYSTRING_H

#include <iostream>
#include <string>
#include <cstring>    // for strchr
#include <cstdlib>
#include <fstream>

// A usefull global function split:
//(Code copied from the stringtok.h file on the gnu libstdc++ page.



namespace {
    inline bool
    isws(char c, const char* wstr="\r\n\t\v\f ")
    {
        return (strchr(wstr,c) != NULL);
    }
}

namespace {
  //
  // Reads through a string. If braces are encountered we do not tokenize inside but try to read them.
  //
  inline void
    skip_until_ws_respect(const std::string &s, std::string::size_type &pos, char const * const ws = "\r\n\t\v\f ")
    {
      char c;

      unsigned count_parent = 0;
      unsigned count_brace  = 0;
      unsigned count_curly  = 0;
      unsigned count_sq     = 0;
      unsigned count_dq     = 0;

      std::string::size_type size = s.size();

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

	if (isws(c, ws) && count_parent==0 && count_brace==0 && count_curly==0 && count_sq==0 && count_dq==0)
	  return;

	++pos;
      }
    }
}

namespace std
{

/*****************************************************************
 * Simplistic and quite Standard, but a bit slow.  This should be
 * templatized on basic_string instead, or on a more generic StringT
 * that just happens to support ::size_type, .substr(), and so on.
 * I had hoped that "whitespace" would be a trait, but it isn't, so
 * the user must supply it.  Enh, this lets them break up strings on
 * different things easier than traits would anyhow.
*/
template <typename Container>
int
split (Container &l, string const &s, char const * const ws = "\r\n\t\v\f ")
{
  l.clear();
    const string::size_type  S = s.size();
          string::size_type  i = 0;

    while (i < S) {
        // eat leading whitespace
        while ((i < S) && (isws(s[i],ws)))  ++i;
        if (i == S)  return l.size();  // nothing left but WS

        // find end of word
        string::size_type  j = i+1;
        while ((j < S) && (!isws(s[j],ws)))  ++j;

        // add word
        l.push_back(s.substr(i,j-i));

        // set up for next loop
        i = j+1;
    }
    return l.size();
}

//
// Tokenizer which does not split braces.
//
template <typename Container>
int
split_respect (Container &l, string const &s, char const * const ws = "\r\n\t\v\f ")
{
  l.clear();
    const string::size_type  S = s.size();
          string::size_type  i = 0;

    while (i < S) {
        // eat leading whitespace
        while ((i < S) && (isws(s[i],ws)))  ++i;
        if (i == S)  return l.size();  // nothing left but WS

        // find end of word
        string::size_type  j = i;  // We can't skip the first char which can be a brace.
	skip_until_ws_respect(s, j, ws);

        // add word
        l.push_back(s.substr(i,j-i));

        // set up for next loop
        i = j+1;
    }
    return l.size();
}

// Multiple successive deliminators result in multiple hits.
template <typename Container>
int
split_strict (Container &l, string const &s, char const * const ws = " \t\n")
{
  l.clear();
  const string::size_type  S = s.size();
        string::size_type  i = 0;

  while (i < S) {
    // find end of word
    string::size_type  j = i;
    while ((j < S) && (!isws(s[j],ws)))  ++j;

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



 class easystring: public string
 {
 public:
   // Default constructor:
   easystring():string(){};

   // Copy constructor:
   easystring(const string& s):string(s){}
   easystring(const string& s, size_type pos, size_type n):string(s,pos,n){}

   // Other constructors:
   easystring(const char* s, size_type n):string(s,n){}
   easystring(size_type n, char c):string(n,c){}
   template <class InputIterator>
     easystring(InputIterator first, InputIterator last):string(first, last){}

   // Type conversion constructors:
   easystring(const char* s):string(s){}
   easystring(const char c):string(1,c){}
   easystring(const int i):string()
   {
     //     std::cout << "Hallo" << std::endl;
     char tmp[21];   // Sufficient for 64 bit numbers + sign
     sprintf(tmp, "%d", i);
     append(tmp);
   }
   easystring(const long i):string()
   {
     //     std::cout << "Hallo" << std::endl;
     char tmp[21];   // Sufficient for 64 bit numbers + sign
     sprintf(tmp, "%ld", i);
     append(tmp);
   }
   easystring(const unsigned i):string()
   {
     //     std::cout << "Hallo" << std::endl;
     char tmp[21];   // Sufficient for 64 bit numbers
     sprintf(tmp, "%u", i);
     append(tmp);
   }
   easystring(const unsigned long i):string()
   {
     //     std::cout << "Hallo" << std::endl;
     char tmp[21];   // Sufficient for 64 bit numbers
     sprintf(tmp, "%lu", i);
     append(tmp);
   }
   easystring(const double x, int pres):string()
   {
     //     std::cout << "Hallo" << std::endl;
     char tmp[25];   // Sufficient for 64 bit numbers
     sprintf(tmp, "%.*f", pres, x);
     append(tmp);
   }


   void                         removeSpacesFront(const char* delims="\r\n\t\v\f ")
   {
           string::size_type  i;
     const string::size_type  n = size();

     for (i=0; i < n && isws((*this)[i], delims); ++i);
     erase(0,i);
   }


   void                         removeSpacesBack(const char* delims="\r\n\t\v\f ")
   {
     unsigned i;
     unsigned n = size();

     if (n > 0)
     {
       for (i=n-1; i != 0 && isws((*this)[i], delims); --i);

       if (i == 0 && isws((*this)[0], delims))
	 erase();
       else
	 erase(i+1);
     }
   }


   void                         ToUpper()
   {
     iterator it, it_end;

     it     = begin();
     it_end = end();

     while (it != it_end)
     {
       *it = toupper(*it);
       ++it;
     }
   }

   void                         ToLower()
   {
     string::iterator it, it_end;

     it     = begin();
     it_end = end();

     while (it != it_end)
     {
       *it = tolower(*it);
       ++it;
     }
   }

   // Convert to unsigned long
   unsigned			ToUnsigned() const
   {
     char *end;
     return strtoul(c_str(), &end, 0);
   }

   // Convert to unsigned long
   unsigned long 		ToUnsignedLong() const
   {
     char *end;
     return strtoul(c_str(), &end, 0);
   }

   // Convert to int
   int				ToInt() const
   {
     char *end;
     return strtol(c_str(), &end, 0);

   }


   // Convert to long
   long				ToLong() const
   {
     char *end;
     return strtol(c_str(), &end, 0);
   }

   long				ToLong(string::size_type  beg_pos,
				       string::size_type& end_pos) const
   {
     const char   *beg = c_str()+beg_pos;
           char   *end;
           long   result;

     result  = strtol(beg, &end, 0);
     end_pos = end-c_str();
     return result;
   }

   unsigned			ToUnsigned(string::size_type  beg_pos,
					   string::size_type& end_pos) const
   {
     const char   *beg = c_str()+beg_pos;
           char   *end;
           long   result;

     result  = strtoul(beg, &end, 0);
     end_pos = end-c_str();
     return result;
   }


   // Convert to double
   double			ToDouble() const
   {
     char *end;
     return strtod(c_str(), &end);
   }

   double			ToDouble(string::size_type  beg_pos,
					 string::size_type& end_pos) const
   {
     const char   *beg = c_str()+beg_pos;
           char   *end;
           double result;

     result  = strtod(beg, &end);
     end_pos = end-c_str();
     return result;
   }

   // not tested
   unsigned countChar(char c)
   {
     unsigned          n = 0;
     string::iterator  b = begin();
     string::iterator  e = end();

     while (b<e)
     {
       if (*b == c)
	 ++n;
       ++b;
     }
     return n;
   }

   // not tested
   void unquote()
   {
     removeSpacesBack();
     removeSpacesFront();
     if (*begin() == '"' && *(end()-1) == '"' )
     {
       this->erase(begin());
       this->erase(end()-1);
     }
     else if (*begin() == '\'' && *(end()-1) == '\'' )
     {
       this->erase(begin());
       this->erase(end()-1);
     }
   }

   // not tested
   void unquote(char s, char e)
   {
     removeSpacesBack();
     removeSpacesFront();
     if (*begin() == s && *(end()-1) == e )
     {
       this->erase(begin());
       this->erase(end()-1);
     }
   }

   // get the next token and remove it from the called object.
   easystring& get_next_token(easystring &str, const char * ws = "\r\n\t\v\f ")
   {
     const string::size_type  n = size();
           string::size_type  i = 0;

     str.clear();

     // eat leading whitespace/deliminators -- this sould be more efficient than
     // calling removeSpacesFront since we only need to erase in the front of the string once. 
     while ( (i < n) && ( isws( (*this)[i], ws)) ) ++i;

     if (i == n)   // nothing left so re remove alle spaces and the string will be empty
     {
       erase();
     }
     else
     {
       // find end of word
       string::size_type  j = i+1;
       while ( (j < n) && (!isws( (*this)[j], ws)) ) ++j;

       str = string(*this, i, j);   // the next token
       erase(0, j);                 // remove token from called object 
     }

     return *this;
   }


   // xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   // Be careful with keywords containing punctuations since they do not
   // count as upper case chars. The upper case loop stops there.
   
   // Possible bug: I think "I" abbreviates "INPUT" which is not OK.
   
   // Use !isalpha() to as equivalent to "isupper() so we prolong the upper region.
   
   bool is_CAPITALIZED_abbreviation_of(const easystring &keyword)
   {
     unsigned i   = 0;
     unsigned max = size();

     // max cannot be greater than keyword.size()
     if ( max > keyword.size() )
       return false;
     
     // We start comparrison in upper case region. Lets see how far we can go.
     while (i < max && toupper( (*this)[i] ) == keyword[i] )
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
     while (i < max && ((*this)[i] == keyword[i] || tolower((*this)[i]) == keyword[i]) )
       ++i;

     // Could we match all letters in str?
     if (i==max)
       return true;
     else
       return false;
   }
   

 };


}


#endif
