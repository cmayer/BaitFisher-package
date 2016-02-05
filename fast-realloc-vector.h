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
//-----------------------------------------------------------------------------
// Author:    Christoph Mayer
//            Gropiusweg 13
//            44801 Bochum
//
// Copyright: Christoph Mayer
//
// Description:
//      A fast vector class for objects which do not require their constructors
//         or destructors to be called.
//      This class uses malloc, realloc and free to allocate and
//         free dynamic memory for the "vector". Especially the efficient
//         realloc has no C++ equivalent - the major drawback of new and
//         delete.
//      Disadvantage: No constructors and destructors are called for
//         objects in the vector. This limits is applicability.
//         The fastvector-nd.h uses new and delete and does not use
//         memcpy. It is therefore well suited for all classes.
//         Unfortunately, its less efficient.
//



#ifndef FASTVECTOR_H
#define FASTVECTOR_H

#include <iostream>
#include <cstdlib>
#include <cstring>
#include <cctype>
// #include <iomanip>

#define _MINBUFFER_CAPACITY 4

template<typename T>
class fastvector
{
  // Possibly even more efficient:
  //    Use pointer _end, _capacity_end instead of _len, _capacity


 private:
  T           *_buf;
  unsigned    _len;
  unsigned    _capacity;   // Capacity of container, i.e number of bytes availible
                           // for writing. If a capacity is requested by the user,
                           // an additional byte is reserved in order to be able to append
                           // the additional \0 which is appended by some function calls.
                           // Consequently, the capacity reported to user is _capacity-1.
 protected:


 public:
  ~fastvector ()
  {
    if (_buf)
      //      delete [] _buf;
      std::free(_buf);
  }

  fastvector ():_buf(NULL),  _len(0), _capacity(0)
  {
  }

  fastvector (const T *v, unsigned len)
  {
    //    std::cout << "New fastvector 2" << std::endl;
    _len      = len;
    _capacity = len;

    if (len)
    {
      _buf      = (T *) malloc(_capacity*sizeof(T)); 
      //_buf      = new T [_capacity]; 
      memcpy(_buf, v, len*sizeof(T));
    }
    else
    {
      _buf = NULL;
    }
  }

  fastvector (const T *v_begin, const T *v_end)
  {
    //    std::cout << "New fastvector 3" << std::endl;
    _len      = (v_end - v_begin);
    _capacity = _len;

    if (_len)
    {
      _buf      = (T *)malloc(_capacity*sizeof(T)); 
      // _buf      = new T [_capacity]; 
      memcpy(_buf, v_begin, _len*sizeof(T));
    }
    else
    {
      _buf = NULL;
    }
  }

  fastvector (const fastvector &a):_len(a._len),_capacity(a._len)
  {
    if (_len)
    {
      _buf = (T *)malloc(_capacity*sizeof(T));
      memcpy(_buf, a._buf, _len*sizeof(T));
    }
    else
    {
      _buf = NULL;
    }      
  }


  void push_back(const T &c)
  {
    if (_capacity == _len)
    {
      if (_buf == NULL)
      {
	_capacity = _MINBUFFER_CAPACITY;
	_buf = (T *)malloc(_capacity*sizeof(T));
      }
      else
      {
	T *tmp;

	_capacity <<= 1;
	tmp = (T *)realloc(_buf, _capacity*sizeof(T));

	if (tmp == NULL)
	{
	  std::cout << "Program aborded due to failed realloc!" << std::endl; 
	  exit(0);
	}
	_buf = tmp;
      }
    }
    _buf[_len] = c;
    ++_len;
  }


  T *begin() const
  {
    return _buf;
  }


  T *end() const
  {
    return _buf+_len;
  }


  T *rbegin() const
  {
    return _buf + (_len - 1);
  }


  T *rend() const
  {
    return _buf - 1;
  }


  void assign(const T *v, int len)
  {
    _len = len;
    reserve(len);
    memcpy(_buf, v, len*sizeof(T));
  }


  void assign(const T *v_begin, const T *v_end)
  {
    _len = v_end - v_begin;
    reserve(_len);
    memcpy(_buf, v_begin, _len*sizeof(T));
  }

  void assign(const fastvector &a)
  {
    _len = a._len;
    reserve(_len);
    memcpy(_buf, a._buf, _len*sizeof(T));
  }


  void reserve(unsigned s)
  {
    //    std::cout << "reserve" << std::endl;
    if (_capacity < s)
    {
      if (_buf == NULL)
      {
	//_buf = new T [s];
	_buf = (T *)malloc(s*sizeof(T));
      }
      else
      {
	T *tmp;

	tmp = (T *)realloc(_buf, s*sizeof(T));
	if (tmp == NULL)
	{
	  std::cout << "Program aborded due to failed realloc!" << std::endl; 
	  exit(0);
	}
	_buf = tmp;
      }
      _capacity = s;
    }
  }

  bool empty() const
  {
    return _len == 0;
  }


  unsigned size() const
  {
    return _len;
  }


  unsigned capacity() const
  {
    return _capacity; // One char is reserved for the \0 at end of string
  }


  void clear()
  {
    //        std::cout << "clear" << std::endl;
    _len = 0;
  }


  void reset()
  {
    //    std::cout << "reset" << std::endl;
    if (_buf)
    {
      // delete [] _buf;
      free(_buf);
      _buf = NULL;
      _capacity = 0;
      _len = 0;
    }
  }


  bool check_pos(unsigned pos)
  {
    return (pos < _len);
  }


  void set_unckecked(unsigned pos, const T &c)
  {
    _buf[pos] = c;
  }


  bool set(unsigned pos, const T &c)
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


  void get_unckecked(unsigned pos, T &c)
  {
    c = _buf[pos];
  }


  T get_unckecked(unsigned pos)
  {
    return _buf[pos];
  }



  bool get(unsigned pos, T &c)
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

  fastvector &operator=(const fastvector &a)
    {
      _len = a._len;
      if (_len)
      {
	reserve(_len);
	memcpy(_buf, a._buf, _len*sizeof(T));
      }
      return *this;
    }

  void swap(fastvector &a)
  {
    T         *_tmpbuf      = _buf;             _buf = a._buf;           a._buf = _tmpbuf; 
    unsigned   _tmplen      = _len;             _len = a._len;           a._len = _tmplen; 
    unsigned   _tmpcapacity = _capacity;   _capacity = a._capacity; a._capacity = _tmpcapacity;     
  }

};


#endif
