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
#ifndef FAST_DYNAMIC_BITSET_H
#define FAST_DYNAMIC_BITSET_H

// 23.7.2014: Since unsigned long has 64 bits on some computers unsinged long cannot be used as
//            the word data type without having to change large parts of the code.

#include <iostream>
#include <fstream>
#include <cstring>
#include <string>

//#define __num_bits_in_word      (( unsigned long)8*sizeof(unsigned long))
//#define __get_word_index(i)     (( unsigned long)i/(__num_bits_in_word))
//#define __get_num_words(size)   (((unsigned long)size+__num_bits_in_word-1ul)/(__num_bits_in_word))
//#define __get_empty_mask(size,numWords)  ( 0xFFFFFFFF >> numWords*__num_bits_in_word - size  )
//#define __get_num_byte(words)   ( words*sizeof(unsigned long))

#define __num_bits_in_word                32
#define __get_word_index(i)               (i >> 5)
#define __get_word_offset(i)              (i & 31)
#define __get_num_words(size)             ((size + 31) >> 5)
#define __get_num_byte(words)             ( words*4 )

// Must be a 32 bit integer type:
//typedef unsigned   bitset_word_type;
typedef unsigned   my_uint_32;

const unsigned char  __num_bits_table[256] =
     {
       0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4, 1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
       1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
       1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
       2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
       1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5, 2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
       2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
       2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6, 3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
       3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7, 4, 5, 5, 6, 5, 6, 6, 7, 5, 6, 6, 7, 6, 7, 7, 8
     };

const my_uint_32 __mask_table[32] =
{
  0x00000001, 0x00000002, 0x00000004, 0x00000008,
  0x00000010, 0x00000020, 0x00000040, 0x00000080,
  0x00000100, 0x00000200, 0x00000400, 0x00000800,
  0x00001000, 0x00002000, 0x00004000, 0x00008000,

  0x00010000, 0x00020000, 0x00040000, 0x00080000,
  0x00100000, 0x00200000, 0x00400000, 0x00800000,
  0x01000000, 0x02000000, 0x04000000, 0x08000000,
  0x10000000, 0x20000000, 0x40000000, 0x80000000
};


const my_uint_32 __empty_mask_table[32] =
{
  0xFFFFFFFF, 0x7FFFFFFF, 0x3FFFFFFF, 0x1FFFFFFF,
  0x0FFFFFFF, 0x07FFFFFF, 0x03FFFFFF, 0x01FFFFFF,
  0x00FFFFFF, 0x007FFFFF, 0x003FFFFF, 0x001FFFFF,
  0x000FFFFF, 0x0007FFFF, 0x0003FFFF, 0x0001FFFF,

  0x0000FFFF, 0x00007FFF, 0x00003FFF, 0x00001FFF,
  0x00000FFF, 0x000007FF, 0x000003FF, 0x000001FF,
  0x000000FF, 0x0000007F, 0x0000003F, 0x0000001F,
  0x0000000F, 0x00000007, 0x00000003, 0x00000001
};


class fast_dynamic_bitset
{
 public:
  class exception {};
  class outOfRange          : public exception {};
  class zeroSizeBitset      : public exception {};
  class UnequalSizedBitsets : public exception {};

 private:
  my_uint_32        *__data;
  my_uint_32         __size;
  my_uint_32         __numWords;
  //  my_uint_32  __empty_mask;

 public:
  fast_dynamic_bitset(my_uint_32 size=0):__size(size),__numWords(__get_num_words(size))
  {
    //    if (size == 0)
    //     throw zeroSizeBitset();

    //    __size     = size;
    __data = new my_uint_32 [__numWords];
    //    __empty_mask = __empty_mask_table[__numWords * __num_bits_in_word - size];
      // 0xFFFFFFFF >> (__numWords * __num_bits_in_word - size);
  }

  fast_dynamic_bitset(const fast_dynamic_bitset &a)
  {
    memcpy((void*)this, (void*) &a, sizeof(fast_dynamic_bitset) );
    __data = new my_uint_32 [__numWords];
    memcpy((void*)__data,(void*)a.__data, __get_num_byte(__numWords) );
  }

  ~fast_dynamic_bitset() {
    delete [] __data;
    //    std::cout << "Dest" << std::endl;
  }

  my_uint_32 size()
  {
    return __size;
  }

  fast_dynamic_bitset& operator=(const fast_dynamic_bitset &a)
  {
    if ( __size != a.__size )
    {
      delete [] __data;
      memcpy((void*)this, (void*) &a, sizeof(fast_dynamic_bitset) );
      __data = new my_uint_32 [__numWords];
      memcpy((void*)__data,(void*)a.__data, __get_num_byte(__numWords) );
    }
    else
    {
      memcpy((void*)__data,(void*)a.__data, __get_num_byte(__numWords) );
    }
    return *this;
  }


  void reset()
  {
    memset((void*)__data, 0, __get_num_byte(__numWords));
  }

  void clear()
  {
    reset();
  }

  void set()
  {
    for(my_uint_32 index=0; index < __numWords-1; ++index)
      __data[index] = ~static_cast<my_uint_32>(0);
    __data[__numWords-1] = __empty_mask_table[__numWords * __num_bits_in_word - __size];
  }

  void flip()
  {
    for(my_uint_32 index=0; index < __numWords; ++index)
      __data[index] = ~__data[index];
    __data[__numWords-1] &= __empty_mask_table[__numWords * __num_bits_in_word - __size];
  }


  void unchecked_set(my_uint_32 i)
  {
    __data[__get_word_index(i)] |= __mask_table[__get_word_offset(i)];
  }

  void set(my_uint_32 i)
  {
    if (i >= __size)
      throw outOfRange();

    __data[__get_word_index(i)] |= __mask_table[__get_word_offset(i)];
  }

  void unchecked_reset(my_uint_32 i)
  {
    __data[__get_word_index(i)] ^= __mask_table[__get_word_offset(i)];
  }

  void reset(my_uint_32 i)
  {
    if (i >= __size)
      throw outOfRange();

    __data[__get_word_index(i)] ^= __mask_table[__get_word_offset(i)];
  }

  void unchecked_set(my_uint_32 i, my_uint_32 x)
  {
    if (x==0)
      unchecked_reset(i);
    else
      unchecked_set(i);
  }

  void set(my_uint_32 i, my_uint_32 x)
  {
    if (x==0)
      reset(i);
    else
      set(i);
  }

  bool unchecked_test(my_uint_32 i) const
  {
    return (__data[__get_word_index(i)] & __mask_table[__get_word_offset(i)]);
  }

  bool test(my_uint_32 i) const
  {
    //    if (i >= __size)
    //      throw outOfRange();

    return (__data[__get_word_index(i)] & __mask_table[__get_word_offset(i)]);
  }

  void unchecked_flip(my_uint_32 i)
  {
    __data[__get_word_index(i)] ^= __mask_table[__get_word_offset(i)];
  }

  void flip(my_uint_32 i)
  {
    //    if (i >= __size)
    //      throw outOfRange();

    __data[__get_word_index(i)] ^= __mask_table[__get_word_offset(i)];
  }

  void operator&=(const fast_dynamic_bitset &a)
  {
    //    if ( a.__numWords != __numWords )
    //      throw UnequalSizedBitsets();

    for(my_uint_32 index=0; index < __numWords; ++index)
      __data[index] &= a.__data[index];
  }

  void operator|=(const fast_dynamic_bitset &a)
  {
    //    if ( a.__numWords != __numWords )
    //     throw UnequalSizedBitsets();

    for(my_uint_32 index=0; index < __numWords; ++index)
      __data[index] |= a.__data[index];
  }

  void operator^=(const fast_dynamic_bitset &a)
  {
    //    if ( a.__numWords != __numWords )
    //      throw UnequalSizedBitsets();

    for(my_uint_32 index=0; index < __numWords; ++index)
      __data[index] ^= a.__data[index];
  }

  fast_dynamic_bitset operator~()
  {
    fast_dynamic_bitset a=*this;
    a.flip();
    return a;
  }

  std::string as_string_highBitsFirst() const
  {
    std::string str;

/*     // Print overhead-empty bit in bitset: */
/*     //------------------------------------- */
/*     my_uint_32 k = __numWords * __num_bits_in_word; */
/*     while (k != __size) */
/*     { */
/*       --k; */
/*       if (unchecked_test(k)) */
/* 	str.push_back('1'); */
/*       else */
/* 	str.push_back('0'); */
/*     } */
/*     //    os << "|"; */

    // Print bitset with high bits first
    //-----------------------------------
    my_uint_32 i = __size;

    while (i != 0)
    {
      --i;
      if (unchecked_test(i))
	str.push_back('1');
      else
	str.push_back('0');
    }
    return str;
  }

  std::string as_string_lowBitsFirst() const
  {
    std::string str;

    // Print bitset with low bits first
    //----------------------------------
    my_uint_32 i = 0;

    while (i < __size )
    {
      if (unchecked_test(i))
	str.push_back('1');
      else
	str.push_back('0');
      ++i;
    }
    return str;
  }


  std::ostream& print_lowBitsFirst(std::ostream& os) const
  {
    // Print bitset with low bits first
    //----------------------------------
    my_uint_32 i = 0;

    while (i < __size )
    {
      if (unchecked_test(i))
	os << "1";
      else
	os << "0";
      ++i;
    }

/*     // Print overhead-empty bit in bitset: */
/*     //------------------------------------- */
/*     os << "|"; */
/*     while ( i < __numWords * __num_bits_in_word) */
/*     { */
/*       if (unchecked_test(i)) */
/* 	os << "1"; */
/*       else */
/* 	os << "0"; */
/*       ++i; */
/*     } */
    return os;
  }

  void binary_out(std::ostream& os)
  {
    os.write((char*)__data,__numWords*sizeof(my_uint_32));
  }

  void binary_in(std::istream& is)
  {
    is.read((char*) __data,__numWords*sizeof(my_uint_32));
  }


  my_uint_32 count() const
  {
    my_uint_32 sum = 0;
    unsigned char *ptr     = (unsigned char*) __data;
    unsigned char *ptr_end = ptr + __get_num_byte(__numWords);

    while ( ptr != ptr_end )
    {
      sum += __num_bits_table[*ptr];
      ++ptr;
    }
    return sum;
  }

  bool any() const
  {
    for(my_uint_32 i=0; i < __numWords; ++i)
    {
      if (__data[i])
	return true;
    }
    return false;
  }

  bool none() const
  {
    return !any();
  }

  void operator++()
  {
    my_uint_32 i=0;

    while ( i < __numWords )
    {
      ++__data[i];
      if (__data[i]) // if __data[i] == 0 we have to carry a bit otherwise we are done.
	break;
      ++i;
    }
  }

  friend bool operator!=(const fast_dynamic_bitset &a, const fast_dynamic_bitset &b)
  {
    //    if ( a.__numWords != b.__numWords )
    //      throw UnequalSizedBitsets();

    for (my_uint_32 i = 0; i < a.__numWords; --i)
    {
      if ( a.__data[i] != b.__data[i] )
	return true;
    }
    return false;
  }

  friend bool operator==(const fast_dynamic_bitset &a, const fast_dynamic_bitset &b)
  {
    return !(a != b);
  }

  friend bool operator<(const fast_dynamic_bitset &a, const fast_dynamic_bitset &b)
  {
    //    if ( a.__numWords != b.__numWords )
    //      throw UnequalSizedBitsets();

    for (my_uint_32 i=a.__numWords-1; i > 0; --i)
    {
      if ( a.__data[i] < b.__data[i] )
	return true;
      else if (a.__data[i] > b.__data[i])
	return false;
    }
    if ( a.__data[0] < b.__data[0] )
      return true;
    else
      return false;
  }

  friend bool operator>(const fast_dynamic_bitset &a, const fast_dynamic_bitset &b)
  {
    //    if ( a.__numWords != b.__numWords )
    //      throw UnequalSizedBitsets();

    for (my_uint_32 i=a.__numWords-1; i > 0; --i)
    {
      if ( a.__data[i] != b.__data[i] )
	return a.__data[i] > b.__data[i];
    }
    return a.__data[0] > b.__data[0];
  }

  friend bool operator<=(const fast_dynamic_bitset &a, const fast_dynamic_bitset &b)
  {
    return !(a>b);
  }

  friend bool operator>=(const fast_dynamic_bitset &a, const fast_dynamic_bitset &b)
  {
    return !(a<b);
  }


};



// Global functions for fast_dynamic_bitsets

inline fast_dynamic_bitset operator&(const fast_dynamic_bitset &a, const fast_dynamic_bitset &b)
{
  fast_dynamic_bitset c(a);
  c &= b;
  return c;
}

inline fast_dynamic_bitset operator|(const fast_dynamic_bitset &a, const fast_dynamic_bitset &b)
{
  fast_dynamic_bitset c(a);
  c |= b;
  return c;
}

inline fast_dynamic_bitset operator^(const fast_dynamic_bitset &a, const fast_dynamic_bitset &b)
{
  fast_dynamic_bitset c(a);
  c ^= b;
  return c;
}

inline  std::ostream& operator<<(std::ostream& os, const fast_dynamic_bitset& b)
{
  b.print_lowBitsFirst(os);
  return os;
}






// To do:
// shift left
// shift right
// readFromString

#endif
