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
#ifndef DNA_ROUTINES_H
#define DNA_ROUTINES_H

/* Todo:  count DNA bases, guess data type, is AA, etc.  */

#include <cstring>
#include <cctype>

static const char DNARNA_symbols[]   = "aAcCgGtTuU";
/////////static const faststring aa_symbols     = "aArRnNdDbBcCeEqQzZgGhHiIlLkKmMfFpPsStTwWyYvVxX";
static const char aa_symbols[]       = "aArRnNdDcCeEqQzZgGhHiIlLkKmMfFpPsStTwWyYvVxX";

static const char aa_symbols_ext[]       = "aArRnNdDcCeEqQzZgGhHiIlLkKmMfFpPsStTwWyYvVxXUuOo*";

static const char non_aa_in_ABC [] = "BJOUZbjouz";

static const char DNARNA_iupac_symbols [] = "aAcCgGtTuURYSWryswKMBDkmbdHVNhvn";
static const char phylip_DNARNA_allowed_symbols [] = "?-nNaAcCgGtTuURYSWryswKMBDkmbdHVhvoOxX";

static const char phylip_neither_aa_dna_symbol [] = "jJzZ";
static const char allowed_non_ABC [] = "-?*";

// O and U are special amino acids that are usually not in the models.
// Chars that are not an aa symbol: B, J, O, U, Z
// X is the amibuity code


inline bool is_DNA_base_upper(char c)
{ 
  if (c == 'A' || c == 'C' || c == 'G' || c == 'T' )
    return true;
  return false;
}

inline bool is_DNA_base_lower(char c)
{ 
  if (c == 'a' || c == 'c' || c == 'g' || c == 't' )
    return true;
  return false;
}

inline bool is_DNA_base(char c)
{
  if (c == 'A' || c == 'C' || c == 'G' || c == 'T' || c == 'a' || c == 'c' || c == 'g' || c == 't' )
    return true;
  return false;
}

inline bool is_DNA_or_RNA_base(char c)
{
  if (c == 'A' || c == 'C' || c == 'G' || c == 'T' || 
      c == 'a' || c == 'c' || c == 'g' || c == 't' || c == 'U' || c == 'u')
    return true;
  return false;
}

inline bool is_allowed_non_ABC(char c)
{
  if (c == '-' || c == '?' || c == '*')
    return true;
  else
    return false;
}

inline bool is_phylip_aa_dna_symbol(char c) // Does not include -? or similar
{
  if (c >='a' && c <= 'z')
  {  
      return true;
  }
  else if (c >='A' && c <= 'Z')
  {
      return true;
  }
  return false;
}


inline bool is_DNA_iupac_ambig(char c)
{
  if (c == 'N'  ||  c == '?'  || c == 'n' ||
      c == 'R'  ||  c == 'Y'  || c == 'S' || c == 'W' ||
      c == 'r'  ||  c == 'y'  || c == 's' || c == 'w' ||
      c == 'K'  ||  c == 'M'  || c == 'B' || c == 'D' ||
      c == 'k'  ||  c == 'm'  || c == 'b' || c == 'd' ||
      c == 'H'  ||  c == 'V'  ||
      c == 'h'  ||  c == 'v')
    return true;
  else
    return false;
}

inline bool is_DNA_iupac_ambig_upper(char c)
{
  if (c == 'N'  ||  c == '?'  ||
      c == 'R'  ||  c == 'Y'  || c == 'S' || c == 'W' ||
      c == 'K'  ||  c == 'M'  || c == 'B' || c == 'D' ||
      c == 'H'  ||  c == 'V')
    return true;
  else
    return false;
}

inline int is_aa_ambig(char c)
{
  if (c == 'b' || c == 'B' || c == '?' || c == 'x' || c == 'X' || 
      c == 'j' || c == 'J' || c == 'z' || c == 'Z')
    return true;
  else
    return false;
}



inline int is_aa_or_aa_ambig(char c) // includes '?' and  B, Z, J, X
{
  if (c >='a' && c <= 'z')
  {
    if (c != 'o' && c != 'u') // ou
      return true;
  }
  else if (c >='A' && c <= 'Z')
  {
    if (c != 'O' && c != 'U') // OU
      return true;
  }

  if (c == '?' || c == '*')
    return true;

  return false;
  //  return (strchr(aa_symbols, c) != NULL);
}

inline int is_aa_or_aa_ambig_extended(char c) // O and U are allowed  // includes '?' and  B, Z, J, X
{
  if (c >='a' && c <= 'z')
  {
      return true;
  }
  else if (c >='A' && c <= 'Z')
  {
      return true;
  }

  if (c == '?' || c == '*')
    return true;

  return false;
  //  return (strchr(aa_symbols, c) != NULL);
}


inline char complementBase_of_upper(char c)
{
  switch (c) {
  case 'A':  c = 'T'; break;
  case 'T':  c = 'A'; break;
  case 'G':  c = 'C'; break;
  case 'C':  c = 'G'; break;
    
  case 'R':  c = 'Y'; break;
  case 'Y':  c = 'R'; break;

  case 'M':  c = 'K'; break;
  case 'K':  c = 'M'; break;

  case 'B':  c = 'V'; break;
  case 'V':  c = 'B'; break;

  case 'D':  c = 'H'; break;
  case 'H':  c = 'D'; break;

  case 'N':  c = 'N'; break;
  case '?':  c = '?'; break;
  case '-':  c = '-'; break;

  default:
    c = '*';
    break;
  }
  return c;
}

inline char complementBase_of_lower(char c)
{
  switch (c) {
  case 'a':  c = 't'; break;
  case 't':  c = 'a'; break;
  case 'g':  c = 'c'; break;
  case 'c':  c = 'g'; break;
    
  case 'r':  c = 'y'; break;
  case 'y':  c = 'r'; break;

  case 'm':  c = 'k'; break;
  case 'k':  c = 'm'; break;

  case 'b':  c = 'v'; break;
  case 'v':  c = 'b'; break;

  case 'd':  c = 'h'; break;
  case 'h':  c = 'd'; break;

  case 'n':  c = 'n'; break;
  case '?':  c = '?'; break;
  case '-':  c = '-'; break;

  default:
    c = '*';
    break;
  }
  return c;
}

inline char complementBase(char c)
{
  if ( isupper(c) )
    return complementBase_of_upper(c);
  else
    return complementBase_of_lower(c);
}


inline void CG_AT_content_in_region(const char *beg, const char *end, unsigned &AT, unsigned &CG)
{
  AT = 0;
  CG = 0;

  char c;

  while (beg != end)
  {
    c = *beg;
    if (c == 'C' || c == 'c' || c == 'G' || c == 'g')
      ++CG;
    else if (c == 'A' || c == 'a' || c == 'T' || c == 't')
      ++AT;

    ++beg;
  }
}

inline void content_in_region_DNA(const char *beg, const char *end, unsigned &AT, unsigned &CG, unsigned &gap, unsigned &ambig, unsigned &unknown)
{
  AT      = 0;
  CG      = 0;
  gap     = 0;
  ambig   = 0;
  unknown = 0;

  char c;

  while (beg != end)
  {
    c = *beg;
    if (c == 'C' || c == 'c' || c == 'G' || c == 'g')
      ++CG;
    else if (c == 'A' || c == 'a' || c == 'T' || c == 't')
      ++AT;
    else if (c == '-')
      ++gap;
    else if (is_DNA_iupac_ambig(c) )
      ++ambig;
    else
      ++unknown;

    ++beg;
  }
}




// Recoding symbols:
// A char can hold 8 bits, enough for 8 characters to be "in or out",
// to be switched on or off. Thus a char can be used to hold
// all bases including ambiguities as follows:
// Say A, C, G, T are the bits 1,2,3,4,
// the ambiguity {A,C} is 1+2 = 3

// A:       65, 97  -> 1
// C:       67, 99  -> 2
// G:       71,103  -> 4
// T:       84,116  -> 8
// R(AG):   82,114  -> 5
// Y(CT):   89,121  -> 10
// S(CG):   83,115  -> 6
// W(TA):   87,119  -> 9
// K(TG):   75,107  -> 12
// M(CA):   77,109  -> 3
// B(CTG):  66,98   -> 14
// D(ATG):  68,100  -> 13
// H(ATC):  72,104  -> 11
// V(ACG):  86,118  -> 7
// N(ACGT): 78,110  -> 15
// -        45      -> 0
// ?        63      -> 15

const char ascii2recode_DNA[128] = {
  0,0,0,0,0,0,0,0,0,0, //  0- 9
  0,0,0,0,0,0,0,0,0,0, // 10-19
  0,0,0,0,0,0,0,0,0,0, // 20-29
  0,0,0,0,0,0,0,0,0,0, // 30-39
  0,0,0,0,0,0,0,0,0,0, // 40-49
  0,0,0,0,0,0,0,0,0,0, // 50-59
  0,0,0,15,0,1,14,2,13,0, // 60-69
  0,4,11,0,0,12,0,3,15,0, // 70-79
  0,0,5,6,8,0,7,9,0,10, // 80-89
  0,0,0,0,0,0,0,1,14,2, // 90-99
  13,0,0,4,11,0,0,0,0,0, //100-109
  0,0,0,0,5,0,8,0,7,9, //110-119
  0,10,0,0,0,0,0,0 //120-127
};

const char recode2ascii_DNA[128] = {
  45,65,67,77,71,82,83,86,84,87, //  0- 9
  89,72,75,68,66,78,0,0,0,0, // 10-19
  0,0,0,0,0,0,0,0,0,0, // 20-29
  0,0,0,0,0,0,0,0,0,0, // 30-39
  0,0,0,0,0,0,0,0,0,0, // 40-49
  0,0,0,0,0,0,0,0,0,0, // 50-59
  0,0,0,0,0,0,0,0,0,0, // 60-69
  0,0,0,0,0,0,0,0,0,0, // 70-79
  0,0,0,0,0,0,0,0,0,0, // 80-89
  0,0,0,0,0,0,0,0,0,0, // 90-99
  0,0,0,0,0,0,0,0,0,0, //100-109
  0,0,0,0,0,0,0,0,0,0, //110-119
  0,0,0,0,0,0,0,0 //120-127
};

const char recodeISDNAamig[128] = {
  0,0,0,1,0,1,1,1,0,1, //  0- 9
  1,1,1,1,1,1,0,0,0,0, // 10-19
  0,0,0,0,0,0,0,0,0,0, // 20-29
  0,0,0,0,0,0,0,0,0,0, // 30-39
  0,0,0,0,0,0,0,0,0,0, // 40-49
  0,0,0,0,0,0,0,0,0,0, // 50-59
  0,0,0,0,0,0,0,0,0,0, // 60-69
  0,0,0,0,0,0,0,0,0,0, // 70-79
  0,0,0,0,0,0,0,0,0,0, // 80-89
  0,0,0,0,0,0,0,0,0,0, // 90-99
  0,0,0,0,0,0,0,0,0,0, //100-109
  0,0,0,0,0,0,0,0,0,0, //110-119
  0,0,0,0,0,0,0,0 //120-127
};

inline unsigned char recode_DNA(unsigned char c)
{
  return ascii2recode_DNA[c];
}

inline unsigned char backrecode_DNA(unsigned char c)
{
  return recode2ascii_DNA[c];
}

inline unsigned char recode_is_DNA_ambig(unsigned char c)
{
  return recodeISDNAamig[c];
}

inline bool is_base_in_recode(unsigned char recode, unsigned char base)
{
  return recode & recode_DNA(base);
}

inline bool is_valid_recode_DNA(unsigned char c)
{
  if (c > 15)
    return false;
  else
    return true;
}

inline bool is_recode_gap(unsigned char c)
{
  return (c == 0);
}



// A	Ala	Alanine         01   65
// R	Arg	Arginine        02   82
// N	Asn	Asparagine      03   78
// D	Asp	Aspartic acid   08   68
// C	Cys	Cysteine        16   67
// Q	Gln	Glutamine       32   81
// E	Glu	Glutamic acid   64   69
// G	Gly	Glycine        128   71
// H	His	Histidine      256   72
// I	Ile	Isoleucine     512   73
// L	Leu	Leucine       1024   76
// K	Lys	Lysine        2048   75 
// M	Met	Methionine    4096   77
// F	Phe	Phenylalanine 8192   70
// P	Pro	Proline      16384   80
// S	Ser	Serine       32768   83
// T	Thr	Threonine    65536   84
// W	Trp	Tryptophan  131072   87
// Y	Tyr	Tyrosine    262144   89
// V	Val	Valine      524288   86
// -        45      -> 0
// ?        63      -> 1048575
// X        88      -> 1048575

const int ascii2recode_aa[128] = {
  0,0,0,0,0,0,0,0,0,0, //  0- 9
  0,0,0,0,0,0,0,0,0,0, // 10-19
  0,0,0,0,0,0,0,0,0,0, // 20-29
  0,0,0,0,0,0,0,0,0,0, // 30-39
  0,0,0,0,0,64,0,0,0,0, // 40-49
  0,0,0,0,0,0,0,0,0,0, // 50-59
  0,0,0,15,0,1,14,2,13,0, // 60-69
  0,4,11,0,0,12,0,3,15,0, // 70-79
  0,0,5,6,8,0,7,9,0,10, // 80-89
  0,0,0,0,0,0,0,1,14,2, // 90-99
  13,0,0,4,11,0,0,0,0,0, //100-109
  0,0,0,0,5,0,8,0,7,9, //110-119
  0,10,0,0,0,0,0,0 //120-127
};


#endif

