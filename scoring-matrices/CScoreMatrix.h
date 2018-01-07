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
#ifndef CSCOREMATRIX_H

#define CSCOREMATRIX_H


#include <iostream>
#include <fstream>
#include "../faststring2.h"
#include <climits>
#include <sstream>

class CScoreMatrix
{
 private:
  faststring     alphabet;
  unsigned char  alphabet_to_index[256];

  unsigned char num_sym;

  int gap_open;
  int gap_ext;


  short **score_mat;

  char type;


 public:
/*   int alphabet_index(char base) */
/*   {     */
/*     return 0; */
/*   } */
  
  CScoreMatrix (const char *filename)
  {
    faststring line;
    std::ifstream ifs;

    ifs.open(filename);

    if (!ifs.fail())
    {
      read_from_file(ifs);
      ifs.close();
    }
    else
    {
      std::cerr << "INFORMATION: Could not find the file with alignment parameters: " << filename << std::endl;
      std::cerr << "The following default values will be assumed as alignment parameters:" << std::endl;

      std::stringstream stringst;
      stringst << ";D standard DNA scoring matrix\n 1 45 80 5 6 80 4\n   -0 -3\n * @ 0 1 2 \n ACGTRYMWSKDHVBN\n 0 1 2 3 0 1 0 0 1 2 0 0 0 1 0\n  1\n -2  1\n -2 -2  1\n -2 -2 -2  1\n  1 -2  1 -2  1\n -2  1 -2  1 -2  1 \n 1  1 -2 -2 -2 -2  1\n  1 -2 -2  1 -2 -2 -2  1\n -2  1  1 -2 -2 -2 -2 -2  1\n -2 -2  1  1 -2 -2 -2 -2 -2  1\n  1 -2  1  1  1 -2 -2  1 -2  1  1\n  1  1 -2  1 -2  1  1  1 -2 -2 -2  1\n  1  1  1 -2  1 -2  1 -2  1 -2 -2 -2  1\n -2  1  1  1 -2  1 -2 -2  1  1 -2 -2 -2  1\n 0  0  0  0  0  0  0  0  0  0  0  0  0  0  0\n";
//        stringst << ";D standard DNA scoring matrix\n 1 45 80 5 6 80 4\n  -12 -4\n * @ 0 1 2 \n ACGTRYMWSKDHVBN\n 0 1 2 3 0 1 0 0 1 2 0 0 0 1 0\n  4\n -3  4\n -3 -3  4\n -3 -3 -3  4\n  1 -1  1 -1  1\n -1  1 -1  1 -3  1\n  1  1 -2 -2  0  0  1\n  1 -2 -2  1  0  0  0  1\n -2  1  1 -2  0  0  0  0  1\n -2 -2  1  1  0  0  0  0  0  1\n  1 -2  1  1  1  0  0  1  0  1  1\n  1  1 -2  1  0  1  1  1  0  0  0  1\n  1  1  1 -2  1  0  1  0  1  0  0  0  1\n -2  1  1  1  0  1  0  0  1  1  0  0  0  1\n 1  1  1  1  0  0  0  0  0  0  0  0  0  0  1\n";
      read_from_file(stringst);

      std::ios_base::fmtflags f = std::cerr.flags();
      std::cerr.flags(std::ios::right);
      std::cerr << "Gap openening penalty:     ";      std::cerr.width(3);      std::cerr << gap_open << std::endl;
      std::cerr << "Gap extension penalty:     ";      std::cerr.width(3);      std::cerr << gap_ext  << std::endl;
      std::cerr << "Match score:               ";      std::cerr.width(3);      std::cerr  << score_compare_symbols('A', 'A') << std::endl;
      std::cerr << "Mismatch score:            ";      std::cerr.width(3);      std::cerr  << score_compare_symbols('A', 'C') << std::endl;
      std::cerr << "Score against N:           ";      std::cerr.width(3);      std::cerr  << score_compare_symbols('A', 'N') << std::endl;
      std::cerr << "Scores against other ambigs: match or mismatch score" << std::endl;
      std::cerr.setf(f);
    }
  }

  ~CScoreMatrix()
  {
    delete [] score_mat;
  }


  int get_Gap_open_score()
  {
    return gap_open;
  }

  int get_Gap_ext_score()
  {
    return gap_ext;
  }

  double get_mean_matchscore_of_first_N_symbols(unsigned N)
  {
    unsigned i;
    double   sum=0;
    for (i=0; i<N; ++i)
    {
      sum += score_mat[i][i];
    }
    return sum/N;
  }

  char get_Type()
  {
    return type;
  }
  
  int get_NumbrOfSymbols()
  {
    return num_sym;
  }

  // Expects an existing file that is openend and reading pointer
  // pointing to the beginning of the file.
  void read_from_file(std::istream &is)
  {
    faststring line;

    // First line
    getline(is, line);
    if (line[0] != ';')
    {
      // Error - must be a ';'
    }
    if (toupper(line[1]) == 'P') // protein
    {
      type = 'P';
    }
    else if (toupper(line[1]) == 'D') // DNA  
    {
      type = 'D';
    }
    else
    {
      // Error, must be either 'D' or 'P'

    }

    // Second line
    getline(is, line); // Unused data
    
    // Third line
    is >> gap_open >> gap_ext;
    getline(is, line); // Read newline of this line.

    // Fourth line
    getline(is, line); // Unused data

    int i;

    for (i=0; i<256; ++i)
      alphabet_to_index[i] = UCHAR_MAX;


    // Fifth line
    unsigned char c = is.get();
    num_sym = 0;
    while (c!= '\n')
    {
      if (!isspace(c))
      {
	alphabet.push_back(c);
	alphabet_to_index[c] = num_sym;
	++num_sym;
      }
      c = is.get();
    }
    
    // Sixth line:
    getline(is, line); // Unused data - Hashkeys

    // Now we read the matrix:
    // Indices:
    // First  index: row
    // Second index: col

    score_mat = new short* [num_sym];
    
    for (i=0; i< num_sym; ++i)
    {
      score_mat[i] = new short[num_sym];
    }
    
    int num_entries = num_sym*(num_sym+1)/2;

    i = 0;
    int ind1 = 0;
    int ind2 = 0;

    short entry;

    while (i < num_entries)
    {
      is >> entry;
      
      score_mat[ind2][ind1] = entry;
      score_mat[ind1][ind2] = entry;

      if (ind1 == ind2)
      {
	++ind2;
	ind1=0;
      }
      else
      {
	++ind1;
      }
      ++i;
    }

    if (!is)
    {
      std::cerr << "Error while reading the score matrix containing the alignment parameters. The file/content seems to be damaged. Please revert to a working version of this file/content." << std::endl;
      exit(10);
    }
  }

  void debug_out(std::ostream &os)
  {
    os << "Data stored in scoring matrix:" << std::endl;
    os << "Type: " << type << std::endl;
    os << "number of symbols: " << (int) num_sym << std::endl;
    os << "Alphabet: " << std::endl;

    int i = 0;
    while (i<num_sym)
    {
      os << i << " " << alphabet.get_unckecked(i) << " " << (int)alphabet_to_index[(unsigned char)alphabet.get_unckecked(i)] << std::endl;
      ++i;
    }
    os << std::endl;
    os << "gap open " << gap_open << " gap_ext " << gap_ext << std::endl;

    int j;

    for (i=0; i< num_sym; ++i) // rows
    {
      for(j=0; j<=i; ++j) // cols
      {
	os << score_mat[i][j] << " ";
      }
      os << std::endl;
    }
  }

  void convert_sym_to_index(unsigned char *beg, unsigned char *end)
  {
    while (beg != end)
    {
      *beg = alphabet_to_index[*beg];
      ++beg;
    }
  }

  void convert_sym_to_index(unsigned char *beg)
  {
    while (*beg)
    {
      *beg = alphabet_to_index[*beg];
      ++beg;
    }
    *beg = UCHAR_MAX; // UCHAR_MAX indicates the end of the string.
  }

  void convert_index_to_sym(unsigned char *beg, unsigned char*end)
  {
    while (beg != end)
    {
      *beg = alphabet.get_unckecked(*beg);
      ++beg;
    }
  }

  void convert_to_index_to_sym(unsigned char *beg)
  {
    while (*beg != UCHAR_MAX)
    {
      *beg = alphabet.get_unckecked(*beg);
      ++beg;
    }
    *beg = 0; // Convert UCHAR_MAX to 0 to indicate the end of the string.
  }



  int score_compare_indizes(unsigned char i, unsigned char j)
  {
    return score_mat[j][i]; // From i to j, col index (i) to row index (j)
  }

  int score_compare_symbols(unsigned char a, unsigned char b)
  {
    return score_mat[alphabet_to_index[b]][alphabet_to_index[a]];
  }
  


};


#endif
