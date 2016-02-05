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

#define  NODEBUG
#ifndef CALIGNER_H

#define CALIGNER_H

// Changes:
// 8.2.2012:   Trailing gaps are never free in the matrix. First I thought this
//             could be handled without much effort in the traceback routine, but
//             things seem to be ambiguous then.
//             One solution would be to look for the maximum along the original
//             back boundary of the dynamic programming matrix.

// 18.9.2014:  TODO: Import coords_of_s1_in_al_string
//                    al_perfection
//             from *_experimental.h

#include <iostream>
#include <fstream>
#include <vector>
#include <string>

#include "../faststring2.h"
#include "../scoring-matrices/CScoreMatrix.h"

#include <limits.h>

//const int DEBUG_LEVEL=5;

#define NDEBUG

#define DEBUG_i_min 15700 
#define DEBUG_i_max 15800



#ifdef  DEBUG
#define DEBUGOUT1(x)           std::cerr << x                << std::endl;
#define DEBUGOUT2(x,y)         std::cerr << x << y           << std::endl;
#define DEBUGOUT3(x,y,z)       std::cerr << x << y << z      << std::endl;
#define DEBUGOUT4(x,y,z,w)     std::cerr << x << y << z << w << std::endl;
 
#define DEBUGOUT4Del_OS(os,x,y,z,w)   os << x << "," << y << "," << z << "," << w << std::endl;

#define DEBUGOUT4Del(x,y,z,w)  std::cerr << x << "," << y << "," << z << "," << w << std::endl;
#define DEBUGCODE(x) {x}

#define IF_IN_DEB_RANGE if(i>=DEBUG_i_min && i<=DEBUG_i_max)
#else
#define DEBUGOUT1(x)
#define DEBUGOUT2(x,y)
#define DEBUGOUT3(x,y,z)
#define DEBUGOUT4(x,y,z,w)
#define DEBUGOUT4Del(x,y,z,w)

#define DEBUGOUT4Del_OS(os,x,y,z,w)

#define DEBUGCODE(x)

#define IF_IN_DEB_RANGE
#endif

 
#define RECENTER_AFTER 3
#define U_ENTRY  INT_MIN + 1000000
#define U_ENTRY2 INT_MIN + 3000000

#define macromax2(a,b)     ((a<b)?b:a)
#define macromin2(a,b)     ((a>b)?b:a)
#define macromax3(a,b,c)   (macromax2(c, macromax2(a,b)))

#define next_odd(x) ((x%2==0)?x+1:x)

// #define MATCHSCORE 1
#define N_SCORE    0
#define NOCOMPILE_NO_OFFSET 1

#define DIRECTION -


#define move_forward_in_seq(pos, b, e, s, l)   {pos += s; while (pos >= e) pos -=l;}
#define move_backward_in_seq(pos, b, e, s, l)  {pos += s; while (pos <= b) pos +=l;}

// #define move_forward_in_seq_pp(pos, b, e, l)   ((++pos >= e)?pos: pos -=l)
// #define move_backward_in_seq_mm(pos, b, e, l)  ((--pos <= b)?pos: pos +=l)



inline void move_forward_in_seq_pp(const char *&pos, const char *b, const char *e, int l)
{
  ++pos;
  if (pos >= e)
    pos -= l;
}

inline void debug_cyclic_checker(const char *pos, const char *s, int l, const char *error_msg)
{
  if (pos < s)
  {
    std::cerr << "CYCLIC Lower bound exceeded: " << error_msg << std::endl;
  }

  if (pos >= s+l)
  {
    std::cerr << "CYCLIC Upper bound exceeded: " << error_msg << std::endl;
  }
}








// typedef int scoretype;

using namespace std;

class align_coord
{
 public:
  int col; 
  int row;

  align_coord(int c=0, int r=0):col(c), row(r)
  {}

  void print(ostream &os)
  {
    os << "(" << col << "," << row << ") ";
  }

  bool valid(int mcol, int mrow)
  {
    return (col >= 0 && col < mcol && row >= 0 && row < mrow);
  }
};


inline int max3(int a, int b, int c)
{
  int tmp = macromax2(a,b);
  return macromax2(tmp,c);
}

inline void print_index(ostream &os, int x, int y)
{
  os << "(" << x << "," << y << ") ";
}

class CAligner
{
  // Convention:
  // s1 has to be envisaged above the first row. s1 index corresponds to first index in the matrices.
  // s2 has to be envisaged before the first column. s2 index corresponds to second index in the matrices.

  // The size of the matrix is l1 x l2.
  //

  ostream &log_file;
  int    **ar_s;                // array of scores, first index will be the column (x), the second index the row (y)
  char   **ar_bto;              // array of bytes indicating the backtrace route. 0: Only initialized, bit1: match , bit2: mismatch, bit3: gap in s2 (score2 und p2), bit4: gap in s1 (score3 und p3)
  // Valid coordinates in these two fields:
  // 0..l1, 0..l2,

  //  int          stopscore;
  //  int          last_col_computed;

  CScoreMatrix &score_mat;

  const char*  s1;  // pointer to sequence 1
  int          l1;  // length of sequence  1

  const char*  s2;  // pointer to sequence 2
  int          l2;  // length of sequence  2

  //  int          l_max;

  int          gap_open; 
  int          gap_ext;
  //  int          mis_pen;
  int          alignment_score;

  short        alignment_mode;  // 0: global_alignment, 1: non-global alignment, 2: cyclic alignment

  //  bool         global_alignment;
  bool free_front_gaps;
  bool free_back_gaps;

 void init_border()
 {
   // init the border
   int i;

   if (free_front_gaps)
   {
     for (i=0; i<=l1; ++i)
     {
       // Fill first row, corresponding to gaps in s2
       ar_s[i][0]   = 0;
       ar_bto[i][0] = 4;
     }
     for (i=0; i<=l2; ++i)
     {
       // Fill first columnn, corresponding to gaps in s1
       ar_s[0][i]   = 0;
       ar_bto[0][i] = 8;
     }
     ar_bto[0][0] = 0; // Has no predesessor
   }
   else
   {
     for (i=0; i<=l1; ++i)
     {
       // Fill first column, corresponding to gaps in s2
       ar_s[i][0]   = -(i*gap_ext + gap_open);
       ar_bto[i][0] = 4;
     }
     for (i=0; i<=l2; ++i)
     {
       // Fill first row, corresponding to gaps in s1
       ar_s[0][i]   = -(i*gap_ext + gap_open);
       ar_bto[0][i] = 8;
     }
     ar_bto[0][0] = 0; // Has no predesessor
   }
 }

 public:
 void print_ar_s(ostream &os)
 {
   int i, j;

   os << "Score matrix:" << endl;
   for (j=0; j <= l2; ++j)
   {
     for (i=0; i <= l1; ++i)
     {
       os.width(3);
       if (ar_s[i][j] < U_ENTRY2)
	 os << "-";
       else
	 os << ar_s[i][j];
     }
     os << endl;
   }
 }

 void print_ar_bto(ostream &os)
 {
   int i, j;

   os << "Traceback matrix:" << endl;
   for (j=0; j <= l2; ++j)
   {
     for (i=0; i <= l1; ++i)
     {
       os.width(3);
       if (ar_s[i][j] < U_ENTRY2)
	 os << "-";
       else
	 os << (int) ar_bto[i][j];
     }
     os << endl;
   }
 }
 

CAligner(const char* param_s1, int pl1, const char* param_s2, int pl2,
	 CScoreMatrix &param_score_mat, short param_alignmentmode,
	 ostream &plog_file=cerr, int pstopscore=INT_MIN):
   log_file(plog_file),
   score_mat(param_score_mat),
   s1(param_s1), l1(pl1), s2(param_s2), l2(pl2),

   gap_open( -param_score_mat.get_Gap_open_score()  ),
   gap_ext(  -param_score_mat.get_Gap_ext_score()   ),
   alignment_mode(param_alignmentmode)
 {
   int i,j;
   ar_s   = new int*[l1+1];         // col pointer vector
   ar_bto = new char*[l1+1];

   for (i=0; i <= l1; ++i)         // for all columns
   {
     ar_s[i]   = new int [l2+1];
     ar_bto[i] = new char[l2+1];
   }

   DEBUGCODE(
	     cerr << "Alignment parameters: " << endl;
	     cerr << "gap opening:   " << -gap_open << endl;
	     cerr << "gap extension: " << -gap_ext << endl << endl;
	     )

   // Init alignment array
   for (i=0; i <= l1; ++i)
   {
     for (j=0; j <= l2; ++j)
     { 
       ar_s[i][j]   = U_ENTRY;
       ar_bto[i][j] = 0;
     }
   }


   DEBUGCODE(cerr << "Mode: (intern) " << alignment_mode << endl;)

   if (alignment_mode == 0)      // global
   { 
     free_front_gaps = false;
     free_back_gaps  = false;
     DEBUGCODE(cerr << "Mode: global mode" << endl;)
   }
   else if (alignment_mode==1) // local
   {
     free_front_gaps = true;
     free_back_gaps  = true;
     DEBUGCODE(cerr << "Mode: local mode" << endl;)
   }
   else // cyclic
   {
     free_front_gaps = false;
     free_back_gaps  = true;
     DEBUGCODE(cerr << "Mode: cyclic mode" << endl;)
   }
 }


 ~CAligner()
 {
   int i;

   for (i=0; i<= l1; ++i)
   {
     delete [] ar_s[i];
     delete [] ar_bto[i];
   }
   delete [] ar_s;
   delete [] ar_bto;
 }

 void print_time_stamp_to_log_file(const char * s)
 {
   log_file << "Time: " << s << " " << time(NULL) << endl;
 }

void align()
 {
   int i,j;

   init_border();

   // print_ar_s(std::cerr);
   // print_ar_bto(std::cerr);

   //   std::cerr << "GAP*******" << gap_open << std::endl;


   int score1;
   int score2;
   int score3;
   char match_mismatch;
   int align_score;

   // TODO: Further improve efficienty: Loop over < l2, <l1.
   //       This removes the if that check for the last row and column.
 
   for (j=1; j <= l2; ++j)        // loop over all rows
   {
     for (i = 1; i <= l1; ++i ) // loop over all cols
     {
       // Caution: Both loops start at 1 since row and column 1 have already been filled.
       // However, row and col 1 correspond to entry 0 of the sequence. So sequence indices
       // are shifted by -1:
 
       DEBUGCODE(
		 cerr << "al: i " << i << " j " << j
		    << " s1[i-1]=" << "s1[" << i-1 << "]=" << s1[i-1] << " s2[j-1]=" <<  "s2[" << j-1 << "]= " << s2[j-1]
		    << " ar_s[i-1][j] "<< ar_s[i-1][j] << " ar_s[i-1][j-1] "<< ar_s[i-1][j-1] << " ar_s[i][j-1] "<< ar_s[i][j-1]
		    << endl;
		 )

       // Do we have a match - what is the match/mismatch score
       align_score = score_mat.score_compare_symbols(s1[i-1], s2[j-1]);
       if (s1[i-1] == s2[j-1])
	 match_mismatch = 1;  // set bit 1
       else
	 match_mismatch = 2;  // set bit 2

       {
	 //                       match                    gap                     gap
	 //	   ar_s[i][j] = macromax3( ar_s[i-1][j-1]+match_score, ar_s[i][j-1]-gap_ext, ar_s[i-1][]-gap_ext);

	 // Compute:
	 // score1, score2, score3:
	 // score1: match or mismatch
	 // score2: gap in s2, horizontal move in matrix
	 // score3: gap in s1, vertical move in matrix

	 score1 = ar_s[i-1][j-1]+align_score;

	 if (!free_back_gaps || (i < l1 && j < l2 ) )
	 {
	   score2 = ar_s[i-1][j] - gap_ext;    // gap in s2
	   score3 = ar_s[i][j-1] - gap_ext;    // gap in s1

	   if ( !(ar_bto[i-1][j]&4) )
	   {
	     score2 -= gap_open;
	   }
	   if ( !(ar_bto[i][j-1]&8) )
	   {
	     score3 -= gap_open;
	   }
	 }
	 else // free_back_gaps and either i or j on the back boundary:
	 {
	   if (i==l1 && j==l2)
	   {
	     score2 = ar_s[i-1][j];    // gap in s2
	     score3 = ar_s[i][j-1];    // gap in s1
	   }
	   else if (i==l1)
	   {
	     score2 = ar_s[i-1][j] - gap_ext;    // gap in s2
	     if ( !(ar_bto[i-1][j]&4) )
	     {
	       score2 -= gap_open;
	     }
	     score3 = ar_s[i][j-1] - 0;          // gap in s1
	   }
	   else // if (j==l2) 
	   {
	     score2 = ar_s[i-1][j] - 0;          // gap in s2
	     score3 = ar_s[i][j-1] - gap_ext;    // gap in s1
	     if ( !(ar_bto[i][j-1]&8) )
	     {
	       score3 -= gap_open;
	     }
	   }
	 }

	 DEBUGOUT4Del("al score1 ",  i-1,  j-1, score1);
	 DEBUGOUT4Del("al score2 ",  i-1,  j,   score2);
	 DEBUGOUT4Del("al score3 ",  i,    j-1, score3);

	 DEBUGOUT4Del("2: al ar_bto[i-1][j]",       i-1,  j,   (int)ar_bto[i-1][j]);
	 DEBUGOUT4Del("3: al ar_bto[i][j-1]",       i,    j-1, (int)ar_bto[i][j-1]);

	 if (score1 >= score2 && score1 >= score3)
	 {
	   ar_s[i][j] = score1;
	   ar_bto[i][j] = match_mismatch;  // set bit 1 or 2
	   if (score1 == score2)
	     ar_bto[i][j] += 4;      // set bit 3 in addition to other bit
	   if (score1 == score3)
	     ar_bto[i][j] += 8;     // set bit 4 in addition to what has been set before
	 }
	 else if (score2 >= score3) // gap in s2
	 {
	   ar_s[i][j]   = score2;
	   ar_bto[i][j] = 4;        // set bit 3
	   if (score2 == score3)
	     ar_bto[i][j] = 12;     // set bit 4 in addition to bit 3
	 }
	 else // score3 > score2 and score1 // gap in s1
	 {
	   ar_s[i][j]   = score3;
	   ar_bto[i][j] = 8;   // set bit 4
	 }

	 DEBUGCODE(
		   log_file << "al: scores: ";
		   print_index(log_file, i, j);
		   log_file << score1 << " " << score2 << " " << score3
		   << " " << ar_s[i][j] << " bto: " << (int)ar_bto[i][j]
		   << endl;
		   );
       } // END Block
     } // END loop over all cols
   } // END loop over all rows


/*    int col_maximum = ar_s[i][find_max_index_of_col(i)]; */

/*      if (col_maximum < stopscore) */
/*        //       break; // This means we leave the for (i=1; i < length; ++i) loop */
/*     { */
/*        DEBUGOUT1("Breaking out of for i loop in align_raw():"); */
/*        // In this case we revert to the previous column to start the traceback */
/*        last_col_computed = i-1; */

/*        DEBUGOUT2("raw: last_col_computed reverted to :", last_col_computed); */

/*        break; // This means we leave the for (i=1; i < length; ++i) loop */
/*      } */


//   print_ar_s(std::cerr);
//   print_ar_bto(std::cerr);

 }  // End align()




 int traceback(faststring &alignment, align_coord &ac_start, align_coord &ac_end)
 { 
   int col, row;

   // Let us determine the starting cooridnates for the backtracking procedure. Later this will be stored in ac_end
   col = l1;
   row = l2;

/*    if (!free_back_gaps) // global alignment */
/*    { */
/*      row = l2; */
/*      if (row <= 0 || row >= width) */
/*      { */
/*        cerr << "Critical out of bounds error: first row index in traceback out of bounds" << endl; */
/*        exit(1); */
/*      } */
/*    } */
/*    else */
/*      row = find_max_index_of_col(col); */

   DEBUGOUT2("Starting in col: ", col);
   DEBUGOUT2("Starting in row: ", row);

   ac_end.col = col;
   ac_end.row = row;

   alignment_score = ar_s[col][row];

   align_coord curr_ac(col, row);

   // Initialise last_ac. Otherwise it might be used uninitialized.
   align_coord p1, p2, p3, last_ac = ac_end;
   bool end_fixed_free_end_gaps     = false;
   bool start_fixed_free_front_gaps = false;

   while (curr_ac.col > 0 || curr_ac.row > 0)
   {
/*      if (curr_ac.col == 1) */
/*      { */
/*        cout << "COORD: "; */
/*        curr_ac.print(cout); */
/*        cout << endl; */
/*      } */

     DEBUGCODE(bool match;) 


     DEBUGCODE(
	       cerr << "curr_ac: ";
	       curr_ac.print(cerr);
	       cerr << "Matrix value: " << ar_s[curr_ac.col][curr_ac.row];
	       cerr << endl;
	       )

     p1.col = curr_ac.col - 1;
     p1.row = curr_ac.row - 1;

     p2.col = curr_ac.col-1;
     p2.row = curr_ac.row;     // horizontal (p2, 3rd bit, gap in s2)

     p3.col = curr_ac.col;     // vertical   (p3, 4th bit, gap in s1)
     p3.row = curr_ac.row-1;

     DEBUGCODE(
	       match = (ar_bto[curr_ac.col][curr_ac.row]&1);

	       cerr << "Potential preds: ";
	       cerr << "p1 ";
	       p1.print(cerr);

	       cerr << "p2 ";
	       p2.print(cerr);

	       cerr << " p3 ";
	       p3.print(cerr);

	       cerr << endl;
 
	       if (ar_bto[curr_ac.col][curr_ac.row]&1)
		 cerr << "Match" << endl;
	       if (ar_bto[curr_ac.col][curr_ac.row]&2)
		 cerr << "Mismatch" << endl;
	       if (ar_bto[curr_ac.col][curr_ac.row]&4)
		 cerr << "gap s2" << endl;
	       if (ar_bto[curr_ac.col][curr_ac.row]&8)
		 cerr << "gap s1" << endl;


	       cerr << "Compare:" << endl;
	       cerr << "Answer: " << ar_s[curr_ac.col][curr_ac.row] << endl;
	       if (p1.valid(l1,l2) )
	       {
		 //		 cerr << "   match? " << ar_s[p1.col][p1.row] << "+" << match_score << endl;
		 //		 cerr << "mismatch? " << ar_s[p1.col][p1.row] << "-" << mis_pen << endl;
	       }
	       else
		 cerr << "   p1 is out of bounds. That is OK if we are in the first row or column." << endl;

	       if (p2.valid(l1,l2) )
		 cerr << "      p2? " << ar_s[p2.col][p2.row] << "-" << gap_ext << endl;
	       else
		 cerr << "   p2 is out of bounds. That is OK if we are in the first row or column." << endl;
	       if (p3.valid(l1,l2) )
		 cerr << "      p3? " << ar_s[p3.col][p3.row] << "-" << gap_ext << endl;
	       else
		 cerr << "   p3 is out of bounds. That is OK if we are in the first row or column." << endl; 
	       );


     // New version for finding the traceback

     DEBUGCODE(if (ar_bto[curr_ac.col][curr_ac.row]==0)
	       {
		 cerr << "Error. Tracebackvalue == 0"<< endl;
		 exit(1);
	       }
	       );

     DEBUGOUT2("ar_bto[curr_ac.col][curr_ac.row] = ",char(ar_bto[curr_ac.col][curr_ac.row]+'0'));

     if (ar_bto[curr_ac.col][curr_ac.row]&1) // If this is true, this was a match
     {
       DEBUGOUT1("Choosing match (p1)");
       alignment.push_back('M');
       curr_ac = p1;
     }
     else if (ar_bto[curr_ac.col][curr_ac.row]&2) 
     {
       DEBUGOUT1("Choosing mismatch (p1)");
       alignment.push_back('m');
       curr_ac = p1;
     }
     else if (ar_bto[curr_ac.col][curr_ac.row]&4)
     {
       DEBUGOUT1("Choosing gap s2 (p2)");     // left <- | gap in s2 | 'i' means: relative insertion in s1.
       alignment.push_back('i');
       curr_ac = p2;
     }
     else if (ar_bto[curr_ac.col][curr_ac.row]&8)
     {
       DEBUGOUT1("Choosing gap s1 (p3)");   //  up ^ | gap in s1 | 'd' means: relative deletion in s1.
       alignment.push_back('d');
       curr_ac = p3;
     }
     else
     {
       cerr << "Critical error in traceback_raw_alignment: No succ found at: ";
       curr_ac.print(cerr);
       cerr << endl;
       exit(0);
     }
     // End: New version for finding the traceback

     // In the case of free end gaps we have to find the end of the
     // alignment dynamically. We enter as long as !end_fixed_free_end_gaps

     if (free_back_gaps && !end_fixed_free_end_gaps)
     {
       if (alignment_score != ar_s[curr_ac.col][curr_ac.row] )
       {
	 ac_end = last_ac;
	 end_fixed_free_end_gaps = true;
       } 
     }
     //     else ac_end is l1, l2

     if (free_front_gaps && !start_fixed_free_front_gaps)
     {
       if (curr_ac.col == 0 || curr_ac.row == 0)
       {
	 ac_start = curr_ac;
	 start_fixed_free_front_gaps = true;
       }
     }
     else if (!free_front_gaps)
     {
       ac_start = curr_ac; // Always the last enty has to be stored.
     }
     last_ac  = curr_ac;
  } 
   return alignment_score;
 }


 // Requires alignment mode 1 (local) for expected results.
 void get_respective_start_end_coordinates(int          &s2_in_s1_start,
					   int          &s2_in_s1_end,
					   int          &s1_in_s2_start,
					   int          &s1_in_s2_end,
					   align_coord  a,
					   align_coord  e,
					   int          l1,
					   int          l2,
					   faststring   &alignment)
 {
   if (a.col > a.row)
   {
     s2_in_s1_start =  a.col;
     s1_in_s2_start = -a.col;
   }
   else
   {
     s2_in_s1_start = -a.row;
     s1_in_s2_start =  a.row;
   }

   if (*(alignment.begin()) == 'i')
   {
     s1_in_s2_end = -(l1 - e.col);
     s2_in_s1_end = e.col;
   }
   else if (*(alignment.begin()) == 'd')
   {
     s1_in_s2_end = e.row;
     s2_in_s1_end = -(l2 - e.row);
   }
   else
   {
     s1_in_s2_end = l2;
     s2_in_s1_end = l1;
   }
 }



 void alignment_strings(faststring &line1, faststring &line2, faststring &line3, faststring alignment)
 {
   line1.clear(); 
   line2.clear();
   line3.clear();

   const char *i_al, *i_al_end;

   const char *ps1 = s1;
   const char *ps2 = s2;

   i_al     = alignment.rbegin();
   i_al_end = alignment.rend();

   // Delete:
   //   int num_m_M = 0;
   //   int s1;

   while (i_al != i_al_end)
   {
     if (*i_al == 'M') // Match
     {
       line1.push_back(*ps1);
       ++ps1;
       line3.push_back('|');
       line2.push_back(*ps2);
       ++ps2;
       //       ++num_m_M;
     }
     else if (*i_al == 'm') // mismatch
     {
       line1.push_back(*ps1);
       ++ps1;
       line3.push_back(' ');
       line2.push_back(*ps2);
       ++ps2;
       //       ++num_m_M;
     }
     else if (*i_al == 'i') // gap in s2, or insertion in s1
     {
       line1.push_back(*ps1);
       ++ps1;
       line3.push_back(' ');
       line2.push_back('-');
     }
     else if (*i_al == 'd') // gap in s1, or insertion in s2
     {
       line1.push_back('-');
       line3.push_back(' ');
       line2.push_back(*ps2);
       ++ps2;
     }
     --i_al;
   }
 }

 int get_alignment_score()
 {
   return alignment_score;
 }

  // NEW
 // Result only verified for mode==1
 void coords_of_s1_in_al_string(faststring &alignment,
				align_coord ac_start,
				align_coord ac_end,
				int len1,
				int len2,
				unsigned &s,
				unsigned &e)
 {
   int          s2_in_s1_start;
   int          s2_in_s1_end;
   int          s1_in_s2_start;
   int          s1_in_s2_end;

   get_respective_start_end_coordinates(s2_in_s1_start,
					s2_in_s1_end,
					s1_in_s2_start,
					s1_in_s2_end,
					ac_start,ac_end,
					len1, len2,
					alignment);

   if (s2_in_s1_start < 0)
     s = 0;
   else
     s = s2_in_s1_start;

   if (s1_in_s2_end < 0) // s1 extends beyond s2 in alignment at its end.
     e = alignment.length() + s1_in_s2_end;
   // If s1 does not extend beyond s2 at its end, s2 goes till the end:
   else // s2_in_s1_end < 0 || s2_in_s1_end > 0 && s1_in_s2_end > 0
     // s2_in_s1_end < 0: s2 extends beyond s1 in alignment at its end.
     // s2_in_s1_end > 0 && s1_in_s2_end > 0: end simultanious.
     e = alignment.length();
 }

 // NEW
 // Usually used for mode==1 only! 
 void al_perfection(faststring &alignment,
		    align_coord ac_start,
		    align_coord ac_end,
		    int len1,
		    int len2,
		    unsigned &len,        // alignment covering s1
		    unsigned &matches,    //
		    unsigned &mismatches, //
		    unsigned &gaps)       //
 {
   unsigned s,e; // Start and end in alignment string

   unsigned al_len = alignment.length();

   coords_of_s1_in_al_string(alignment,
			     ac_start,
			     ac_end,
			     len1,
			     len2,
			     s,
			     e);
 
   if (s > e)
   {
     cerr << "Internal error: s>e in method al_perfection is not possible." << endl;
     exit(-33);
   }

   if (e > al_len)
   {
     cerr << "Internal error: e>alignment.length() in method al_perfection is not possible." << endl;
     exit(-33);
   }


   len        = e - s;
   matches    = 0;
   mismatches = 0;
   gaps       = 0;

   unsigned i, N;
   char c;

   i = al_len-e;
   N = al_len-s;

   for (; i<N; ++i)
   {
     c = alignment[i];
     if (c == 'M')
       ++matches;
     else if (c == 'm')
       ++mismatches;
     else
       ++gaps;
   }
 }

};

///////


inline void pretty_print_alignment(FILE *of, faststring &line1, faststring &line2, faststring &line3, int char_per_line=50)
{ 
  //        0     .    :    .    :    .    :    .    :    .    :
  faststring top_str = "    .    :    .    :    .    :    .    :    .    :";

  int i=50;
  while (i<char_per_line)
  {
    i += char_per_line;
    top_str += "    .    :    .    :    .    :    .    :    .    :";
  }

  unsigned len = line1.size();
  unsigned pos = 0;
  unsigned back_pos = pos + char_per_line;

  const char* l1 = line1.c_str();
  const char* l2 = line2.c_str();
  const char* l3 = line3.c_str();

  unsigned n_print_now;
  while (pos < len)
  {
    back_pos = pos + char_per_line;

    if (back_pos > len)
    {
      back_pos    = len;
    }
    n_print_now = (unsigned)(back_pos - pos);

    //    fprintf(of, "debug: %d\n", n_print_now);

    fprintf(of, "%7d %.*s\n",     pos+1, n_print_now, top_str.data());
    fprintf(of, "        %.*s\n", n_print_now, l1+pos);
    fprintf(of, "        %.*s\n", n_print_now, l3+pos);
    fprintf(of, "        %.*s\n\n", n_print_now, l2+pos);

    pos += char_per_line;

  }
}

inline void pretty_print_alignment_base0(FILE *of, faststring &line1, faststring &line2, faststring &line3, int char_per_line=50)
{
  //        0     .    :    .    :    .    :    .    :    .    :
  faststring top_str = "     .    :    .    :    .    :    .    :    .    :";

  int i=50;
  while (i<char_per_line)
  {
    i += char_per_line;
    top_str += "    .    :    .    :    .    :    .    :    .    :";
  }

  unsigned len = line1.size();
  unsigned pos = 0;
  unsigned back_pos = pos + char_per_line;

  const char* l1 = line1.c_str();
  const char* l2 = line2.c_str();
  const char* l3 = line3.c_str();

  unsigned n_print_now;
  while (pos < len)
  { 
    back_pos = pos + char_per_line;

    if (back_pos > len)
    {
      back_pos    = len;
    }
    n_print_now = (unsigned)(back_pos - pos);

    //    fprintf(of, "debug: %d\n", n_print_now);

    fprintf(of, "%7d %.*s\n",     pos, n_print_now, top_str.data());
    fprintf(of, "        %.*s\n", n_print_now, l1+pos);
    fprintf(of, "        %.*s\n", n_print_now, l3+pos);
    fprintf(of, "        %.*s\n\n", n_print_now, l2+pos);

    pos += char_per_line;

  }
}

#endif
