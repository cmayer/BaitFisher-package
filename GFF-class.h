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
#ifndef GFF_RECORD_H
#define GFF_RECORD_H

#include "easystring.h"
#include <vector>
#include "range_functions.h"
#include "cstdio"               // for FILE
#include "CSequence_Mol2_1.h"

//***********************************************************************************
// In the long run, we intend to follow the conventions specified in: http://www.bioperl.org/wiki/GFF3
// In particular the handling of special characters does not conform to this standard yet.

// Convention: fstart and fend are the direct gff field entries 4 and 5. Minimum value is 1. The end value is the last position in the feature.
//             This convention is not identical to the convention in CgenomeFeature class, where all values are 0 based.
//             In this class the philisoph is to be as close as possible to the value in the file.
//             Example: First three bases of seq.: 1,3 
//
//             Convention has been checked.
//***********************************************************************************

using namespace std;

static const easystring empty_easy_string;

// Forward declaration of friend functions of the classes in this file. This avoids the friend-injection problem:
class GFF_record;
bool position_lessThan_GFF_record(GFF_record &, GFF_record &);
bool position_lessThan_GFF_record_pointer(GFF_record *, GFF_record *);


class GFF_record {

 private:
  easystring                   parse_error;

  easystring                   seqname;     // gff field 1
  easystring                   source;      // gff field 2
  easystring                   feature;     // gff field 3
  unsigned                     fstart;      // gff field 4
  unsigned                     fend;        // gff field 5 // According to the gff file conventions this is an inclusive end.
  easystring                   score;       // gff field 6
  char                         strand;      // gff field 7
  char                         phase;       // gff field 8
  std::vector<easystring>      att_names;
  std::vector<easystring>      att_names_lower_case;
  std::vector<easystring>      att_values;


 public:
  GFF_record (const std::string &str)
  {
    parse_string(str);
  }

  GFF_record (const easystring &str)
  {
    parse_string(str);
  }

  GFF_record (easystring &p_seqname,
	      easystring &p_source,
	      easystring &p_feature,
	      unsigned   p_fstart,
	      unsigned   p_fend,
	      easystring &p_score,
	      char       p_strand,
              char       p_phase,
	      easystring &p_att
	      ):
  seqname(p_seqname), source(p_source), feature(p_feature), fstart(p_fstart),
    fend(p_fend), score(p_score), strand(p_strand), phase(p_phase)
  {
    parse_attributes(p_att);
  }


  unsigned good() const
  {
    return parse_error.empty();
  }

  std::string getErrorReason() const
  {
    return parse_error;
  }


 private:

  bool parse_attributes(const easystring &tmp_att)
  {
    int i, att_num;
    std::vector<std::easystring> sv;

    att_num = split_respect(sv, tmp_att, ";");
    
    //    std::cerr << "+++ Number of attributes found: " << att_num << std::endl;

    easystring   key;
    easystring   val;
    std::size_t  index;
    //    unsigned     num_quotes1, num_quotes2;

    for (i=0; i<att_num; ++i)
    {
      sv[i].removeSpacesFront();
      sv[i].removeSpacesBack();
      
      index = sv[i].find_first_of(" =");
      if (index == std::string::npos)
      {
	parse_error = "Unrecognized attribute: \"" + sv[i] + "\"";
	return false;
      }
      key = easystring(sv[i].begin(),     sv[i].begin()+index);
      //      cerr << "new key: " << key << endl;
      while (index < sv[i].size() &&
	     (sv[i][index] == ' ' || sv[i][index] == '=') )
	++index;
      if (index == sv[i].size() )
      {
	parse_error = "Unrecognized attribute: \"" + sv[i] + "\"";
	return false;	
      }
      val = easystring(sv[i].begin()+index, sv[i].end() );
      //      cerr << "new val: " << val << endl;

      // Drop the following checks since quotes can occur as singletons:
      // Example: "3' UTR"
/*       // More error checks concerning val could be useful */
/*       num_quotes1 = val.countChar('\"'); */
/*       num_quotes2 = val.countChar('\''); */

/*       if (num_quotes1 % 2 != 0 || num_quotes2 % 2 != 0) */
/*       { */
/* 	parse_error = "unmatched quotes on in attribute field. Attribute: " + key;  */
/* 	return false;	 */
/*       } */
      att_names.push_back(key);
      key.ToLower();
      att_names_lower_case.push_back(key);
      att_values.push_back(val);
    }
    return true;
  }

  bool parse_string(const easystring &str)
  {
    std::vector<std::easystring> sv;
    easystring                   tmp_att;
    int                          i, split_num;

    // We use split_strict since successive tabs are not allowed
    split_num = split_strict(sv, str, "\t");
    
    if (split_num != 9)
    {
      parse_error = "wrong number of tab deliminated fields";
      return false;
    }

    for (i=0; i<9; ++i)
    {
      sv[i].removeSpacesFront();
      sv[i].removeSpacesBack();

      if ( sv[i].empty() )
      {
	parse_error = "on the tab deliminated line the field with position" + easystring(char(i+49)) + " is empty";
	return false;
      }
    }

    seqname = sv[0];
    source  = sv[1];
    feature = sv[2];
    fstart  = sv[3].ToUnsigned();
    fend    = sv[4].ToUnsigned();
    score   = sv[5];
    strand  = sv[6][0];          // Result must be checked.
    phase   = sv[7][0];          // Result must be checked.
    tmp_att = sv[8];

    if (fstart > fend)
    {
      parse_error = "start > end";
      return false;
    }
    
    if (strchr("+-?.0",strand) ==  NULL)
    {
      parse_error = easystring("wrong strand symbol \"") + strand + "\". Must be one of \"+-?.0\"";
      return false;
    }
    
    if (strchr("012+-.", phase) ==  NULL)
    {
      parse_error = easystring("wrong phase symbol \"")  + phase + "\". Must be one of \"012+-.\"";
      return false;
    }

    return parse_attributes(tmp_att);
  }

 public:

  easystring         getSeqid   () { return seqname; }
  easystring         getSource  () { return source;  }
  easystring         getFeature () { return feature; }
  unsigned           getStart   () { return fstart;  }
  unsigned           getEnd     () { return fend;    }

  unsigned           getStartIndex   () { return fstart-1;  }  // Lowest index is 0
  unsigned           getEndIndex     () { return fend;      }  // Is index after last index that lies in the feature.
                                                               // This is equal to the last position in feature. fend has lowest value 1.
                                                               // Example: Start codon at beginning of sequence: gff positions 1,3. Sequence indices: 0,3

  easystring         getScoreAsString() { return score; }

/*   std::pair<bool, double> getScore()  */
/*   { */
/*     if (score.size() == 1 && score[0] == '.') */
/*       return make_pair(false, 0); */
/*     else */
/*       return make_pair(true, score.ToDouble());  */
/*   } */

  char               getStrand () { return strand; }
  char               getPhase  () { return phase;  }

  unsigned           length()     { return fend-fstart; }

  const vector<easystring>& getAttributeNames()
  {
    return att_names;
  }

  const vector<easystring>& getAttributeNamesLowerCase()
  {
    return att_names_lower_case;
  }

  int getIndexOfAttribute(const easystring &s, bool lower_case_comparison = false)
  {
    int i = 0, n = att_names.size();

    if (!lower_case_comparison)
      while (i < n && att_names[i] != s)
      {
	++i;
      }
    else
      while (i < n && att_names_lower_case[i] != s)
      {
	++i;
      }

    if (i==n)
      return -1;
    else
      return i;
  }


  unsigned countValuesOfAttribute(const easystring &s, bool lower_case_comparison = false)
  {
    int i = getIndexOfAttribute(s, lower_case_comparison);

    if (i == -1)
      return 0;
    else
      return att_values[i].countChar(',')+1;
  }

  unsigned countValuesOfAttribute(unsigned i)
  {
    if (i >= att_names.size())
      return 0;
    else
      return att_values[i].countChar(',')+1;
  }

  const easystring &getValueOfAttribute(const easystring &s, bool lower_case_comparison = false)
  {
    int i = getIndexOfAttribute(s, lower_case_comparison);

    if (i == -1)
      return empty_easy_string;
    else
      return att_values[i];
  }

  const easystring &getValueOfAttribute(unsigned i)
  {
    if (i >= att_names.size())
      return empty_easy_string;
    else
      return att_values[i];
  }

  bool isFeature(const easystring &str)
  {
    return (str == feature);
  }


  void print(ostream &os, int flag=0)
  {
    unsigned i;

    if ( !good() )
    {
      os << "Invalid gff recored." << endl;
      os << "Reason: " << parse_error << endl;
      return;
    }
    if (flag == 0) // gff output format
    {
      os << seqname << "\t"
	 << source  << "\t"
	 << feature << "\t"
	 << fstart  << "\t"
	 << fend    << "\t"
	 << score   << "\t"
	 << strand  << "\t"
	 << phase   << "\t";

      for (i=0; i < att_names.size()-1; ++i)
      {
	os << att_names[i] << "=" << getValueOfAttribute(att_names[i]) << "; ";
      }
      os << att_names[i] << "=" << getValueOfAttribute(att_names[i]);
    }
    else if (flag == 1) // long multiline format
    {
      os << "seqname: " << seqname << endl
	 << "source:  " << source  << endl
	 << "feature: " << feature << endl
	 << "start:   " << fstart  << endl
	 << "end:     " << fend    << endl
	 << "score:   " << score   << endl
	 << "strand:  " << strand  << endl
	 << "phase:   " << phase   << endl;

      
      for (i=0; i < att_names.size(); ++i)
      {
	os << "attribute: (" << att_names[i]
	   << ") number of values: " << countValuesOfAttribute(att_names[i])
	   << ", value: " << getValueOfAttribute(att_names[i]) << endl;
      }
    }
  }

/*   bool is_in_range(unsigned a, unsigned e) */
/*   { */
/*     return is_r2_Subrange_of_r1(a, e, fstart, fend); */
/*   } */

/*   unsigned overlap_with_range(unsigned a, unsigned e) */
/*   { */
/*     return overlap(a, e, fstart, fend); */
/*   } */

/*   void union_with_range(unsigned a, unsigned e, unsigned &c1, unsigned &c2) */
/*   { */
/*     range_union(a, e, fstart, fend, c1, c2); */
/*   } */

/*   void intersection_with_range(unsigned a, unsigned e, unsigned &c1, unsigned &c2) */
/*   { */
/*     range_intersection(a, e, fstart, fend, c1, c2); */
/*   } */

  bool print_sequence(FILE *of, easystring format, CSequence_Mol &seq, int char_per_line,
		      bool ref_comp_if_neg_direction = true )
  {
    // the name of the sequence we print:
    easystring  seq_name;
    string      strand_desc;
    unsigned    i;

    if (seq.getName() != seqname)
      return false;

    // %I  seqID
    // %S  source
    // %F  feature
    // %b  begin
    // %e  end
    // %s  score
    // %d  strand (Direction) short
    // %D  strand (Direction) long
    // %p  phase
    // %A  attribute 0

    if ( format.empty() )
    {
      format = "%I#%S#%F#%b-%e#%p#%A%D";
    }

    unsigned N = format.size();

    for (i=0; i < N; ++i)
    {
      if (format[i] == '\\' )
      {
	++i;
	if (i == N)
	  break;
	seq_name += format[i];
      }
      else if (format[i] == '%')
      {
	++i;
	if (i == N)
	  break;
	switch(format[i])
	{
	case 'I': seq_name += seqname;              break;
	case 'S': seq_name += source;               break;
	case 'F': seq_name += feature;              break;
	case 'b': seq_name += easystring(fstart);   break;
	case 'e': seq_name += easystring(fend);      break;
	case 's': seq_name += score;                break;
	case 'd': seq_name += strand;               break;
	case 'D':
	  if ( strand == '+') {
	    strand_desc = "on+strand_is+strand";
	  }
	  else if ( strand == '-' ) {
	    if (ref_comp_if_neg_direction) {
	      seq_name += "on-strand_is-strand";
	    }
	    else {
	      seq_name += "on-strand_is+strand";
	    }
	  }
	  else {
	    seq_name += "not-stranded";
	  }
	  break;

	case 'p': seq_name += easystring( phase );               break;
	case 'A': seq_name += getValueOfAttribute(0);            break;
	default:  seq_name += easystring(format[i]);             break;
	}
      }
      else
	seq_name += easystring(format[i]);
    }

/*     seq_name = */
/*       seqname + */
/*       string("|") + */
/*       source + */
/*       string("|") + */
/*       feature + */
/*       string("|") + */
/*       easystring( fstart ) + */
/*       string(":") + */
/*       easystring( fend ) + */
/*       string("|") + */
/*       score + */
/*       string("|") + */
/*       easystring( strand ) + */
/*       string("|") + */
/*       easystring( phase ) + */
/*       string("|") + */
/*       getValueOfAttribute(0) + */
/*       string("|") + */
/*       strand_desc; */

      bool revComp = ref_comp_if_neg_direction && strand == '-';

      // fend is position at which feature ends.
      // The following write function expects 0 based indizes, 
      // Example: print first three bases, pass 0 and 3 in both cases, we have 1 and 3

      if (revComp)
	seq.writeSequence_revComp(of, fstart-1, fend, seq_name.c_str(), char_per_line);
      else
 	seq.writeSequence(of, fstart-1, fend, seq_name.c_str(), char_per_line);
      return true;
  }


  // friend functions:

  friend unsigned overlap(GFF_record &a, GFF_record &b)
  {
    if (a.seqname != b.seqname)
      return 0;
    return overlap(a.fstart, a.fend, b.fstart, b.fend);
  }

  friend unsigned overlap_in_same_direction(GFF_record &a, GFF_record &b)
  {
    if (a.strand != b.strand)
      return 0;

    if (a.seqname != b.seqname)
      return 0;

    return overlap(a.fstart, a.fend, b.fstart, b.fend);
  }

  friend bool position_lessThan_GFF_record(GFF_record &a, GFF_record &b)
  {
    if (a.seqname != b.seqname)
      return a.seqname < b.seqname;

    if (a.fstart < b.fstart)
      return true;

    if (a.fstart > b.fstart)
      return false;

    // We have a.fstart == b.fstart now!!!
    if ( equalRange(a, b) )
    {
      if ( a.feature == "exon" )
	return true;
      else
	return false;
    }

    if ( b.fend < a.fend )
      return true;
    else
      return false;
  }

  friend bool position_lessThan_GFF_record_pointer(GFF_record *a, GFF_record *b)
  {
    if (a->seqname != b->seqname)
      return a->seqname < b->seqname;

    if (a->fstart < b->fstart)
      return true;

    if (a->fstart > b->fstart)
      return false;

    // We have a.fstart == b.fstart now!!!
    if ( equalRange(*a, *b) )
    {
      if ( a->feature == "exon" )
	return true;
      else
	return false;
    }

    if ( b->fend < a->fend )
      return true;
    else
      return false;
  }



  friend bool operator==(GFF_record &a, GFF_record &b)
  {
    if (a.fstart  != b.fstart || a.fend != b.fend)
      return false;
    if (a.feature != b.feature)
      return false;
    if (a.seqname != b.seqname)
      return false;
    if (a.source  != b.source)
      return false;

    if (a.score   != b.score)
      return false;
    if (a.strand  != b.strand)
      return false;
    if (a.phase   != b.phase)
      return false;

    if (a.att_names != b.att_names)
      return false;
    if (a.att_values != b.att_values)
      return false;

    return true;
  }

  friend bool equalRange(GFF_record &a, GFF_record &b)
  {
    return equalRange(a.fstart, a.fend, b.fstart, b.fend);
  }

  friend bool isSubfeature_Arg2_of_Arg1(GFF_record &a, GFF_record &b)
  {
    return is_r2_Subrange_of_r1(a.fstart, a.fend, b.fstart, b.fend);
  }


  

};




#endif
