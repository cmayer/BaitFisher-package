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
#ifndef GFF_COLLECTION_H
#define GFF_COLLECTION_H

#include "cstdio"
#include "easystring.h"
#include <vector>
#include <list>
#include <set>
#include <map>
#include "range_functions.h"
#include "GFF-class.h"


using namespace std;

static int ____verbosity = 1;


// Currently, GFF_record_list has no destructor.
// This can lead to memory leaks, since this class stores
// pointers to GFF_record, which will never be destructed explicitly.
// TODO!!!!

class GFF_record_list : public list<GFF_record*>
{
 public:
  void add(GFF_record *pf)
  {
    //    cout << "Adding: ";  pf->print(cout); cout << endl; 
    push_back(pf);
  }

  // Keep only entries with a feature value f.
  void filter_by_feature(easystring &f)
  {
     list<GFF_record*>::iterator it, it_end;
     it     = begin();
     it_end = end();

     while (it != it_end)
     {
       if ((**it).getFeature() != f )
	 it = erase(it);
       else
	 ++it;
     }
  }

  // Find all attribute values for a given set of attribute names:
  set<easystring> getSetOfAttributeValues(vector<easystring> att_name_vec, bool case_sensitive = true)
  {
    set<easystring> s;
/*     list<GFF_record*>::iterator it, it_end; */

/*     it     = begin(); */
/*     it_end = end(); */

/*     while (it != it_end) */
/*     { */
/*       s.insert((**it).getGeneName()); */
/*       ++it; */
/*     } */
    return s;
  }

  void getListOfGeneLocations(GFF_record_list &logl)
  {
/*     list<GeneLocation*>  logl; */
/*     GeneLocation         *gl; */
/*     set<string> s = getSetOfGeneNames(); */
/*     set<string>::iterator  s_it     = s.begin(); */
/*     set<string>::iterator  s_it_end = s.end(); */
/*     string    geneName; */
/*     CfeatureList::iterator l_it; */
/*     CfeatureList::iterator l_it_end; */
/*     unsigned range_start, range_end; */
/*     int      direction = 0; */
/*     bool     found; */

/*     while(s_it != s_it_end) */
/*     { */
/*       geneName = *s_it; */
/*       l_it     = begin(); */
/*       l_it_end = end(); */
/*       found = false; */

/*       // Find first occurence: */
/*       while(l_it != l_it_end) */
/*       { */
/* 	if ((**l_it).getGeneName() == geneName) */
/* 	{ */
/* 	  range_start = (**l_it).getStart(); */
/* 	  range_end   = (**l_it).getEnd(); */
/* 	  direction   = (**l_it).getDirection(); */
/* 	  found = true; */
/* 	  break; */
/* 	} */
/* 	++l_it; */
/*       } */

/*       // Find other occurences: */
/*       while(l_it != l_it_end) */
/*       { */
/* 	if ((**l_it).getGeneName() == geneName) */
/* 	{ */
/* 	  range_union( range_start,          range_end, */
/* 		       (**l_it).getStart(),  (**l_it).getEnd(), */
/* 		       range_start,          range_end); */
/* 	} */
/* 	++l_it; */
/*       } */

/*       if (!found) */
/*       { */
/* 	cerr << "Fatal error in element function getListOfGeneLocations. Gene with name " */
/* 	     << geneName << " not found." << endl; */
/* 	exit(0); */
/*       } */

/*       gl = new GeneLocation(geneName, range_start, range_end, direction); */
/*       logl.push_back(gl); */
/*       ++s_it; */
/*     } */
/*     return logl; */
  }



  void getVectorOfFeatures_withGivenName(const easystring &the_geneName,
					 vector<GFF_record*>  &vogffr )
  {
/*     CfeatureList::iterator it, it_end; */

/*     vogf.clear(); */

/*     it     = begin(); */
/*     it_end = end(); */

/*     while (it != it_end) */
/*     { */
/*       if ( (**it).getGeneName() == the_geneName ) */
/*       { */
/* 	vogf.push_back(*it); */
/* 	//	list_of_features.push_back(*it); */
/*       } */
/*       ++it; */
/*     } */
/*     stable_sort(vogf.begin(), vogf.end(), lessThan_genomeFeature); */
  }

/*   CfeatureList::iterator find_equal_feature(CgenomeFeature &gf) */
/*   { */
/*     CfeatureList::iterator  it      = begin(); */
/*     CfeatureList::iterator  it_end  = end(); */

/*     it     = lower_bound(it, it_end, &gf, lessThan_genomeFeature); */
/*     it_end = upper_bound(it, it_end, &gf, lessThan_genomeFeature); */

/*     while (it != it_end) */
/*     { */
/*       if ( **it == gf ) */
/* 	return it; */
/*       ++it; */
/*     } */
/*     return end(); */
/*   } */

  void print(ostream &os, int flag = 0)
  {
    list<GFF_record*>::iterator it, it_end;
    
    it     = begin();
    it_end = end();

    while (it != it_end)
    {
      (**it).print(os, flag);
      os << endl;
      ++it;
    }
  }

  unsigned print_sequences(FILE *of, easystring &format, CSequence_Mol &seq, int char_per_line,
			   bool ref_comp_if_neg_direction = true)
  {
    GFF_record_list::iterator  it, it_end;
    easystring                 seq_name = seq.getName();
    unsigned                   count    = 0;

    it     = begin();
    it_end = end();

    while (it != it_end)
    {
      //      if (it->getSeqid == seq_name) // This is also checked in the print_sequence member function of the GFF_record.
      {
	if ( (**it).print_sequence(of, format, seq, char_per_line, ref_comp_if_neg_direction) )
	  ++count;
      }
      ++it;
    }
    return count;
  }

}; // End class GFF_record_list : public list<GFF_record*>







class GFF_collection
{
  map<easystring, int>   sequences;  // counts the GFF-records associated to sequences
  map<easystring, int>   features;   // counts the GFF-records associated to features
  unsigned               total_number;
  unsigned               max_seq_name_length;

  // Type reminder: class GFF_record_list : public list<GFF_record*>

  GFF_record_list        gff_record_list;  // not inherited since then we can lose control
                                           // over features added to the list.
  
 public:
 GFF_collection(int v=0):total_number(0), max_seq_name_length(0)
  {}

  void reset()
  {
    sequences.clear();
    features.clear();
    gff_record_list.clear();

    total_number = 0;
    max_seq_name_length = 0;
  }


  void add(GFF_record* a)
  {
    map<easystring, int>::iterator  it;

    easystring sid = a->getSeqid();
    easystring fea = a->getFeature();

    ++total_number;
    if ( max_seq_name_length < sid.length() )
      max_seq_name_length = sid.length();

    // Search for sequnece name in map
    it  = sequences.find(sid);
    if (it != sequences.end() )
    {
      ++(it->second);  // Count the feature for this list
    }
    else // Add seq to map and set counter to 1 
    {
      sequences[sid] = 1;
    }

    // Search for feature name in map
    it  = features.find(fea);
    if (it != features.end() )
    {
      ++(it->second);  // Count the feature for this list
    }
    else // Add seq to map and set counter to 1 
    {
      features[fea] = 1;
    }

    gff_record_list.add(a);
  }

  unsigned size()
  {
    return gff_record_list.size();
  }


  void add(GFF_collection &col)
  {
    GFF_record_list::iterator it, it_end;
    
    it     = col.gff_record_list.begin();
    it_end = col.gff_record_list.end();
    
    while (it != it_end)
    {
      add(*it);
      ++it;
    }

  }

  static void set_verbosity(int v)
  {
    ____verbosity = v;
  }


  const GFF_record_list& getGFFrecordList()
  {
    return gff_record_list;
  }

  unsigned getNumGFFrecords()
  {
    return total_number;
  }

  unsigned getMaxSeqidLength()
  {
    return max_seq_name_length;
  }

  void print(ostream &os, int flag = 0)
  {
    gff_record_list.print(os, flag);
  }

  void sort()
  { 
    gff_record_list.sort(position_lessThan_GFF_record_pointer);
  }

  void get_List_of_gff_records_of_a_sequence(const easystring  &seq,
					     GFF_record_list   &list_of_records)
  {
    GFF_record_list::iterator it, it_end;

    it     = gff_record_list.begin();
    it_end = gff_record_list.end();

    while (it != it_end)
    {
      if ( (**it).getSeqid() == seq )
      {
	list_of_records.push_back(*it);
      }
      ++it;
    }
  }


  // Worries: Perfomrance.
  // Note: attribute_key must be lower case, since we only do a lower case comparison.
  void get_List_of_gff_records_for_attribute(const easystring  &attribute_key,
					     const easystring  &attribute_value,
					     GFF_record_list   &list_of_records)
  {
    GFF_record_list::iterator it, it_end;

    it     = gff_record_list.begin();
    it_end = gff_record_list.end();

    list_of_records.clear();

    while (it != it_end)
    {
      easystring val = (**it).getValueOfAttribute(attribute_key, true);

      if ( val == attribute_value )
      {
	list_of_records.push_back(*it);
      }
      ++it;
    }
  }


  vector<easystring> getVectorOfSeqid()
  {
    vector<easystring>  v;
    map<easystring, int>::iterator it, it_end;

    it     = sequences.begin();
    it_end = sequences.end();
    while (it != it_end)
    {
      v.push_back( it->first );
      ++it;
    }
    return v;
  }

  void get_min_max_number_of_features_for_Seqid(int &mi, int &ma)
  {
    mi = ma = 0; // Only used when collection is empty
    map<easystring, int>::iterator it, it_end;

    it     = sequences.begin();
    it_end = sequences.end();

    if (it != it_end)
      mi = ma = it->second;
    while (it != it_end)
    {
      if (it->second > ma)
	ma = it->second;
      if (it->second < mi)
	mi = it->second;
      ++it;
    }
  }


  void addSeqidsToSet(set<easystring> &s)
  {
    map<easystring, int>::iterator it, it_end;

    it     = sequences.begin();
    it_end = sequences.end();

    while (it != it_end)
    {
      s.insert( it->first );
      ++it;
    }
  }

  
  void readGFFfile(easystring filename,
		   unsigned   &good_lines,
		   unsigned   &bad_lines,
		   bool       &file_error )
  {
    ifstream    is(filename.c_str() );
    easystring  line;
    GFF_record  *gffr;
    bool        match_OK;

    file_error = false;
    good_lines = 0;
    bad_lines  = 0;
    
    if ( is.fail() )
    {
      file_error = true;
      return;
    }

    getline(is, line);
    while (is)
    {
      line.removeSpacesBack();
      if (line[0] == '#' || line.empty() )
      {
	getline(is, line);
	continue;
      }

      gffr     = new GFF_record(line);
      match_OK = gffr->good();

      if (match_OK)
      {
	add(gffr);
	++good_lines;
      }
      else
      {
	++bad_lines;
	std::cerr << endl << "Problematic line: " << line << endl; 
	std::cerr << "Reason for error: " << gffr->getErrorReason() << endl; 
	delete gffr;
      }
      getline(is, line);
    }
    is.close();

    //    gff_record_list.sort(lessThan_genomeFeature);
  }

  void print_sequences_all_features(FILE *of, easystring &format, CFile *seq_file, int char_per_line,
				    bool ref_comp_if_neg_direction,
				    unsigned &count_seq, unsigned &count_printed_features)
  {
    // Variables to read sequence file:
    CSequence_Mol   seq(CSequence_Mol::dna);

    count_seq=0;
    count_printed_features=0;

    //std::ios::sync_with_stdio(false); 

      // For each sequence: analyse the features
    while (!seq_file->eof())
    {
      seq.readRawFastaSequence(*seq_file);
      ++count_seq;

      //      cerr << count_seq << " " << seq.getName() << endl;

      count_printed_features += gff_record_list.print_sequences(of, format, seq, char_per_line, ref_comp_if_neg_direction);
    }
  }


  void getVectorOfFeatures(vector<easystring> &the_features)
  {
    map<easystring, int>::iterator it, it_end;
    it     = features.begin(); 
    it_end = features.end(); 

    the_features.clear();

    while (it != it_end)
    {
      the_features.push_back(it->first);
      ++it;
    }
  }


  void getRangeList(CRangeList &rl, const easystring &seq)
  {
    GFF_record_list::iterator  it, it_end;
    it     = gff_record_list.begin();
    it_end = gff_record_list.end();

    rl.clear();

    while (it != it_end)
    {
      if ( (**it).getSeqid() == seq )
      {
	// Remember: Ranges should be 0 based and the end should be the point after the last element in the range.
	// Coordinates in the gff file are 1 based. 
	rl.add( (**it).getStart(),
		(**it).getEnd()+1);  // the range class expects
      }
      ++it;
    }
  }


  void getConsolidatedRangeListOfGFFCollection_for_Sequence(CRangeList &crl, easystring &seq)
  {
    CRangeList rl;
    getRangeList(rl, seq);
    rl.get_consolidated_range_list(crl);
  }


  unsigned get_coverage_loop(std::ostream &os_error)
  {
    map<easystring, int>::iterator it, it_end;

    it     = sequences.begin();
    it_end = sequences.end();

    unsigned all_len=0;
    unsigned n;

    while (it != it_end)
    {
      // Do something with it->first:

      os_error << "Analysing sequence coverage:" << it->first << endl;
      CRangeList rl1;

      getRangeList(rl1, it->first);
      os_error << "Elements: " << rl1.size() << endl;
      rl1.sort();
      n = rl1.coverage_len_loop();
      os_error << "Contribution: " << n << endl;
      all_len += n;
      ++it;
    }
    return all_len;
  }


  unsigned get_coverage(std::ostream &os_error)
  {
    map<easystring, int>::iterator it, it_end;

    it     = sequences.begin();
    it_end = sequences.end();

    unsigned all_len=0;
    unsigned n;

    while (it != it_end)
    {
      // Do something with it->first:

      //   os_error << "Analysing sequence coverage:" << it->first << endl;
      CRangeList rl1;
      getRangeList(rl1, it->first);
      //   os_error << "Elements: " << rl1.size() << endl;
      n = rl1.coverage_len();
      //   os_error << "Contribution: " << n << endl;
      all_len += n;
      ++it;
    }
    return all_len;
  }



  friend void get_coverage_info_loop(GFF_collection &gff1, GFF_collection &gff2,
				     unsigned &len_gff1,  unsigned &len_gff2,
				     unsigned &len_union, unsigned &len_intersection,
				     unsigned &only_gff1, unsigned &only_gff2)
  {
    set<easystring> s;

    gff1.addSeqidsToSet(s);
    gff2.addSeqidsToSet(s);

    Crange limit(0, -1u);
    bool   err1, err2;
    CRangeList rl1;
    CRangeList rl2;

    unsigned len_gff1_this_seq;
    unsigned len_gff2_this_seq;
    unsigned len_union_this_seq;
    unsigned len_intersection_this_seq;
    unsigned limit_len_dummy;

    len_gff1         = 0;
    len_gff2         = 0;
    len_union        = 0;
    len_intersection = 0;
    only_gff1        = 0;
    only_gff2        = 0;

    set<easystring>::iterator it, it_end;
    map<easystring, int>::iterator it_find;

    it     = s.begin();
    it_end = s.end();

    while (it != it_end)
    {
      err1 = false;
      err2 = false;

      it_find = gff1.sequences.find(*it);
      if (it_find == gff1.sequences.end() )
      {
	if (____verbosity > 1)
	{
	  cerr << "Note: Sequence " << *it
	       << " is present in gff file 2 but not in file 1." << endl;
	}
	err1 = true;
      }
      it_find = gff2.sequences.find(*it);
      if (it_find == gff2.sequences.end() )
      {
	if (____verbosity > 1)
	{
	  cerr << "Note: Sequence " << *it
	       << " is present in gff file 1 but not in file 2." << endl;
	}
	err2 = true;
      }

      if (!(err1 || err2))
      {
	if (____verbosity > 0)
	{
	  cerr << "Analysing sequence " << *it << endl;
	}

	gff1.getRangeList(rl1, *it);
	gff2.getRangeList(rl2, *it);

	Crange b1 = rl1.bounds();
	Crange b2 = rl2.bounds();

	Crange limit = range_span(b1, b2);

	rl1.sort();
	rl2.sort();
	coverage_values_loop(limit, rl1, rl2, limit_len_dummy,
			     len_gff1_this_seq,
			     len_gff2_this_seq,
			     len_union_this_seq,
			     len_intersection_this_seq);
	
	//	cerr << len_gff1_this_seq         << endl;
	//	cerr << len_gff2_this_seq         << endl;
	//	cerr << len_union_this_seq        << endl;
	//	cerr << len_intersection_this_seq << endl;

	len_gff1         += len_gff1_this_seq;
	len_gff2         += len_gff2_this_seq;
	len_union        += len_union_this_seq;
	len_intersection += len_intersection_this_seq;
      }
      else // we still want to update the individual lengths as well as the union length:
      {
	if (!err1)
	{
	  gff1.getRangeList(rl1, *it);
	  rl1.sort();
	  unsigned n = rl1.coverage_len_loop();
	  len_gff1  += n;
	  len_union += n;
	  ++only_gff1;
	}
	if (!err2)
	{
	  gff2.getRangeList(rl2, *it);
	  rl2.sort();
	  unsigned n = rl2.coverage_len_loop();
	  len_gff2  += n;
	  len_union += n;
	  ++only_gff2;
	}
      }
      ++it;
    }
  }


  friend void get_coverage_info(GFF_collection &gff1, GFF_collection &gff2,
				unsigned &len_gff1,  unsigned &len_gff2,
				unsigned &len_union, unsigned &len_intersection,
				unsigned &only_gff1, unsigned &only_gff2)
  {
    set<easystring> s;

    gff1.addSeqidsToSet(s);
    gff2.addSeqidsToSet(s);

    Crange limit(0, -1u);
    bool   err1, err2;
    CRangeList rl1;
    CRangeList rl2;

    unsigned len_gff1_this_seq;
    unsigned len_gff2_this_seq;
    unsigned len_union_this_seq;
    unsigned len_intersection_this_seq;
    unsigned limit_len_dummy;

    len_gff1         = 0;
    len_gff2         = 0;
    len_union        = 0;
    len_intersection = 0;
    only_gff1        = 0;
    only_gff2        = 0;

    set<easystring>::iterator it, it_end;
    map<easystring, int>::iterator it_find;

    it     = s.begin();
    it_end = s.end();

    while (it != it_end)
    {
      err1 = false;
      err2 = false;
      it_find = gff1.sequences.find(*it);
      if (it_find == gff1.sequences.end() )
      {
	if (____verbosity > 1)
	{
	  cerr << "Note: Sequence " << *it
	       << " is present in gff file 1 but not in file 2." << endl;
	}
	err1 = true;
      }
      it_find = gff2.sequences.find(*it);
      if (it_find == gff2.sequences.end() )
      {
	if (____verbosity > 1)
	{
	  cerr << "Note: Sequence " << *it
	       << " is present in gff file 2 but not in file 1." << endl;
	}
	err2 = true;
      }
      
      if (!( err1 || err2 ))
      {
	if (____verbosity > 0)
	{
	  cerr << "Analysing sequence " << *it << endl;
	}

	gff1.getRangeList(rl1, *it);
	gff2.getRangeList(rl2, *it);

	Crange b1 = rl1.bounds();
	Crange b2 = rl2.bounds();

	Crange limit = range_span(b1, b2);

	rl1.sort();
	rl2.sort();
	coverage_values_loop(limit, rl1, rl2, limit_len_dummy,
			     len_gff1_this_seq,
			     len_gff2_this_seq,
			     len_union_this_seq,
			     len_intersection_this_seq);
	
	//	cerr << len_gff1_this_seq         << endl;
	//	cerr << len_gff2_this_seq         << endl;
	//	cerr << len_union_this_seq        << endl;
	//	cerr << len_intersection_this_seq << endl;

	len_gff1         += len_gff1_this_seq;
	len_gff2         += len_gff2_this_seq;
	len_union        += len_union_this_seq;
	len_intersection += len_intersection_this_seq;
      }
      else // we still want to update the individual lengths as well as the union length:
      {
	if (!err1)
	{
	  gff1.getRangeList(rl1, *it);
	  unsigned n = rl1.coverage_len();
	  len_gff1  += n;
	  len_union += n;
	  ++only_gff1;
	}
	if (!err2)
	{
	  gff2.getRangeList(rl2, *it);
	  unsigned n = rl2.coverage_len();
	  len_gff2  += n;
	  len_union += n;
	  ++only_gff2;
	}
      }
      ++it;
    } // while
  } // END   friend void get_coverage_info


  // Main purpose: Find in1in2not3not4 = (1^2)\(3 U 4)
  friend void get_coverage_info_in1in2not3no4_loop
    (GFF_collection &gff1,     GFF_collection &gff2,
     GFF_collection &gff3,     GFF_collection &gff4,
     unsigned &len_gff1,       unsigned &len_gff2,
     unsigned &len_gff3,       unsigned &len_gff4,
     unsigned &in1not3not4,    unsigned &in2not3not4,
     unsigned &in1in2not3not4, unsigned &in3in4not1not2
     )
  {
    len_gff1       = 0;
    len_gff2       = 0;
    len_gff3       = 0;
    len_gff4       = 0;
    in1not3not4    = 0;
    in2not3not4    = 0;
    in1in2not3not4 = 0;
    in3in4not1not2 = 0;

    set<easystring> s;

    gff1.addSeqidsToSet(s);
    gff2.addSeqidsToSet(s);
    gff3.addSeqidsToSet(s);
    gff4.addSeqidsToSet(s);

    //    bool   noseq1, noseq2, noseq3, noseq4;
    CRangeList rl1;
    CRangeList rl2;
    CRangeList rl3;
    CRangeList rl4;

    unsigned len_gff1_this_seq;
    unsigned len_gff2_this_seq;
    unsigned len_gff3_this_seq;
    unsigned len_gff4_this_seq;
    unsigned in1not3not4_this_seq;
    unsigned in2not3not4_this_seq;
    unsigned in1in2not3not4_this_seq;
    unsigned in3in4not1not2_this_seq;

    set<easystring>::iterator it, it_end;
    map<easystring, int>::iterator it_find;

    it     = s.begin();
    it_end = s.end();

    while (it != it_end)
    {
      if (____verbosity > 0)
      {
	cerr << "Analysing sequence: " << *it << endl;
      }
/*       noseq1 = false; */
/*       noseq2 = false; */
/*       noseq3 = false; */
/*       noseq4 = false; */

      it_find = gff1.sequences.find(*it);
      if (it_find == gff1.sequences.end() )
      {
	if (____verbosity > 1)
	{
	  cerr << "Note: Sequence " << *it
	       << " is not present in gff file 1." << endl;
	}
	//	noseq1 = true;
	rl1 = CRangeList();
      }
      else
      {
	gff1.getRangeList(rl1, *it);
      }

      it_find = gff2.sequences.find(*it);
      if (it_find == gff2.sequences.end() )
      {
	if (____verbosity > 1)
	{
	  cerr << "Note: Sequence " << *it
	       << " is not present in gff file 2." << endl;
	}
	//	noseq2 = true;
	rl2 = CRangeList();
      }
      else
      {
	gff2.getRangeList(rl2, *it);
      }

      it_find = gff3.sequences.find(*it);
      if (it_find == gff3.sequences.end() )
      {
	if (____verbosity > 1)
	{
	  cerr << "Note: Sequence " << *it
	       << " is not present in gff file 3." << endl;
	}
	//	noseq3 = true;
	rl3 = CRangeList();
      }
      else
      {
	gff3.getRangeList(rl3, *it);
      }

      it_find = gff4.sequences.find(*it);
      if (it_find == gff4.sequences.end() )
      {
	if (____verbosity > 1)
	{
	  cerr << "Note: Sequence " << *it
	       << " is not present in gff file 4." << endl;
	}
	//	noseq4 = true;
	rl4 = CRangeList();
      }
      else
      {
	gff4.getRangeList(rl4, *it);
      }


	
      Crange b1 = rl1.bounds();
      Crange b2 = rl2.bounds();
      Crange b3 = rl2.bounds();
      Crange b4 = rl2.bounds();

      Crange limit = range_span(b1, b2);
             limit = range_span(limit, b3);
             limit = range_span(limit, b4);

      rl1.sort();
      rl2.sort();
      rl3.sort();
      rl4.sort();

      coverage_values_loop(limit,
			   rl1, rl2, rl3, rl4,
			   len_gff1_this_seq,
			   len_gff2_this_seq,
			   len_gff3_this_seq,
			   len_gff4_this_seq,
			   in1not3not4_this_seq,
			   in2not3not4_this_seq,
			   in1in2not3not4_this_seq,
			   in3in4not1not2_this_seq
			   );

/*       cerr << 	"len1: " << len_gff1_this_seq         << endl; */
/*       cerr << 	"len2: " << len_gff2_this_seq         << endl; */
/*       cerr << 	"len3: " << len_gff3_this_seq         << endl; */
/*       cerr << 	"len4: " << len_gff4_this_seq         << endl; */

/*       cerr << 	"in1n34: " << in1not3not4_this_seq      << endl; */
/*       cerr << 	"in2n34: " << in2not3not4_this_seq      << endl; */
/*       cerr << 	"in12n34: " << in1in2not3not4_this_seq   << endl; */
/*       cerr << 	"in34n12: " << in3in4not1not2_this_seq   << endl; */

      len_gff1 += len_gff1_this_seq;
      len_gff2 += len_gff2_this_seq;
      len_gff3 += len_gff3_this_seq;
      len_gff4 += len_gff4_this_seq;

/*       cerr << 	"int-res: len1: " << len_gff1         << endl; */
/*       cerr << 	"int-res: len2: " << len_gff2         << endl; */
/*       cerr << 	"int-res: len3: " << len_gff3         << endl; */
/*       cerr << 	"int-res: len4: " << len_gff4         << endl; */


      in1not3not4 += in1not3not4_this_seq;
      in2not3not4 += in2not3not4_this_seq;
      in1in2not3not4 += in1in2not3not4_this_seq;
      in3in4not1not2 += in3in4not1not2_this_seq;

      ++it;
    }
  } // END friend void get_coverage_info_in1in2not3no4_loop


  void set_to_complement_in_universe(GFF_collection &universe, GFF_collection &gff,
				    easystring creator, easystring feature)
  {
    int mi, ma;
    universe.get_min_max_number_of_features_for_Seqid(mi, ma);

    if (mi != 1 || ma != 1)
    {
      std::cerr << "Error in set_to_complement_in_inverse: In the universe collection each sequence id can only occurr once."
		<< std::endl;
      return;
    }

    vector<easystring> vec_seqIds = universe.getVectorOfSeqid();

    unsigned    i, N=vec_seqIds.size();
    CRangeList  rl1, rl2, rl2res;
    Crange      limit;
    GFF_record  *p;
    easystring  name;

    for (i=0; i<N; ++i)
    {
      std::cerr << "Compute complement in " << vec_seqIds[i] << endl;
      universe.getRangeList(rl1, vec_seqIds[i]);
      gff.getRangeList(rl2, vec_seqIds[i]);
      limit = *(rl1.begin());

      rl2res.set_to_complement_of(rl2, limit);

      //      std::cerr << "SIZE: " << rl2res.size() << endl;

      CRangeList::CRangeList_iterator_const rl_it, rl_it_end;

      rl_it     = rl2res.begin();
      rl_it_end = rl2res.end();

      easystring n100 = ".";

      while (rl_it != rl_it_end)
      {
	// Remember the conventions here:
	// The ranges in rl2res are non inclusive for the end of the range.
	// But GFF_records store ranges such that the end coordinate is
	// the last coordinate in the range. So we have to subract 1.
	name = easystring("name=") + vec_seqIds[i] + "_" + easystring(rl_it->begin()) + "_" + easystring(rl_it->end()-1); 
	p = new GFF_record(vec_seqIds[i], creator, feature,
			   rl_it->begin(),  rl_it->end()-1, n100, '.', '.', name);
	add(p);

	//	p->print(cout);
	//	std::cout << endl;

	//	delete p;

	++rl_it;
      }

    }
  }


  void set_to_union_of(GFF_collection &gff1, GFF_collection &gff2,
		       easystring creator, easystring feature)
  {
    // Determine set of all Seqids:
    set<easystring> s;

    gff1.addSeqidsToSet(s);
    gff2.addSeqidsToSet(s);

    set<easystring>::iterator it, it_end;
    map<easystring, int>::iterator it_find;

    it     = s.begin();
    it_end = s.end();

    //    bool present1, present2;

    CRangeList rl1;
    CRangeList rl2;
    CRangeList rlres;

    easystring  name;

    // For all Seqids 
    while (it != it_end)
    {
/*       present1 = false; */
/*       present2 = false; */
      rl1.clear();
      rl2.clear();
      rlres.clear();
 
      it_find = gff1.sequences.find(*it);
      if (it_find == gff1.sequences.end() )
      {
	if (____verbosity > 1)
	{
	  cerr << "Note: Sequence " << *it
	       << " is present in gff file 2 but not in file 1." << endl;
	}
      }
      else
      {
	//	present1 = true;
	gff1.getRangeList(rl1, *it);
      }
      
      it_find = gff2.sequences.find(*it);
      if (it_find == gff2.sequences.end() )
      {
	if (____verbosity > 1)
	{
	  cerr << "Note: Sequence " << *it
	       << " is present in gff file 1 but not in file 2." << endl;
	}
      }
      else
      {
	//	present2 = true;
	gff2.getRangeList(rl2, *it);
      }

      rlres.set_to_range_list_union(rl1, rl2);

      easystring n100 = ".";
      easystring seq_name = *it;
      GFF_record  *p;
 
      CRangeList::CRangeList_iterator_const rl_it, rl_it_end;
      rl_it     = rlres.begin();
      rl_it_end = rlres.end();
    
      while (rl_it != rl_it_end)
      {
	// Remember the conventions here:
	// The ranges in rlres are non inclusive for the end of the range.
	// But GFF_records store ranges such that the end coordinate is
	// the last coordinate in the range. So we have to subract 1.
	name = easystring("name=") + seq_name + "_" + easystring(rl_it->begin()) + "_" + easystring(rl_it->end()-1); 
	p = new GFF_record(seq_name, creator, feature,
			   rl_it->begin(),  rl_it->end()-1, n100, '.', '.', name);
	add(p);
	//	p->print(cout);
	//	std::cout << endl;
	//	delete p;
	++rl_it;
      } // while (rl_it != rl_it_end)
      ++it;
    } // while (it != it_end) 
    
  }

  void set_to_intersection_of(GFF_collection &gff1, GFF_collection &gff2,
			      easystring creator, easystring feature)
  {
    // Determine set of all Seqids:
    set<easystring> s;

    gff1.addSeqidsToSet(s);
    gff2.addSeqidsToSet(s);

    set<easystring>::iterator it, it_end;
    map<easystring, int>::iterator it_find;

    it     = s.begin();
    it_end = s.end();

    //    bool present1, present2;

    CRangeList rl1;
    CRangeList rl2;
    CRangeList rlres;

    easystring  name;

    // For all Seqids 
    while (it != it_end)
    {
      //      present1 = false;
      //      present2 = false;
      rl1.clear();
      rl2.clear();
      rlres.clear();
 
      it_find = gff1.sequences.find(*it);
      if (it_find == gff1.sequences.end() )
      {
	if (____verbosity > 1)
	{
	  cerr << "Note: Sequence " << *it
	       << " is present in gff file 2 but not in file 1." << endl;
	}
      }
      else
      {
	//	present1 = true;
	gff1.getRangeList(rl1, *it);
      }
      
      it_find = gff2.sequences.find(*it);
      if (it_find == gff2.sequences.end() )
      {
	if (____verbosity > 1)
	{
	  cerr << "Note: Sequence " << *it
	       << " is present in gff file 1 but not in file 2." << endl;
	}
      }
      else
      {
	//	present2 = true;
	gff2.getRangeList(rl2, *it);
      }

      rlres.set_to_range_list_intersection(rl1, rl2);

      //            std::cerr << "rl1" << endl;
      //            rl1. print_bracketed_list(std::cerr);
      //            std::cerr << "rl2" << endl;
      //            rl2. print_bracketed_list(std::cerr);
      //            std::cerr << "rlres" << endl;
      //            rlres. print_bracketed_list(std::cerr);

      easystring n100 = ".";
      easystring seq_name = *it;
      GFF_record  *p;
 
      CRangeList::CRangeList_iterator_const rl_it, rl_it_end;
      rl_it     = rlres.begin();
      rl_it_end = rlres.end();
    
      while (rl_it != rl_it_end)
      {
	// Remember the conventions here:
	// The ranges in rlres are non inclusive for the end of the range.
	// But GFF_records store ranges such that the end coordinate is
	// the last coordinate in the range. So we have to subract 1.
	name = easystring("name=") + seq_name + "_" + easystring(rl_it->begin()) + "_" + easystring(rl_it->end()-1); 
	p = new GFF_record(seq_name, creator, feature,
			   rl_it->begin(),  rl_it->end()-1, n100, '.', '.', name);
	add(p);
	//	p->print(cout);
	//	std::cout << endl;
	//	delete p;
	++rl_it;
      } // while (rl_it != rl_it_end)
      ++it;
    } // while (it != it_end) 

  }





}; // End class GFF_collection



#endif
