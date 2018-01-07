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
// Some parts of this header file could be made more efficient.
// Much of the functionality relies on a sorted list of Blast hits.
// The question arises whether a list is the best data structure for such 
// a purpose. It could be more efficient to use a data structure such as
// a set, which is sorted automatically and which should in principle be
// more efficient.

// GI number will not be supported by NCBI in the future! So far we have not removed the code.
// This should be done soon.

#ifndef CBLASTPARSER_H
#define CBLASTPARSER_H

// The following two defines allow us to use fopen() in "CSeqNameList.h"
// without error message when using the Visual Studio 2015 compiler.
#define _CRT_SECURE_NO_DEPRECATE
#define _CRT_SECURE_NO_WARNINGS

#include <iostream>
#include "CFile/CFile2_1.h"
#include "faststring2.h"
#include <vector>
#include <set>
#include <list>
#include "CSeqNameList.h"

#include "range_functions.h"
#include "CHistogram.h"

#include "typedefs.h"
#include <climits>

class CTaxonomy;
class CBlastHit;

// Forward declaration of friend functions in order to make them globally visible (friend injection problem):
bool less_than_CBlastHit_query(const CBlastHit &a, const CBlastHit &b);
bool less_than_CBlastHit_query_quality(const CBlastHit &a, const CBlastHit &b);
bool less_than_CBlastHit_target(const CBlastHit &a, const CBlastHit &b);
bool operator<(const CBlastHit &a, const CBlastHit &b);
bool less_than_CBlastHit_pointer_query(const CBlastHit *a, const CBlastHit *b);
bool less_than_CBlastHit_pointer_query_quality(const CBlastHit *a, const CBlastHit *b);
bool less_than_CBlastHit_pointer_target(const CBlastHit *a, const CBlastHit *b);


signed char CB_is_gi_of_type(CTaxonomy *tax, unsigned gi, unsigned type);

struct struct_less_than_CBlastHit_pointer_query_quality
{
  bool operator()(const CBlastHit *a, const CBlastHit *b)
  {
    return less_than_CBlastHit_pointer_query_quality(a,b);
  }
};


static void __add_or_count(std::map<unsigned, unsigned> &m, unsigned x)
{
  std::map<unsigned, unsigned>::iterator it;

  it = m.find(x);
  if (it == m.end() )
  {
    m[x] = 1;
  }
  else
  {
    ++it->second;
  }
}

/* //Not used in most applications. Not tested yet. Needs to be tested before it is used.
static void majority_rule(std::vector<unsigned> v, int &val, double &prop)
{
  CHistogram h(v, -10);

  double value;

  h.dominant_category(value, prop);

  val = (int)(value+0.5); // Adding 0.5 ensures we are not returning a wrong int due to
                          // rounding errors.
}
*/

class CBlastHit
{
 private:
  unsigned internal_number;

  faststring target_seq_name;
  faststring query_seq_name;

  double perfection;
  unsigned length;
  unsigned mismatches;
  unsigned gaps;

  bool reverse_query;
  bool reverse_target;
  
  unsigned start_query;
  unsigned end_query;
  unsigned start_target;
  unsigned end_target;

  double   e_value;
  double   score;

 public:

 CBlastHit(const faststring& line, unsigned in):internal_number(in)
  {
    std::vector<faststring> v;

    split(v, line);

    //    std::cout << v.size() << std::endl;

    if (v.size() == 12)
    {
      
      query_seq_name  = v[0];
      target_seq_name = v[1];

      perfection      = v[2].ToDouble();
      length          = v[3].ToUnsigned();
      mismatches      = v[4].ToUnsigned();
      gaps            = v[5].ToUnsigned();
  
      start_query     = v[6].ToUnsigned();
      end_query       = v[7].ToUnsigned();

      //      std::cerr << "Identfied start_query: " <<  start_query << std::endl;
      //      std::cerr << "Identfied end_query:   " <<  end_query << std::endl;

      start_target    = v[8].ToUnsigned();
      end_target      = v[9].ToUnsigned();
      
      e_value         = v[10].ToDouble();
      score           = v[11].ToDouble();

      if (start_query > end_query)
      {
	reverse_query= true;
	unsigned tmp;
	tmp         = start_query;
	start_query = end_query;
	end_query   = tmp;
      }
      else
	reverse_query= false;

      if (start_target > end_target)
      {
	reverse_target = true;
	unsigned tmp;
	tmp          = start_target;
	start_target = end_target;
	end_target   = tmp;
      }
      else
	reverse_target= false; 
    }
    else
    {
      std::cerr << "Unrecognised line of input passed to CBlastHit constructor: " << line << "\n";
    }
  }

  double get_perfection()
  {
    return perfection;
  }

  unsigned get_length()
  {
    return length;
  }

  unsigned get_query_span()
  {
    if (end_query >= start_query)
      return end_query - start_query+1;
    else
      return start_query - end_query+1;
  }

  // Searches for gi|xxxxxxxx| in the query sequence and returns the xxx string. 
  void extract_GI_number_query(faststring &gi)
  {
    faststring::size_t pos1 = query_seq_name.find("gi|");
    faststring::size_t pos2;
    if (pos1 != faststring::npos)
    {
      pos2 = query_seq_name.find("|", pos1+3);
      if (pos2 == faststring::npos) // No proper gi number found so we return empty string; 
      {
	gi.clear();
	return;
      }
      gi   = query_seq_name.substr(pos1+3, pos2-pos1-3);
      return;
    }
    else // No proper gi number found so we return empty string;
    {
      gi.clear();
      return;
    }
  }

  void extract_GI_number_target(faststring &gi)
  {
    faststring::size_t pos1 = target_seq_name.find("gi|");
    faststring::size_t pos2;
    if (pos1 != faststring::npos)
    {
      pos2 = target_seq_name.find("|", pos1+3);
      if (pos2 == faststring::npos) // No proper gi number found so we return empty string; 
      {
	gi.clear();
	return;
      }
      gi   = target_seq_name.substr(pos1+3, pos2-pos1-3);
      return;
    }
    else // No proper gi number found so we return empty string;
    {
      gi.clear();
      return;
    }
  }


  double get_evalue()
  {
    return e_value;
  }

  unsigned get_start_query()
  {
    return start_query;
  }

  unsigned get_end_query()
  {
    return end_query;
  }

  unsigned get_start_target()
  {
    return start_target;
  }

  unsigned get_end_target()
  {
    return end_target;
  }

  const faststring &get_query_seq_name()
  {
    return query_seq_name;
  }

  const faststring &get_target_seq_name()
  {
    return target_seq_name;
  }


  friend bool less_than_CBlastHit_query(const CBlastHit &a, const CBlastHit &b)
  {
    // sort first by -query name, -query start, -query end, - target name, target start, target end, -length, -perfection

    //    a.print(std::cerr);
    //    b.print(cerr);

    // Primary criterion: Query seq name:
    if (a.query_seq_name != b.query_seq_name)
      return a.query_seq_name < b.query_seq_name;

    if (a.start_query != b.start_query)
      return a.start_query < b.start_query;
    if (a.end_query != b.end_query)
      return a.end_query < b.end_query;

    if (a.target_seq_name != b.target_seq_name)
      return a.target_seq_name < b.target_seq_name;

    if (a.start_target != b.start_target)
      return a.start_target < b.start_target;

    if (a.end_target != b.end_target)
      return a.end_target < b.end_target;

    // They should be irrelevant since we would have multiple identical hits
    // if we come here.
    if (a.length != b.length)
      return a.length < b.length;

    //    if (a.perfection != b.perfection)   // Including this line is a bug. This was a bug.
      return a.perfection < b.perfection;
  }

  friend bool less_than_CBlastHit_query_quality(const CBlastHit &a, const CBlastHit &b)
  {
    // First criterion: name of query sequence: ascending (smaller is "better")
    // Second criterion: e_value ascending (smaller is better)
    // Third criterion: score descending (higher is better)
    // Fourth criterion: perfection descending (higher is better)
    // Fifth criterion: sequence positions
    // Sixth criterion: target sequence name and positions

    // sort first by -query name, -query start, -query end, - target name, target start, target end, -length, -perfection

    //    a.print(std::cerr);
    //    b.print(cerr);

    // Primary criterion: Query seq name:
    if (a.query_seq_name != b.query_seq_name)
      return a.query_seq_name < b.query_seq_name;

    // Secondary criterion: e_value
    if (a.e_value != b.e_value)
      return a.e_value < b.e_value;

    if (a.score != b.score)
      return a.score > b.score;

    if (a.perfection != b.perfection)
      return a.perfection > b.perfection;

    if (a.start_query != b.start_query)
      return a.start_query < b.start_query;

    if (a.end_query != b.end_query)
      return a.end_query < b.end_query;

    if (a.target_seq_name != b.target_seq_name)
      return a.target_seq_name < b.target_seq_name;

    if (a.start_target != b.start_target)
      return a.start_target < b.start_target;

    if (a.end_target != b.end_target)
      return a.end_target < b.end_target;

    // They should be irrelevant since we would have multiple identical hits
    // if we come here.
    //    if (a.length != b.length)   // Including this line is a bug. This was a bug.
      return a.length < b.length;
  }




  friend bool less_than_CBlastHit_target(const CBlastHit &a, const CBlastHit &b)
  {
    // sort first by -query name, -query start, -query end, - target name, target start, target end, -length, -perfection

    //    a.print(std::cerr);
    //    b.print(cerr);

    // Primary criterion: Target seq name:


    if (a.target_seq_name != b.target_seq_name)
      return a.target_seq_name < b.target_seq_name;

    if (a.start_target != b.start_target)
      return a.start_target < b.start_target;

    if (a.end_target != b.end_target)
      return a.end_target < b.end_target;

    if (a.query_seq_name != b.query_seq_name)
      return a.query_seq_name < b.query_seq_name;

    if (a.start_query != b.start_query)
      return a.start_query < b.start_query;

    if (a.end_query != b.end_query)
      return a.end_query < b.end_query;

    // They should be irrelevant since we would have multiple identical hits
    // if we come here.
    if (a.length != b.length)
      return a.length < b.length;

    //    if (a.perfection != b.perfection)  // Including this line is a bug. This was a bug.
      return a.perfection < b.perfection;
  }

  friend bool operator<(const CBlastHit &a, const CBlastHit &b)
  {
    //    std::cout << "---------------" << std::endl;
    return less_than_CBlastHit_query(a,b);
  }

  friend bool less_than_CBlastHit_pointer_query(const CBlastHit *a, const CBlastHit *b)
  {
    return (*a < *b);
  }

  friend bool less_than_CBlastHit_pointer_query_quality(const CBlastHit *a, const CBlastHit *b)
  {
    return (less_than_CBlastHit_query_quality(*a, *b));
  }

  friend bool less_than_CBlastHit_pointer_target(const CBlastHit *a, const CBlastHit *b)
  {
    return (less_than_CBlastHit_target(*a, *b));
  }

  void print(std::ostream &os, char delim = '\t') const
  {
    os << internal_number << delim;
    os << query_seq_name << delim;
    os << target_seq_name << delim;

    os << perfection << delim;
    os << length << delim;
    os << mismatches << delim;
    os << gaps << delim;
  
    os << start_query << delim;
    os << end_query << delim;
    os << start_target << delim;
    os << end_target << delim;
    
    os << e_value << delim;
    os << score;
    os << '\n';
  }

  void print(std::ostream &os,
	     bool pinternal_number,
	     bool pqueryname,
	     bool ptargetname,
	     bool pperfection,
	     bool plength,
	     bool pmismatches,
	     bool pgaps,
	     bool pstart_query,
	     bool pend_query,
	     bool pstart_target, 
	     bool pend_target,
	     bool pe_value, 
	     bool pscore,
	     char delim = '\t') const
  {
    if (pinternal_number)
      os << internal_number << delim;

    if (pqueryname)
      os << query_seq_name << delim;

    if (ptargetname)
      os << target_seq_name << delim;

    if (pperfection)
      os << perfection << delim;

    if (plength)
      os << length << delim;
    
    if (pmismatches)
      os << mismatches << delim;

    if (pgaps)
      os << gaps << delim;
  
    if (pstart_query)
      os << start_query << delim;
    
    if (pend_query)
      os << end_query << delim;

    if (pstart_target)
      os << start_target << delim;

    if (pend_target)
      os << end_target << delim;
    
    if (pe_value)
      os << e_value << delim;

    if (pscore)
      os << score;

    os << '\n';
  }

  void print_essential_1(std::ostream &os, char delim = '\t') const
  {
    os << query_seq_name << delim;
    os << "E: "   << e_value << delim;
    os << "len: " << length << delim;
    os << "%: "   << perfection << delim;
    os << '\n';
  }



};


class CBlast_parser
{
  std::list<CBlastHit*> lbh;
  // Previous versions of this header file used a map in the declaration of bh_map.
  // This map is/was only used in the get_max_hit_query_span__matching_beginning2 function.
  // In a map in which elements are inserted with the insert function we only inserted the first (and best entry) for
  // each query. This was fine for all cases for which we used this map so far.
  // However, a multimap has the same functionality and is more secure, since it will store all entries.

  std::multimap<faststring, CBlastHit*>  bh_map;

  // Stores coverage in bp. If query seq. length is known it can be used to determine the relative sequence length.
  std::map<faststring, unsigned> total_coverage_map; 

  bool query_name_quality_sorted;
  //TODO: add delete in destructor concept to improve performance in some constructors by avoiding to copy the data.

  void add(CBlastHit* p)
  {
    lbh.push_back(p);
    faststring tmp = p->get_query_seq_name();
    bh_map.insert(bh_map.end(), std::make_pair(tmp, p));
  }

 public:
 CBlast_parser(const char *filename):query_name_quality_sorted(false)
  {
    CFile f;
    faststring line;
    unsigned counter=0;

    f.ffopen(filename);

#ifdef FINANA
    std::cerr << "We are in the FINANA mode.\n";
#endif

    f.getline(line);
    while (!f.fail() )
    {
      line.removeSpacesFront();
      line.removeSpacesBack();
      
      if (!line.empty() && line[0] != '#') // if not a comment - allows to read the -m 9 in addition to the -m 8 output format
      {
	CBlastHit* p = new CBlastHit(line, counter);

	// Debugging code:
	//	if (counter > 1)
	//	  break;

	++counter;
	add(p);
      }
      f.getline(line);
    }
  }

  // Default constructor - rarely used.
 CBlast_parser():query_name_quality_sorted(false)
  {}


  ~CBlast_parser()
  {
    std::list<CBlastHit*>::iterator it, it_end;
    it     = lbh.begin();
    it_end = lbh.end();

    while (it != it_end)
    {
      //      (**it).print(std::cout);
      delete *it;
      ++it;
    }
  }

  // See also keep_only_n_top_hits
  // Copy constructor that can be used to filter the n first blast hits for a query.
  // Will be sorted by query name, since the source will be sorted.
 CBlast_parser(CBlast_parser &b, unsigned n=UINT_MAX):query_name_quality_sorted(true)
  {
    b.sort_by_query_name_quality();

    std::list<CBlastHit*>::iterator it, it_end;
    it     = b.lbh.begin();
    it_end = b.lbh.end();

    faststring last_qname;
    unsigned count = 0;

    while (it != it_end)
    {
      if ((**it).get_query_seq_name() == last_qname)
      {
	++count;
      }
      else
      {
	count = 1;
	last_qname = (**it).get_query_seq_name();
      }

      CBlastHit* tmp_copy;
      if (count <= n)
      {
	tmp_copy = new CBlastHit(**it);
	add(tmp_copy);
      }

      ++it;
    }
  }

  // Copy constructor which filters hits that have a specific taxonomy
  // of the target.
  CBlast_parser(CBlast_parser &b, CTaxonomy *tax, unsigned type,
		signed char (*func)(CTaxonomy *, unsigned, unsigned)):query_name_quality_sorted(b.query_name_quality_sorted)
  {
    std::list<CBlastHit*>::iterator it, it_end;
    it     = b.lbh.begin();
    it_end = b.lbh.end();

    unsigned    target_gi;
    faststring  tmp;
    CBlastHit   *tmp_copy;
    signed char res;

    while (it != it_end)
    {
      (**it).extract_GI_number_target(tmp);
      target_gi = tmp.ToUnsigned();

      res = func(tax,target_gi, type);
      if (res < 0)
      {
	std::cerr << "Warning: In CBlast_parser. res < 0.\n";
	std::cerr << "GI number seems to be unknown: String: " << tmp << " unsigned: " << target_gi << '\n';
	//	exit(0);
      }
      
      if (res > 0)
      {
	tmp_copy = new CBlastHit(**it);
	add(tmp_copy);
      }
      ++it;
    }
  }

  // Copy constructor which filters hits that have a match in the title field
  CBlast_parser(CBlast_parser &b, CTaxonomy *tax, std::vector<faststring> &query_strings,
		signed char (*func)(CTaxonomy *, unsigned, std::vector<faststring> &query_strings)):query_name_quality_sorted(b.query_name_quality_sorted)
  {
    std::list<CBlastHit*>::iterator it, it_end;
    it     = b.lbh.begin();
    it_end = b.lbh.end();

    unsigned     target_gi;
    faststring   tmp;
    CBlastHit    *tmp_copy;
    signed char  res;

    while (it != it_end)
    {
      (**it).extract_GI_number_target(tmp);
      target_gi = tmp.ToUnsigned();

      res = func(tax,target_gi, query_strings);
      if (res < 0)
      {
	std::cerr << "Error in CBlast_parser. res < 0.\n\n";
	exit(0);
      }
      
      if (res > 0)
      {
	tmp_copy = new CBlastHit(**it);
	add(tmp_copy);
      }
      ++it;
    }
  }

  void clear()
  {
    lbh.clear();
    query_name_quality_sorted = false;
  }


  unsigned get_hit_count()
  {
    return lbh.size();
  }


  void keep_only_n_top_hits(unsigned n=1)
  {
    sort_by_query_name_quality();
    //    lbh.sort(less_than_CBlastHit_pointer_query_quality);

    std::list<CBlastHit*>::iterator it, it2, it_end;
    it     = lbh.begin();
    it_end = lbh.end();

    faststring last_name;
    unsigned i=1;

    while (it != it_end)
    {
      if ( (**it).get_query_seq_name() == last_name )
      {
	++i; // Count this hit
	if (i>n)
	{
	  // Remove *it:
	  delete *it;
	  it = lbh.erase(it);
	} // it now points to the hit behind the hit that has been removed.
	  // it should not be incremented any more.
      }
      else
      {
	last_name = (**it).get_query_seq_name();
	i=1;
	++it;
      }
    }
  }

  // Remove all hits for func does not return true when given gi:
  void keep_only_hits_for_which_func_returns_true_for_parameter_gi(CTaxonomy *tax,
								   bool (*func)(CTaxonomy *, unsigned gi),
								   const char *s)
  {
    std::list<CBlastHit*>::iterator it, it2, it_end;
    it     = lbh.begin();
    it_end = lbh.end();

    unsigned gi;
    faststring tmp;
    unsigned count = 0;

    while (it != it_end)
    {
      (**it).extract_GI_number_target(tmp);
      gi = tmp.ToUnsigned();

      if (!func(tax, gi))
      {
	std::cerr << s << (**it).get_target_seq_name() << " Gi: " << gi << ". Blast hit will be removed.\n";
	++count;
	
	// Remove *it:
	delete *it;
	it = lbh.erase(it);
      }
      else
      {
	++it;
      }
    }

    std::cerr << "WARNING: Number of blast hits removed since no taxid could be found for gi: "
	      << count
	      << '\n'; 
  }


  void keep_only_e_select(double max_e, double factor_max_decrease)
  {
    // Sort by query name, then by e-value.
    // Filter all those with e < max_e and e < best-hit-this-query * factor_max_decrease

    sort_by_query_name_quality();
    //    lbh.sort(less_than_CBlastHit_pointer_query_quality);

    std::list<CBlastHit*>::iterator it, it_end;
    it     = lbh.begin();
    it_end = lbh.end();

    faststring last_hit_query = "";
    double     best_e_value_this_query;
    bool       erase_this;


    while (it != it_end)
    {
      faststring  this_hit_query   = (*it)->get_query_seq_name();
      double      this_hit_e       = (*it)->get_evalue();
      erase_this                   = false;

      if (this_hit_e < max_e)
      {
	if (this_hit_query == last_hit_query) // Still the same query
	{
	  // 
	  if (this_hit_e > best_e_value_this_query*factor_max_decrease)
	  {
	    erase_this = true;
	  }
	  // else // if "<=" then keep
	  
	  if (this_hit_e < best_e_value_this_query)
	  {
	    // Should never happen since this is sorted.
	    std::cerr << "Error better Evalue found: ";
	    (*it)->print(std::cerr);
	    std::cerr << '\n';
	  }
	}
	else  // A new query sequence
	{
	  // Must be best hit this query
	  last_hit_query          = this_hit_query;
	  best_e_value_this_query = this_hit_e; 
	  // keep
	}
      } // if (this_hit_e < max_e)
      else // erase
      {
	erase_this = true;
      }

      // Simple logic to remove this hit.
      if (erase_this)
      {
	  delete *it;
	  it = lbh.erase(it);
      }
      else
      {
	++it;
      }
    }
  }





  void print_all(std::ostream &os)
  {
    std::list<CBlastHit*>::iterator it, it_end;
    it     = lbh.begin();
    it_end = lbh.end();

    while (it != it_end)
    {
      (*it)->print(os);
      ++it;
    }
  }

  void print_all(std::ostream &os,
	     bool pinternal_number,
	     bool pqueryname,
	     bool ptargetname,
	     bool pperfection,
	     bool plength,
	     bool pmismatches,
	     bool pgaps,
	     bool pstart_query,
	     bool pend_query,
	     bool pstart_target, 
	     bool pend_target,
	     bool pe_value, 
	     bool pscore,
	     char delim = '\t')
  {
    std::list<CBlastHit*>::iterator it, it_end;
    it     = lbh.begin();
    it_end = lbh.end();

    while (it != it_end)
    {
      (*it)->print(os, 
		   pinternal_number,
		   pqueryname,
		   ptargetname,
		   pperfection,
		   plength,
		   pmismatches,
		   pgaps,
		   pstart_query,
		   pend_query,
		   pstart_target, 
		   pend_target,
		   pe_value, 
		   pscore,
		   delim = '\t');
      ++it;
    }
  }

  void print_with_title(std::ostream &os, const std::map<unsigned, faststring> &gi2title_map)
  {
    std::list<CBlastHit*>::iterator it, it_end;
    faststring  tmp;
    unsigned    target_gi;

    it     = lbh.begin();
    it_end = lbh.end();

    while (it != it_end)
    {
      (*it)->print_essential_1(os);
      os << std::endl;
      (**it).extract_GI_number_target(tmp);
      target_gi = tmp.ToUnsigned();

      const_iterator_U_Faststring_map find_it;
      find_it = gi2title_map.find(target_gi);
      
      if (find_it != gi2title_map.end() )
      {
	os << find_it->second << std::endl;
      }
      else
      {
	os << "*** No title available. for gi" << target_gi << " ***" << std::endl;
      }
      ++it;	
    }


  }

  void sort_by_query_name_quality()
  {
    //    std::cerr << "RESORT" << std::endl;
    if (!query_name_quality_sorted)
    {
      lbh.sort(less_than_CBlastHit_pointer_query_quality);
      query_name_quality_sorted = true;
    }
  }

  // Calling this function is discouraged. We recommend to call:
  // sort_by_query_name_quality(), since this is more efficient
  // in the case of multiple calls.
  void sort_query()
  {
    lbh.sort(less_than_CBlastHit_pointer_query);
    query_name_quality_sorted = false;
  }


  void sort_target()
  {
    lbh.sort(less_than_CBlastHit_pointer_target);
    query_name_quality_sorted = false;
  }


  std::list<CBlastHit*>::iterator begin()
  {
    return lbh.begin();
  } 


  std::list<CBlastHit*>::iterator end()
  {
    return lbh.end();
  }

  /*
  void erase(std::list<CBlastHit*>::iterator it)
  {
    lbh.erase(it);
  }
  */

  void filter_CSeqNameList_with_major_blast_hit(CSeqNameList &orig, CSeqNameList &result_with,
						double min_query_prop, double min_perfection, double max_evalue,
						bool min_query_prop_with_respect_to_all_hits_this_query,
						const char * accepted_hits_Log_file_name=NULL,
						CSeqNameList *target_list=NULL, const char * search_expr=NULL)
  {
    if ( orig.is_only_names() )
    {
      std::cerr << "Object orig was created with the option <only_names>." << std::endl;
      std::cerr << "I do not know how to proceed." << std::endl;
      exit(0);
    }



    if (min_query_prop_with_respect_to_all_hits_this_query)
      compute_coverage_map(min_perfection, max_evalue);

    bool check_target_seq_match = target_list!= NULL && search_expr != NULL;

    std::list<CBlastHit*>::iterator it, it_end;

    sort_by_query_name_quality();
    //    sort_query(); // This guarantees that hits are ordered by query name.
                        // This will be used here.

    it     = lbh.begin();
    it_end = lbh.end();

    // List of seq names:
    unsigned id_of_hit;

    //    std::map<faststring *, unsigned, less_than_pointer_to_faststring_struct>::iterator it_find;

    CBlastHit* tmp;

    faststring sname; // = tmp->get_query_seq_name();
    faststring lastname;

    FILE *accepted_hits_log;
    if (accepted_hits_Log_file_name != NULL)
    {
      accepted_hits_log = fopen(accepted_hits_Log_file_name, "w");
    }

    // We move through all Blast hits:
    while (it != it_end)
    {
      tmp = *it;

      unsigned seq_len;
      unsigned hit_len;
      double   hit_coverage;
      double   hit_perfection;
      double   hit_evalue;

      // Get name of query of blast hit:
      sname = tmp->get_query_seq_name();
      
      //      std::cerr << "Working on hit with sname " << sname << std::endl;

      if (sname != lastname) // We never want to add the same sequence twice.
	                     // If it has been added (see below) we won't check hits with the same query name
	                     // If the last hit has not been added lastname
	                     // should not have been updated to this name.
      {
	// Look up the name of the query in the list of sequences in the fasta file.
	// Remember: blast tells us the short name only.
	id_of_hit = orig.get_id_short_name(sname); 
 
	//	std::cerr << "sname " << sname << " id-of-hit "  << id_of_hit << std::endl;
	
	if (id_of_hit == UINT_MAX)
	{
	  std::cerr << "Error: Could not find sequence name: " 
		    << tmp->get_query_seq_name()
		    << " for which a blast hit has been found in fasta file. I guess the files are incompatible or target and query are interchanged!!"
		    << std::endl;
	}
	else // Sequence of blast hit has been found
	{
	  seq_len        = orig.get_seq_length(id_of_hit);

	  //	  std::cerr << "sname " << sname << " "  << id_of_hit << " " << seq_len << std::endl;

	  if (min_query_prop_with_respect_to_all_hits_this_query)
	  {
	    std::map<faststring, unsigned>::iterator it_find;
	    it_find = total_coverage_map.find(sname);
	    if (it_find != total_coverage_map.end())
	      hit_len = it_find->second;
	    else
	    {
	      std::cerr << "Warning: Could not find sequence "<< sname << " in total_coverage_map." << std::endl;
	      hit_len = 0;
	    }
	  }
	  else
	  {
	    //	    hit_len        = tmp->get_length();
	    // length might be wrong for protein / nucleotide balsts
	    hit_len = tmp->get_end_query() - tmp->get_start_query() +1;  // We need to add 1 since the end is the last index in hit range.
	  }
	  hit_perfection = tmp->get_perfection();
	  hit_evalue     = tmp->get_evalue();
	  
	  hit_coverage = (double)hit_len/seq_len;

	  //	  
	  //	  std::cerr << "seq-len hit_len coverage "<< sname << " " << seq_len << " " << hit_len  << " " << hit_coverage << std::endl;
	  //

	  if (hit_coverage >= min_query_prop && hit_perfection >= min_perfection && hit_evalue <= max_evalue)
	  {
	    // Query test OK!!
	    if (check_target_seq_match) // Do we check the target seq. full name for a match??
	    {
	      unsigned  t_id;
	      faststring target_short_name = tmp->get_target_seq_name();
	      t_id = target_list->get_id_short_name(target_short_name);
	      faststring target_full_name = target_list->get_name(t_id);

	      // Strange effect: hit_coverage can 1.01 - is this a bug?
	      if (hit_coverage > 1)
	      {
		std::cerr << "Detected hit_coverage>1: " << hit_coverage << "#" << hit_len << "/" << seq_len << "#" << sname << "#" << target_full_name << std::endl;
	      }
	      faststring::size_t find_pos = target_full_name.find(search_expr);
	      if (find_pos != faststring::npos) // search_expr has been found - we have a match
	      {
		// This is an accepted hit:

		if (accepted_hits_Log_file_name != NULL)
		{
		  fprintf(accepted_hits_log, "Accepted hit (qT and tT passed):\t%s\t%s\thc:\t%.2lf\tperf:\t%.3lf\te-val:\t%8lf\n",
			  sname.c_str(), target_full_name.c_str(), hit_coverage, hit_perfection, hit_evalue);
		}
		result_with.add_non_redundant(orig,id_of_hit);
		lastname = sname;  // In case we add the hit, we will ignore other blast hits with the same query name.
		                   // In other words: A query with one or more hits is only inserted once.
		//	  std::cerr << "Would add: " << id_of_hit << std::endl;
	      }
	      else
	      {
		if (accepted_hits_Log_file_name != NULL)
		{
		  fprintf(accepted_hits_log, "Hit not accepted (qT OK but tT failed): %s hc: %.2lf perf: %.3lf e-val: %8lf\n",
			  sname.c_str(), hit_coverage, hit_perfection, hit_evalue);
		}
		
		std::cerr << "--Hit not accepted (qery-test OK, target does not match). Target seq: " << target_full_name << " ";
		tmp->print(std::cerr);
	      }
	    }
	    else // Do not check target seq full name - so we have an accepted hit
	    {
	      //	      unsigned  t_id;
	      faststring target_short_name = tmp->get_target_seq_name();

	      // This is an accepted hit:
	      result_with.add_non_redundant(orig, id_of_hit);
	      
	      if (accepted_hits_Log_file_name != NULL)
	      {
		fprintf(accepted_hits_log, "Accepted hit (qT OK, no tT):\t%s\t%s\thc:\t%.2lf\tperf:\t%.3lf\te-val:\t%8lf\n",
//		fprintf(accepted_hits_log, "Accepted hit (qT OK, no tT):  %s      hc:  %.2lf perf:   %.3lf e-val:   %8lf\n",
		sname.c_str(), target_short_name.c_str(), hit_coverage, hit_perfection, hit_evalue);
	      }

	      lastname = sname;  // In case we add the hit, we will ignore other blast hits with the same query name.
	      //	  std::cerr << "Would add: " << id_of_hit << std::endl;
	    }
	  }
	  else // Does not pass the query test. coverage, perfection, Evalue
	  {
	    if (accepted_hits_Log_file_name != NULL)
	    {
	      fprintf(accepted_hits_log, "Hit not accepted (query-test not passed): %s hc: %.2lf perf: %.3lf e-val: %8lf\n",
		      sname.c_str(), hit_coverage, hit_perfection, hit_evalue);
	    }

#ifdef FINANA
	    std::cerr << "--Hit not accepted (qery-test failed). ";
	    tmp->print(std::cerr);
#endif
	  }
	}
	
      } // End if (sname != lastname)
	
      ++it;
    } // END while (it != it_end)

    if (accepted_hits_Log_file_name != NULL)
    {
      fclose(accepted_hits_log);
    }
  }


  // Typical application:
  // For a given query, the best hit has an Evalue of e-12,
  // the second best hit has an Evalue of e-10.
  // This indicates two independent good hits for the same query.
  // These hits shall be saved in the double_hits data structure.

  // Mode 0:
  //     First  hit has fixed Evalue threshold.
  //     Second hit has a floating threshold depending on the the Evalue of the first hit.
  // Mode 1:
  //     First and second hits have fixed Evalue thresholds.
  //     This mode makes more sense. Example: First hit: E=10^-90
  //     Second hits needs e.g. 10^-87, while 10^-10  is already more than significant.

  // Depricated function name detect_2_target_hits
  void detect_2_target_hits(double max_e_1, double max_e_2, CBlast_parser &double_hits, int mode=1)
  {
    detect_2_good_hits_per_query(max_e_1, max_e_2, double_hits, mode);
  }

  void detect_2_good_hits_per_query(double max_e_1, double max_e_2, CBlast_parser &double_hits, int mode=1)
  {
    // Sort by query name, then by e-value.
    // Filter all those with e < max_e and e < best-hit-this-query * factor_max_decrease

    sort_by_query_name_quality();
    //    lbh.sort(less_than_CBlastHit_pointer_query_quality);

    std::list<CBlastHit*>::iterator it, it_end, tmp;
    it     = lbh.begin();
    it_end = lbh.end();

    faststring last_hit_query = "";
    double     best_e_value_this_query = 100; // Initialised to avoid a warning. The code below should always assign a value to it before it is used.
                                              // This should be the case, since it is only used if (this_hit_query == last_hit_query) which should only be
                                              // the case if it has been assigned a value. 
    bool       still_increment_it;

    // mode 0:
    // First hit threshold:  max_e_1
    // Second hit threshold: factor_max_decrease*best_e_value_this_query (factor_max_decrease := max_e_2)

    if (mode == 0) 
    {
      while (it != it_end)
      {
	faststring  this_hit_query      = (*it)->get_query_seq_name();
	double      this_hit_e          = (*it)->get_evalue();
	double      factor_max_decrease = max_e_2;

	still_increment_it = true;
	
	if (this_hit_query == last_hit_query) // Still the same query
	{
	  // 
	  if (best_e_value_this_query < max_e_1 && 
	      this_hit_e < best_e_value_this_query*factor_max_decrease)
	  {
	    double_hits.add(*it);
	    it = lbh.erase(it);  // This line does an implicit ++it
	    still_increment_it = false;
	  }

	  if (this_hit_e < best_e_value_this_query)
	  {
	    // Should never happen since this list is sorted.
	    std::cerr << "Error better Evalue found: ";
	    (*it)->print(std::cerr);
	    std::cerr << std::endl;
	    exit(-1);
	  }
	}
	else // New query. We store the new reference data.
	{
	  // Must be best hit this query
	  last_hit_query          = this_hit_query;
	  best_e_value_this_query = this_hit_e; 
	}
	if (still_increment_it)
	  ++it;
      }
    }
    // mode 1:
    // First hit threshold:  max_e_1
    // Second hit threshold: max_e_2

    else if (mode == 1)
    {

      while (it != it_end)
      {
	faststring  this_hit_query   = (*it)->get_query_seq_name();
	double      this_hit_e       = (*it)->get_evalue();
	
	still_increment_it = true;
	
	if (this_hit_query == last_hit_query) // Still the same query
	{
	  // 
	  if ( best_e_value_this_query < max_e_1 && 
	       this_hit_e              < max_e_2   )
	  {
	    double_hits.add(*it);
	    it = lbh.erase(it);  // This line does an implicit ++it
	    still_increment_it = false;
	  }

	  if (this_hit_e < best_e_value_this_query)
	  {
	    // Should never happen since this list is sorted.
	    std::cerr << "Error better Evalue found: ";
	    (*it)->print(std::cerr);
	    std::cerr << std::endl;
	    exit(-1);
	  }
	}
	else // New query. We store the new reference data.
	{
	  // Must be best hit this query
	  last_hit_query          = this_hit_query;
	  best_e_value_this_query = this_hit_e; 
	}
	if (still_increment_it)
	  ++it;
      }
    }
    else
    {
      std::cerr << "Intern error: Unknown mode in detect_2_target_hits" << std::endl;
      exit(-2);
    }
  }


  // The reverse mode needs to be revised. CSeqNameList does not allow to add redundant elements since
  // redundant entries are not handled in a save way at the moment.


  /*
  void filter_CSeqNameList_with_major_reverse_blast_hit(CSeqNameList &orig, CSeqNameList &result_with,
							double min_query_prop, double min_perfection, double max_evalue, bool multi=true)
  // The multi paramter determines whether multiple entries are made if the target sequence has multiple hits.
  // Multiple hits of the query are often not of interest, but multiple hits to the same target should usually be considered as independent.
  {
    std::cerr << "Called filter_CSeqNameList_with_major_reverse_blast_hit and multi=" << multi << std::endl; 

    //    compute_coverage_map();

    std::list<CBlastHit*>::iterator it, it_end;

    sort_target(); // This guarantees that hits are ordered by query name.
                   // This will be used here.

    it     = lbh.begin();
    it_end = lbh.end();

    // List of seq names:
    unsigned id_of_hit;

    //    std::map<faststring *, unsigned, less_than_pointer_to_faststring_struct>::iterator it_find;

    CBlastHit* tmp;

    faststring sname;
    faststring lastname;

    // We move through all Blast hits:
    while (it != it_end)
    {
      tmp = *it;

      unsigned seq_len;
      unsigned hit_len;
      double   hit_coverage;
      double   hit_perfection;
      double   hit_evalue;

      // Get name of target:
      sname = tmp->get_target_seq_name();

      //      std::cerr << "Target sequence name: " << sname << std::endl;

      if (multi || sname != lastname) // In contrast to the forward direction we usually want to keep multiple hits against the target.
	                              // If it has been added (see below) we won't check hits with the same query name
	                              // Here we need that the sequences are ordered - hits to the same sequence are successors.
      {
	id_of_hit = orig.get_id_short_name(sname); // Look up the name in the list of sequences in the fasta file.
 
	//      std::cerr << "sname " << sname << " "  << id_of_hit << std::endl;
	
	if (id_of_hit == UINT_MAX)
	{
	  std::cerr << "Error: Could not find sequence name: " 
		    << sname
		    << " for which a blast hit (target) has been found in fasta file. I guess the files are incompatible or target and query are interchanged!!"
		    << std::endl;
	}
	else // Sequence of blast hit has been found
	{
	  seq_len        = orig.get_seq_length(id_of_hit);

	  //	  hit_len        = tmp->get_length();
	  // length might be wrong for protein / nucleotide balsts

	  // We need to add 1 since the end is the last index in hit range.

	  hit_len = tmp->get_end_target() - tmp->get_start_target() +1;


	  hit_perfection = tmp->get_perfection();
	  hit_evalue     = tmp->get_evalue();
	  
	  hit_coverage = (double)hit_len/seq_len;
	  
	  if (hit_coverage >= min_query_prop && hit_perfection >= min_perfection && hit_evalue <= max_evalue)
	  {
	    // This is an accepted hit:
	    //	    result_with.add_non_redundant(orig,id_of_hit);
	    result_with.add(orig,id_of_hit); // multi=true allows redundant entries.
	    if (!multi)
	      lastname = sname;  // In case we add the hit, we will ignore other blast hits with the same target name.
	    //	  std::cerr << "Would add: " << id_of_hit << std::endl;
	  }
	  else
	  {
	    std::cerr << "--Hit not accepted (target test). ID: " << id_of_hit << std::endl; 
	  }
	}
	
      } // End if (sname != lastname)
	
      ++it;
    } // END while (it != it_end)
      
  } // END filter_CSeqNameList_with_major_reverse_blast_hit
  */

  // Input:  CSeqNameList of all query names
  // Putput: CSeqNameList of all query names that have a blast hit in *this.

  void filter_CSeqNameList_has_simple_hit(CSeqNameList &orig, CSeqNameList &result_with)
  {
    std::list<CBlastHit*>::iterator it, it_end;

    sort_by_query_name_quality();
    //    sort_query(); // This guarantees that hits are ordered by query name.
                  // This simplifies the further algorithm.
                  // 

    it     = lbh.begin();
    it_end = lbh.end();

    // List of seq names:
    unsigned id_of_hit;

    CBlastHit* tmp;

    faststring query_name;
    faststring last_query_name;

    // We move through all Blast hits:
    //   - we fetch the sequence name and copy the data to the new CSeqNameList.

    while (it != it_end)
    {
      tmp = *it;

      // Get name of query of blast hit:
      query_name = tmp->get_query_seq_name();

      if (query_name != last_query_name) // We never want to add the same sequence twice.
      {
	// Look up the name of the query in the list of sequences in the fasta file.
	// Remember: blast tells us the short name only.
	id_of_hit = orig.get_id_short_name(query_name); 
 
	//      std::cerr << "query_name " << query_name << " "  << id_of_hit << std::endl;
	
	if (id_of_hit == UINT_MAX)
	{
	  std::cerr << "Error: Could not find sequence name: " 
		    << query_name
		    << " for which a blast hit has been found in fasta file. I guess the files are incompatible or target and query are interchanged!!"
		    << std::endl;
	}
	else // Sequence of blast hit has been found
	{
	  // This is an accepted hit:
	  result_with.add_non_redundant(orig, id_of_hit);
	  last_query_name = query_name;
	}
	
      } // End if (sname != lastname)
	
      ++it;
    } // END while (it != it_end)
  }



  void compute_coverage_map(double min_perfection, double max_e_value)
  {
    std::cerr << "Started to compute coverage values" << std::endl;

    sort_by_query_name_quality();
    //    sort_query(); // Sort by query name and coordinates

    std::list<CBlastHit*>::iterator it, it2, it_end;
    it     = lbh.begin();
    it_end = lbh.end();
   
    const faststring *current_query, *last_query;

    last_query    = NULL;
    current_query = NULL;

    if (it == it_end)
    {
      return;
    }

    CRangeList rl_current_query;
    Crange     r_all(0, UINT_MAX);
    unsigned   cov;

    // Move through all hits:
    while (it != it_end)
    {
      // Get sequence name of query of this hit.
      current_query = &((**it).get_query_seq_name());

      if (last_query && *current_query != *last_query) // next query - not for first!
      {
	// Save coverage of last query
	cov = rl_current_query.coverage_len(r_all);

	// current_query is the next query so last_query contains the correct name for this:
	std::cerr << *last_query << " cov: " << cov << std::endl;

	total_coverage_map[*last_query] = cov;

	// Clear range list:
	rl_current_query.clear();
      }

      // For all queries, We have to add all hits that meet the criteria to the list of ranges: rl_current_query
      // We add this hit:
      // Remember: The blast result contains the last position in the range so we have to add one.

      // Let us check the criteria before we add the hit:
      if ((**it).get_evalue() <= max_e_value && (**it).get_perfection() >= min_perfection)
	{
	  rl_current_query.add((**it).get_start_query(), (**it).get_end_query()+1 );
	  
	  std::cerr << "Adding to range of " << *current_query << " "
		    << (**it).get_start_query() << " " << (**it).get_end_query()+1 << std::endl;
	}
      else
	{
	  std::cerr << "Not adding to range of " << *current_query << " "
		    << (**it).get_start_query() << " " << (**it).get_end_query()+1 << std::endl;
	}

      ++it;
      last_query = current_query;
    }
    // After exiting we have to save the last query-sequence:
    cov = rl_current_query.coverage_len(r_all);
    std::cerr << *last_query << " cov: " << cov << std::endl;
    total_coverage_map[*last_query] = cov;

    //    std::cerr << "Coverage map computed- " << std::endl;

    std::cerr << "Finished computing coverage values" << std::endl;
  }

  /*
  std::list<CBlastHit*>::iterator blc_lower_bound(faststring &value)
  {
    
    std::list<CBlastHit*>::iterator it, first;
    unsigned count, step;
    count = std::distance(lbh.begin(), lbh.end());
    first = lbh.begin();

    while (count > 0)
    {
        it = first; 
        step = count / 2; 
        std::advance(it, step);
        if ((**it).get_query_seq_name() < value)
	{
            first = ++it; 
            count -= step + 1; 
        }
        else
	{
            count = step;
	}
    }
    return first;
  }

  // 
  // Search all blast hits that have a query name that begins identical to the query_beginning_match_str
  // and determine the maximum hit length in this set.
  // More formally: Search for all hits with a query that is identical to query_beginning_match_str
  // when restricted to the first query_beginning_match_str.size() characters.
  int get_max_hit_query_span__matching_beginning1(faststring &query_beginning_match_str)
  {
    int max_len;

    std::cout << "Entering: get_max_hit_query_span__matching_beginning1 with " << query_beginning_match_str << std::endl;

    sort_by_query_name_quality();

    // Find first hit which start with the search string:
    std::list<CBlastHit*>::iterator it, it_end;
    it      = blc_lower_bound(query_beginning_match_str);
    it_end  = lbh.end();

    if (it == it_end)
    {
      std::cout << "No hit found (p1) for partial name of query: " << query_beginning_match_str
		<< std::endl; 
      return -1;
    }

    // If there is no Blast hit for the query with this name,
    // it can be but in most cases will not be equal to it_end
    // The definitive check is that it points to a hit that matches

    if (fstrncmp(query_beginning_match_str, (**it).get_query_seq_name(),query_beginning_match_str.size() ) != 0)
    {
      std::cout << "No hit found (p2) for partial name of query: " << query_beginning_match_str
		<< std::endl; 
      return -2;
    }

    max_len = (**it).get_query_span();
    ++it;

    while (it != it_end)
    {
      // If the strings to not match over the length of query_beginning_match_str.size(), we break
      if (fstrncmp(query_beginning_match_str, (**it).get_query_seq_name(),query_beginning_match_str.size() ) != 0)
	break;
      // If they are still equal over the length of query_beginning_match_str.size(),
      if ( (int)(**it).get_query_span() > max_len)
	max_len = (int)(**it).get_query_span();
      ++it;
    }

    std::cout << "Found partial name of query: " << query_beginning_match_str
              << " " << max_len << std::endl; 

    return max_len;
  }
  */

  
  int get_max_hit_query_span__matching_beginning2(faststring &query_beginning_match_str, unsigned verbosity)
  {
    int max_len;

    if (verbosity > 100)
      std::cout << "Entering: get_max_hit_query_span__matching_beginning2 with " << query_beginning_match_str << std::endl;

    // Find first hit which start with the search string:
    std::multimap<faststring, CBlastHit*>::iterator it, it_end;
    it      = bh_map.lower_bound(query_beginning_match_str);
    it_end  = bh_map.end();

    if (it == it_end)
    {
      if (verbosity > 100)
      {
	std::cout << "No hit found (p1) for partial name of query: " << query_beginning_match_str
		  << std::endl; 
      }
      return -1;
    }

    // If there is no Blast hit for the query with this name,
    // the iterator it can be but in most cases will not be equal to it_end
    // The definitive check is that it points to a hit that matches

    if (fstrncmp(query_beginning_match_str, (it->second)->get_query_seq_name(),query_beginning_match_str.size() ) != 0)
    {
      if (verbosity > 100)
      {    
	std::cout << "No hit found (p2) for partial name of query: " << query_beginning_match_str
		  << std::endl; 
      }
      return -2;
    }

    max_len = (it->second)->get_query_span();
    ++it;

    while (it != it_end)
    {
      // If the strings to not match over the length of query_beginning_match_str.size(), we break
      if (fstrncmp(query_beginning_match_str, (it->second)->get_query_seq_name(),query_beginning_match_str.size() ) != 0)
	break;
      // If they are still equal over the length of query_beginning_match_str.size(),
      if ( (int)(it->second)->get_query_span() > max_len)
	max_len = (int)(it->second)->get_query_span();
      ++it;
    }

    if (verbosity > 100)
    {
      std::cout << "Found partial name of query: " << query_beginning_match_str
		<< " " << max_len << std::endl; 
    }

    return max_len;
  }
  


  void set_of_gi_numbers_target(std::set<unsigned> &gi_set)
  {
    std::list<CBlastHit*>::iterator it, it_end;
    it     = lbh.begin();
    it_end = lbh.end();

    faststring tmp_gi;

    while (it != it_end)
    {
      (**it).extract_GI_number_target(tmp_gi);
      gi_set.insert(tmp_gi.ToUnsigned());
      ++it;
    }
  }

  void set_of_query_names(std::set<faststring> &s)
  {
    std::list<CBlastHit*>::iterator it, it_end;
    it     = lbh.begin();
    it_end = lbh.end();

    while (it != it_end)
    {
      s.insert( (**it).get_query_seq_name() );
      ++it;
    }
  }

/*   void get_max_query_span_map_for_matching_queries(std::map<faststring, int> m, char stop_char, int num_stop_char) */
/*   { */

/*   } */

  void write_set_of_gi_numbers_to_file(const char *fname)
  {
    std::set<unsigned> gi_set;
    set_of_gi_numbers_target(gi_set);
    std::set<unsigned>::iterator it,it_end;
    it      = gi_set.begin();
    it_end  = gi_set.end();
    
    std::ofstream os(fname);
    while (it != it_end)
    {
      os << *it << std::endl;
      ++it;
    }
    os.close();
  } // END write_set_of_gi_numbers_to_file




  // fill a map with the following information:
  // For all gi numbers of the target sequences, the map contains the number of occurrences of this target.
  void get_gi_occurence_map(std::map<unsigned, unsigned> &gim)
  {
    std::list<CBlastHit*>::iterator it, it_end;
    it     = lbh.begin();
    it_end = lbh.end();

    faststring tmp_gi;

    while (it != it_end)
    {
      (**it).extract_GI_number_target(tmp_gi);
      __add_or_count(gim, tmp_gi.ToUnsigned());
      ++it;
    }
  
  }


};

#endif
