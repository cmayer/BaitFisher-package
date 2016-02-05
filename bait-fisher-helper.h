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
#include "faststring2.h"
#include "CSequences2.h"
#include "CSequence_Mol2_1.h"



faststring extract_taxon_name(faststring, char, unsigned short);
faststring extract_gene_name(faststring, char, unsigned short);
void read_alignment_sequence_of_species2(faststring, faststring,
					 CSequence_Mol*&, CSequences2*&,
					 unsigned&, char, unsigned short);

extern int DEBUG_COUNT_NEW_CALLS_FEATURE_ALIGNMENT;
extern int DEBUG_COUNT_NEW_CALLS_CBAIT_REGION_INFORMATION;
extern int DEBUG_COUNT_NEW_CALLS_CSEQUENCE_LOCI_CLUSTER_COLLECTION;
extern int DEBUG_COUNT_NEW_CALLS_MSA;
extern int DEBUG_COUNT_NEW_ALIGNMENT_FROM_FILE;

extern int DEBUG_COUNT_DELETE_CALLS_FEATURE_ALIGNMENT;
extern int DEBUG_COUNT_DELETE_CALLS_CBAIT_REGION_INFORMATION;
extern int DEBUG_COUNT_DELETE_CALLS_CSEQUENCE_LOCI_CLUSTER_COLLECTION;
extern int DEBUG_COUNT_DELETE_CALLS_MSA;
extern int DEBUG_COUNT_DELETE_ALIGNMENT_FROM_FILE;


