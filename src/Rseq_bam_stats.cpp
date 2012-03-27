/*****************************************************************************
  
  Bam_stats_RNA-seq.cpp
  1) basic stats of the alignment of RNA-seq reads against ref.genome
  (current mapper: GSNAP)
  this will only work for a bam file or files sorted by read name
  2) also calculate paired-end mapping to generate "best" alignment of a fragment
  (based on flag containing "Primary" Record), write the uniquely mapped reads into a new bam file
  
  solely based on paired-end sequencing (fragment level)

  (c) 2011 - Sun Ruping
  Dept. Vingron (Computational Mol. Bio.)
  Max-Planck-Institute for Molecular Genetics
  Ihnestr. 73, D-14195, Berlin, Germany   

  ruping@molgen.mpg.de

g++ 
-I /scratch/ngsvin2/RNA-seq/ruping/Tools/bamtools/include/ 
   /bips/lib/packages/NGS_tools/bamtools-1.0.2/include/
-L /scratch/ngsvin2/RNA-seq/ruping/Tools/bamtools/lib/
   /bips/lib/packages/NGS_tools/bamtools-1.0.2/lib/ 
-lbamtools -lz -static
******************************************************************************/

#include <api/BamReader.h>
#include <api/BamWriter.h>
#include <api/BamMultiReader.h>

using namespace BamTools;

#include <iostream>
#include <fstream>
#include <cstdlib>
#include <vector>
#include <deque>
#include <set>
#include <string>
#include <cstring>
#include <sstream>
#include "Rseq_bam_stats.h"
using namespace std;


struct RseqSTATS {
  unsigned int num_Reads;
  unsigned int num_Duplicates;
  unsigned int num_FailedQC;
  unsigned int num_Mapped;
  unsigned int num_Unique;
  unsigned int num_spliced;
  unsigned int num_Singletons;
  unsigned int num_ProperPair;
  unsigned int num_WrongPair;
  unsigned int num_Multi;
  unsigned int num_Unmapped;
  unsigned int num_UniqueHalf;
};

struct Alignment {
  string chr1;
  unsigned int start1;
  unsigned int end1;
  string chr2;
  unsigned int start2;
  unsigned int end2;
  unsigned int cate; // 1:unmapped; 2:multi; 3:singleton; 4:unique;  (UN 5-9) 5:one_end_mapped; 6:1uniq; 7:2uniq; 8:1multi; 9:2multi; 10:saved_one_end_unique
  bool junction;
  BamAlignment mate1;
  BamAlignment mate2;
};


inline void ParseCigar(const vector<CigarOp> &cigar, vector<int> &blockStarts, vector<int> &blockEnds, unsigned int &alignmentEnd);
inline void splitstring(const string &str, vector<string> &elements, const string &delimiter);
inline string int2str(unsigned int &i);
inline void print_stats(struct RseqSTATS &rstats);

int main (int argc, char *argv[]) {
 
  struct parameters *param = 0;
  param = interface(param, argc, argv);

  //-------------------------------------------------------------------------------------------------------+
  // BAM input (file or filenames?)                                                                        |
  //-------------------------------------------------------------------------------------------------------+
  char *fof = param->mapping_f;
  FILE *IN=NULL;
  char line[5000];
  int filecount=0;
  vector <string> fnames;

  if (strchr(fof,' ')!=NULL) {
    char *ptr;
    ptr=strtok(fof," ");
    while (ptr!=NULL) {
      fnames.push_back(ptr);
      filecount++;
      ptr=strtok(NULL," ");
    }
  } else {
    IN=fopen(fof,"rt");
    if (IN!=NULL) {
      long linecount=0;
      while (fgets(line,5000-1,IN)!=NULL) {
        linecount++;
        if (line[0]!='#' && line[0]!='\n') {
          char *ptr=strchr(line,'\n');
          if (ptr!=NULL && ptr[0]=='\n') {
            ptr[0]='\0';
          }
          FILE *dummy=NULL;
          dummy=fopen(line,"rt");
          if (dummy!=NULL) {     // seems to be a file of filenames...
            fclose(dummy);
            fnames.push_back(line);
            filecount++;
          } else if (filecount==0 || linecount>=1000-1) {  // seems to be a single file
            fnames.push_back(fof);
            filecount++;
            break;
          }
        }
      }
      fclose(IN);
    }
  } //file or file name decided and stored in vector "fnames"

  cerr << "the input mapping files are:" << endl;
  vector <string>::iterator fit = fnames.begin();
  for(; fit != fnames.end(); fit++){
    cerr << *fit << endl;
  }
  //-------------------------------------------------------------------------------------------------------+
  // end of file or filenames                                                                              |
  //-------------------------------------------------------------------------------------------------------+

  //bam input and generate index if not yet
  BamMultiReader reader;
  reader.Open(fnames);   // the mapping bam file is opened 

  // get header & reference information
  string header = reader.GetHeaderText();
  RefVector refs = reader.GetReferenceData();

  // attempt to open BamWriter
  BamWriter writer;
  if ( !writer.Open(param->writer, header, refs) ) {
    cerr << "Could not open output BAM file" << endl;
    exit(0);
  }

  // attempt to write arp reads
  ofstream arp_f;
  string arp = param->arp;
  if ( arp != "" ) {
    arp_f.open(param->arp);
  }

  // statistics
  struct RseqSTATS BAMSTATS = {0,0,0,0,0,0,0,0,0,0,0,0};

  map <string, struct Alignment> fragment; // map for fragment

  // type == "s" or type == "p" ?
  string type = param->type;
  
  string old_frag = "SRP";

  BamAlignment bam;
  while ( reader.GetNextAlignment(bam) ) {

    BamAlignment cBAM;

    unsigned int unique = 0;
    string XS = "SRP";
    bool jc = false;
    string chrom = "SRP";
    string strand = "+";
    unsigned int alignmentStart = 0;
    unsigned int alignmentEnd = 0;
    //unsigned int mapq = 0;
    
    if ( bam.IsMapped() == true) {
      bam.GetTag("NH", unique);                     // uniqueness
      bam.GetTag("XS", XS);                         // juction reads?
      //if (XS == "+" || XS == "-") jc = true;
      if (XS != "SRP") jc = true;
      chrom  = refs.at(bam.RefID).RefName;          // chromosome
      if (bam.IsReverseStrand()) strand = "-";      // strand -
      alignmentStart = bam.Position+1;              // start
      alignmentEnd   = bam.GetEndPosition();        // end
    }
    unsigned int mate = 1;
    if ( bam.IsFirstMate() == false ) mate = 2;           // second mate

    if ( bam.Name != old_frag ) {  // new frag

      if (old_frag != "SRP") {
        fragment.erase(old_frag);    // remove old_frag
      }

      if ( bam.IsMapped() == false && bam.IsMateMapped() == false ) {  // unmapped
        ++BAMSTATS.num_Reads;
        ++BAMSTATS.num_Unmapped;
        struct Alignment tmp = {"UM", 0, 0, "UM", 0, 0, 1, jc, cBAM, cBAM};
        fragment.insert( pair<string, struct Alignment>(bam.Name, tmp) );
      } // unmapped
      else if ( bam.IsMapped() == false && bam.IsMateMapped() == true ) {  // one end is not mappable
        if (mate == 1){
          struct Alignment tmp = {"UM", 0, 0, "SRP", 0, 0, 5, jc, bam, cBAM};
          fragment.insert( pair<string, struct Alignment>(bam.Name, tmp) );
        }
        else {
          struct Alignment tmp = {"SRP", 0, 0, "UM", 0, 0, 5, jc, cBAM, bam};
          fragment.insert( pair<string, struct Alignment>(bam.Name, tmp) );
        }
      } // undecided one end not mappable
      else if ( bam.IsMapped() == true && bam.IsMateMapped() == false ) {  // one other end is not mappable

        if (unique > 1) {  // one end multiple mapped, the other end not mappable
          ++BAMSTATS.num_Reads;
          ++BAMSTATS.num_Mapped;
          ++BAMSTATS.num_Multi;
          if (mate == 1) {
            struct Alignment tmp = {"MM", 0, 0, "UM", 0, 0, 2, jc, cBAM, cBAM};
            fragment.insert( pair<string, struct Alignment>(bam.Name, tmp) );
          }
          else {
            struct Alignment tmp = {"UM", 0, 0, "MM", 0, 0, 2, jc, cBAM, cBAM};
            fragment.insert( pair<string, struct Alignment>(bam.Name, tmp) );
          }
        } // one end multiple mapped, the other end not mappable
        else { // Singletons (it should be output here, since TopHat does not output the alignment of other mate)
           ++BAMSTATS.num_Reads;
           ++BAMSTATS.num_Mapped;
           ++BAMSTATS.num_Unique;
           ++BAMSTATS.num_Singletons;
           ++BAMSTATS.num_WrongPair;
           if (jc == true) ++BAMSTATS.num_spliced;
           if (mate == 1) {
             struct Alignment tmp = {chrom, alignmentStart, alignmentEnd, "UM", 0, 0, 3, jc, bam, cBAM};
             fragment.insert( pair<string, struct Alignment>(bam.Name, tmp) );
           }
           else {
             struct Alignment tmp = {"UM", 0, 0, chrom, alignmentStart, alignmentEnd, 3, jc, cBAM, bam};
             fragment.insert( pair<string, struct Alignment>(bam.Name, tmp) );
           }
           writer.SaveAlignment(bam);                        // write
           if ( arp != "" ) arp_f << bam.Name << endl;       // write arp
        } // Singletons

      } //one another end is not mappable

      else {  // both ends mapped
        if ( unique == 1 ) {  // current end is uniquelly mapped
          if (mate == 1) {
            struct Alignment tmp = {chrom, alignmentStart, alignmentEnd, refs.at(bam.MateRefID).RefName, (bam.MatePosition+1), 0, 6, jc, bam, cBAM};
            fragment.insert( pair<string, struct Alignment>(bam.Name, tmp) );
          }
          else {
            struct Alignment tmp = {refs.at(bam.MateRefID).RefName, (bam.MatePosition+1), 0, chrom, alignmentStart, alignmentEnd, 7, jc, cBAM, bam};
            fragment.insert( pair<string, struct Alignment>(bam.Name, tmp) );
          }
        } // current unique
        else {  // current end is not unique
          if (mate == 1) {
            if ( bam.IsPrimaryAlignment() == true ){
               struct Alignment tmp = {chrom, alignmentStart, alignmentEnd, refs.at(bam.MateRefID).RefName, (bam.MatePosition+1), 0, 8, jc, bam, cBAM};
               fragment.insert( pair<string, struct Alignment>(bam.Name, tmp) );
            }
            else {
               struct Alignment tmp = {"SRP", 0, 0, "SRP", 0, 0, 8, jc, cBAM, cBAM};
               fragment.insert( pair<string, struct Alignment>(bam.Name, tmp) );
            }
          }
          else { // mate 2
            if ( bam.IsPrimaryAlignment() == true ){
               struct Alignment tmp = {refs.at(bam.MateRefID).RefName, (bam.MatePosition+1), 0, chrom, alignmentStart, alignmentEnd, 9, jc, cBAM, bam};
               fragment.insert( pair<string, struct Alignment>(bam.Name, tmp) );
            }
            else {
               struct Alignment tmp = {"SRP", 0, 0, "SRP", 0, 0, 9, jc, cBAM, cBAM};
               fragment.insert( pair<string, struct Alignment>(bam.Name, tmp) );
            }
          }
        } // current multi
      } // both ends mapped;

      old_frag = bam.Name; // reset old frag

    } // a new frag;

    else {  // IT IS AN OLD FRAGMENT!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if (fragment[bam.Name].cate == 5) { // one end mapped the other end not, but the mapped end is not decided;
         if ( unique > 1 ) {  // one end multiple mapped, the other end not mappable
          ++BAMSTATS.num_Reads;
          ++BAMSTATS.num_Mapped;
          ++BAMSTATS.num_Multi;
          if (mate == 1){
            fragment[bam.Name].chr1 = "MM";          
          }
          else {
            fragment[bam.Name].chr2 = "MM";
          }
          fragment[bam.Name].cate = 2;
        } // one end multiple mapped, the other end not mappable
        else { // Singletons
           ++BAMSTATS.num_Reads;
           ++BAMSTATS.num_Mapped;
           ++BAMSTATS.num_Unique;
           ++BAMSTATS.num_Singletons;
           ++BAMSTATS.num_WrongPair;
           if (jc == true) ++BAMSTATS.num_spliced;
           if (mate == 1){
             fragment[bam.Name].chr1 = chrom;
             fragment[bam.Name].start1 = alignmentStart;
             fragment[bam.Name].end1 = alignmentEnd;
             fragment[bam.Name].junction = jc;
             fragment[bam.Name].mate1 = bam;
             writer.SaveAlignment(bam);                        // write
             writer.SaveAlignment(fragment[bam.Name].mate2);   // write
             if ( arp != "" ) arp_f << bam.Name << endl;       // write arp
           }
           else {
             fragment[bam.Name].chr2 = chrom;
             fragment[bam.Name].start2 = alignmentStart;
             fragment[bam.Name].end2 = alignmentEnd;
             fragment[bam.Name].junction = jc;
             fragment[bam.Name].mate2 = bam;
             writer.SaveAlignment(fragment[bam.Name].mate1);   // write
             writer.SaveAlignment(bam);                        // write
             if ( arp != "" ) arp_f << bam.Name << endl;       // write arp
           }
           fragment[bam.Name].cate = 3;
        } // Singletons
      } //cate == 5

      else if (fragment[bam.Name].cate == 6) { // mate 1 is unique
        if (mate == 1) {cerr << "mate1 unique inconsistency, exit\n"; exit(0);} 
        if (unique == 1) { // both ends are unique VERY GOOD
           ++BAMSTATS.num_Reads;
           ++BAMSTATS.num_Mapped;
           ++BAMSTATS.num_Unique;
           if ( bam.IsDuplicate() ) ++BAMSTATS.num_Duplicates;
           if (  bam.IsFailedQC() ) ++BAMSTATS.num_FailedQC;
           if ( bam.IsProperPair()) ++BAMSTATS.num_ProperPair;
           else                     ++BAMSTATS.num_WrongPair;
           if (jc == true || fragment[bam.Name].junction == true) ++BAMSTATS.num_spliced;
           fragment[bam.Name].end2 = alignmentEnd;
           fragment[bam.Name].cate = 4;  //unique
           fragment[bam.Name].mate2 = bam;
           writer.SaveAlignment(fragment[bam.Name].mate1);     // write
           writer.SaveAlignment(bam);                          // write
           if ( arp != "" ) {
             int dis = fragment[bam.Name].start1 - fragment[bam.Name].start2;
             if ( (fragment[bam.Name].chr1 != fragment[bam.Name].chr2) || (abs(dis) > 100000) ) 
               arp_f << bam.Name << endl;
             //if ( fragment[bam.Name].chr1 != fragment[bam.Name].chr2 ) cerr << bam.Name << endl;
           }

        } // both ends are unique
        else { // the mate 2 is multi, try to figure out the "primary" record
          if ( bam.IsPrimaryAlignment() == true ){ // if this is a primary result
            ++BAMSTATS.num_Reads;
            ++BAMSTATS.num_Mapped;
            ++BAMSTATS.num_UniqueHalf;
            fragment[bam.Name].chr2   = chrom;
            fragment[bam.Name].start2 = alignmentStart;
            fragment[bam.Name].end2   = alignmentEnd;
            fragment[bam.Name].cate   = 10;
            fragment[bam.Name].mate2 = bam;
            writer.SaveAlignment(fragment[bam.Name].mate1);    // write
            writer.SaveAlignment(bam);                         // write
            if ( arp != "" ) {
              int dis = fragment[bam.Name].start1 - fragment[bam.Name].start2;
              if ( (fragment[bam.Name].chr1 != fragment[bam.Name].chr2) || (abs(dis) > 100000) ) 
                arp_f << bam.Name << endl;
              //if ( fragment[bam.Name].chr1 != fragment[bam.Name].chr2 ) cerr << bam.Name << endl;
            }
          }
        }
      } // cate == 6         
      else if (fragment[bam.Name].cate == 7) { // mate 2 is unique
        if (mate == 2) {cerr << "mate2 unique inconsistency, exit\n"; exit(0);}
        if (unique == 1) { // both ends are unique VERY GOOD
          ++BAMSTATS.num_Reads;
          ++BAMSTATS.num_Mapped;
          ++BAMSTATS.num_Unique;
          if ( bam.IsDuplicate() ) ++BAMSTATS.num_Duplicates;
          if (  bam.IsFailedQC() ) ++BAMSTATS.num_FailedQC;
          if ( bam.IsProperPair()) ++BAMSTATS.num_ProperPair;
          else                     ++BAMSTATS.num_WrongPair;
          if (jc == true || fragment[bam.Name].junction == true) ++BAMSTATS.num_spliced;
          fragment[bam.Name].end1 = alignmentEnd;
          fragment[bam.Name].cate = 4;  //unique
          fragment[bam.Name].mate1 = bam;
          writer.SaveAlignment(bam);                           // write
          writer.SaveAlignment(fragment[bam.Name].mate2);      // write
          if ( arp != "" ) {
            int dis = fragment[bam.Name].start1 - fragment[bam.Name].start2;
            if ( (fragment[bam.Name].chr1 != fragment[bam.Name].chr2) || (abs(dis) > 100000) ) 
              arp_f << bam.Name << endl;
            //if ( fragment[bam.Name].chr1 != fragment[bam.Name].chr2 ) cerr << bam.Name << endl;
          }
        } // both ends are unique
        else { // the mate 1 is multi
          if ( bam.IsPrimaryAlignment() == true ) { // if this is a primary result
            ++BAMSTATS.num_Reads;
            ++BAMSTATS.num_Mapped;
            ++BAMSTATS.num_UniqueHalf;
            fragment[bam.Name].chr1   = chrom;
            fragment[bam.Name].start1 = alignmentStart;
            fragment[bam.Name].end1   = alignmentEnd;
            fragment[bam.Name].cate   = 10;
            fragment[bam.Name].mate1 = bam;
            writer.SaveAlignment(bam);                         // write
            writer.SaveAlignment(fragment[bam.Name].mate2);    // write
            if ( arp != "" ) {
              int dis = fragment[bam.Name].start1 - fragment[bam.Name].start2;
              if ( (fragment[bam.Name].chr1 != fragment[bam.Name].chr2) || (abs(dis) > 100000) ) 
                arp_f << bam.Name << endl;
              //if ( fragment[bam.Name].chr1 != fragment[bam.Name].chr2 ) cerr << bam.Name << endl;
            }
          }
        }
      } // cate == 7
      else if (fragment[bam.Name].cate == 8) {   // mate 1 is multi
        if ( fragment[bam.Name].chr1 == "SRP" ){  // meaning it is multi, but the primary is not decided
          if ( mate == 1 ){
            if ( bam.IsPrimaryAlignment() == true ){
              fragment[bam.Name].chr1   = chrom;
              fragment[bam.Name].start1 = alignmentStart;
              fragment[bam.Name].end1   = alignmentEnd;
              fragment[bam.Name].mate1  = bam;
            }
          } 
        }  // the primary is not decided

        if ( mate == 2 ) {     // check mate 2
          if ( unique == 1 ) { // if mate 2 is unique
            ++BAMSTATS.num_Reads;
            ++BAMSTATS.num_Mapped;
            ++BAMSTATS.num_UniqueHalf;
            fragment[bam.Name].chr2   = chrom;
            fragment[bam.Name].start2 = alignmentStart;
            fragment[bam.Name].end2   = alignmentEnd;
            fragment[bam.Name].cate   = 10;
            fragment[bam.Name].mate2  = bam;
            writer.SaveAlignment(fragment[bam.Name].mate1);    // write
            writer.SaveAlignment(bam);                         // write
            if ( arp != "" ) {
              int dis = fragment[bam.Name].start1 - fragment[bam.Name].start2;
              if ( (fragment[bam.Name].chr1 != fragment[bam.Name].chr2) || (abs(dis) > 100000) ) 
                arp_f << bam.Name << endl;
              //if ( fragment[bam.Name].chr1 != fragment[bam.Name].chr2 ) cerr << bam.Name << endl;
            }
          }
          else { // if mate 2 is multi
            ++BAMSTATS.num_Reads;
            ++BAMSTATS.num_Mapped;
            ++BAMSTATS.num_Multi;
            fragment[bam.Name].cate = 2;
          }
        } // check mate 2  
      } // cate == 8
      else if (fragment[bam.Name].cate == 9) { // mate 2 is multi
        if ( fragment[bam.Name].chr2 == "SRP" ) {  // meaning it is multi, but the primary is not decided
          if ( mate == 2 ){
            if ( bam.IsPrimaryAlignment() == true ){
              fragment[bam.Name].chr2   = chrom;
              fragment[bam.Name].start2 = alignmentStart;
              fragment[bam.Name].end2   = alignmentEnd;
              fragment[bam.Name].mate2  = bam;
            }
          } 
        }  // the primary is not decided

        if ( mate == 1 ) {     // check mate 1
          if ( unique == 1 ) { // if mate 1 is unique
            ++BAMSTATS.num_Reads;
            ++BAMSTATS.num_Mapped;
            ++BAMSTATS.num_UniqueHalf;
            fragment[bam.Name].chr1   = chrom;
            fragment[bam.Name].start1 = alignmentStart;
            fragment[bam.Name].end1   = alignmentEnd;
            fragment[bam.Name].cate   = 10;
            fragment[bam.Name].mate1  = bam;
            writer.SaveAlignment(bam);                         // write
            writer.SaveAlignment(fragment[bam.Name].mate2);    // write
            if ( arp != "" ) {
              int dis = fragment[bam.Name].start1 - fragment[bam.Name].start2;
              if ( (fragment[bam.Name].chr1 != fragment[bam.Name].chr2) || (abs(dis) > 100000) ) 
                arp_f << bam.Name << endl;
              //if ( fragment[bam.Name].chr1 != fragment[bam.Name].chr2 ) cerr << bam.Name << endl;
            }
          }
          else { // if mate 1 is multi
            ++BAMSTATS.num_Reads;
            ++BAMSTATS.num_Mapped;
            ++BAMSTATS.num_Multi;
            fragment[bam.Name].cate = 2;
          }
        } // check mate 1 

      } // cate == 9

    } // old fragment

  }  //  read a bam
      
  reader.Close();
  writer.Close();
  arp_f.close();

  print_stats(BAMSTATS);

  return 0;
} //main

inline string int2str(unsigned int &i){
  string s;
  stringstream ss(s);
  ss << i;
  return ss.str();
}

inline void splitstring(const string &str, vector<string> &elements, const string &delimiter) {
  string::size_type lastPos = str.find_first_not_of(delimiter, 0);
  string::size_type pos     = str.find_first_of(delimiter, lastPos);

  while (string::npos != pos || string::npos != lastPos) {
    elements.push_back(str.substr(lastPos, pos - lastPos));
    lastPos = str.find_first_not_of(delimiter, pos);
    pos = str.find_first_of(delimiter, lastPos);
  }
}


inline void ParseCigar(const vector<CigarOp> &cigar, vector<int> &blockStarts, vector<int> &blockLengths, unsigned int &alignmentEnd) {

  int currPosition = 0;
  int blockLength  = 0;

  //  Rip through the CIGAR ops and figure out if there is more
  //  than one block for this alignment
  vector<CigarOp>::const_iterator cigItr = cigar.begin();
  vector<CigarOp>::const_iterator cigEnd = cigar.end();
  for (; cigItr != cigEnd; ++cigItr) {
    switch (cigItr->Type) {
    case ('M') :
      blockLength  += cigItr->Length;
      currPosition += cigItr->Length;
    case ('I') : break;
    case ('S') : break;
    case ('D') : break;
      blockLength  += cigItr->Length;
      currPosition += cigItr->Length;
    case ('P') : break;
    case ('N') :
      blockStarts.push_back(currPosition + cigItr->Length);
      blockLengths.push_back(blockLength);
      currPosition += cigItr->Length;
      blockLength = 0;
    case ('H') : break;                             // for 'H' - do nothing, move to next op
    default    :
      printf("ERROR: Invalid Cigar op type\n");   // shouldn't get here
      exit(1);
    }
  }
  // add the kast block and set the
  // alignment end (i.e., relative to the start)
  blockLengths.push_back(blockLength);
  alignmentEnd = currPosition;
}


inline void print_stats(struct RseqSTATS &rstats){
  cout << "Reads:      " << rstats.num_Reads      << endl;
  cout << "Mapped:     " << rstats.num_Mapped     << endl;
  cout << "Unmapped:   " << rstats.num_Unmapped   << endl;
  cout << "Unique:     " << rstats.num_Unique     << endl;
  cout << "Uniquehalf: " << rstats.num_UniqueHalf << endl;
  cout << "duplicates: " << rstats.num_Duplicates << endl;
  cout << "failed_QC:  " << rstats.num_FailedQC   << endl;
  cout << "singletons: " << rstats.num_Singletons << endl;
  cout << "ProperPair: " << rstats.num_ProperPair << endl;
  cout << "WrongPair:  " << rstats.num_WrongPair  << endl;
  cout << "Spliced:   "  << rstats.num_spliced    << endl;
  cout << "MultiMap:   " << rstats.num_Multi      << endl;
}
