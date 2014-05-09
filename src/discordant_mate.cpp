/*****************************************************************************
  
  get the discordant mate information of paired-end reads

  (c) 2011 - Sun Ruping
  Dept. Vingron (Computational Mol. Bio.)
  Max-Planck-Institute for Molecular Genetics
  Ihnestr. 73, D-14195, Berlin, Germany   

  current: Department of Systems Biology, Columbia University, NY, USA
  rs3412@c2b2.columbia.edu


g++ 
-I /scratch/ngsvin2/RNA-seq/ruping/Tools/bamtools/include/ 
   /bips/lib/packages/NGS_tools/bamtools-1.0.2/include/
-L /scratch/ngsvin2/RNA-seq/ruping/Tools/bamtools/lib/
   /bips/lib/packages/NGS_tools/bamtools-1.0.2/lib/ 
-lbamtools -lz -static
******************************************************************************/

#include <api/BamReader.h>
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
#include "discordant_mate.h"
using namespace std;

struct region {  // a txt file containing breakpoints
  unsigned int id;
  string type;
  string chr;
  unsigned int coor;
  unsigned int support_no;
  unsigned int pw;              //proper or wrong
  unsigned int ct;              //clip type
  string rep;
  unsigned int start;
  unsigned int end;
  map < string, set<unsigned int> > tags;
};

inline void ParseCigar(const vector<CigarOp> &cigar, vector<int> &blockStarts, vector<int> &blockEnds, unsigned int &alignmentEnd, bool &keep);
inline void splitstring(const string &str, vector<string> &elements, const string &delimiter);
inline string int2str(unsigned int &i);
inline string float2str(float &f);
inline void output_processing(struct region &region, unsigned int &consisCountTh);

int main ( int argc, char *argv[] ) { 

  struct parameters *param = 0;
  param = interface(param, argc, argv);

  unsigned int idstart = (param->idstart);

  //bam input and generate index if not yet
  //-------------------------------------------------------------------------------------------------------+
  // BAM input (file or filenames?)                                                                        |
  //-------------------------------------------------------------------------------------------------------+
  char *fof = param->mapping_f;
  FILE *IN=NULL;
  char linefof[5000];
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
      while (fgets(linefof,5000-1,IN)!=NULL) {
        linecount++;
        if (linefof[0]!='#' && linefof[0]!='\n') {
          char *ptr=strchr(linefof,'\n');
          if (ptr!=NULL && ptr[0]=='\n') {
            ptr[0]='\0';
          }
          FILE *dummy=NULL;
          dummy=fopen(linefof,"rt");
          if (dummy!=NULL) {     // seems to be a file of filenames...
            fclose(dummy);
            fnames.push_back(linefof);
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
  }  //file or file name decided and stored in vector "fnames"

  cerr << "the input mapping files are:" << endl;
  vector <string>::iterator fit = fnames.begin();
  for(; fit != fnames.end(); fit++) {
    cerr << *fit << endl;
  }

  //-------------------------------------------------------------------------------------------------------+
  // end of file or filenames                                                                              |
  //-------------------------------------------------------------------------------------------------------+

  // open the BAM file(s)
  BamMultiReader reader;
  reader.Open(fnames);

  // get header & reference information
  string header = reader.GetHeaderText();
  RefVector refs = reader.GetReferenceData();

  if ( ! reader.LocateIndexes() )     // opens any existing index files that match our BAM files
    reader.CreateIndexes();         // creates index files for BAM files that still lack one


  string type = param->type;

  unsigned int consisCountTh = param->consisCount;
  if (consisCountTh == 0){
    cerr << "error: consistent Count threshold not set!" << endl;
    exit(1);
  } else {
    cerr << "consistent Count threshold: " << consisCountTh << endl;
  }

  struct region discoRegion;  //for discordant calculation
  unsigned int discoID = idstart;

  string oldChrom = "SRP";

  BamAlignment bam;

  discoRegion.id = 0;
  
  while (reader.GetNextAlignment(bam)) {

    if ( bam.IsMapped() == false ) continue;   // skip unaligned reads

    unsigned int unique;
    bam.GetTag("NH", unique);
    if (param->unique == 1) {
      if (unique != 1) {                         // skipe uniquelly mapped reads
        continue;
      }
    }

    //cout << bam.Name << endl;
    string chrom = refs.at(bam.RefID).RefName;
    string strand = "+";
    string strandMate = "+";
    if (bam.IsReverseStrand()) strand = "-";
    if (bam.IsMateReverseStrand()) strandMate = "-";

    unsigned int alignmentStart =  bam.Position+1;
    unsigned int alignmentEnd = bam.GetEndPosition();
    unsigned int threeEndPos = alignmentEnd;
    if (strand == "-") threeEndPos = alignmentStart;
    unsigned int cigarEnd;
    bool keep = false;
    vector <int> blockLengths;
    vector <int> blockStarts;
    blockStarts.push_back(0);
    ParseCigar(bam.CigarData, blockStarts, blockLengths, cigarEnd, keep);
    string mateChr = "SRP";
    unsigned int matePos = 0;

    if ( bam.IsMateMapped() == true) {
      mateChr = refs.at(bam.MateRefID).RefName;
      matePos = bam.MatePosition;
      int mateDistance = matePos-alignmentStart;
      if (mateChr != chrom || abs(mateDistance) > 230000)
        keep = true;          
    }

    if (keep == false) continue;   // skip normal reads

    if (chrom != oldChrom){
      if (oldChrom != "SRP") {
        output_processing(discoRegion, consisCountTh);
        discoRegion.id = 0;
        (discoRegion.tags).clear();
      }
      oldChrom = chrom;
    }
     
    if (discoRegion.id == 0){ //initialize
      discoID++;
      discoRegion.id    = discoID;
      discoRegion.type  = "D";
      discoRegion.chr   = chrom;
      discoRegion.coor  = threeEndPos;
      discoRegion.start = threeEndPos;
      discoRegion.end   = threeEndPos;
      discoRegion.support_no = 1;
      discoRegion.pw = 1;
      discoRegion.ct = 0;
      discoRegion.rep = "U";
      set<unsigned int> mateposTemp;
      mateposTemp.insert(matePos);
      (discoRegion.tags).insert( pair< string, set<unsigned int> >(mateChr, mateposTemp) );
    } else { //compare and add
      if (chrom == discoRegion.chr){
        int disDis = threeEndPos - discoRegion.coor;
        if (abs(disDis) < 430) {
          discoRegion.support_no++;
          discoRegion.pw++;
          if ( (discoRegion.tags).count(mateChr) > 0 ){
            (discoRegion.tags)[mateChr].insert(matePos);
          } else {
            set<unsigned int> mateposTemp;
            mateposTemp.insert(matePos);
            (discoRegion.tags).insert( pair< string, set<unsigned int> >(mateChr, mateposTemp) );
          }
        } else { //process the old one and initialize new one
          //process the old one
          output_processing(discoRegion, consisCountTh);
          //re-initialized
          discoRegion.id = 0;
          (discoRegion.tags).clear();
          discoID++;
          discoRegion.id    = discoID;
          discoRegion.type  = "D";
          discoRegion.chr   = chrom;
          discoRegion.coor  = threeEndPos;
          discoRegion.start = threeEndPos;
          discoRegion.end   = threeEndPos;
          discoRegion.support_no = 1;
          discoRegion.pw = 1;
          discoRegion.ct = 0;
          discoRegion.rep = "U";
          set<unsigned int> mateposTemp;
          mateposTemp.insert(matePos);
          (discoRegion.tags).insert( pair< string, set<unsigned int> >(mateChr, mateposTemp) );             
        }
      } else {
        cerr << "error: chromosme names is un-equal in the disco processing part!" << endl;
        exit(1);
      }
    } //compare and add

  }  // read a bam
 
  output_processing(discoRegion, consisCountTh);
  (discoRegion.tags).clear();
  reader.Close();
  return 0;

} //main


inline string int2str(unsigned int &i){
  string s;
  stringstream ss(s);
  ss << i;
  return ss.str();
}


inline string float2str(float &f){
  string s;
  stringstream ss(s);
  ss << f;
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


inline void ParseCigar(const vector<CigarOp> &cigar, vector<int> &blockStarts, vector<int> &blockLengths, unsigned int &alignmentEnd, bool &keep) {

  int currPosition = 0;
  int blockLength  = 0;

  //  Rip through the CIGAR ops and figure out if there is more
  //  than one block for this alignment
  vector<CigarOp>::const_iterator cigItr = cigar.begin();
  vector<CigarOp>::const_iterator cigEnd = cigar.end();
  for (; cigItr != cigEnd; ++cigItr) {
    switch (cigItr->Type) {
    case ('M') :                           // matching
      blockLength  += cigItr->Length;
      currPosition += cigItr->Length;
      break;
    case ('I') : break;                    // insertion
    case ('S') :                           // soft-clipping
      break;
    case ('D') :                           // deletion
      blockLength  += cigItr->Length;
      currPosition += cigItr->Length;
      break;
    case ('P') : break;                    // padding
    case ('N') :                           // skipped region
      blockStarts.push_back(currPosition + cigItr->Length);
      blockLengths.push_back(blockLength);
      currPosition += cigItr->Length;
      blockLength = 0;                     // a new block
      break;
    case ('H') :                           // for 'H' - do nothing, move to next op
      break;
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


inline void output_processing (struct region &region, unsigned int &consisCountTh) {

  string current_bp = int2str(region.id)+"\t"+region.type+"\t"+region.chr+"\t"+int2str(region.coor)+"\t"+int2str(region.support_no)+"\t"+int2str(region.pw)+"\t"+int2str(region.ct)+"\t"+region.rep;
  vector<unsigned int> consisCounts;
  unsigned int consisCount = 0;

  map < string, set<unsigned int> >::iterator tit = region.tags.begin();
  for(; tit != region.tags.end(); tit++){ //each new chr
    unsigned int lastPos = 0;
    string currentChr = tit->first;
    set <unsigned int>::iterator pt = (region.tags[currentChr]).begin(); //the set for position
    for (; pt != (region.tags[currentChr]).end(); pt++){
       if (lastPos == 0){
         consisCounts.push_back(1);
         lastPos = *pt;
       } else {
         int cdis = *pt - lastPos;
         if (abs(cdis) > 430) {
           consisCounts.push_back(1);
           lastPos = *pt;
         } else {
           int csize = consisCounts.size();
           consisCounts[csize-1] += 1;
         }
      }
    } //new Pos
  } //new chr  

  vector<unsigned int>::iterator cit = consisCounts.begin();
  for (; cit != consisCounts.end(); cit++){
    if (consisCount < *cit)
      consisCount = *cit;
  }

  if (consisCount >= consisCountTh){
    cout << current_bp << "\t" << consisCount << endl;
  }
   
  return;

}


