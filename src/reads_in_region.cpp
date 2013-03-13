/*****************************************************************************
  
  get reads in a region from a sorted bam file

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
#include "reads_in_region.h"
using namespace std;

struct region {  // a txt file containing breakpoints
  unsigned int id;
  string type;
  string chr;
  unsigned int coor;
  unsigned int support_no;
  unsigned int pw;              //proper wrong pair counts
  unsigned int start;
  unsigned int end;
  set <string> tags;
};

map < unsigned int, vector<string> > rembp;
map < unsigned int, set<string> > remtag;

inline void ParseCigar(const vector<CigarOp> &cigar, vector<int> &blockStarts, vector<int> &blockEnds, unsigned int &alignmentEnd, bool &keep);
inline void splitstring(const string &str, vector<string> &elements, const string &delimiter);
inline void eatline(const string &str, deque <struct region> &region_ref);
inline string int2str(unsigned int &i);
inline string float2str(float &f);
inline void output_processing(struct region &region);

int main ( int argc, char *argv[] ) { 

  struct parameters *param = 0;
  param = interface(param, argc, argv);

  //region file input (the region file should be sorted as the same way as the bam file)
  ifstream region_f;
  region_f.open(param->region_f, ios_base::in);  // the region file is opened

  //bam input and generate index if not yet
  BamReader reader;
  reader.Open(param->mapping_f);   // the mapping bam file is opened 

  // get header & reference information
  string header = reader.GetHeaderText();
  RefVector refs = reader.GetReferenceData();

  if ( !reader.LocateIndex() )  //create index incase 
     reader.CreateIndex();

  //should decide which chromosome
  string line;
  string old_chr = "SRP";
  string type = param->type;

  //regions for the input of region file
  deque <struct region> regions;

  getline(region_f, line); //get the first line
  eatline(line,regions);
  
  deque <struct region>::iterator it = regions.begin();

  while ( it->chr != old_chr ) {

    old_chr = it->chr;  // set the current chr as old chr

    int chr_id  = reader.GetReferenceID(it->chr);

    if ( chr_id == -1 ) {  //reference not found

      for (; it != regions.end() && it->chr == old_chr; ) {
        output_processing(*it);           // print the old region info
        it = regions.erase(it);         // erase the current region
      }
  
      while ( regions.empty() ) {    
        getline(region_f, line);
        if ( region_f.eof() ){
          cerr << "finished: end of region file, zone 0" << endl;
          break;
        }
        eatline(line, regions);
        it = regions.begin();
        if (it->chr == old_chr){  
          output_processing(*it);      
          regions.clear();
          continue;
        }
      }
      continue;
    }

    int chr_len = refs.at(chr_id).RefLength;

    if ( !reader.SetRegion(chr_id, 1, chr_id, chr_len) ) // here set region
      {
        cerr << "bamtools count ERROR: Jump region failed " << it->chr << endl;
        reader.Close();
        exit(1);
      }


    BamAlignment bam;
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
      if (bam.IsReverseStrand()) strand = "-";

      unsigned int alignmentStart =  bam.Position+1;
      unsigned int alignmentEnd = bam.GetEndPosition();
      unsigned int cigarEnd;
      bool keep = false;
      vector <int> blockLengths;
      vector <int> blockStarts;
      blockStarts.push_back(0);
      ParseCigar(bam.CigarData, blockStarts, blockLengths, cigarEnd, keep);

      //if (bam.Name == "Run0009Lane4Tile83x12044y19781Multi0"){
      //  if (keep == false) 
      //    cerr << "yes" << endl;
      // else {
      //   cerr << "shit" << endl;
      // }
      //}

      if (keep == false) continue;   // skip normal reads

      if ( regions.empty() ) {
        getline(region_f, line);
        if ( region_f.eof() ){
          cerr << "finished: end of region file, zone 2.5" << endl;
          exit(0);
        }
      }
      
      deque <struct region>::iterator iter = regions.begin();

      if ( iter->start > alignmentEnd ) continue;  // skip reads not overlapping with the first region
  
      while ( iter->chr == old_chr && iter->start <= alignmentEnd && iter != regions.end() ) {

        if (iter->end < alignmentStart) {            // the region end is beyond the alignmentStart

          output_processing(*iter);            // processing
          iter = regions.erase(iter);                // this region should be removed
          if ( regions.empty() ) { 
            getline(region_f, line);                        // get a line of region file
            if ( ! region_f.eof() ) {
              eatline(line, regions);                         // eat a line and put it into the duque
              iter = regions.begin();
            }
            else {  // it's reaching the end of the region file
              cerr << "finished: end of region file, zone 3" << endl;
              break;
            }
          }
          continue;
        }

        if (iter->end >= alignmentStart && iter->start <= alignmentEnd) {  //overlapping, should take action

          (iter->tags).insert(bam.Name);

        }  // overlapping take action!

        if ( (iter+1) != regions.end() )
          iter++;                                           // if this region is not the last element in the deque
        else {                                              // the last element
          getline(region_f, line);                          // get a line of region file
          if ( ! region_f.eof() ){
            eatline(line, regions);                         // eat a line and put it into the duque
            iter = regions.end();
            iter--;
          }
          else {  //it's reaching the end of the region file
            cerr << "finished: end of region file, zone 4" << endl;
            break;
          }
        }

      } //while

    }  // read a bam
 
    //somehow to loop back
    it = regions.begin();                   //reset to begin
    for (; it != regions.end() && it->chr == old_chr; ) {
      output_processing(*it);              // print the old region info
      it = regions.erase(it);             // erase the current region
    }

    while ( regions.empty() ) {    
      getline(region_f, line);
      if ( region_f.eof() ){
        cerr << "finished: end of region file, zone 5" << endl;
        exit(0);
      }
      eatline(line, regions);
      it = regions.begin();
      if (it->chr == old_chr){
        output_processing(*it);      
        regions.clear();
        continue;
      }
    }

  } // region chr != old chr
      
  regions.clear();
  reader.Close();
  region_f.close();
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


inline void eatline(const string &str, deque <struct region> &region_ref) {
  
   vector <string> line_content;
   //split line and then put it into a deque

   splitstring(str, line_content, "\t");
   vector <string>::iterator iter = line_content.begin();
   unsigned int i;

   struct region tmp;
 
   for(i = 1; iter != line_content.end(); iter++, i++){
     switch (i) {
     case 1:  // id
       tmp.id = atoi((*iter).c_str());
       continue;
     case 2:  // type
       tmp.type = *iter;
       continue;
     case 3:  // chr
       tmp.chr = *iter;
       continue;
     case 4:  // coor
       tmp.coor = atoi((*iter).c_str());
       tmp.start = tmp.coor - 200;
       tmp.end   = tmp.coor + 200;
       continue;
     case 5:  // support
       tmp.support_no = atof((*iter).c_str());
       continue;
     case 6:  // pw
       tmp.pw = atof((*iter).c_str());
       continue;
     default:
       break;
     }
   }

   region_ref.push_back(tmp);
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
    case ('I') : break;                    // insertion
    case ('S') :                           // soft-clipping
      if (cigItr->Length >= 8)
        keep = true;
      break;
    case ('D') : break;                    // deletion
      blockLength  += cigItr->Length;
      currPosition += cigItr->Length;
    case ('P') : break;                    // padding
    case ('N') :                           // skipped region
      blockStarts.push_back(currPosition + cigItr->Length);
      blockLengths.push_back(blockLength);
      currPosition += cigItr->Length;
      blockLength = 0;                     // a new block
      break;
    case ('H') :                           // for 'H' - do nothing, move to next op
      if (cigItr->Length >= 8)    
        keep = true;
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


inline void output_processing (struct region &region) {

  string current_bp = int2str(region.id)+"\t"+region.chr+"\t"+int2str(region.coor)+"\t"+int2str(region.support_no)+"\t"+int2str(region.pw)+"\t"+region.type;

  if ( region.type == "s" ){
    cout << current_bp << endl;
    set <string>::iterator tagit = (region.tags).begin();
    for (; tagit != (region.tags).end(); tagit++){
      cout << *tagit << endl;
    }
    return;
  }

  //rest is for the type of "p"

  if ( rembp.count(region.id) > 0 ){   // the id is already there
    if ( rembp[region.id].size() == 1 ) {  // should add the other
      rembp[region.id].push_back(current_bp);
      set <string>::iterator tagit = (region.tags).begin();
      for (; tagit != (region.tags).end(); tagit++){
        if ( remtag[region.id].count(*tagit) == 0 ){
          remtag[region.id].insert(*tagit);
        }
      }
    } //add for the paired    
    if ( rembp[region.id].size() == 2 ) {  // should print out

      string bp1 = (rembp[region.id])[0];
      string bp2 = (rembp[region.id])[1];
      cout << bp1 << endl;
      cout << bp2 << endl;
      set <string>::iterator tagit2 = (remtag[region.id]).begin();
      for (; tagit2 != (remtag[region.id]).end(); tagit2++){
        cout << *tagit2 << endl;
      }

      map < unsigned int, vector <string> >::iterator erait_bp;
      map < unsigned int, set <string> >::iterator erait_tag;
      erait_bp = rembp.find(region.id);
      erait_tag = remtag.find(region.id);
      rembp.erase(erait_bp);     
      remtag.erase(erait_tag);
    } // print out and delete
  }
  else {  // an new id
    remtag.insert( pair< unsigned int, set<string> > (region.id, region.tags));
    vector <string> tmp_bp;
    tmp_bp.push_back(current_bp);
    rembp.insert( pair< unsigned int, vector<string> > (region.id, tmp_bp));
  }

}


