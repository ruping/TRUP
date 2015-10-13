/*****************************************************************************
  
  Rseq_bam_reads2expr.cpp (improved c++ version of reads2expr.pl from Ho)
  1) calculate gene expr based on an bed file with exon information
  2) calculate the locus bias for each gene

  (c) 2011 - Sun Ruping
  Dept. Vingron (Computational Mol. Bio.)
  Max-Planck-Institute for Molecular Genetics
  Ihnestr. 73, D-14195, Berlin, Germany   

  current: Department of Systems Biology, Columbia University, NY, USA
  email: rs3412@c2b2.columbia.edu

g++ Rseq_bam_reads2expr.cpp
-I/ifs/home/c2b2/ac_lab/rs3412/tools/bamtools/include/
-I/ifs/home/c2b2/ac_lab/rs3412/tools/zlib-1.2.8/include/
-L/ifs/home/c2b2/ac_lab/rs3412/tools/bamtools/lib/
-L/ifs/home/c2b2/ac_lab/rs3412/tools/zlib-1.2.8/lib/
-lbamtools -lz -static -o 
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
#include "Rseq_bam_reads2expr.h"
using namespace std;

struct region {  // a bed file containing gene annotations
  string chr;
  unsigned int start;
  unsigned int end;
  string ensemble_id;
  float score;
  string strand;
  unsigned int thickStart;
  unsigned int thickEnd;
  string itemRgb;
  unsigned int blockCount;
  string blockSizes;
  string blockStarts;
  // results storing here
  map <int, float> tags;  // storing position and count of tags relative to a gene
};

struct lb {  // for locus bias
  float ps;
  unsigned int hs;
  float pe;
  unsigned int he;
};

unsigned int read_length = 0;

inline void ParseCigar(const vector<CigarOp> &cigar, vector<int> &blockStarts, vector<int> &blockEnds, unsigned int &alignmentEnd);
inline void splitstring(const string &str, vector<string> &elements, const string &delimiter);
inline void eatline(const string &str, deque <struct region> &region_ref, bool &noChr);
inline string int2str(unsigned int &i);
inline string float2str(float &f);
inline void gene_processing(struct region &gene, vector <struct lb> &locusb);

int main ( int argc, char *argv[] ) { 

  struct parameters *param = 0;
  param = interface(param, argc, argv);

  //region file input (the region file should be sorted as the same way as the bam file)
  ifstream region_f;
  region_f.open(param->region_f, ios_base::in);  // the region file is opened

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


  // locus bias
  struct lb empty_profile = {0,0,0,0};
  vector <struct lb> locus_b(1000, empty_profile);
  // output locus bias file
  string locus_bias_set = param->lbias;
  ofstream locus_bias;
  if ( locus_bias_set != "" ) {
    locus_bias.open(param->lbias);
    if ( !locus_bias ) {
      cerr << "can not open locus_bias file.\n";
      exit(0);
    }
  }

  //should decide which chromosome
  string line;
  string old_chr = "SRP";
  string type = param->type;

  //whether do some position-level pile-up stuff
  bool posc = false;
  ofstream posc_f;
  ofstream chrmap_f;
  string poscset = param->posc;
  if ( poscset != "" ) {
    posc = true;
    posc_f.open(param->posc);
    chrmap_f.open(param->chrmap);
  }

  bool noChr;
  if ( param->nochr == 1 ){
    noChr = true;
  } else {
    noChr = false;
  }

  //regions for the input of region file
  deque <struct region> regions;

  getline(region_f, line); //get the first line
  eatline(line,regions,noChr);
  
  deque <struct region>::iterator it = regions.begin();

  while ( it->chr != old_chr ) {

    old_chr = it->chr;  // set the current chr as old chr

    int chr_id  = reader.GetReferenceID(it->chr);

    if ( chr_id == -1 ) {  //reference not found

      for (; it != regions.end() && it->chr == old_chr; ) {
        gene_processing(*it,locus_b);           // print the old region info
        it = regions.erase(it);         // erase the current region
      }
  
      while ( regions.empty() ) {    
        getline(region_f, line);
        if ( region_f.eof() ){
          cerr << "finished: end of region file, zone 0" << endl;
          break;
        }
        eatline(line, regions,noChr);
        it = regions.begin();
        if (it->chr == old_chr){  
          gene_processing(*it,locus_b);      
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

    //pile-up pos stats
    set <string> fragment;
    map <string, unsigned int> pileup;
    bool isposPileup = false;
    unsigned int old_start   = 0;
    unsigned int total_tags  = 0;
    unsigned int total_pos   = 0;
    unsigned int pileup_pos  = 0;


    BamAlignment bam;
    while (reader.GetNextAlignment(bam)) {

      if ( bam.IsMapped() == false ) continue;   // skip unaligned reads

      unsigned int unique;
      if ( bam.HasTag("NH") ) {
        bam.GetTag("NH", unique);                   // rnaseq aligners, such as gsnap and tophat                                                                                                                                                                                                                             
      } else {
        if (bam.MapQuality > 10) {                  // other aligner                                                                                                                                                                                                                                                         
          unique = 1;
        }
      }

      if (param->unique == 1) {
        if (unique != 1) {                       // skipe multi mapped reads
          continue;
        }
      }

      if (read_length == 0){
        read_length = bam.Length;
      }

      //cout << bam.Name << endl;
      string chrom = refs.at(bam.RefID).RefName;
      string strand = "+";
      if (bam.IsReverseStrand()) strand = "-";

      unsigned int alignmentStart =  bam.Position+1;
      unsigned int mateStart;
      if (type == "p") mateStart = bam.MatePosition+1;
      unsigned int alignmentEnd = bam.GetEndPosition();
      unsigned int cigarEnd;
      vector <int> blockLengths;
      vector <int> blockStarts;
      blockStarts.push_back(0);
      ParseCigar(bam.CigarData, blockStarts, blockLengths, cigarEnd);


      // position check for unique mapped reads (because is paired-end reads, shoule base on fragment level for paired end reads)
      if (posc == true && unique == 1) {

        if (type == "p" && fragment.count(bam.Name) > 0) 
          fragment.erase(bam.Name);

        else {

          total_tags++;
          if (type == "p"){
            fragment.insert(bam.Name);
          }
          string alignSum;
          if (type == "p") {
             alignSum = int2str(alignmentStart) + "\t" + int2str(mateStart) + "\t.\t" + strand;
          } else {
             alignSum = int2str(alignmentStart) + "\t" + int2str(alignmentEnd) + "\t.\t" + strand;
          }

          if ( alignmentStart != old_start ) {
            isposPileup = false;
            map <string, unsigned int>::iterator pit = pileup.begin();            
            for (; pit != pileup.end(); pit++) {
              posc_f << chrom << "\truping\tpileup\t" << pit->first << "\t.\t" << "Pileup=" << pit->second << endl;     //print pileup
            }
            pileup.clear();           //clear pileup set
            pileup.insert( pair <string, unsigned int> (alignSum, 1) );  //insert the new read
            total_pos++;
          }

          else if ( alignmentStart == old_start ) { // same starts
            if ( pileup.count(alignSum) > 0 ) {  // pileup
              if ( pileup[alignSum] == 1 && isposPileup == false ) { 
                pileup_pos++; isposPileup = true;
              }
              pileup[alignSum]++;
            }
            else {
              pileup.insert( pair <string, unsigned int> (alignSum, 1) );
            }
          } //same starts

        }   //new fragment

        old_start = alignmentStart;
      } // do pos check



      float incre = 1.;
      if (blockStarts.size() > 1) incre = 0.5;     // incre half for junction reads
      incre /= static_cast < float >(unique);        // for multi aligned reads

      deque <struct region>::iterator iter = regions.begin();

      if ( iter->start > alignmentEnd ) continue;  // skip reads not overlapping with the first region

      while ( iter->chr == old_chr && iter->start <= alignmentEnd && iter != regions.end() ) {

        if (iter->end < alignmentStart) {            // the region end is beyond the alignmentStart

          gene_processing(*iter,locus_b);            // processing
          iter = regions.erase(iter);                // this region should be removed
          if ( regions.empty() ) { 
            getline(region_f, line);                        // get a line of region file
            if ( ! region_f.eof() ) {
              eatline(line, regions, noChr);                         // eat a line and put it into the duque
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

          vector <int>::iterator cigit = blockStarts.begin();
          for (; cigit != blockStarts.end(); cigit++) {
            unsigned int current_start = *cigit + alignmentStart;
            int current_pos = current_start - (iter->start);
            //cout << iter->chr << "\t" << iter->start << "\t" << iter->end << "\t" << current_start << endl;
            if ( (iter->tags).count(current_pos) > 0 ) {
              (iter->tags)[current_pos] += incre;
            }
            else
              (iter->tags).insert( pair<int, float>(current_pos, incre) );  
          }

        }  // overlapping take action!

        if ( (iter+1) != regions.end() )
          iter++;                                           // if this region is not the last element in the deque
        else {                                              // the last element
          getline(region_f, line);                          // get a line of region file
          if ( ! region_f.eof() ){
            eatline(line, regions, noChr);                         // eat a line and put it into the duque
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


    // print chr map
    if (posc == true) {
      chrmap_f << old_chr << "\t" << total_tags << "\t" << total_pos << "\t" << pileup_pos << endl;
    } 
 
    //somehow to loop back
    it = regions.begin();                   //reset to begin
    for (; it != regions.end() && it->chr == old_chr; ) {
      gene_processing(*it,locus_b);              // print the old region info
      it = regions.erase(it);             // erase the current region
    }
  
    while ( regions.empty() ) {    

      getline(region_f, line);
      if ( region_f.eof() ){
        cerr << "finished: end of region file, zone 5" << endl;
        //print locus bias
        for (unsigned int l = 0; l < 1000; l++){
	  locus_bias << l << "\t" << locus_b[l].ps << "\t" << locus_b[l].hs << "\t" << locus_b[l].pe << "\t" << locus_b[l].he << endl;
	}
        exit(0);
      }
      eatline(line, regions, noChr);
      it = regions.begin();
      if (it->chr == old_chr){
        gene_processing(*it, locus_b);      
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


inline void eatline(const string &str, deque <struct region> &region_ref, bool &noChr) {
  
   vector <string> line_content;
   //split line and then put it into a deque

   splitstring(str, line_content, "\t");
   vector <string>::iterator iter = line_content.begin();
   unsigned int i;

   struct region tmp;
 
   for(i = 1; iter != line_content.end(); iter++, i++){
     switch (i) {
     case 1:  // chr
       tmp.chr = *iter;
       if (noChr == true){
         tmp.chr = (*iter).substr(3);
         if (tmp.chr == "M"){
           tmp.chr = "MT";
         }
       }
       continue;
     case 2:  // start
       tmp.start = atoi((*iter).c_str()) + 1;
       continue;
     case 3:  // end
       tmp.end = atoi((*iter).c_str());
       continue;
     case 4:  // ensemble_id
       tmp.ensemble_id = *iter;
       continue;
     case 5:  // score
       tmp.score = atof((*iter).c_str());
       continue;
     case 6:  // strand
       tmp.strand = *iter;
       continue;
     case 7:  // thickStart
       tmp.thickStart = atoi((*iter).c_str());
       continue;
     case 8:  // thickEnd
       tmp.thickEnd = atoi((*iter).c_str());
       continue;
     case 9: // itemRgb
       tmp.itemRgb = *iter;
       continue;
     case 10: // blockCount
       tmp.blockCount = atoi((*iter).c_str());
       continue;
     case 11: // blockSizes
       tmp.blockSizes = *iter;
       continue;
     case 12: // blockStarts
       tmp.blockStarts = *iter;
       continue;
     default:
       break;
     }
   }

   region_ref.push_back(tmp);
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
    case ('M') :                           // matching
      blockLength  += cigItr->Length;
      currPosition += cigItr->Length;
      break;
    case ('I') : break;                    // insertion
    case ('S') : break;                    // soft-clipping
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
    case ('H') : break;                    // for 'H' - do nothing, move to next op
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


inline void gene_processing(struct region &gene, vector <struct lb> &locusb) {

  vector <string> exonL;
  vector <string> exonS;
  unsigned int idx = 0;
  vector <float> exonT(gene.blockCount, 0.);
  float totalT = 0.;
  float totalP = 0.;

  splitstring(gene.blockSizes, exonL, ",");
  splitstring(gene.blockStarts, exonS, ",");
  unsigned int total_length = 0;
  vector <string>::iterator lit = exonL.begin();
  for (; lit != exonL.end(); lit++) {
    total_length += atoi((*lit).c_str());
  }
  vector <float> exonPos2tags(total_length,0.);
  vector <unsigned int> exonPos2hit(total_length,0);
  vector <unsigned int> exonPos2depth(total_length,0);

  unsigned int current_exon_start = 0;

  for (unsigned int i = 0; i < gene.blockCount; i++) { // exon i

    unsigned int exon_length = atoi(exonL[i].c_str());
    
    for (unsigned int j = 0; j < exon_length; j++) {    // each position of exon i
      unsigned int pos = atoi(exonS[i].c_str()) + j;
 
      if ( gene.tags.count(pos) > 0 ) {  // current position having tags
	totalT += (gene.tags)[pos];
        totalP += 1;
        exonT[i] += (gene.tags)[pos];
        exonPos2tags[idx] += (gene.tags)[pos];
        exonPos2hit[idx] = 1;

        // for base depth
        
        for (unsigned int k = idx; k < (current_exon_start + exon_length) && k < (idx + read_length); k++){
           exonPos2depth[k] += 1;
        }
        
      }
      idx++;
    } // j
    current_exon_start += exon_length;
    exonT[i] /= exon_length;
  } // i

  if (idx > 0) {

     gene.score = totalT/idx;
     float poscore = totalP/idx;

     vector <float>::iterator eit = exonT.begin();
     string exontags;
     for (; eit != exonT.end(); eit++ ){
       if ((eit+1) != exonT.end()) exontags += float2str(*eit)+",";
       else exontags += float2str(*eit);
     }

     float depthNonempty = 0.;
     vector <unsigned int>::iterator depit = exonPos2depth.begin();
     for (; depit != exonPos2depth.end(); depit++){
       if ( (*depit) > 0 )
         depthNonempty += 1; 
     }
     float depthScore = depthNonempty/idx;

     cout << gene.chr << "\t" << gene.start << "\t" << gene.end << "\t" << gene.ensemble_id << "\t" << gene.score << "\t" << poscore << "\t" << depthScore << "\t" << totalT << "\t" << idx << "\t" << gene.strand << "\t" << gene.blockCount << "\t" << gene.blockSizes << "\t" << gene.blockStarts << "\t" << exontags << endl;

     //here for locus bias
     unsigned int maxL = 1000;
     if (maxL > idx) maxL = idx;

     if (gene.strand == "+") {
       for (unsigned int m = 0; m < maxL; m++){
         locusb[m].ps += exonPos2tags[m];
         locusb[m].hs++;
         locusb[m].pe += exonPos2tags[idx-1-m];
         locusb[m].he++;
       }
     } // + strand
     else {
       for (unsigned int n = 0; n < maxL; n++){
	 locusb[n].pe += exonPos2tags[n];
         locusb[n].he++;
         locusb[n].ps += exonPos2tags[idx-1-n];
         locusb[n].hs++;
       }
     } // - strand

  } // idx > 0
}


