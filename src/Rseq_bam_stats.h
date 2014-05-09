#include <cstdio>
#include <getopt.h>
#include <cstdlib>
#include <cstring>


struct parameters {
  char* mapping_f;
  char* pileup;
  char* type;
  char* writer;
  char* unmapped;
  char* arp;
  char* breakpoint;
  unsigned int readlength;
};

struct parameters* interface(struct parameters* param,int argc, char *argv[]);
void delete_param(struct parameters* param);
void usage(void);

const char* program_name;

struct parameters* interface(struct parameters* param, int argc, char *argv[]){

  program_name = argv[0];
  int c;     // the next argument
  int help = 0;

  if (argc < 2){
    usage();
    exit(0);
  }

  param = new struct parameters;
  param->pileup = new char;
  param->mapping_f = new char;
  param->type = new char;
  param->writer = new char;
  param->unmapped = new char; 
  param->arp = new char;
  param->breakpoint = new char;

  const struct option long_options[] ={
    {"mapping",1,0,'m'},
    {"pileup",1,0,'p'},
    {"type",1,0,'t'},
    {"writer",1,0,'w'},
    {"unmapped",1,0,'u'},
    {"arp",1,0,'a'},
    {"breakpoint",1,0,'b'},
    {"readlength",1,0,'l'},
    {"help",0,0,'h'},
    {0, 0, 0, 0}
  };


  while (1){

    int option_index = 0;
    c = getopt_long_only (argc, argv,"hm:t:p:w:u:a:b:l:",long_options, &option_index);

    if (c == -1){
      break;
    }

    switch(c) {
    case 0:
      break;
    case 'm':
      param->mapping_f = optarg;
      break;
    case 't':
      param->type = optarg;
      break;
    case 'p':
      param->pileup = optarg;
      break;
    case 'w':
      param->writer = optarg;
      break;
    case 'u':
      param->unmapped = optarg;
      break;
    case 'a':
      param->arp = optarg;
      break;
    case 'b':
      param->breakpoint = optarg;
      break;
    case 'l':
      param->readlength = atoi(optarg);
      break;
    case 'h':
      help = 1;
      break;
    case '?':
      help = 1;
      break;
    default:
      help = 1;
      break;
    }
  }

  if(help){
    usage();
    delete_param(param);
    exit(0);
  }

  return param;
}

void usage()
{
  fprintf(stdout, "\nRseq_bam_stats, Copyright (C) 2011 Sun Ruping <ruping@molgen.mpg.de>\n");
  fprintf(stdout, "\n");
  fprintf(stdout, "Usage: %s options [inputfile] \n\n", program_name);
  fprintf(stdout, "-h --help        print the help message\n");
  fprintf(stdout, "-m --mapping     mapping_file (bam file)\n");
  fprintf(stdout, "-p --pileup      <forget it, currently it is no use> yes: allow pileup, no: skip redundant reads \n");
  fprintf(stdout, "-w --writer      the bam output name (for mapped ailgnments).\n");
  fprintf(stdout, "-u --unmapped    write the unmapped tag names into this file.\n");
  fprintf(stdout, "-a --arp         filename of the arp read name (for fusion assembly use, default not write).\n");
  fprintf(stdout, "-b --breakpoint  the file for output of potential breakpoint.\n");
  fprintf(stdout, "-l --readlength  the length of the reads.\n");
  fprintf(stdout, "-t --type        (p)aired-end or (s)ingle-end.\n");
  fprintf(stdout, "\n");
}


void delete_param(struct parameters* param)
{
  delete(param->mapping_f);
  delete(param->type);
  delete(param->pileup);
  delete(param->writer);
  delete(param->unmapped);
  delete(param->arp);
  delete(param->breakpoint);
  delete(param);
}
