#include <cstdio>
#include <getopt.h>
#include <cstdlib>
#include <cstring>


struct parameters {
  char* mapping_f;
  char* type;
  unsigned int idstart;
  unsigned int consisCount;
  unsigned int unique;
};

struct parameters* interface(struct parameters* param,int argc, char *argv[]);
void delete_param(struct parameters* param);
void usage(void);

const char* program_name;

struct parameters* interface(struct parameters* param, int argc, char *argv[]){

  program_name = argv[0];
  int c;     // the next argument
  int help = 0;

  if (argc < 2) {
    usage();
    exit(0);
  }

  param = new struct parameters;
  param->mapping_f = new char;
  param->type = new char;
 
  const struct option long_options[] ={
    {"mapping",1,0,'m'},
    {"type",1,0,'t'},
    {"unique",0,0,'u'},
    {"idstart",1,0,'i'},
    {"consisCount",1,0,'c'},
    {"help",0,0,'h'},
    {0, 0, 0, 0}
  };

  while (1){

    int option_index = 0;
    c = getopt_long_only (argc, argv,"hum:t:i:c:",long_options, &option_index);

    if (c == -1) {
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
    case 'u':
      param->unique = 1;
      break;
    case 'i':
      param->idstart = atoi(optarg);
      break;
    case 'c':
      param->consisCount = atoi(optarg);
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

  if(help) {
    usage();
    delete_param(param);
    exit(0);
  }

  return param;
}

void usage()
{
  fprintf(stdout, "\ndiscordant_mate, Copyright (C) 2011 Sun Ruping <ruping@molgen.mpg.de>\n");
  fprintf(stdout, "\n");
  fprintf(stdout, "Usage: %s options [inputfile] \n\n", program_name);
  fprintf(stdout, "-h --help    print the help message\n");
  fprintf(stdout, "-m --mapping <filename>  mapping_file (RNA-seq bam file, chromosomes and coordinates sorted also)\n");
  fprintf(stdout, "-m --idstart INT  starting id\n");
  fprintf(stdout, "-m --consisCount INT  threshold for consistent read pairs with discordant mapping\n");
  fprintf(stdout, "-q --unique              only calculate for uniquely mapped reads (you may not set this when the bam files only contain unique reads).\n");
  fprintf(stdout, "-t --type    <p/s>       current no use\n");
  fprintf(stdout, "\n");
}


void delete_param(struct parameters* param)
{
  delete(param->mapping_f);
  delete(param->type);
  delete(param);
}
