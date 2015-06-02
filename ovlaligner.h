
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "klib/kvec.h" // C dynamic vector

typedef kvec_t(unsigned char) charvec;

struct aln_result {
  int score;
  int qstart;
  int qend;
  int tstart;
  int tend;
};
typedef struct aln_result result;

unsigned char MATCH = 0, INS = 1, DEL = 2, MISMATCH = 3;
int SCORES[] = {1, -1, -1, -1}; // corresponds to events above
int LOW = -2000000000; // almost the lowest 32-bit integer

#include "aligner_utils.c" // C dynamic vector
