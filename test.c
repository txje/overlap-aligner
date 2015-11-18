
#include "ovlaligner.h"
#include "full_matrix_aligner.c"
#include "diagonal_matrix_aligner.c"
#include "diag_funnel_aligner.c"

void funnel_test(char* query, char* target, int qlen, int tlen) {
  printf("\nDiagonal funnel alignment:\n");
  charvec path;
  kv_init(path);
  int window_size = 10;
  int funnel_size = 20;
  int funnel_len = 5;
  result res = align_diagonal_funnel(query, target, qlen, tlen, window_size, funnel_size, funnel_len, &path);
  printf("Score: %d\n", res.score);
  printf("Query: %d - %d\n", res.qstart, res.qend);
  printf("Target: %d - %d\n", res.tstart, res.tend);
  int path_size = kv_size(path);
  printf("Path length: %d\n", path_size);
  unsigned char cigar[path_size];
  int i;
  for(i = 0; i < path_size; i++) {
    cigar[i] = kv_pop(path);
    printf("%d, ", cigar[i]);
  }
  printf("\n");
  float acc = cigar_accuracy(cigar, path_size);
  printf("Accuracy: %f\n", acc);
  kv_destroy(path);
}

void diag_test(char* query, char* target, int qlen, int tlen, int semilocal) {
  printf("\nDiagonal alignment:\n");
  charvec path;
  kv_init(path);
  int window_size = 10;
  result res = align_diagonal(query, target, qlen, tlen, window_size, &path, semilocal);
  printf("Score: %d\n", res.score);
  printf("Query: %d - %d\n", res.qstart, res.qend);
  printf("Target: %d - %d\n", res.tstart, res.tend);
  int path_size = kv_size(path);
  printf("Path length: %d\n", path_size);
  unsigned char cigar[path_size];
  int i;
  for(i = 0; i < path_size; i++) {
    cigar[i] = kv_pop(path);
    printf("%d, ", cigar[i]);
  }
  printf("\n");
  float acc = cigar_accuracy(cigar, path_size);
  printf("Accuracy: %f\n", acc);
  kv_destroy(path);
}

void full_test(char* query, char* target, int qlen, int tlen, int semilocal) {
  printf("\nFull matrix alignment:\n");
  charvec path;
  kv_init(path);
  result res = align_full_matrix(query, target, qlen, tlen, &path, semilocal);
  printf("Score: %d\n", res.score);
  printf("Query: %d - %d\n", res.qstart, res.qend);
  printf("Target: %d - %d\n", res.tstart, res.tend);
  int path_size = kv_size(path);
  printf("Path length: %d\n", path_size);
  unsigned char cigar[path_size];
  int i;
  for(i = 0; i < path_size; i++) {
    cigar[i] = kv_pop(path);
    printf("%d, ", cigar[i]);
  }
  printf("\n");
  float acc = cigar_accuracy(cigar, path_size);
  printf("Accuracy: %f\n", acc);
  kv_destroy(path);
}

int main(int argc, char *argv[]) {
  char* query = "ACCGATAAATTGATCCGATAGAAATTATACGATGGTAGACTAAGAT";
  char* target = "AGGTGTTGACTGAGGACCTGTGAGTACTGATTGATAACCGATGATAGTCCTGATA";
  printf("\nquery: %s\n", query);
  printf("target: %s\n", target);
  int qlen = strlen(query);
  int tlen = strlen(target);
  full_test(query, target, qlen, tlen, 0);
  full_test(query, target, qlen, tlen, 1);
  diag_test(query, target, qlen, tlen, 0);
  diag_test(query, target, qlen, tlen, 1);
  funnel_test(query, target, qlen, tlen);

  query = "TTTAAGTAGTACGCCGATGACCGTAGATGGGAAGATGAGATCGCGTGAGTCGATG";
  target = "TTTAAGTAGTACGCCGATGACCGTAGATGGTGTGAGATTGTTTTGAGGTAGCGTAG";
  printf("\nquery: %s\n", query);
  printf("target: %s\n", target);
  qlen = strlen(query);
  tlen = strlen(target);
  full_test(query, target, qlen, tlen, 0);
  full_test(query, target, qlen, tlen, 1);
  diag_test(query, target, qlen, tlen, 0);
  diag_test(query, target, qlen, tlen, 1);
  funnel_test(query, target, qlen, tlen);
}
