float cigar_accuracy(unsigned char cigar[], int cigar_len) {
  int matches = 0, total = 0;
  int i;
  for(i = 0; i < cigar_len; i++) {
    if(cigar[i] == MATCH) {
      matches++;
    }
    total++;
  }
  return ((float)matches) / total;
}
