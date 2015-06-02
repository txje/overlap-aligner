
// query along y-axis, target along x
result align_full_matrix(char* query, char* target, int qlen, int tlen, charvec *path) {

  int score_matrix[qlen][tlen];
  unsigned char direction_matrix[qlen][tlen];

  int x, y;

  for(y = 0; y < qlen; y++) {
    for(x = 0; x < tlen; x++) {
      int match_score, ins_score, del_score;

      // match
      if(query[y] == target[x]) {
        match_score = SCORES[MATCH];
      } else {
        match_score = SCORES[MISMATCH];
      }
      if(y > 0 && x > 0) {
        match_score = score_matrix[y-1][x-1] + match_score;
      }

      // ins
      if(y == 0) {
        ins_score = SCORES[INS];
      } else {
        ins_score = score_matrix[y-1][x] + SCORES[INS];
      }

      // del
      if(x == 0) {
        del_score = SCORES[DEL];
      } else {
        del_score = score_matrix[y][x-1] + SCORES[DEL];
      }

      // compare
      if(match_score >= ins_score && match_score >= del_score) {
        score_matrix[y][x] = match_score;
        if(query[y] == target[x]) {
          direction_matrix[y][x] = MATCH;
        } else {
          direction_matrix[y][x] = MISMATCH;
        }
      } else if(ins_score >= del_score) {
        score_matrix[y][x] = ins_score;
        direction_matrix[y][x] = INS;
      } else {
        score_matrix[y][x] = del_score;
        direction_matrix[y][x] = DEL;
      }
    }
  }

  // compute maximum score position
  int max_x = 0;
  int max_y = qlen - 1;
  for(x = 1; x < tlen; x++) { // check last row
    if(score_matrix[qlen - 1][x] > score_matrix[max_y][max_x]) {
      max_y = qlen - 1;
      max_x = x;
    }
  }
  for(y = 0; y < qlen; y++) { // check last column
    if(score_matrix[y][tlen - 1] > score_matrix[max_y][max_x]) {
      max_y = y;
      max_x = tlen - 1;
    }
  }

  x = max_x;
  y = max_y;
  int lastx, lasty;
  while(x >= 0 && y >= 0) {
    kv_push(unsigned char, *path, direction_matrix[y][x]);
    lastx = x;
    lasty = y;
    if(direction_matrix[y][x] == MATCH || direction_matrix[y][x] == MISMATCH) {
      x--;
      y--;
    } else if(direction_matrix[y][x] == INS) {
      y--;
    } else if(direction_matrix[y][x] == DEL) {
      x--;
    }
  }

  result res;
  res.score = score_matrix[max_y][max_x];
  res.qstart = lasty;
  res.qend = max_y;
  res.tstart = lastx;
  res.tend = max_x;

  return res;
}
