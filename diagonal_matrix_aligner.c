
result align_diagonal(char* query, char* target, int qlen, int tlen, int window_size, charvec *path) {

  int rows = tlen + qlen - 1;
  int score_matrix[rows][window_size];
  int direction_matrix[rows][window_size];

  // init matrices
  // I don't believe this should make a difference, but it doesn't work if we don't do it
  int k, l;
  for(k = 0; k < rows; k++) {
    for(l = 0; l < window_size; l++) {
      score_matrix[k][l] = 0;
      direction_matrix[k][l] = 0;
    }
  }

  int xoffset[rows];
  xoffset[0] = 0;
  // yoffset is implicitly the complement to the xoffset

  // fill in the matrix
  int y, x;
  int row_max, xlim;
  int qpos, tpos;
  int match, insertion, deletion;
  int xdiff;

  for(y = 0; y < rows; y++) {
    row_max = 0;
    xlim = y + 1 - xoffset[y]; // index of the first cell which goes outside the matrix
    if(window_size < xlim) xlim = window_size;
    if(tlen - xoffset[y] < xlim) xlim = tlen - xoffset[y];
    for(x = 0; x < xlim; x++) { // do not evaluate parts of the row in the overflow/padding at the top of the matrix
      qpos = y - xoffset[y] - x;
      tpos = x + xoffset[y];

      if(qpos >= qlen)
        continue;

      if(y > 0) {
        xdiff = xoffset[y] - xoffset[y-1];
        // test for out of bounds for deletion
        if(x-1+xdiff < 0 || x-1+xdiff >= xlim)
          deletion = LOW;
        else
          deletion = score_matrix[y-1][x-1+xdiff] + SCORES[DEL];
        // test for out of bounds for insertion
        if(x+xdiff < 0 || x+xdiff >= xlim)
          insertion = LOW;
        else
          insertion = score_matrix[y-1][x+xdiff] + SCORES[INS];
      } else {
        deletion = SCORES[DEL];
        insertion = SCORES[INS];
      }

      if(y > 1) {
        xdiff = xoffset[y] - xoffset[y-2];
        // test for out of bounds for match
        if(x-1+xdiff < 0 || x-1+xdiff >= xlim)
          match = LOW;
        else if(query[qpos] == target[tpos])
          match = score_matrix[y-2][x-1+xdiff] + SCORES[MATCH];
        else
          match = score_matrix[y-2][x-1+xdiff] + SCORES[MISMATCH];
      } else if(query[qpos] == target[tpos])
        match = SCORES[MATCH];
      else
        match = SCORES[MISMATCH];

      // set matrix value and direction
      if(match >= insertion && match >= deletion) {
        if(query[qpos] == target[tpos]) {
          score_matrix[y][x] = match;
          direction_matrix[y][x] = MATCH;
        } else {
          score_matrix[y][x] = match;
          direction_matrix[y][x] = MISMATCH;
        }
      } else if(insertion >= deletion) {
        score_matrix[y][x] = insertion;
        direction_matrix[y][x] = INS;
      } else {
        score_matrix[y][x] = deletion;
        direction_matrix[y][x] = DEL;
      }

      // keep track of maximum in this row
      if(score_matrix[y][x] > score_matrix[y][row_max])
        row_max = x;
    }

    if(y < rows - 1) {

      // shift offset by 1 toward the row max
      // bounded by both ends by edges of the target sequence

      xoffset[y+1] = xoffset[y];
      if(row_max >= window_size/2)
        xoffset[y+1]++;
      if(xoffset[y+1] < 0) xoffset[y+1] = 0;
      if(xoffset[y+1] > tlen - window_size) xoffset[y+1] = tlen - window_size;
    }
  }

  /*
   traceback
  */
  int qstart, tstart, qend, tend;

  // find the start position - this is tricky because our last "row" is not the set of terminal positions
  int max_x = 0, max_y = 0;
  for(y = rows-1; y > -1; y--) {
    for(x = 0; x < window_size; x++) {
      qpos = y - xoffset[y] - x;
      tpos = x + xoffset[y];
      if((tpos >= 0 && qpos == qlen - 1) || (qpos > 0 && tpos == tlen - 1)) {
        if((max_y == 0 && max_x == 0) || score_matrix[y][x] > score_matrix[max_y][max_x]) {
          max_y = y;
          max_x = x;
          qend = qpos;
          tend = tpos;
        }
      }
    }
  }

  // go backwards through matrix
  y = max_y;
  x = max_x;
  while(y >= 0 && y - xoffset[y] - x >= 0) {

    tpos = x + xoffset[y];
    qpos = y - xoffset[y] - x;

    kv_push(unsigned char, *path, direction_matrix[y][x]);
    qstart = qpos;
    tstart = tpos;

    if(direction_matrix[y][x] == MATCH || direction_matrix[y][x] == MISMATCH) {
      y = y - 2;
      if(y >= 0) { // if it is not, it won't matter what x becomes
        xdiff = xoffset[y+2] - xoffset[y];
        x = x - 1 + xdiff;
      }
    } else if(direction_matrix[y][x] == INS) {
      y--;
      if(y >= 0) {
        xdiff = xoffset[y+1] - xoffset[y];
        x = x + xdiff;
      }
    } else {
      y--;
      if(y >= 0) {
        xdiff = xoffset[y+1] - xoffset[y];
        x = x - 1 + xdiff;
      }
    }
  }

  result res;
  res.score = score_matrix[max_y][max_x];
  res.qstart = qstart;
  res.qend = qend;
  res.tstart = tstart;
  res.tend = tend;
  return res;
}
