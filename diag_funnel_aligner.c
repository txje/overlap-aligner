
/*
 * establish out-of-bounds rules for the x position:
 */
int get_xlim(int y, int xoffset, int row_size, int tlen) {
  int xlim = y + 1 - xoffset; // index of the first cell which goes outside the matrix
  if(row_size < xlim) xlim = row_size;
  if(tlen - xoffset < xlim) xlim = tlen - xoffset; // remaining length of target
  return xlim;
}

result align_diagonal_funnel(char* query, char* target, int qlen, int tlen, int window_size, int funnel_size, int funnel_len, charvec *path) {

  int rows = tlen + qlen - 1;
  if(rows < funnel_len) {
    funnel_len = rows;
  }
  int upper_score_matrix[funnel_len][funnel_size];
  int lower_score_matrix[rows - funnel_len][window_size];
  int upper_direction_matrix[funnel_len][funnel_size];
  int lower_direction_matrix[rows - funnel_len][window_size];

  // init matrices
  // I don't believe this should make a difference, but it doesn't work if we don't do it
  int k, l;
  for(k = 0; k < rows; k++) {
    if(k < funnel_len) {
      for(l = 0; l < funnel_size; l++) {
        upper_score_matrix[k][l] = 0;
        upper_direction_matrix[k][l] = 0;
      }
    } else {
      for(l = 0; l < window_size; l++) {
        lower_score_matrix[k-funnel_len][l] = 0;
        lower_direction_matrix[k-funnel_len][l] = 0;
      }
    }
  }

  int xoffset[rows];
  xoffset[0] = 0;
  // yoffset is implicitly the complement to the xoffset

  // fill in the matrix
  int y, x;
  int row_max, xlim, prev_xlim, pprev_xlim; // prev: y-1, pprev: y-2
  int qpos, tpos;
  int match, insertion, deletion;
  int xdiff;

  for(y = 0; y < rows; y++) {
    row_max = 0;

    xlim = get_xlim(y, xoffset[y], (y < funnel_len ? funnel_size : window_size), tlen);

    /*
     * loop across diagonal row
     */
    for(x = 0; x < xlim; x++) { // do not evaluate parts of the row in the overflow/padding at the top of the matrix
      qpos = y - xoffset[y] - x;
      tpos = x + xoffset[y];

      if(qpos >= qlen)
        continue;

      if(y > 0) {
        xdiff = xoffset[y] - xoffset[y-1];
        // test for out of bounds for deletion
        prev_xlim = get_xlim(y-1, xoffset[y-1], (y-1 < funnel_len ? funnel_size : window_size), tlen);
        if(x-1+xdiff < 0 || x-1+xdiff >= prev_xlim)
          deletion = LOW;
        else {
          if(y-1 < funnel_len) {
            deletion = upper_score_matrix[y-1][x-1+xdiff] + SCORES[DEL];
          } else {
            deletion = lower_score_matrix[y-1-funnel_len][x-1+xdiff] + SCORES[DEL];
          }
        }
        // test for out of bounds for insertion
        if(x+xdiff < 0 || x+xdiff >= prev_xlim)
          insertion = LOW;
        else {
          if(y-1 < funnel_len) {
            insertion = upper_score_matrix[y-1][x+xdiff] + SCORES[INS];
          } else {
            insertion = lower_score_matrix[y-1-funnel_len][x+xdiff] + SCORES[INS];
          }
        }
      } else {
        deletion = SCORES[DEL];
        insertion = SCORES[INS];
      }

      if(y > 1) {
        xdiff = xoffset[y] - xoffset[y-2];
        // test for out of bounds for match
        pprev_xlim = get_xlim(y-2, xoffset[y-2], (y-2 < funnel_len ? funnel_size : window_size), tlen);
        if(x-1+xdiff < 0 || x-1+xdiff >= pprev_xlim)
          match = LOW;
        else {
          if(y-2 < funnel_len) {
            match = upper_score_matrix[y-2][x-1+xdiff] + (query[qpos] == target[tpos] ? SCORES[MATCH] : SCORES[MISMATCH]);
          } else {
            match = lower_score_matrix[y-2-funnel_len][x-1+xdiff] + (query[qpos] == target[tpos] ? SCORES[MATCH] : SCORES[MISMATCH]);
          }
        }
      } else if(query[qpos] == target[tpos]) {
        match = SCORES[MATCH];
      } else {
        match = SCORES[MISMATCH];
      }

      // set matrix value and direction
      if(match >= insertion && match >= deletion) {
        if(y < funnel_len) {
          upper_score_matrix[y][x] = match;
          upper_direction_matrix[y][x] = (query[qpos] == target[tpos] ? MATCH : MISMATCH);
        } else {
          lower_score_matrix[y-funnel_len][x] = match;
          lower_direction_matrix[y-funnel_len][x] = (query[qpos] == target[tpos] ? MATCH : MISMATCH);
        }
      } else if(insertion >= deletion) {
        if(y < funnel_len) {
          upper_score_matrix[y][x] = insertion;
          upper_direction_matrix[y][x] = INS;
        } else {
          lower_score_matrix[y-funnel_len][x] = insertion;
          lower_direction_matrix[y-funnel_len][x] = INS;
        }
      } else {
        if(y < funnel_len) {
          upper_score_matrix[y][x] = deletion;
          upper_direction_matrix[y][x] = DEL;
        } else {
          lower_score_matrix[y-funnel_len][x] = deletion;
          lower_direction_matrix[y-funnel_len][x] = DEL;
        }
      }

      // keep track of maximum in this row
      if(y < funnel_len) {
        if(upper_score_matrix[y][x] > upper_score_matrix[y][row_max])
          row_max = x;
      } else {
        if(lower_score_matrix[y-funnel_len][x] > lower_score_matrix[y-funnel_len][row_max])
          row_max = x;
      }
    }

    if(y < rows - 1) {

      // shift offset by 1 toward the row max
      // bounded by both ends by edges of the target sequence

      xoffset[y+1] = (y+1 == funnel_len ? xoffset[y]+(funnel_size-window_size)/2 : xoffset[y]);
      if(row_max >= (y<funnel_len-1 ? funnel_size : window_size)/2)
        xoffset[y+1]++;
      if(xoffset[y+1] < 0)
        xoffset[y+1] = 0;
      if(xoffset[y+1] > tlen - (y<funnel_len-1 ? funnel_size : window_size))
        xoffset[y+1] = tlen - (y<funnel_len-1 ? funnel_size : window_size);
    }
  }

  /*
   traceback
  */
  int qstart, tstart, qend, tend;

  // find the start position - this is tricky because our last "row" is not the set of terminal positions
  int max_x = 0, max_y = 0;
  for(y = rows-1; y > -1; y--) {
    for(x = 0; x < (y < funnel_len ? funnel_size : window_size); x++) {
      qpos = y - xoffset[y] - x;
      tpos = x + xoffset[y];
      if((tpos >= 0 && qpos == qlen - 1) || (qpos > 0 && tpos == tlen - 1)) {
        if((max_y == 0 && max_x == 0) || (y < funnel_len ? upper_score_matrix : lower_score_matrix)[y][x] > (max_y < funnel_len ? upper_score_matrix : lower_score_matrix)[max_y][max_x]) {
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

    if(y < funnel_len)
      kv_push(unsigned char, *path, upper_direction_matrix[y][x]);
    else
      kv_push(unsigned char, *path, lower_direction_matrix[y-funnel_len][x]);
    qstart = qpos;
    tstart = tpos;

    if((y<funnel_len ? upper_direction_matrix[y] : lower_direction_matrix[y-funnel_len])[x] == MATCH || (y<funnel_len ? upper_direction_matrix[y] : lower_direction_matrix[y-funnel_len])[x] == MISMATCH) {
      y = y - 2;
      if(y >= 0) { // if it is not, it won't matter what x becomes
        xdiff = xoffset[y+2] - xoffset[y];
        x = x - 1 + xdiff;
      }
    } else if((y<funnel_len ? upper_direction_matrix[y] : lower_direction_matrix[y-funnel_len])[x] == INS) {
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
  res.score = (max_y < funnel_len ? upper_score_matrix : lower_score_matrix)[max_y][max_x];
  res.qstart = qstart;
  res.qend = qend;
  res.tstart = tstart;
  res.tend = tend;
  return res;
}
