Tools for semi-global/overlap/edge alignment in C
=================================================

Included are two variants of Smith-Waterman/Needleman-Wunsch-type dynamic programming alignment. Alignments are restricted to overlaps between the query and target sequences - either end-to-end edge overlap or with one sequence contained within the other.

  1. full_matrix_aligner: computes the entire DP matrix, guaranteed optimal

  2. diagonal_matrix_aligner: computes only a fixed-size diagonal "window" which may shift up and down to track the current highest-scoring path 


For other SW/NW implementations which vary in their features and capabilities (notably, greater code quality and speed), see:


Global, local, and semi-global vectorized (fast) alignment (scores and stats only, no traceback)

https://github.com/jeffdaily/parasail


Vectorized (fast) local alignment, with traceback

https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library


C API
-----

See test.c:

  gcc test.c -o test
  ./test


Python API
----------

  python setup.py install

  import ovlalign
  query = "..."
  target = "..."
  res = ovlalign.full_align(query, target, len(query), len(target))
  res = ovlalign.diag_align(query, target, len(query), len(target), window_size=10)
  print "Score:", res[0]
  print "Query range: %i - %i" % (res[1], res[2])
  print "Target range: %i - %i" % (res[3], res[4])
  print "Path:", path
  print "0: Match"
  print "1: Insertion"
  print "2: Deletion"
  print "3: Mismatch"


License
-------

MIT
