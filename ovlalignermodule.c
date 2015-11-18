#include "ovlaligner.h"
#include "full_matrix_aligner.c"
#include "diagonal_matrix_aligner.c"
#include <Python.h>

static PyObject* diag_wrapper(PyObject * self, PyObject * args) {
  char *query, *target;
  int qlen, tlen;
  int semilocal; // boolean flag
  int window_size;
  if (!PyArg_ParseTuple(args, "ssiiii", &query, &target, &qlen, &tlen, &window_size, &semilocal)) {
    Py_RETURN_NONE;
  }

  charvec path;
  kv_init(path);
  result res = align_diagonal(query, target, qlen, tlen, window_size, &path, semilocal);

  if(res.failed) { // could not align, probably because qlen or tlen was 0
    Py_RETURN_NONE;
  }

  int path_size = kv_size(path);
  PyObject* cigar = PyList_New(path_size);
  //unsigned char cigar[path_size];
  int i;
  for(i = 0; i < path_size; i++) {
    //cigar[i] = kv_pop(path);
    PyList_SetItem(cigar, i, (PyObject*)PyInt_FromLong((int)kv_pop(path)));
  }
  kv_destroy(path);

  PyObject* ret = PyList_New(6);
  PyList_SetItem(ret, 0, (PyObject*)PyInt_FromLong(res.score));
  PyList_SetItem(ret, 1, (PyObject*)PyInt_FromLong(res.qstart));
  PyList_SetItem(ret, 2, (PyObject*)PyInt_FromLong(res.qend));
  PyList_SetItem(ret, 3, (PyObject*)PyInt_FromLong(res.tstart));
  PyList_SetItem(ret, 4, (PyObject*)PyInt_FromLong(res.tend));
  //PyList_SetItem(ret, 5, (PyObject*)PyByteArray_FromStringAndSize(cigar, path_size));
  PyList_SetItem(ret, 5, cigar);

  return ret;
}

static PyObject* full_wrapper(PyObject * self, PyObject * args) {
  char *query, *target;
  int qlen, tlen;
  int semilocal; // boolean flag
  if (!PyArg_ParseTuple(args, "ssiii", &query, &target, &qlen, &tlen, &semilocal)) {
    Py_RETURN_NONE;
  }

  charvec path;
  kv_init(path);
  result res = align_full_matrix(query, target, qlen, tlen, &path, semilocal);

  if(res.failed) { // could not align, probably because qlen or tlen was 0
    Py_RETURN_NONE;
  }

  int path_size = kv_size(path);
  PyObject* cigar = PyList_New(path_size);
  //unsigned char cigar[path_size];
  int i;
  for(i = 0; i < path_size; i++) {
    //cigar[i] = kv_pop(path);
    PyList_SetItem(cigar, i, (PyObject*)PyInt_FromLong((int)kv_pop(path)));
  }
  kv_destroy(path);

  PyObject* ret = PyList_New(6);
  PyList_SetItem(ret, 0, (PyObject*)PyInt_FromLong(res.score));
  PyList_SetItem(ret, 1, (PyObject*)PyInt_FromLong(res.qstart));
  PyList_SetItem(ret, 2, (PyObject*)PyInt_FromLong(res.qend));
  PyList_SetItem(ret, 3, (PyObject*)PyInt_FromLong(res.tstart));
  PyList_SetItem(ret, 4, (PyObject*)PyInt_FromLong(res.tend));
  PyList_SetItem(ret, 5, cigar);

  return ret;
}

static PyMethodDef OvlAlignerMethods[] = {
  { "diag_align", diag_wrapper, METH_VARARGS, "DiagAlign" },
  { "full_align", full_wrapper, METH_VARARGS, "FullAlign" },
  { NULL, NULL, 0, NULL }
};

DL_EXPORT(void) initovlalign(void) {
  Py_InitModule("ovlalign", OvlAlignerMethods);
}
