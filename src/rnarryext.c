#define PY_SSIZE_T_CLEAN
#include "Python.h"
#include <numpy/arrayobject.h>

static PyObject *ErrorObject;

#define NUMBASES        5


static double
shannon_entropy(const unsigned int *counts)
{
    ssize_t i;
    unsigned int maxreads, totalreads;
    double entropy;
    double totalreadsf;
    maxreads = totalreads = 0;

    for (i = 0; i < NUMBASES; i++) {
        totalreads += counts[i];
        if (maxreads < counts[i])
            maxreads = counts[i];
    }  

    if (maxreads == totalreads)
        return 0.;

    totalreadsf = (double)totalreads;
    entropy = 0.;
    for (i = 0; i < NUMBASES; i++) {
        double p = counts[i] / totalreadsf;
        if (p > 0.)
            entropy -= log(p) * p;
    }  

    return entropy;
}

static PyObject *
Py_shannon_entropy(PyObject *self, PyObject *args)
{
    PyObject *reads;
    double entropy;

    if (!PyArg_ParseTuple(args, "O:shannon_entropy", &reads))
        return NULL;

    if (!PyArray_Check(reads)) {
        PyErr_SetString(PyExc_TypeError, "needs a numpy array.");
        return NULL;
    }

    if (PyArray_NDIM(reads) != 1 || *PyArray_DIMS(reads) != NUMBASES ||
            PyArray_DESCR(reads) != PyArray_DescrFromType(NPY_UINT)) {
        PyErr_Format(PyExc_TypeError, "array size must be %d.", NUMBASES);
        return NULL;
    }

    entropy = shannon_entropy(PyArray_DATA(reads));
    return PyFloat_FromDouble(entropy);
}

/* List of functions defined in the module */

static PyMethodDef rnarryext_methods[] = {
    {"shannon_entropy", Py_shannon_entropy, METH_VARARGS,
        PyDoc_STR("shannon_entropy([A, C, G, T, D]) -> float")},
    {NULL, NULL}           /* sentinel */
};

PyDoc_STRVAR(module_doc,
"This is a template module just for instruction.");

/* Initialization function for the module (*must* be called initxx) */

PyMODINIT_FUNC
initrnarryext(void)
{
    PyObject *m;

    import_array();

    /* Create the module and add the functions */
    m = Py_InitModule3("rnarryext", rnarryext_methods, module_doc);
    if (m == NULL)
        return;

    /* Add some symbolic constants to the module */
    if (ErrorObject == NULL) {
        ErrorObject = PyErr_NewException("rnarryext.error", NULL, NULL);
        if (ErrorObject == NULL)
            return;
    }
    Py_INCREF(ErrorObject);
    PyModule_AddObject(m, "error", ErrorObject);
}
