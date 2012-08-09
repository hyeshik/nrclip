#define PY_SSIZE_T_CLEAN
#include "Python.h"
#include <numpy/arrayobject.h>

static PyObject *ErrorObject;

#define MAXREADLENGTH   200
#define NUMBASES        5
#define BASEDELETION    4
static const int nucleobase2int[256] = {
    9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,
    9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,
    9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 4, 9, 9,
    9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,
    9, 0, 9, 1, 9, 9, 9, 2, 9, 9, 9, 9, 9, 9, 9, 9,
    9, 9, 9, 9, 3, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,
    9, 0, 9, 1, 9, 9, 9, 2, 9, 9, 9, 9, 9, 9, 9, 9,
    9, 9, 9, 9, 3, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,
    9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,
    9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,
    9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,
    9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,
    9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,
    9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,
    9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,
    9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9, 9,
};
static const char *nucleobases = "ACGT-";

typedef struct {
    PyObject_HEAD
    unsigned char *src[MAXREADLENGTH][NUMBASES];
                                    /* Randomized source string of read */
    const unsigned char *srcend[MAXREADLENGTH][NUMBASES];
                                /* The last element of src (not inclusive) */
    const unsigned char *srccur[MAXREADLENGTH][NUMBASES];
                                    /* The element to be read next */
} RandomReadCounterObject;

static PyTypeObject RandomReadCounter_Type;

#define RandomReadCounterObject_Check(v)      (Py_TYPE(v) == &RandomReadCounter_Type)

static void
RRC_dealloc(RandomReadCounterObject *self)
{
    int i, j;

    for (i = 0; i < MAXREADLENGTH; i++)
        for (j = 0; j < NUMBASES; j++)
            if (self->src[i][j] != NULL) {
                PyMem_Free(self->src[i][j]);
                self->src[i][j] = NULL;
            }

    Py_TYPE(self)->tp_free(self);
}

static PyObject *
RRC_new(PyTypeObject *type, PyObject *args, PyObject *kwds)
{
    static char *kwargs[] = {"src", NULL};
    RandomReadCounterObject *self;
    PyObject *srcdict;
    size_t arrsize;
    int i, j;

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "O:randomreadcounter", kwargs,
                                     &srcdict))
        return NULL;

    if (!PyDict_Check(srcdict)) {
        PyErr_SetString(PyExc_TypeError, "arg 1 must be a dictionary");
        return NULL;
    }

    self = (RandomReadCounterObject *)type->tp_alloc(type, 0);
    if (self == NULL)
        return NULL;

    arrsize = sizeof(unsigned char *) * MAXREADLENGTH * NUMBASES;
    memset(self->src, 0, arrsize);
    memset(self->srcend, 0, arrsize);
    memset(self->srccur, 0, arrsize);

    for (i = 0; i < MAXREADLENGTH; i++) {
        for (j = 0; j < NUMBASES; j++) {
            PyObject *key, *srcstr;

            key = Py_BuildValue("(ic)", i, nucleobases[j]);
            srcstr = PyDict_GetItem(srcdict, key);
            Py_DECREF(key);

            if (srcstr == NULL) { /* assign uniform reads for undefined  */
                unsigned char *uniform = (unsigned char *)PyMem_Malloc(128);

                self->srccur[i][j] = self->src[i][j] = uniform;
                self->srcend[i][j] = uniform + 127;
                memset(uniform, j, 127);
                uniform[127] = 0;
            }
            else {
                Py_ssize_t poolsize, k;
                unsigned char *converted, *origsrcstr;

                if (!PyString_Check(srcstr)) {
                    PyErr_SetString(PyExc_ValueError, "Value must be strings");
                    goto onError;
                }

                origsrcstr = (unsigned char *)PyString_AS_STRING(srcstr);
                poolsize = PyString_GET_SIZE(srcstr);

                self->src[i][j] = converted = (
                                    (unsigned char *)PyMem_Malloc(poolsize));
                self->srcend[i][j] = converted + poolsize;
                self->srccur[i][j] = converted;

                for (k = 0; k < poolsize; k++) {
                    converted[k] = nucleobase2int[origsrcstr[k]];
                    if (converted[k] >= NUMBASES) {
                        PyErr_Format(PyExc_ValueError,
                                "Illegal base '%c' included in (%d, '%c').",
                                origsrcstr[k], i, nucleobases[j]);
                        goto onError;
                    }
                }
            }
        }
    }

    return (PyObject *)self;

  onError:
    Py_DECREF(self);
    return NULL;
}

static PyObject *
RRC_simulate(RandomReadCounterObject *self, PyObject *args)
{
    Py_ssize_t nreads, refseqsize, i, j, k;
    PyObject *ret;
    char *refseq, r;
    unsigned int *retdata;
    npy_intp readreturndims[2];

    if (!PyArg_ParseTuple(args, "z#n:read", &refseq, &refseqsize, &nreads))
        return NULL;

    if (refseqsize > MAXREADLENGTH) {
        PyErr_SetString(PyExc_ValueError,
            "The length of reference sequence exceeds random pool.");
        return NULL;
    }

    readreturndims[0] = refseqsize;
    readreturndims[1] = NUMBASES;

    refseq = strndup(refseq, refseqsize); /* convert to number-coded sequence in-place */
    for (i = 0; i < refseqsize; i++)
        if ((r = nucleobase2int[(unsigned char)refseq[i]])
                >= NUMBASES) {
            PyErr_Format(PyExc_ValueError,
                         "Illegal base '%c' included in arg 1", refseq[i]);
            free(refseq);
            return NULL;
        }
        else
            refseq[i] = r;

    ret = PyArray_ZEROS(2, readreturndims, NPY_UINT, 0);
    retdata = (unsigned int *)PyArray_DATA(ret);

    for (i = 0; i < nreads; i++) {
        /* Simulate a read iteration */
        for (j = k = 0; j < refseqsize; j++, k++) {
            /* j for refseq position, k for read position */
            int refbase = refseq[j];
            int readgen = *(self->srccur[k][refbase]++);

            if (self->srccur[k][refbase] >= self->srcend[k][refbase])
                self->srccur[k][refbase] = self->src[k][refbase];

            retdata[j * NUMBASES + readgen]++;

            /* lagging refseq index for deletions */
            if (readgen == BASEDELETION)
                k--;
        }
    }

    free(refseq);
    return ret;
}

static PyMethodDef RRC_methods[] = {
    {"simulate",    (PyCFunction)RRC_simulate,  METH_VARARGS,
                    PyDoc_STR("simulate() -> None")},
    {NULL,          NULL}           /* sentinel */
};

static PyTypeObject RandomReadCounter_Type = {
    /* The ob_type field must be initialized in the module init function
     * to be portable to Windows without using C++. */
    PyVarObject_HEAD_INIT(NULL, 0)
    "xxmodule.rnarryext",             /*tp_name*/
    sizeof(RandomReadCounterObject), /*tp_basicsize*/
    0,                          /*tp_itemsize*/
    /* methods */
    (destructor)RRC_dealloc,    /*tp_dealloc*/
    0,                          /*tp_print*/
    0,                          /*tp_getattr*/
    0,                          /*tp_setattr*/
    0,                          /*tp_compare*/
    0,                          /*tp_repr*/
    0,                          /*tp_as_number*/
    0,                          /*tp_as_sequence*/
    0,                          /*tp_as_mapping*/
    0,                          /*tp_hash*/
    0,                          /*tp_call*/
    0,                          /*tp_str*/
    PyObject_GenericGetAttr,    /*tp_getattro*/
    0,                          /*tp_setattro*/
    0,                          /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT,         /*tp_flags*/
    0,                          /*tp_doc*/
    0,                          /*tp_traverse*/
    0,                          /*tp_clear*/
    0,                          /*tp_richcompare*/
    0,                          /*tp_weaklistoffset*/
    0,                          /*tp_iter*/
    0,                          /*tp_iternext*/
    RRC_methods,                /*tp_methods*/
    0,                          /*tp_members*/
    0,                          /*tp_getset*/
    0,                          /*tp_base*/
    0,                          /*tp_dict*/
    0,                          /*tp_descr_get*/
    0,                          /*tp_descr_set*/
    0,                          /*tp_dictoffset*/
    0,                          /*tp_init*/
    PyType_GenericAlloc,        /*tp_alloc*/
    RRC_new,                    /*tp_new*/
    PyObject_Del,               /*tp_free*/
    0,                          /*tp_is_gc*/
};
/* --------------------------------------------------------------------- */

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

    /* Finalize the type object including setting type of the new type
     * object; doing it here is required for portability, too. */
    if (PyType_Ready(&RandomReadCounter_Type) < 0)
        return;

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

    Py_INCREF(&RandomReadCounter_Type);
    PyModule_AddObject(m, "RandomReadCounter",
                       (PyObject *)&RandomReadCounter_Type);
}
