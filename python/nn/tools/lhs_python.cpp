#include <stdlib.h>
#include <stdio.h>
#include <errno.h>
#include <time.h>
#include <stdint.h>
#undef NDEBUG
#include <assert.h>
#include <omp.h>

#include <Python.h>
#include <numpy/ndarrayobject.h>

typedef uint32_t u32;
typedef uint64_t u64;
typedef float f32;

static void init(u32 *arr, size_t D, size_t N);
static void init_f32(f32 *arr, size_t D, size_t N);
template<typename T> static void shuffle(T *arr, size_t D, size_t N);
static PyObject *method_lhs(PyObject *self, PyObject *args);
static PyObject *method_lhs_float(PyObject *self, PyObject *args);

static PyMethodDef lhsMethods[] = {
    {"lhs", method_lhs, METH_VARARGS, "latin hypercube sampling (integer values)"},
    {"lhs_float", method_lhs_float, METH_VARARGS, "latin hypercube sampling (float values)"},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef lhsmodule = {
    PyModuleDef_HEAD_INIT,
    "lhs",
    "latin hypercube sampling",
    -1,
    lhsMethods
};

PyMODINIT_FUNC PyInit_lhs() {
    import_array();
    return PyModule_Create(&lhsmodule);
}

static PyObject *method_lhs(PyObject *self, PyObject *args) {
    size_t D, N;

    if (!PyArg_ParseTuple(args, "nn", &D, &N)) {
        return NULL;
    }

    u32 *arr = (u32*) malloc(D * N * sizeof(u32));
    if (arr == NULL) {
        return PyErr_NoMemory();
    }

    init(arr, D, N);
    shuffle(arr, D, N);

    const npy_intp shape[2] = {D, N};
    PyObject *result = PyArray_SimpleNewFromData(2, shape, NPY_INT32, (void*) arr);
    PyArray_ENABLEFLAGS((PyArrayObject*)result, NPY_ARRAY_OWNDATA);

    return result;
}

static PyObject *method_lhs_float(PyObject *self, PyObject *args) {
    size_t D, N;

    if (!PyArg_ParseTuple(args, "nn", &D, &N)) {
        return NULL;
    }

    f32 *arr = (f32*) malloc(D * N * sizeof(f32));
    if (arr == NULL) {
        return PyErr_NoMemory();
    }

    init_f32(arr, D, N);
    shuffle(arr, D, N);

    const npy_intp shape[2] = {D, N};
    PyObject *result = PyArray_SimpleNewFromData(2, shape, NPY_FLOAT32, (void*) arr);
    PyArray_ENABLEFLAGS((PyArrayObject*)result, NPY_ARRAY_OWNDATA);

    return result;
}

static void init(u32 *arr, size_t D, size_t N) {
    u32 *cur = arr;
    for (size_t d = 0; d < D; ++d) {
        for (size_t i = 0; i < N; ++i) {
            *cur++ = i;
        }
    }
}

static void init_f32(f32 *arr, size_t D, size_t N) {
    f32 *cur = arr;
    const f32 delta = 1.0f / (N - 1);
    for (size_t d = 0; d < D; ++d) {
        for (size_t i = 0; i < N; ++i) {
            *cur++ = i * delta;
        }
    }
}

template<typename T>
static void shuffle(T *arr, size_t D, size_t N) {
    T tmp;
    size_t i, j;
#pragma omp parallel for private(i, j, tmp) schedule(static)
    for (size_t d = 0; d < D; ++d) {
        uint64_t x64 = 88172645463325252ull ^ d ^ time(NULL);
        T *row = arr + d * N;
        for (i = N - 1; i > 0; i--) {
            x64 ^= x64 << 13;
            x64 ^= x64 >> 17;
            x64 ^= x64 << 5;

            j = x64 % (i + 1);

            tmp = row[j];
            row[j] = row[i];
            row[i] = tmp;
        }
    }
}
