#define PY_SSIZE_T_CLEAN
#include "Python.h"
#include <string.h>

#if PY_MAJOR_VERSION < 3
inline PyObject* PyModule_GetFilenameObject(PyObject* module)
{
    return PyString_FromString(PyModule_GetFilename(module));
}

inline const char* PyUnicode_AsUTF8AndSize(PyObject* s, void* p)
{
    return PyString_AsString(s);
}
#endif

const char* get_classy_dir()
{
    PyObject* os_path = NULL;
    PyObject* os_path_dirname = NULL;
    PyObject* os_path_join = NULL;
    PyObject* os_path_exists = NULL;
    PyObject* module = NULL;
    PyObject* file = NULL;
    PyObject* rootdir = NULL;
    PyObject* subdir = NULL;
    PyObject* subexists = NULL;
    static char buf[1024] = {0};

    if(strlen(buf))
        return buf;

    module = PyImport_ImportModule("classy");

    if(module)
    {
        file = PyModule_GetFilenameObject(module);

        if(file)
        {
            os_path = PyImport_ImportModule("os.path");
            os_path_dirname = PyObject_GetAttrString(os_path, "dirname");
            os_path_join = PyObject_GetAttrString(os_path, "join");
            os_path_exists = PyObject_GetAttrString(os_path, "exists");

            rootdir = PyObject_CallFunctionObjArgs(os_path_dirname, file, NULL);
            subdir = PyObject_CallFunction(os_path_join, "Os", rootdir, "classy");
            subexists = PyObject_CallFunctionObjArgs(os_path_exists, subdir, NULL);

            if(PyObject_IsTrue(subexists))
                strncpy(buf, PyUnicode_AsUTF8AndSize(subdir, NULL), (sizeof buf)-1);
            else
                strncpy(buf, PyUnicode_AsUTF8AndSize(rootdir, NULL), (sizeof buf)-1);

            Py_XDECREF(subexists);
            Py_XDECREF(subdir);
            Py_XDECREF(rootdir);
            Py_XDECREF(os_path_exists);
            Py_XDECREF(os_path_join);
            Py_XDECREF(os_path_dirname);
            Py_XDECREF(os_path);
        }

        Py_XDECREF(file);
    }

    Py_XDECREF(module);

    if(!strlen(buf))
        buf[0] = '.';

    return buf;
}
