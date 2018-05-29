/****************************************************************************
 > File Name: linalgmodule.cpp
 > Author: Netcan
 > Blog: http://www.netcan666.com/
 > Mail: netcan1996@gmail.com
 > Created Time: 2018-05-28 -- 10:09
 ****************************************************************************/
#include <string>
#include "../include/linalgmodule.h"
using std::string;
using std::move;


// 模块方法
static PyMethodDef LinAlgMethods[] = {
    {NULL, NULL, 0, NULL}
};

// 模块定义
static PyModuleDef linalgmodule = {
    PyModuleDef_HEAD_INIT,
    "linalg",   /* name of module */
    NULL, /* module documentation, may be NULL */
    -1,       /* size of per-interpreter state of the module,
                 or -1 if the module keeps state in global variables. */
    LinAlgMethods
};


// 向量Vector类 {
static PyMemberDef PyVector_members[] = {
    {NULL}
};

static PyMethodDef PyVector_methods[] = {
    {"T", (PyCFunction)PyVector_T, METH_NOARGS, "change vector type"},
    {"copy", (PyCFunction)PyVector_Copy, METH_NOARGS, "deep copy of vector"},
    {NULL}  /* Sentinel */
};

static PySequenceMethods PyVectorSeq_methods = {
    .sq_length = (lenfunc)PyVector_len,
    .sq_item = (ssizeargfunc)PyVector_GetItem,
    .sq_ass_item = (ssizeobjargproc)PyVector_SetItem,
};

static PyNumberMethods PyVectorNum_methods = {
    .nb_inplace_multiply = (binaryfunc)PyVector_imul,
    .nb_inplace_add = (binaryfunc)PyVector_iadd,
    .nb_inplace_subtract = (binaryfunc)PyVector_isub,
    .nb_multiply = (binaryfunc)PyVector_mul,
    .nb_add = (binaryfunc)PyVector_add,
    .nb_subtract = (binaryfunc)PyVector_sub,
};

static PyTypeObject PyVectorType = {
    PyObject_HEAD_INIT(NULL)
    .tp_name = "linalg.Vector",
    .tp_doc = "Vector objects",
    .tp_basicsize = sizeof(PyVectorObject),
    .tp_itemsize = 0,
    .tp_dealloc = (destructor) PyLinAlg_dealloc,
    .tp_new = PyType_GenericNew,
    .tp_init = (initproc) PyVector_init,
    .tp_members = PyVector_members,
    .tp_methods = PyVector_methods,
    .tp_flags = Py_TPFLAGS_DEFAULT,
    .tp_str = (reprfunc)PyVector_str,
    .tp_repr = (reprfunc)PyVector_str,
    .tp_as_sequence = &PyVectorSeq_methods,
    .tp_as_number = &PyVectorNum_methods,
};

// 函数定义
static bool
isNumber(PyObject *o, bool throwErr) {
    bool ret = (PyLong_Check(o) || PyFloat_Check(o));
    if(throwErr && !ret) PyErr_SetString(PyExc_TypeError, "value should be number");
    return ret;
}

static double
getNumber(PyObject *o) {
    double ret = 0.0;
    if(PyLong_Check(o)) ret = PyLong_AsDouble(o);
    else if(PyFloat_Check(o)) ret = PyFloat_AsDouble(o);
    return ret;
}

// 深复制
static PyVectorObject *
PyVector_Copy(PyVectorObject *self) {
    PyVectorObject *tmp = (PyVectorObject*)PyType_GenericAlloc(&PyVectorType, 0);
    tmp->ob_vector = self->ob_vector;
    return tmp;
}

// *=
static PyObject *
PyVector_imul(PyVectorObject *self, PyObject *arg) {
    if(!arg || ! isNumber(arg)) return NULL;
    Py_XINCREF(self); // 不加这个会崩溃...因为会被释放掉
    self->ob_vector *= getNumber(arg);
    return (PyObject*)self;
}

// +=
static PyObject *
PyVector_iadd(PyVectorObject *self, PyVectorObject *arg) {
    if(!arg || !PyObject_TypeCheck(arg, &PyVectorType)) {
        PyErr_SetString(PyExc_TypeError, "vector add vector");
        return NULL;
    }
	if(self->ob_vector.getSize() != arg->ob_vector.getSize()) {
        PyErr_SetString(PyExc_TypeError, "two vectors' size mismatch");
        return NULL;
    }
    Py_XINCREF(self);
    self->ob_vector += arg->ob_vector;
    return (PyObject*)self;
}

// -=
static PyObject *
PyVector_isub(PyVectorObject *self, PyVectorObject *arg) {
    if(!arg || !PyObject_TypeCheck(arg, &PyVectorType)) {
        PyErr_SetString(PyExc_TypeError, "vector sub vector");
        return NULL;
    }
	if(self->ob_vector.getSize() != arg->ob_vector.getSize()) {
        PyErr_SetString(PyExc_TypeError, "two vectors' size mismatch");
        return NULL;
    }
    Py_XINCREF(self);
    self->ob_vector -= arg->ob_vector;
    return (PyObject*)self;
}

// 乘法，数乘或者内积
static PyVectorObject *
PyVector_mul(PyVectorObject *self, PyObject *arg) {
    PyVectorObject *tmp = PyVector_Copy(self); // refcnt += 1
    if(isNumber(arg, false)) { // 数乘
        PyVector_imul(tmp, arg);
    } else {
        Py_XDECREF(tmp); return NULL;
    }
    Py_XDECREF(tmp); // 记得减去多余的1
    return tmp;
}

// 向量加法
static PyObject *
PyVector_add(PyVectorObject *self, PyVectorObject *arg) {
    PyVectorObject *tmp = PyVector_Copy(self);
    if(! PyVector_iadd(tmp, arg)) {
        Py_XDECREF(tmp); return NULL;
    }
    Py_XDECREF(tmp); // 记得减去多余的1
    return (PyObject*)tmp;
}

// 向量减法
static PyObject *
PyVector_sub(PyVectorObject *self, PyVectorObject *arg) {
    PyVectorObject *tmp = PyVector_Copy(self);
    if(! PyVector_isub(tmp, arg)) {
        Py_XDECREF(tmp); return NULL;
    }
    Py_XDECREF(tmp); // 记得减去多余的1
    return (PyObject*)tmp;
}

static PyObject *
PyVector_GetItem(PyVectorObject *self, Py_ssize_t i) {
    if((size_t)i >= self->ob_vector.getSize()) {
        PyErr_SetString(PyExc_IndexError, "vector index out of range");
        return NULL;
    }
    return PyFloat_FromDouble(self->ob_vector[i]);
}

static int
PyVector_SetItem(PyVectorObject *self, Py_ssize_t i, PyObject *v) {
    if((size_t)i >= self->ob_vector.getSize()) {
        PyErr_SetString(PyExc_IndexError, "vector index out of range");
        return -1;
    }
	if(!isNumber(v)) return -1;
    self->ob_vector[i] = getNumber(v);
    return 0;
}

// 析构函数
static void
PyLinAlg_dealloc(PyObject *self) {
    Py_TYPE(self)->tp_free((PyObject *)self);
}

// 构造函数
static int
PyVector_init(PyVectorObject *self, PyObject *args, PyObject *kwds) {
    PyObject *pList, *pItem;
    if (!PyArg_ParseTuple(args, "O!", &PyList_Type, &pList)) {
        PyErr_SetString(PyExc_TypeError, "parameter must be a list.");
        return -1;
    }
    Py_ssize_t n = PyList_Size(pList);
	if(n <= 0) {
        PyErr_SetString(PyExc_TypeError, "size of list must be greater than 0");
        return -1;
	}
    vector<double> data;
    for(Py_ssize_t i = 0; i < n; ++i) {
        pItem = PyList_GetItem(pList, i);
        if(! isNumber(pItem))  return -1;
		else data.push_back(getNumber(pItem));
    }
    self->ob_vector = Vector(data);
    return 0;
}

// 打印用
static PyObject *
PyVector_str(PyVectorObject *self) {
    static char str[10 << 10], buffer[1<<9];
    char *ps = str, *pb = buffer;

    *ps++ = self->ob_vector.getType() == VecType::Row ? 'R':'C';
    *ps++ = '(';

    for(size_t i = 0; i < self->ob_vector.getSize(); ++i) {
        sprintf(buffer, "%f%c", self->ob_vector[i], i == self->ob_vector.getSize() - 1 ? ')':' ');
        pb = buffer; while(*pb && (*ps++ = *pb++));
    }
    *ps = '\0';
    return Py_BuildValue("s", str);
}

// 向量长度
static Py_ssize_t
PyVector_len(PyVectorObject *self) {
    return self->ob_vector.getSize();
}

// 改变向量行列属性
static PyObject *
PyVector_T(PyVectorObject *self) {
    PyVectorObject *t = PyVector_Copy(self);
    t->ob_vector = self->ob_vector;
    t->ob_vector.setType(t->ob_vector.getType() == VecType::Row ?
                        VecType::Col : VecType::Row);
    return (PyObject*)t;
}

// }

// 矩阵Matrix类 {
static PyMemberDef PyMatrix_members[] = {
    {NULL}
};

static PyMethodDef PyMatrix_methods[] = {
	{"Jacobi", (PyCFunction)PyMatrix_Jacobi, METH_NOARGS, "jacobi method for symmetric matrix's eigenvals"},
	{"inv", (PyCFunction)PyMatrix_inv, METH_NOARGS, "invertible matrix"},
	{"det", (PyCFunction)PyMatrix_det, METH_NOARGS, "determinant of matrix"},
	{"LUsolve", (PyCFunction)PyMatrix_LUsolve, METH_O, "LU solve Ax = b"},
	{"LUdecomp", (PyCFunction)PyMatrix_LUdecomp, METH_NOARGS, "LU decompose"},
    {"copy", (PyCFunction)PyMatrix_Copy, METH_NOARGS, "deep copy of matrix"},
	{"getRowSize", (PyCFunction)PyMatrix_getRowSize, METH_NOARGS, "get matrix row size"},
	{"getColSize", (PyCFunction)PyMatrix_getColSize, METH_NOARGS, "get matrix col size"},
	{"getNRowVec", (PyCFunction)PyMatrix_getNRowVec, METH_O, "get nth row vector of matrix"},
	{"getNColVec", (PyCFunction)PyMatrix_getNColVec, METH_O, "get nth col vector of matrix"},
    {NULL}  /* Sentinel */
};

static PyTypeObject PyMatrixType = {
    PyObject_HEAD_INIT(NULL)
    .tp_name = "linalg.Matrix",
    .tp_doc = "Matrix objects",
    .tp_basicsize = sizeof(PyMatrixObject),
    .tp_itemsize = 0,
    .tp_dealloc = (destructor) PyLinAlg_dealloc,
    .tp_new = PyType_GenericNew,
    .tp_init = (initproc) PyMatrix_init,
    .tp_members = PyMatrix_members,
    .tp_methods = PyMatrix_methods,
    .tp_flags = Py_TPFLAGS_DEFAULT,
    .tp_str = (reprfunc)PyMatrix_str,
    .tp_repr = (reprfunc)PyMatrix_str,
};

static PyVectorObject *
PyMatrix_Jacobi(PyMatrixObject *self) {
	if(! self->ob_matrix.isSquare() || ! self->ob_matrix.isSymmetric()) {
		PyErr_SetString(PyExc_ValueError, "matrix is not square or is not symmetric");
		return NULL;
	}
	PyMatrixObject * tmp = PyMatrix_Copy(self);
    PyVectorObject *ret = (PyVectorObject*)PyType_GenericAlloc(&PyVectorType, 0);
	ret->ob_vector = move(tmp->ob_matrix.Jacobi());
	Py_XDECREF(tmp);
	return ret;
}

static PyMatrixObject *
PyMatrix_inv(PyMatrixObject *self) {
	if(! self->ob_matrix.isSquare() || isZero(self->ob_matrix.det())) {
		PyErr_SetString(PyExc_ValueError, "matrix is not square or is singular matrix");
		return NULL;
	}
	PyMatrixObject *ret = PyMatrix_Copy(self);
	ret->ob_matrix = move(ret->ob_matrix.inv());
	return ret;
}

static PyObject *
PyMatrix_det(PyMatrixObject *self) {
	if(! self->ob_matrix.isSquare()) {
		PyErr_SetString(PyExc_ValueError, "matrix is not square");
		return NULL;
	}
	double det = self->ob_matrix.det();
	return PyFloat_FromDouble(isZero(det)?0.0:det); // 美化0 = =
}

static PyVectorObject *
PyMatrix_LUsolve(PyMatrixObject *self, PyVectorObject *y) {
    if(!PyObject_TypeCheck(y, &PyVectorType)) {
		PyErr_SetString(PyExc_TypeError, "y must be linalg.Vector");
		return NULL;
	}
	if(y->ob_vector.getSize() != self->ob_matrix.getRowSize()) {
		PyErr_SetString(PyExc_ValueError, "matrix row size must be equal y size");
		return NULL;
	}
	if(! self->ob_matrix.isSquare() || isZero(self->ob_matrix.det())) {
		PyErr_SetString(PyExc_ValueError, "matrix is not square or is singular matrix");
		return NULL;
	}
    PyVectorObject *ret = (PyVectorObject*)PyType_GenericAlloc(&PyVectorType, 0);
	ret->ob_vector = move(self->ob_matrix.LUsolve(y->ob_vector));
	return ret;
}


static PyMatrixObject *
PyMatrix_LUdecomp(PyMatrixObject *self) {
	PyMatrixObject *ret = PyMatrix_Copy(self);
	ret->ob_matrix = move(ret->ob_matrix.LUdecomp());
	return ret;
}

static PyMatrixObject *
PyMatrix_Copy(PyMatrixObject *self) {
    PyMatrixObject *tmp = (PyMatrixObject*)PyType_GenericAlloc(&PyMatrixType, 0);
    tmp->ob_matrix = self->ob_matrix;
    return tmp;
}

static PyObject*
PyMatrix_getRowSize(PyMatrixObject *self) {
	return PyLong_FromLong(self->ob_matrix.getRowSize());
}

static PyObject*
PyMatrix_getColSize(PyMatrixObject *self) {
	return PyLong_FromLong(self->ob_matrix.getColSize());
}

static PyVectorObject*
PyMatrix_getNRowVec(PyMatrixObject *self, PyObject* i) {
	if(! PyLong_Check(i)) {
		PyErr_SetString(PyExc_TypeError, "n must be number");
		return NULL;
	}
	size_t ii = getNumber(i);
	if(ii >= self->ob_matrix.getRowSize()) {
		PyErr_SetString(PyExc_IndexError, "i must be less matrix row size");
		return NULL;
	}

    PyVectorObject *vec = (PyVectorObject*)PyType_GenericAlloc(&PyVectorType, 0);
	vec->ob_vector = move(self->ob_matrix.getNRowVec(ii));
	return vec;
}

static PyVectorObject*
PyMatrix_getNColVec(PyMatrixObject *self, PyObject* i) {
	if(! PyLong_Check(i)) {
		PyErr_SetString(PyExc_TypeError, "n must be number");
		return NULL;
	}
	size_t ii = getNumber(i);
	if(ii >= self->ob_matrix.getColSize()) {
		PyErr_SetString(PyExc_IndexError, "i must be less matrix row size");
		return NULL;
	}

    PyVectorObject *vec = (PyVectorObject*)PyType_GenericAlloc(&PyVectorType, 0);
	vec->ob_vector = move(self->ob_matrix.getNColVec(ii));
	return vec;
}

static int
PyMatrix_init(PyMatrixObject *self, PyObject *args, PyObject *kwds) {
    PyObject *pList, *pItem, *pVal;
    if (!PyArg_ParseTuple(args, "O!", &PyList_Type, &pList)) {
        PyErr_SetString(PyExc_TypeError, "parameter must be a list of list.");
        return -1;
    }
    Py_ssize_t m = PyList_Size(pList);
	if(m <= 0) {
        PyErr_SetString(PyExc_TypeError, "size of list must be greater than 0");
        return -1;
	}
	pItem = PyList_GetItem(pList, 0);
	if(! PyObject_TypeCheck(pItem, &PyList_Type)) {
        PyErr_SetString(PyExc_TypeError, "parameter must be a list of list.");
        return -1;
	}

    Py_ssize_t n = PyList_Size(pItem);
	if(n <= 0) {
        PyErr_SetString(PyExc_TypeError, "size of list must be greater than 0");
        return -1;
	}

    vector<vector<double>> data(m, vector<double>(n, 0));
	for(Py_ssize_t i = 0; i < m; ++i) {
		pItem = PyList_GetItem(pList, i);
		if(! PyObject_TypeCheck(pItem, &PyList_Type)) {
			PyErr_SetString(PyExc_TypeError, "parameter must be a list of list.");
			return -1;
		}
		if(PyList_Size(pItem) != n) {
			PyErr_SetString(PyExc_TypeError, "list of list must be matrix");
			return -1;
		}

		for(Py_ssize_t j = 0; j < n; ++j) {
			pVal = PyList_GetItem(pItem, j);
			if(! isNumber(pVal)) return -1;
			data[i][j] = getNumber(pVal);
		}
	}
	self->ob_matrix = Matrix(data);
	return 0;
}

static PyObject *
PyMatrix_str(PyMatrixObject *self) {
    static char str[100 << 10], buffer[1<<9];
    char *ps = str, *pb = buffer;
	for(size_t i = 0; i < self->ob_matrix.getRowSize(); ++i) {
		for(size_t j = 0; j < self->ob_matrix.getColSize(); ++j) {
			sprintf(buffer, "%10.6f ", self->ob_matrix.at(i, j));
			pb = buffer; while(*pb && (*ps++ = *pb++));
		}
		*ps++ = '\n';
	}
	*ps = '\0';
    return Py_BuildValue("s", str);
}

// }


// 初始化模块
PyMODINIT_FUNC
PyInit_linalg(void) {
    PyObject *m;

    m = PyModule_Create(&linalgmodule);
    if (m == NULL) return NULL;


    // 添加Vector/Matrix类
    if (PyType_Ready(&PyVectorType) < 0) return NULL;
    if (PyType_Ready(&PyMatrixType) < 0) return NULL;
    Py_INCREF(&PyVectorType);
    Py_INCREF(&PyMatrixType);
    PyModule_AddObject(m, "Vector", (PyObject*)&PyVectorType);
    PyModule_AddObject(m, "Matrix", (PyObject*)&PyMatrixType);

    return m;
}

