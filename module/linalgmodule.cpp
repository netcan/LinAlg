/****************************************************************************
 > File Name: linalgmodule.cpp
 > Author: Netcan
 > Blog: http://www.netcan666.com/
 > Mail: netcan1996@gmail.com
 > Created Time: 2018-05-28 -- 10:09
 ****************************************************************************/

#include <string>
#include "linalgmodule.h"
using std::string;

static PyObject *LinAlgError;

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
	.nb_multiply = (binaryfunc)PyVector_mul,
	.nb_add = (binaryfunc)PyVector_add,
	.nb_subtract = (binaryfunc)PyVector_sub,
};

static PyTypeObject PyVectorType = {
    PyVarObject_HEAD_INIT(NULL, 0)
    .tp_name = "linalg.Vector",
    .tp_doc = "Vector objects",
    .tp_basicsize = sizeof(PyVectorObject),
    .tp_itemsize = 0,
	.tp_dealloc = (destructor) PyVector_dealloc,
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
isNumber(PyObject *o) {
	bool ret = (PyLong_Check(o) || PyFloat_Check(o));
	if(!ret) PyErr_SetString(PyExc_TypeError, "value should be number");
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
	PyVectorObject *tmp = (PyVectorObject*)PyType_GenericAlloc(&PyVectorType, 1);
	tmp->ob_vector = self->ob_vector;
	return tmp;
}

// *=
static PyObject *
PyVector_imul(PyVectorObject *self, PyObject *arg) {
	if(! isNumber(arg)) {
		return NULL;
	}
	self->ob_vector *= getNumber(arg);
	return (PyObject*)self;
}

// 乘法
static PyVectorObject *
PyVector_mul(PyVectorObject *self, PyObject *arg) {
	PyVectorObject *tmp = PyVector_Copy(self);
	PyVector_imul(tmp, arg);
	return tmp;
}

// 向量加法
static PyObject *
PyVector_add(PyVectorObject *self, PyVectorObject *arg) {
	if(!PyObject_TypeCheck(arg, &PyVectorType)) {
		PyErr_SetString(PyExc_IndexError, "vector add vector");
		return NULL;
	} else if(self->ob_vector.getSize() != arg->ob_vector.getSize()) {
		PyErr_SetString(PyExc_IndexError, "two vectors' size mismatch");
		return NULL;
	}
	PyVectorObject *tmp = PyVector_Copy(self);
	tmp->ob_vector += arg->ob_vector;
	return (PyObject*)tmp;
}
// 向量减法
static PyObject *
PyVector_sub(PyVectorObject *self, PyVectorObject *arg) {
	return PyVector_add(self, PyVector_mul(arg, Py_BuildValue("i", -1)));
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
	} else if(!isNumber(v)) {
		return -1;
	}
	self->ob_vector[i] = getNumber(v);
	return 0;
}

// 析构函数
static void
PyVector_dealloc(PyVectorObject *self) {
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
	vector<double> data;
	for(Py_ssize_t i = 0; i < n; ++i) {
		pItem = PyList_GetItem(pList, i);
		if(! isNumber(pItem)) {
			PyErr_SetString(PyExc_TypeError, "list items must be numbers.");
			return -1;
		} else data.push_back(getNumber(pItem));
	}
	self->ob_vector = Vector(data);
	return 0;
}

// 打印用
static PyObject *
PyVector_str(PyVectorObject *self) {
	char str[10 << 10], buffer[1<<10];
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

// 初始化模块
PyMODINIT_FUNC
PyInit_linalg(void) {
    PyObject *m;

    m = PyModule_Create(&linalgmodule);
    if (m == NULL) return NULL;

    LinAlgError = PyErr_NewException("linalg.error", NULL, NULL);
    Py_INCREF(LinAlgError);

	// 添加Vector类
	if (PyType_Ready(&PyVectorType) < 0) return NULL;
	Py_INCREF(&PyVectorType);
	PyModule_AddObject(m, "Vector", (PyObject*)&PyVectorType);

    PyModule_AddObject(m, "error", LinAlgError);
    return m;
}

