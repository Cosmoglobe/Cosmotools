#include "Python.h"
#include "arrayobject.h"

static PyObject *accumulate_slice1(PyObject *self, PyObject *args);
static PyObject *accumulate_slice2(PyObject *self, PyObject *args);
static PyObject *get_alpha_arr(PyObject *self, PyObject *args);
static PyObject *get_azimuth_arr(PyObject *self, PyObject *args);

#define jz_likely(x)       __builtin_expect((x),1)
#define jz_unlikely(x)     __builtin_expect((x),0)

extern double  fortran_getalpha_(double*,double*,double*);
extern double  fortran_getazimuth_(double*,double*,double*);


static PyMethodDef _quietPyExtMethods[] = {
    {"accumulate_slice1", accumulate_slice1, METH_VARARGS},
	    {"accumulate_slice2", accumulate_slice2, METH_VARARGS},
	    {"get_alpha_arr", get_alpha_arr, METH_VARARGS},
	    {"get_azimuth_arr", get_azimuth_arr, METH_VARARGS},
      {NULL, NULL}     /* Sentinel - marks the end of this structure */
};

void init_quietPyExt()  {
    (void) Py_InitModule("_quietPyExt", _quietPyExtMethods);
    import_array();  // Must be present for NumPy.  Called first after above line.
}




double *pyvector_to_double_ptr(PyArrayObject *arrayin)  {
    int i,n;    
    n=arrayin->dimensions[0];
    return (double *) arrayin->data;  /* pointer to arrayin data as double */
}


int64_t *pyvector_to_int64_ptr(PyArrayObject *arrayin)  {
    int i,n;    
    n=arrayin->dimensions[0];
	return (int64_t *) arrayin->data;  /* pointer to arrayin data as double */
}




/* ==== Check that PyArrayObject is a double (Float) type and a vector ==============
    return 1 if an error and raise exception */ 
int  not_doublevector(PyArrayObject *vec)  {
    if (vec->descr->type_num != NPY_DOUBLE || vec->nd != 1)  {
        PyErr_SetString(PyExc_ValueError,
            "In not_doublevector: array must be of type Float and 1 dimensional (n).");
        return 1;  }
    return 0;
}

/* ==== Check that PyArrayObject is a double (Float) type and a vector ==============
    return 1 if an error and raise exception */ 
int  not_int64vector(PyArrayObject *vec)  {
    if (vec->descr->type_num != NPY_INT64 || vec->nd != 1)  {
        PyErr_SetString(PyExc_ValueError,
            "In not_doublevector: array must be of type int64 and 1 dimensional (n).");
        return 1;  }
    return 0;
}
static PyObject *accumulate_slice1(PyObject *self, PyObject *args){
	PyArrayObject *indexVector, *sourceVector, *targetVector;
	if (!PyArg_ParseTuple(args, "O!O!O!", &PyArray_Type, &targetVector,&PyArray_Type, &indexVector, &PyArray_Type, &sourceVector))  return NULL;
    if (not_doublevector(sourceVector)) return NULL;
    if (not_doublevector(targetVector)) return NULL;
    if (not_int64vector(indexVector)) return NULL;

	npy_intp n = indexVector->dimensions[0];
	npy_intp tmax = targetVector->dimensions[0];

	if (!sourceVector->dimensions[0]==n){
        PyErr_SetString(PyExc_ValueError,
            "Source vector and index vector must have same length.");
		return NULL;
		}

	double *target = pyvector_to_double_ptr(targetVector);
	double *source = pyvector_to_double_ptr(sourceVector);
	int64_t * index = pyvector_to_int64_ptr(indexVector);

	int64_t i,ind;
	for(i=0;i<n;i++){
//		printf("%lld\n",i);
		ind=index[i];
		if (jz_unlikely(ind<0)) ind=tmax+ind;
		if (jz_unlikely(ind<0||ind>=tmax)){
	        PyErr_SetString(PyExc_ValueError,
	            "Index out of range.  Also your array is now corrupted.  Nice one.");
			return NULL;
			}
		target[ind]+=source[i];
	}
	Py_INCREF(Py_None);
	return Py_None;
}

static PyObject *accumulate_slice2(PyObject *self, PyObject *args){
	PyArrayObject *indexVector,*targetVector;
	double source;
	if (!PyArg_ParseTuple(args, "O!O!d", &PyArray_Type, &targetVector,&PyArray_Type, &indexVector,&source))  return NULL;
    if (not_doublevector(targetVector)) return NULL;
    if (not_int64vector(indexVector)) return NULL;

	npy_intp n = indexVector->dimensions[0];
	npy_intp tmax = targetVector->dimensions[0];


	double *target = pyvector_to_double_ptr(targetVector);
	int64_t * index = pyvector_to_int64_ptr(indexVector);

	int64_t i,ind;
	for(i=0;i<n;i++){
		ind=index[i];
		if (jz_unlikely(ind<0)) ind=tmax+ind;
		if (jz_unlikely(ind<0||ind>tmax)){
	        PyErr_SetString(PyExc_ValueError,
	            "Index out of range.  Also your array is now corrupted.  Nice one.");
			return NULL;
			}
		target[ind]+=source;
	}
	
}




static PyObject *get_alpha_arr(PyObject *self, PyObject *args){
	PyArrayObject *psiVector,*thetaVector,*phiVector,*alphaVector;
	int one=1;
	int galactic=0;
	if (!PyArg_ParseTuple(args, "O!O!O!", &PyArray_Type, &phiVector,&PyArray_Type, &thetaVector, &PyArray_Type, &psiVector))  return NULL;


    if (not_doublevector(phiVector)) return NULL;
    if (not_doublevector(thetaVector)) return NULL;
    if (not_doublevector(psiVector)) return NULL;

	npy_intp np = phiVector->dimensions[0];
	int n = (int) np;
	if ((thetaVector->dimensions[0]!=n)||(psiVector->dimensions[0]!=n)){
        PyErr_SetString(PyExc_ValueError,
            "All vectors must have same length.");
		return NULL;
	}
	alphaVector = (PyArrayObject *)PyArray_SimpleNew(one,&np,PyArray_FLOAT64);

	double *phi = pyvector_to_double_ptr(phiVector);
	double *theta = pyvector_to_double_ptr(thetaVector);
	double *psi = pyvector_to_double_ptr(psiVector);
	double *alpha = pyvector_to_double_ptr(alphaVector);


	npy_intp phiStride =  phiVector->strides[0]/sizeof(double);
	npy_intp thetaStride =  thetaVector->strides[0]/sizeof(double);
	npy_intp psiStride =  psiVector->strides[0]/sizeof(double);
	npy_intp alphaStride =  alphaVector->strides[0]/sizeof(double);
	int i;
	double phi_i,theta_i,psi_i;
	double r;
	for(i=0;i<n;i++){
		phi_i = phi[i*phiStride];
		theta_i = theta[i*thetaStride];
		psi_i = psi[i*psiStride];
		alpha[i*alphaStride] =  fortran_getalpha_(&phi_i,&theta_i,&psi_i);
	}
	
	Py_INCREF(alphaVector);
	return alphaVector;
}



static PyObject *get_azimuth_arr(PyObject *self, PyObject *args){
	PyArrayObject *phiVector,*thetaVector,*timeVector,*azimuthVector;
	int one=1;
	int galactic=0;
	if (!PyArg_ParseTuple(args, "O!O!O!", &PyArray_Type, &phiVector,&PyArray_Type, &thetaVector, &PyArray_Type, &timeVector))  return NULL;

    if (not_doublevector(phiVector)) return NULL;
    if (not_doublevector(thetaVector)) return NULL;
    if (not_doublevector(timeVector)) return NULL;

	npy_intp np = phiVector->dimensions[0];
	int n = (int) np;
	if ((thetaVector->dimensions[0]!=n)||(timeVector->dimensions[0]!=n)){
        PyErr_SetString(PyExc_ValueError,
            "All vectors must have same length.");
		return NULL;
	}
	azimuthVector = (PyArrayObject *)PyArray_SimpleNew(one,&np,PyArray_FLOAT64);

	double *phi = pyvector_to_double_ptr(phiVector);
	double *theta = pyvector_to_double_ptr(thetaVector);
	double *mjd = pyvector_to_double_ptr(timeVector);
	double *azimuth = pyvector_to_double_ptr(azimuthVector);


	npy_intp phiStride =  phiVector->strides[0]/sizeof(double);
	npy_intp thetaStride =  thetaVector->strides[0]/sizeof(double);
	npy_intp timeStride =  timeVector->strides[0]/sizeof(double);
	npy_intp azimuthStride =  azimuthVector->strides[0]/sizeof(double);
	int i;
	double phi_i,theta_i,mjd_i;
	double r;
	for(i=0;i<n;i++){
		phi_i = phi[i*phiStride];
		theta_i = theta[i*thetaStride];
		mjd_i = mjd[i*timeStride];
		azimuth[i*azimuthStride] =  fortran_getazimuth_(&phi_i,&theta_i,&mjd_i);
	}
	Py_INCREF(azimuthVector);
	return azimuthVector;
}
