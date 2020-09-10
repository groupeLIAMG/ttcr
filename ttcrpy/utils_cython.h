

#include <cassert>
#include <vector>

#include "Python.h"
#include <numpy/arrayobject.h>

#include "ttcr_t.h"

#ifndef _utils_cython_h_
#define _utils_cython_h_

template<typename T>
int build_matrix_siv(const size_t M,
                     const size_t N,
                     const std::vector<std::vector<ttcr::siv<T>>>& m_data,
                     PyObject* m_tuple) {

    assert( M == m_data.size() );
    import_array();

    size_t nnz = 0;
    for ( size_t i=0; i<m_data.size(); ++i ) {
        nnz += m_data[i].size();// number of non-zero elements per row.
    }
    npy_intp dims[] = {static_cast<npy_intp>(nnz)};
    double* data_p = new double[nnz];

    PyObject* data = PyArray_SimpleNewFromData(1, dims, NPY_DOUBLE, data_p);

    int64_t* indices_p = new int64_t[nnz];
    PyObject* indices = PyArray_SimpleNewFromData(1, dims, NPY_INT64, indices_p);

    dims[0] = M+1;
    int64_t* indptr_p = new int64_t[M+1];
    PyObject* indptr = PyArray_SimpleNewFromData(1, dims, NPY_INT64, indptr_p);

    size_t k = 0;
    for ( size_t i=0; i<M; ++i ) {
        indptr_p[i] = k;
        for ( size_t j=0; j<N; ++j ) {
            for ( size_t nn=0; nn<m_data[i].size(); ++nn) {
                if ( m_data[i][nn].i == j ) {
                    indices_p[k] = j;
                    data_p[k] = m_data[i][nn].v;
                    k++;
                }
            }
        }
    }
    indptr_p[M] = k;

    PyTuple_SetItem(m_tuple, 0, data);
    PyTuple_SetItem(m_tuple, 1, indices);
    PyTuple_SetItem(m_tuple, 2, indptr);

    return 0;
}

#endif
