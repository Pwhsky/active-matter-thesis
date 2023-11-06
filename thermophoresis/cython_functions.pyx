# cython_functions.pyx
import numpy as np
cimport numpy as np


def histogram2d_cython(np.ndarray[np.float64_t, ndim=1] x,
                           np.ndarray[np.float64_t, ndim=1] z,
                           np.ndarray[np.float64_t, ndim=1] grad,
                           np.ndarray[np.float64_t, ndim=1] x_bins,
                           np.ndarray[np.float64_t, ndim=1] y_bins):
    cdef np.ndarray[np.float64_t, ndim=2] H
    cdef np.ndarray[np.float64_t, ndim=2] H_counts
    cdef int i, j, k

    H, xedges, yedges = np.histogram2d(x, z, bins=[x_bins, y_bins], weights=grad)
    H_counts, _, _ = np.histogram2d(x, z, bins=[x_bins, y_bins])

    for j in range(len(x_bins) - 1):
        for k in range(len(y_bins) - 1):
            if H_counts[j, k] > 0:
                H[j, k] /= H_counts[j, k]

    return H, xedges, yedges
    
def gradient_cython(np.ndarray[double, ndim=2] arr, double dx, double dy):
    cdef np.ndarray[double, ndim=2] grad_x = np.empty_like(arr)
    cdef np.ndarray[double, ndim=2] grad_y = np.empty_like(arr)

    cdef int nx = arr.shape[0]
    cdef int ny = arr.shape[1]

    cdef int i, j
    for i in range(1, nx - 1):
        for j in range(1, ny - 1):
            grad_x[i, j] = (arr[i + 1, j] - arr[i - 1, j]) / (2 * dx)
            grad_y[i, j] = (arr[i, j + 1] - arr[i, j - 1]) / (2 * dy)

    return grad_x, grad_y
    
    
def transpose_cython(np.ndarray[double, ndim=2] arr):
    return np.ascontiguousarray(arr.T)
