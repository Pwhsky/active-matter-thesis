# cython_functions.pyx
import numpy as np
cimport numpy as np


def accelerate_histogram2d(np.ndarray[np.float64_t, ndim=1] x,
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
