from libc.math cimport fabs
import numpy as np

def c_emd(double[:] p, double[:] q, unsigned int size):
    cdef Py_ssize_t i

    cdef double to_move = 0, cost = 0
    for i in range(size):
        to_move += p[i] - q[i]
        cost += fabs(to_move)

    return cost

def match_graphs(double[:, :] p, double[:, :] q, double[:, :] p1,
                 double[:, :] q1, unsigned int size):
    cdef Py_ssize_t i, j

    for i in range(size):
        for j in range(size):
            if p[i][j] == 0:
                q1[i][j] = 0

            if q[i][j] == 0:
                p1[i][j] = 0
