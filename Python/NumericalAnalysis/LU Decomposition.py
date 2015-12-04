'''Computes LU Decomposition of Matrix using Numpy.  - Might be Crout or Doolittle dunno'''
from __future__ import division
import numpy as np


def row_switch(A, indices):
    '''takes a matrix and 2 row indices and switches those cols.
    Also returns associated permutation matrix'''
    P = np.eye(len(np.diag(A)))
    # Permute the Matrix
    first_row = np.copy(A[indices[0], :])
    A[indices[0], :] = A[indices[1], :]
    A[indices[1], :] = first_row
    # Permute the Identity
    id_first_row = np.copy(P[indices[0], :])
    P[indices[0], :] = P[indices[1], :]
    P[indices[1], :] = id_first_row
    return [P, A]


def test_row_switch(noisy=False):
    A = np.array([[1., 3., 5.], [2., 4., 7.], [1., 1., 0.]])
    old_A = np.copy(A)
    [P, new_A] = row_switch(A, [0, 1])
    if noisy:
        print new_A, '\nfrom\n', old_A, '\nthrough\n', P
    assert not (old_A - P.dot(new_A)).any()
    print 'row switch OK'


def partial_pivot(A, noisy=False):
    '''runs through a matrix and switches largest elements of
    rows onto the diagonal'''
    P = np.eye(len(np.diag(A)))
    for col in enumerate(A):
        maxelem = max(col[1])
        for elem in enumerate(col[1]):
            if elem[1] == maxelem:
                if noisy:
                    print 'switching %d to (%d, %d), maxelem: %d'%(
                        elem[1], col[0], elem[0], maxelem)
                    print col[1], '-> '
                [inc_P, A] = row_switch(A, [col[0], elem[0]])
                P = inc_P.dot(P)
                if noisy:
                    print col[1]
    if noisy:
        print 'P: \n', P, '\nA: \n', A
    return [P, A]


def test_partial_pivot(noisy=False):
    A = np.array([[1, 3, 5], [2, 4, 7], [1, 1, 0]])
    B = np.array([[3, 2, 5], [2, 3, 6], [10, 2, 1]])
    C = np.random.rand(5, 5)
    old_A = np.copy(A)
    old_B = np.copy(B)
    old_C = np.copy(C)
    [P1, new_A] = partial_pivot(A)
    [P2, new_B] = partial_pivot(B)
    [P3, new_C] = partial_pivot(C)
    if noisy:
        print new_A.dot(P1)
        print '\nPivoted:\n', old_A, '\n=>\n', new_A, '\nWith:\n', P1, '\n'
    assert not (new_A-P1.dot(old_A)).any()
    assert not (new_B-P2.dot(old_B)).any()
    assert not (new_C-P3.dot(old_C)).any()
    print 'Partial pivot OK'


def PLU_decomposition(A, noisy=False):
    '''performs LU decomposition with partial pivoting'''
    L = np.identity(len(np.diag(A)))
    U = np.zeros_like(A)
    [P, A] = partial_pivot(A, False)
    if noisy:
        print A
    for row in enumerate(A):
        if noisy:
            print 'row: ', row[1]
        # row has indices in first place and row in second
        for col in enumerate(row[1]):
            # col has column index in first place and element in second
            i = row[0]
            j = col[0]
            if noisy:
                print 'elem: ', col[1], ' at ', (i, j)
            if i <= j:
                U[i, j] = A[i, j] - np.sum(
                    U[k, j] * L[i, k] for k in range(i))
            if i > j:
                L[i, j] = (1 / U[j, j]) * (
                    A[i, j] - np.sum(
                        U[k, j] * L[i, k] for k in range(j)))
            if noisy:
                print "U\n", U
                print "L\n", L
    return [P, L, U]


def test_PLU_decomposition(noisy=False):
    A = np.array([[1., 3., 5.], [2., 4., 7.], [1., 1., 0.]])
    old_A = np.copy(A)
    [P1, L1, U1] = PLU_decomposition(A)
    B = np.array(
        [[11., 9., 24., 2.],
         [1., 5., 2., 6.],
         [3., 17., 18., 1.],
         [2., 5., 7., 1.]])
    old_B = np.copy(B)
    [P2, L2, U2] = PLU_decomposition(B)
    C = np.random.rand(5, 5)
    old_C = np.copy(C)
    [P3, L3, U3] = PLU_decomposition(C)
    if noisy:
        print "P=\n", P1, "\nL=\n", L1, '\nU=\n', U1
        if (P1.T.dot(L1).dot(U1)-old_A).any():
            print "failed\n", P1.T.dot(L1).dot(U1), '\ndoesnt equal\n', old_A
        else:
            print "succeeded: P^T L U = \n", P1.T.dot(L1).dot(U1), '\nis equal to A\n', old_A
    print
    assert np.allclose(P2.T.dot(L2).dot(U2), old_B)
    assert np.allclose(P2.T.dot(L2).dot(U2), old_B)
    assert np.allclose(P3.T.dot(L3).dot(U3), old_C)
    print 'PLU decomposition OK'


def run_tests(noisy=False):
    test_PLU_decomposition(noisy)
    test_partial_pivot(noisy)
    test_row_switch(noisy)
    print 'All tests OK\n'

run_tests(True)

