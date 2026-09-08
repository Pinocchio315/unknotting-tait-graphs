"""Exact finite-group algebra for Section 3, ``Linking pairings''.

A nonsingular integral symmetric matrix presents a finite abelian group. Its
inverse presents a linking pairing, up to the boundary-orientation convention.
Lickorish's obstruction tests both signs, so that convention does not affect
the conclusion. Passing this necessary condition does not prove u = 1.
"""
import math

from sympy import Matrix, factorint


def cyclic_generator(matrix, order=None):
    """Return a vector generating coker(matrix), or None if it is noncyclic.

    The order of [x] is the least common denominator of matrix**(-1) x.
    The exponent of the cokernel is therefore the lcm of the column orders.
    It equals the group order exactly when the group is cyclic. CRT
    idempotents combine generators of its primary components; unlike a
    search over a few short vectors this construction is exhaustive.
    """
    matrix = Matrix(matrix)
    if matrix.rows != matrix.cols or matrix != matrix.T:
        raise ValueError('an integral symmetric square matrix is required')
    if any(x.q != 1 for x in matrix):
        raise ValueError('the presentation must be integral')
    determinant = abs(int(matrix.det()))
    if determinant == 0:
        raise ValueError('a finite discriminant group requires nonzero determinant')
    if order is not None and order != determinant:
        raise ValueError('group order does not match the determinant')
    if determinant == 1:
        return Matrix.zeros(matrix.rows, 1)
    inverse = matrix.inv()
    orders = [math.lcm(*(int(t.q) for t in inverse[:, j]))
              for j in range(matrix.cols)]
    if math.lcm(*orders) != determinant:
        return None
    vector = Matrix.zeros(matrix.rows, 1)
    for p, exponent in factorint(determinant).items():
        primary = int(p ** exponent)
        j = next(j for j, d in enumerate(orders) if d % primary == 0)
        complement = determinant // primary
        vector[j] += complement * pow(complement, -1, primary)
    return vector.applyfunc(lambda x: x % determinant)


def generator_self_linking(matrix, order):
    """Return a with lambda(g,g)=a/order using +matrix**(-1), or None."""
    matrix = Matrix(matrix)
    vector = cyclic_generator(matrix, order)
    if vector is None:
        return None
    value = (vector.T * matrix.inv() * vector)[0]
    numerator = value * order
    if numerator.q != 1:
        raise ArithmeticError('self-linking denominator does not divide the group order')
    return int(numerator) % order
