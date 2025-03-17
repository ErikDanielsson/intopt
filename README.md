# intopt

A parallel integer program solver for that won the course competition in the Efficient C course, EDAG01, at LTH, Lund.

To solve an integer program of the form

    maximize    c^T x
    subject to  Ax <= b
                x >= 0
where A is a m x n matrix, b is a m-vector, and c is an n-vector,
give the following input to the solver:

    m n
    c1 c2 ... cn 
    a11 a12 ... a1n
    a21 a22 ... a2n
    ... ... ... ...
    am1 am2 ... amn
    b1 b2 ... bm
Example input is available in **examples/**.
