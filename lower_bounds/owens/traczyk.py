"""Traczyk's e^{i pi/3} criterion (Contemp. Math. 233 (1999); used by Owens for 9_35), in the exact form
needed to complement Owens' Theorem 1 for |sigma| = 2 knots.

Facts (Lickorish–Millett, Jones): V_K(e^{i pi/3}) = eps * (i sqrt3)^d with eps in {+-1} and
d = dim_{Z/3} H_1(Sigma_2(K); Z/3).  Working mod the minimal polynomial t^2 - t + 1 (so t = e^{i pi/3},
i sqrt3 = 2t - 1, (2t-1)^2 = -3): V == s * 3^k (2t-1)^{d0} with d0 in {0,1}, d = 2k + d0, eps = s(-1)^k.

Skein relation at t = e^{i pi/3}:  t^{-1} V(K+) - t V(K-) = i V(K0)  forces, for a single crossing change
(with V(K±) = eps± (i sqrt3)^{d±}):
    d- = d+ + 1  =>  eps+ = -eps-         d+ = d- + 1  =>  eps+ = eps-          d+ = d-  =>  no constraint.
Reading a change of a positive crossing as K+ -> K- and of a negative crossing as K- -> K+:
    d decreases by 1: eps factor +1 for a positive change, -1 for a negative change;
    d increases by 1: eps factor -1 for a positive change, +1 for a negative change;
    d unchanged: no information.
For an unknotting sequence (end value d = 0, eps = +1) in which EVERY change must lower d by one
(i.e. d(K) equals the number of changes), this gives  eps(K) = (-1)^{#negative changes}.

Application to u = 2, sigma(K) = +2 (mirror-normalised): the sign patterns compatible with the signature
are (p, n) = (1, 1) and (0, 2).  If d(K) = 2, both changes lower d, so pattern (0, 2) forces eps(K) = +1.
Hence if d(K) = 2 and eps(K) = -1, the (0, 2) pattern is impossible and only (1, 1) remains — which is
exactly the case covered by Owens' Theorem 1 with n = 1.
"""
from __future__ import annotations


def jones_mod_t2t1(jones: dict[int, int]) -> tuple[int, int]:
    """Reduce V = sum c_e t^e modulo t^2 - t + 1 (exact integers): returns (alpha, beta) with V == alpha+beta t.
    Uses t^2 == t - 1, i.e. the sequence t^e == A_e + B_e t with period 6; t^{-1} == 1 - t."""
    a, b = 0, 0
    for e, c in jones.items():
        # compute t^e mod (t^2 - t + 1): iterate from 0
        A, B = 1, 0
        if e >= 0:
            for _ in range(e % 6):
                A, B = -B, A + B          # (A + Bt) t = At + Bt^2 = At + B(t-1) = -B + (A+B) t
        else:
            for _ in range((-e) % 6):
                A, B = A + B, -A          # (A + Bt) t^{-1} = A(1-t) + B = (A+B) - A t
        a += c * A
        b += c * B
    return a, b


def eps_and_dim3(jones: dict[int, int]) -> tuple[int, int]:
    """(eps, d) with V(e^{i pi/3}) = eps (i sqrt3)^d.  jones = {exponent: coefficient} of the knot."""
    a, b = jones_mod_t2t1(jones)
    # match  s 3^k          (d0 = 0):  (a, b) = (s 3^k, 0)
    # match  s 3^k (2t-1)   (d0 = 1):  (a, b) = (-s 3^k, 2 s 3^k)
    if b == 0:
        s3 = a
        d0 = 0
    elif a * 2 == -b:
        s3 = -a
        d0 = 1
    else:
        raise ValueError(f'V mod t^2-t+1 = {a} + {b} t is not of the form +-3^k (2t-1)^d')
    s = 1 if s3 > 0 else -1
    m = abs(s3)
    if m == 0:
        # The knot evaluation cannot vanish. Without this check, repeatedly
        # dividing zero by three never terminates on malformed input.
        raise ValueError('the Jones evaluation at exp(i*pi/3) must be nonzero')
    k = 0
    while m % 3 == 0:
        m //= 3
        k += 1
    if m != 1:
        raise ValueError(f'|value| {abs(s3)} is not a power of 3')
    return s * (-1) ** k, 2 * k + d0


def mirror_jones(jones: dict[int, int]) -> dict[int, int]:
    return {-e: c for e, c in jones.items()}


def traczyk_excludes_00_2(jones_sigma_pos: dict[int, int]) -> tuple[bool, int, int]:
    """For the sigma = +2 representative: does the criterion exclude the (p, n) = (0, 2) pattern?
    Returns (excluded, eps, d): excluded iff d == 2 and eps == -1."""
    eps, d = eps_and_dim3(jones_sigma_pos)
    return (d == 2 and eps == -1), eps, d
