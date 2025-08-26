
from __future__ import annotations
import math
from typing import Iterable, List, Set, Tuple, Dict, Iterator
from .wheel import MOD, COPRIME_RESIDUES_210, COPRIME_RESIDUES_210_ORDERED
from .ptp import ptp_residues

# --------- helpers ---------
def _simple_sieve(n: int) -> List[int]:
    if n < 2:
        return []
    sieve = bytearray(b"\x01") * (n + 1)
    sieve[0:2] = b"\x00\x00"
    p = 2
    while p * p <= n:
        if sieve[p]:
            start = p * p
            step = p
            sieve[start:n+1:step] = b"\x00" * (((n - start)//step) + 1)
        p += 1
    return [i for i in range(n + 1) if sieve[i]]

def _is_square(n: int) -> bool:
    r = int(n**0.5)
    return r*r == n

def _miller_rabin_base(n: int, a: int) -> bool:
    # one base MR
    d = n - 1
    s = 0
    while d % 2 == 0:
        d //= 2
        s += 1
    x = pow(a, d, n)
    if x == 1 or x == n-1:
        return True
    for _ in range(s - 1):
        x = (x * x) % n
        if x == n - 1:
            return True
    return False

def _jacobi(a: int, n: int) -> int:
    if n <= 0 or n % 2 == 0:
        return 0
    a = a % n
    result = 1
    while a != 0:
        while a % 2 == 0:
            a //= 2
            r = n % 8
            if r in (3,5):
                result = -result
        a, n = n, a
        if a % 4 == 3 and n % 4 == 3:
            result = -result
        a %= n
    return result if n == 1 else 0

def _lucas_selfridge_params(n: int) -> Tuple[int,int,int]:
    # find D such that Jacobi(D/n) = -1, with D = 5, -7, 9, ...
    D = 5
    while True:
        j = _jacobi(D, n)
        if j == -1:
            break
        if j == 0 and abs(D) < n:
            return (0,0,0)
        D = -D + 2 if D > 0 else -D + 2
    P = 1
    Q = (1 - D) // 4
    return (D, P, Q)

def _lucas_prp(n: int) -> bool:
    # Strong Lucas probable prime test (Selfridge method)
    if n % 2 == 0:
        return n == 2
    D, P, Q = _lucas_selfridge_params(n)
    if D == 0:
        return False
    # n - Jacobi(D/n) must be > 0 even
    # write n+1 = d*2^s
    d = n + 1
    s = 0
    while d % 2 == 0:
        d //= 2
        s += 1
    # Lucas sequences
    U = 0
    V = 2
    Qk = 1
    # compute U_d, V_d
    bits = bin(d)[2:]
    for bit in bits:
        U, V = (U * V) % n, (V * V - 2 * Qk) % n
        Qk = (Qk * Qk) % n
        if bit == '1':
            U, V = ( (P*U + V) % n, ( (D*U + P*V) ) % n )
            Qk = (Qk * Q) % n
    if U % n == 0 or V % n == 0:
        return True
    for _ in range(s - 1):
        V = (V * V - 2 * Qk) % n
        Qk = (Qk * Qk) % n
        if V % n == 0:
            return True
    return False

def is_probable_prime_bpsw(n: int) -> bool:
    """Baillie–PSW probable prime test (no known counterexample)."""
    if n < 2:
        return False
    small_primes = [2,3,5,7,11,13,17,19,23,29,31,37]
    for p in small_primes:
        if n == p:
            return True
        if n % p == 0:
            return False
    if _is_square(n):
        return False
    if not _miller_rabin_base(n, 2):
        return False
    return _lucas_prp(n)

# --------- segmented 210-wheel sieve ---------
def _eligible_residue_set(seeds: int | None, cover_all: bool) -> Set[int]:
    if cover_all:
        return set(COPRIME_RESIDUES_210_ORDERED)
    if seeds is None or seeds <= 0:
        return set()
    # dedup PTP residues
    rs = ptp_residues(seeds, filter_coprime=True)
    seen = set()
    out = []
    for r in rs:
        if r not in seen:
            seen.add(r)
            out.append(r)
    return set(out)

def primes_segmented(limit: int, *, seeds: int = 0, cover_all: bool = False,
                     segment_size: int = 1_000_000) -> List[int]:
    """
    Segmented sieve over [2, limit], using a 210-wheel prefilter.
    - If cover_all=True: returns ALL primes <= limit (proper sieve).
    - Else: returns only primes in the residue classes selected by PTP seeds.
    """
    if limit < 2:
        return []
    base = _simple_sieve(int(limit**0.5) + 1)
    # we'll include 2,3,5,7 unconditionally
    primes = []
    # residue selection
    scan_residues = _eligible_residue_set(seeds, cover_all)
    # process in segments
    low = 2
    while low <= limit:
        high = min(low + segment_size - 1, limit)
        size = high - low + 1
        mark = bytearray(b"\x01") * size  # True means "potentially prime"
        # presieve by smallest primes: cross off multiples starting at p*p
        for p in base:
            start = max(p*p, ((low + p - 1)//p)*p)
            if start > high:
                continue
            step = p
            # mark composites
            for m in range(start, high+1, step):
                mark[m - low] = 0
        # collect primes in this segment
        for i in range(size):
            n = low + i
            if n < 2:
                continue
            if mark[i]:
                if cover_all:
                    primes.append(n)
                else:
                    r = n % MOD
                    if r in scan_residues:
                        primes.append(n)
        low = high + 1
    # ensure small primes are included even if residues omitted
    for p in (2,3,5,7):
        if p <= limit and p not in primes and (cover_all or p % MOD in scan_residues or p in (2,3,5,7)):
            primes.insert(0, p)
    primes.sort()
    return primes

def primes_segmented_bpsw(limit: int, *, seeds: int = 0, cover_all: bool = False,
                          segment_size: int = 1_000_000) -> List[int]:
    """
    Variant that presieves and then verifies survivors with BPSW (useful if we later add
    partial presieving only). For now this just calls primes_segmented; included for API parity.
    """
    return primes_segmented(limit, seeds=seeds, cover_all=cover_all, segment_size=segment_size)
