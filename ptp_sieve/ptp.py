
"""
PTP-210 sieve: prime generation guided by a symbolic residue generator mod 210.
Implements the formula from the prompt image:

    a(n) = floor( ((zeta(3) + phi + kappa3(n)) ** n) / log(n + pi) ) (mod 210)

where
  - zeta(3) ~ 1.2020569  (Apéry's constant)
  - phi = (1 + sqrt(5)) / 2
  - kappa3(n) = arg(n + i * sin(n ** phi))

We use these residues as seeds r (coprime to 210), then scan arithmetic
progressions r + 210k up to a user-specified limit and test primality.
"""
from __future__ import annotations

import math
from typing import Iterable, Iterator, List, Set, Tuple, Dict

from .wheel import MOD, COPRIME_RESIDUES_210, COPRIME_RESIDUES_210_ORDERED

ZETA3 = 1.2020569
PHI = (1 + 5 ** 0.5) / 2.0

# ------------------------ symbolic sequence ------------------------

def kappa3(n: int) -> float:
    """kappa^3(n) = arg(n + i * sin(n ** phi))."""
    # arg(x + i y) = atan2(y, x)
    return math.atan2(math.sin(n ** PHI), float(n))

def a_n(n: int) -> int:
    """Compute a(n) as defined, then reduce mod 210."""
    x = (ZETA3 + PHI + kappa3(n)) ** n
    denom = math.log(n + math.pi)
    if denom == 0.0:
        denom = 1e-12
    val = int(x / denom)
    return val % MOD

def ptp_residues(seed_count: int, *, filter_coprime: bool = True) -> List[int]:
    """
    Produce the first `seed_count` residues a(n) mod 210.
    If `filter_coprime` is True, keep only residues coprime to 210
    (the only classes that can contain primes > 7). Duplicates are kept
    to preserve the original distribution, but you can de-dup downstream.
    """
    out: List[int] = []
    n = 2  # start from 2
    while len(out) < seed_count:
        r = a_n(n)
        if (not filter_coprime) or (r in COPRIME_RESIDUES_210):
            out.append(r)
        n += 1
    return out

# ------------------------ primality testing ------------------------

def _small_prime_trial(n: int) -> bool:
    """Quick trial division by small primes. Returns False if composite, True if inconclusive."""
    if n < 2:
        return False
    small_primes = [2,3,5,7,11,13,17,19,23,29,31,37]
    for p in small_primes:
        if n == p:
            return True
        if n % p == 0:
            return False
    return True

def _miller_rabin(n: int, bases: Iterable[int]) -> bool:
    """Deterministic Miller–Rabin for 64-bit using a conservative base set."""
    if n < 2:
        return False
    # even?
    if n % 2 == 0:
        return n == 2
    # write n-1 = d * 2^s with d odd
    d = n - 1
    s = 0
    while d % 2 == 0:
        d //= 2
        s += 1
    for a in bases:
        if a % n == 0:
            continue
        x = pow(a, d, n)
        if x == 1 or x == n - 1:
            continue
        for _ in range(s - 1):
            x = pow(x, 2, n)
            if x == n - 1:
                break
        else:
            return False
    return True

# A conservative base set that is deterministic for all 64-bit integers
# (superset of minimal known sets; fine for performance here).
_MR_BASES_64 = (2,3,5,7,11,13,17,19,23,29,31,37)

def is_prime(n: int) -> bool:
    if not _small_prime_trial(n):
        return False
    # for small n, after small trial we can finish quickly
    if n < 10_000_000:
        # trial divide by odd numbers up to sqrt
        r = int(n**0.5)
        d = 41
        while d <= r:
            if n % d == 0:
                return False
            d += 2
        return True
    # Otherwise use Miller-Rabin with conservative base set
    return _miller_rabin(n, _MR_BASES_64)


# ------------------------ main generation ------------------------

def primes_via_ptp(
    limit: int,
    *,
    seeds: int = 64,
    dedup_residues: bool = True,
    log_csv: str | None = None,
    log_json: str | None = None,
    print_table: bool = False
) -> List[int]:
    """
    Generate primes <= limit using arithmetic progressions guided by the PTP residues.
    Parameters
    ----------
    limit : int
        Upper bound for primes (inclusive).
    seeds : int
        Number of residues a(n) to use.
    dedup_residues : bool
        Use each residue at most once (first occurrence).
    log_csv : str | None
        If provided, write per-residue statistics to this CSV file.
    log_json : str | None
        If provided, write full summary (including top residues) to this JSON file.
    print_table : bool
        If True, print a compact per-residue summary table to stdout.
    """
    primes = []
    # Include small primes
    for p in (2,3,5,7):
        if p <= limit:
            primes.append(p)

    residues = ptp_residues(seeds, filter_coprime=True)
    if dedup_residues:
        seen = set()
        deduped = []
        for r in residues:
            if r not in seen:
                seen.add(r)
                deduped.append(r)
        residues = deduped

    stats = []  # per-residue stats

    for r in residues:
        if math.gcd(r, MOD) != 1:
            continue
        if r > limit:
            candidates = 0
            k_max = -1
        else:
            k_max = (limit - r) // MOD
            candidates = k_max + 1

        prime_hits = 0
        first_hit = None
        last_hit = None

        for k in range(k_max + 1):
            m = r + MOD * k
            if is_prime(m):
                prime_hits += 1
                last_hit = m
                if first_hit is None:
                    first_hit = m
                primes.append(m)

        density = (prime_hits / candidates) if candidates > 0 else 0.0
        stats.append(
            {
                "residue": r,
                "candidates": candidates,
                "primes_found": prime_hits,
                "density": density,
                "first_prime": first_hit,
                "last_prime": last_hit,
            }
        )

    primes = sorted(set(primes))

    # Optional logging outputs
    if log_csv:
        os.makedirs(os.path.dirname(log_csv) or ".", exist_ok=True)
        with open(log_csv, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=["residue","candidates","primes_found","density","first_prime","last_prime"])
            w.writeheader()
            for row in sorted(stats, key=lambda x: (-x["primes_found"], -x["density"], x["residue"])):
                w.writerow(row)

    if log_json:
        os.makedirs(os.path.dirname(log_json) or ".", exist_ok=True)
        summary = {
            "limit": limit,
            "seeds": seeds,
            "dedup_residues": dedup_residues,
            "total_primes": len(primes),
            "top_residues": sorted(stats, key=lambda x: (-x["primes_found"], -x["density"]))[:10],
            "all_stats": stats,
        }
        with open(log_json, "w") as f:
            json.dump(summary, f, indent=2)

    if print_table:
        # Pretty console table
        header = f"{'res':>4}  {'cand':>6}  {'hits':>5}  {'dens':>6}  {'first':>8}  {'last':>8}"
        print(header)
        print("-"*len(header))
        for s in sorted(stats, key=lambda x: (-x["primes_found"], -x["density"], x["residue"]))[:20]:
            print(f"{s['residue']:>4}  {s['candidates']:>6}  {s['primes_found']:>5}  {s['density']:>6.3f}  {str(s['first_prime'] or ''):>8}  {str(s['last_prime'] or ''):>8}")

    return primes


# ------------------------ diagnostics ------------------------


def recall_vs_truth(limit: int, *, seeds: int = 64) -> Tuple[float, int, int]:
    """
    Compare primes found by PTP-guided search against ground truth primes up to `limit`
    (computed by a classic sieve). Returns (recall, found, truth_count).
    """
    truth = _eratosthenes(limit)
    found = set(primes_via_ptp(limit, seeds=seeds, dedup_residues=True))
    # exclude 1
    truth_primes = {p for p in truth if p >= 2}
    recall = len(found & truth_primes) / len(truth_primes) if truth_primes else 1.0
    return (recall, len(found), len(truth_primes))



# ------------------------ per-class coverage report ------------------------

def per_residue_coverage(limit: int, *, seeds: int = 64, dedup_residues: bool = True,
                         print_table: bool = True, log_csv: str | None = None, log_json: str | None = None) -> Dict:
    """
    Compute how many primes (in %) we capture per residue class mod 210 and in total, and
    estimate the asymptotic ceiling based on how many coprime classes we scan.
    Returns a dict with summary + per-residue rows.
    """
    # Ground truth primes and their per-residue counts (exclude <=7)
    truth = _eratosthenes(limit)
    truth_primes = [p for p in truth if p > 7]
    truth_total = len(truth_primes)
    truth_counts: Dict[int, int] = {}
    for p in truth_primes:
        r = p % MOD
        if math.gcd(r, MOD) == 1:
            truth_counts[r] = truth_counts.get(r, 0) + 1

    # Generate with PTP and also gather per-residue stats from primes_via_ptp's internal logic
    # We'll reuse some parts: compute the residue set we actually scan
    residues = ptp_residues(seeds, filter_coprime=True)
    if dedup_residues:
        seen = set(); residues_d = []
        for r in residues:
            if r not in seen:
                seen.add(r); residues_d.append(r)
        residues = residues_d
    scanned_set = {r for r in residues if math.gcd(r, MOD) == 1}

    # Find primes via our method
    found_primes = set(primes_via_ptp(limit, seeds=seeds, dedup_residues=dedup_residues))
    found_primes = {p for p in found_primes if p > 7}
    found_total = len(found_primes)

    found_counts: Dict[int, int] = {}
    for p in found_primes:
        r = p % MOD
        if math.gcd(r, MOD) == 1:
            found_counts[r] = found_counts.get(r, 0) + 1

    rows = []
    for r in sorted(COPRIME_RESIDUES_210_ORDERED):
        truth_r = truth_counts.get(r, 0)
        found_r = found_counts.get(r, 0)
        recall_r = (found_r / truth_r * 100.0) if truth_r else 0.0
        share_total = (found_r / truth_total * 100.0) if truth_total else 0.0
        scanned = r in scanned_set
        rows.append({
            "residue": r,
            "truth_primes": truth_r,
            "found_primes": found_r,
            "recall_pct": recall_r,
            "share_of_total_pct": share_total,
            "scanned": scanned
        })

    overall_recall_pct = (found_total / truth_total * 100.0) if truth_total else 100.0
    asymptotic_ceiling_pct = (len(scanned_set) / len(COPRIME_RESIDUES_210) * 100.0) if COPRIME_RESIDUES_210 else 100.0

    summary = {
        "limit": limit,
        "seeds": seeds,
        "distinct_residues_scanned": len(scanned_set),
        "phi_210": len(COPRIME_RESIDUES_210),
        "overall_recall_pct": overall_recall_pct,
        "asymptotic_ceiling_pct": asymptotic_ceiling_pct,
        "note": "Asymptotic ceiling assumes primes equidistribute among coprime classes (Dirichlet)."
    }

    if print_table:
        print(f"PTP-210 Coverage up to {limit} (seeds={seeds})")
        print(f"Overall recall: {overall_recall_pct:.2f}% | Asymptotic ceiling (with {len(scanned_set)}/{len(COPRIME_RESIDUES_210)} classes): {asymptotic_ceiling_pct:.2f}%")
        header = f"{'res':>4}  {'truth':>6}  {'found':>5}  {'recall%':>7}  {'share%':>7}  {'scan':>4}"
        print(header); print("-"*len(header))
        for row in rows:
            mark = "✓" if row["scanned"] else ""
            print(f"{row['residue']:>4}  {row['truth_primes']:>6}  {row['found_primes']:>5}  {row['recall_pct']:>7.2f}  {row['share_of_total_pct']:>7.2f}  {mark:>4}")

    if log_csv:
        import csv, os
        os.makedirs(os.path.dirname(log_csv) or ".", exist_ok=True)
        with open(log_csv, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=["residue","truth_primes","found_primes","recall_pct","share_of_total_pct","scanned"])
            w.writeheader()
            for row in rows:
                w.writerow(row)

    if log_json:
        import json, os
        os.makedirs(os.path.dirname(log_json) or ".", exist_ok=True)
        with open(log_json, "w") as f:
            json.dump({"summary": summary, "rows": rows}, f, indent=2)

    return {"summary": summary, "rows": rows}
def _eratosthenes(n: int) -> List[int]:
    """Ground-truth sieve for diagnostics."""
    if n < 2:
        return []
    sieve = bytearray(b"\x01") * (n + 1)
    sieve[0:2] = b"\x00\x00"
    p = 2
    while p * p <= n:
        if sieve[p]:
            step = p
            start = p * p
            sieve[start:n+1:step] = b"\x00" * (((n - start) // step) + 1)
        p += 1
    return [i for i in range(n + 1) if sieve[i]]
