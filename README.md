
# PTP‑210 Prime Generator

> _Primes via symbolic residue classes modulo 210 (2·3·5·7) guided by the PTP a(n) sequence._

[![CI](https://img.shields.io/github/actions/workflow/status/USER/REPO/ci.yml?branch=main)](https://github.com/USER/REPO/actions)

This project turns the PTP idea into a practical prime generator that
**only scans the 210‑wheel coprime residue classes suggested by a(n)** and tests numbers of the form
$$\(a(n) + 210k\)$$. It includes a CLI, library API, **per‑residue logging**, and a small test.

## 📐 The sequence

$$\[ a(n) = \left\lfloor \frac{(\zeta(3)+\varphi+\kappa^3(n))^n}{\log(n+\pi)} \right\rfloor \bmod 210 \]$$

- $$\(\zeta(3) \approx 1.2020569\)$$ (Apéry's constant)  
- $$\(\varphi = \tfrac{1+\sqrt 5}2\)$$  
- $$\(\kappa^3(n) = \arg(n + i\cdot \sin(n^\varphi))\)$$  
- Reduce mod 210 (product of first four primes). We **filter** to the 48 classes coprime to 210.

## 🚀 Quick start

```bash
python -m venv .venv && . .venv/bin/activate
pip install -e .

# run the CLI
python ptp210.py --limit 200000 --seeds 96 --show-residues --diagnostics --print-table                  --log-csv logs/ptp_200k.csv --log-json logs/ptp_200k.json
```

**Example output (table):**
```
 res    cand   hits   dens     first      last
----------------------------------------------
  11     952     16  0.017       431    199931
 119     952     15  0.016      1199    197579
 ...
```

- `cand` = number of tested candidates in that residue class (≤ limit)  
- `hits` = primes found in that class  
- `dens` = hits / cand  
- `first`, `last` = first/last prime hit in that class

## 🔧 CLI options

- `--limit` upper bound for primes (inclusive)  
- `--seeds` how many `a(n)` residues to use (we de‑dup by default)  
- `--show-residues` print the residue list actually used  
- `--diagnostics` compare against a ground‑truth sieve (recall)  
- `--print-table` print a top‑20 residue table with densities  
- `--log-csv PATH` write per‑residue stats to CSV  
- `--log-json PATH` write full summary + stats to JSON  
- `--json` output primes as a JSON array

## 📚 Library usage

```python
from ptp_sieve.ptp import primes_via_ptp

primes = primes_via_ptp(
    limit=1_000_000, seeds=128,
    log_csv="logs/ptp_1M.csv", log_json="logs/ptp_1M.json", print_table=True
)
```

## 🧪 CI (GitHub Actions)

A simple CI runs tests on Python 3.10–3.12. After pushing to GitHub, replace the badge `USER/REPO` with your repo path.

```
.github/workflows/ci.yml
```

## ⚠️ Scope & caveats

- This is **exploratory**: the PTP sequence is a heuristic prioritization across the 210‑wheel classes.
- Deterministic Miller–Rabin is used for large numbers; for small limits we finish by trial division.
- For massive ranges you’ll want segmentation and/or parallelization.

## 🗺️ Roadmap

- Adaptive seeding (multi‑armed bandits over residue performance)  
- Segmented scanning over [A, B] windows (distributed searches)  
- Multiprocessing for residue partitions  
- Precision upgrades for constants and angle function variants

## License

MIT


## 📊 Coverage & asymptotics

Use `--coverage` to see **per-residue capture %** and overall coverage, plus an **asymptotic ceiling** estimate.

```bash
python ptp210.py --limit 300000 --seeds 96 --coverage \
  --log-csv logs/ptp_300k.csv --log-json logs/ptp_300k.json
```

- **Overall recall** = % of all true primes ≤ limit that we captured.  
- **Per-class recall** = % of true primes in that residue class we captured.  
- **Asymptotic ceiling** ≈ `(# distinct scanned coprime classes) / 48` — assuming primes equidistribute across the 48 coprime classes mod 210 (Dirichlet). This is the **max %** we can ever capture with our current class set, regardless of limit.


## ⚡ Segmented 210-wheel sieve (fast + complete)

You can now run a **proper segmented sieve** that enumerates primes up to `--limit` efficiently.
Use PTP residues to **prioritize** classes, or add `--cover-all` to get a **complete** prime list.

```bash
# fast, PTP-driven subset of classes
python ptp210.py --sieve --limit 2000000 --seeds 96 --segment-size 1000000

# full correctness (all primes), still fast
python ptp210.py --sieve --cover-all --limit 2000000 --segment-size 1000000
```

- Precomputes base primes up to √N once.
- Segments `[low, high)` to keep memory flat.
- When `--cover-all` is on, this is a **true sieve** (no MR/BPSW needed) and matches a classic sieve.
- Without `--cover-all`, it returns only primes in the chosen residue classes (good for targeted hunts).

