
from ptp_sieve.ptp import primes_via_ptp, ptp_residues, recall_vs_truth

LIMIT = 20000
SEEDS = 64

print("PTP residues (dedup, first 20):", list(dict.fromkeys(ptp_residues(SEEDS)))[:20])

primes = primes_via_ptp(
    LIMIT, seeds=SEEDS,
    log_csv="logs/demo_20k.csv",
    log_json="logs/demo_20k.json",
    print_table=True
)
print(f"Found {len(primes)} primes <= {LIMIT}. First 30:", primes[:30])

recall, found, truth = recall_vs_truth(LIMIT, seeds=SEEDS)
print(f"Recall vs ground truth up to {LIMIT}: {recall:.2%}  (found {found} / truth {truth})")
