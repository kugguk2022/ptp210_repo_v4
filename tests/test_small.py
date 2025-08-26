
from ptp_sieve.ptp import primes_via_ptp, recall_vs_truth

def test_small_recall():
    recall, found, truth = recall_vs_truth(10000, seeds=64)
    # We don't enforce a particular recall, just ensure we find > 800 primes up to 10k
    # (there are 1229 true primes <= 10k).
    assert found > 800
