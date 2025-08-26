
from ptp_sieve.segmented import primes_segmented
from ptp_sieve.ptp import _eratosthenes

def test_segmented_complete_matches_truth():
    limit = 20000
    truth = set(_eratosthenes(limit))
    got = set(primes_segmented(limit, cover_all=True, segment_size=5000))
    assert truth == got
