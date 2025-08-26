
import math
from typing import List, Set

MOD = 210
COPRIME_RESIDUES_210: Set[int] = {r for r in range(MOD) if math.gcd(r, MOD) == 1}
# Canonical ordering (sorted) – useful for iterating deterministically
COPRIME_RESIDUES_210_ORDERED: List[int] = sorted(COPRIME_RESIDUES_210)

def is_coprime_to_210(x: int) -> bool:
    return math.gcd(x, MOD) == 1
