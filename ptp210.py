
#!/usr/bin/env python3
import argparse, json
from ptp_sieve.ptp import primes_via_ptp, recall_vs_truth, ptp_residues, COPRIME_RESIDUES_210_ORDERED, per_residue_coverage
from ptp_sieve.segmented import primes_segmented

def main():
    ap = argparse.ArgumentParser(description="PTP-210 Sieve: prime generator guided by symbolic residues mod 210.")
    ap.add_argument("--limit", type=int, default=100000, help="Upper bound (inclusive) for generated primes.")
    ap.add_argument("--seeds", type=int, default=64, help="Number of residues a(n) to use.")
    ap.add_argument("--json", action="store_true", help="Print primes as JSON array.")
    ap.add_argument("--log-csv", type=str, default=None, help="Write per-residue stats to CSV file.")
    ap.add_argument("--log-json", type=str, default=None, help="Write summary and per-residue stats to JSON file.")
    ap.add_argument("--print-table", action="store_true", help="Print top-20 residues table.")
    ap.add_argument("--sieve", action="store_true", help="Use segmented 210-wheel sieve (fast).")
    ap.add_argument("--cover-all", action="store_true", help="When used with --sieve, cover all 48 residue classes (complete prime list).")
    ap.add_argument("--segment-size", type=int, default=1000000, help="Segment size for the segmented sieve.")

    ap.add_argument("--show-residues", action="store_true", help="Print the residues used (mod 210).")
    ap.add_argument("--diagnostics", action="store_true", help="Print recall vs ground truth.")
    ap.add_argument("--coverage", action="store_true", help="Per-residue coverage table + asymptotic ceiling.")
    args = ap.parse_args()

    if args.show_residues:
        residues = ptp_residues(args.seeds, filter_coprime=True)
        # Deduplicate while maintaining order
        seen = set(); residues_dedup = []
        for r in residues:
            if r not in seen:
                seen.add(r); residues_dedup.append(r)
        print("Residues used (mod 210):", residues_dedup)

    if args.sieve:
        primes = primes_segmented(args.limit, seeds=args.seeds, cover_all=args.cover_all, segment_size=args.segment_size)
    else:
        primes = primes_via_ptp(args.limit, seeds=args.seeds, dedup_residues=True,
                         log_csv=args.log_csv, log_json=args.log_json, print_table=args.print_table)
    if args.json:
        print(json.dumps(primes))
    else:
        print(f"Primes <= {args.limit} (found {len(primes)}):")
        print(primes[:50], "..." if len(primes) > 50 else "")

    if args.diagnostics:
        recall, found, truth = recall_vs_truth(args.limit, seeds=args.seeds)
        print(f"Recall vs ground truth: {recall:.3%}  (found {found} / truth {truth})")
    if args.coverage:
        per_residue_coverage(args.limit, seeds=args.seeds, dedup_residues=True,
                          print_table=True,
                          log_csv=args.log_csv.replace(".csv","_coverage.csv") if args.log_csv else None,
                          log_json=args.log_json.replace(".json","_coverage.json") if args.log_json else None)

if __name__ == "__main__":
    main()
