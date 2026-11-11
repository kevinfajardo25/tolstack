"""Command line: python -m tolstack examples/hinge_gap.csv --lsl 0.10 --usl 0.60"""
import argparse

from .stack import load_csv


def main() -> None:
    ap = argparse.ArgumentParser(prog="tolstack", description="1D tolerance stack-up analysis")
    ap.add_argument("csv", help="stack definition (name, nominal, plus, minus, direction, dist)")
    ap.add_argument("--lsl", type=float, help="lower spec limit for the gap")
    ap.add_argument("--usl", type=float, help="upper spec limit for the gap")
    ap.add_argument("-n", type=int, default=200_000, help="Monte Carlo samples")
    args = ap.parse_args()

    stack = load_csv(args.csv)
    print(f"{len(stack.dims)} dimensions, nominal gap {stack.nominal:.4f}, mean gap {stack.mean:.4f}\n")
    print(stack.worst_case())
    print(stack.rss(args.lsl, args.usl))
    print(stack.monte_carlo(args.n, lsl=args.lsl, usl=args.usl))

    wc = stack.worst_case()
    if args.lsl is not None and wc.low < args.lsl or args.usl is not None and wc.high > args.usl:
        print("\nWorst case violates the spec: the design relies on statistical assumptions.")
    print("\nVariance contribution:")
    for name, pct in stack.contributions():
        print(f"  {name:<24} {pct:5.1f}%  {'#' * round(pct / 2)}")


if __name__ == "__main__":
    main()
