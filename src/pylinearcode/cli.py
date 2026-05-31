"""Small command-line front end."""

from __future__ import annotations

import argparse

from .constructors import bch_code, hamming_code, reed_muller_code, reed_solomon_code
from .fields import GF


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="pylinearcode")
    sub = parser.add_subparsers(dest="command", required=True)

    hamming = sub.add_parser("hamming", help="construct a q-ary Hamming code")
    hamming.add_argument("--q", type=int, required=True)
    hamming.add_argument("--r", type=int, required=True)

    rs = sub.add_parser("rs", help="construct a Reed-Solomon code")
    rs.add_argument("--q", type=int, required=True)
    rs.add_argument("--n", type=int, required=True)
    rs.add_argument("--k", type=int, required=True)

    bch = sub.add_parser("bch", help="construct a primitive narrow-sense BCH code")
    bch.add_argument("--q", type=int, required=True)
    bch.add_argument("--m", type=int, required=True, help="extension degree")
    bch.add_argument("--delta", type=int, required=True, help="designed distance")

    rm = sub.add_parser("rm", help="construct a binary Reed-Muller code")
    rm.add_argument("--r", type=int, required=True, help="order")
    rm.add_argument("--m", type=int, required=True, help="number of variables")

    args = parser.parse_args(argv)
    if args.command == "hamming":
        field = GF(args.q)
        code = hamming_code(field, args.r)
    elif args.command == "rs":
        field = GF(args.q)
        code = reed_solomon_code(field, length=args.n, dimension=args.k)
    elif args.command == "bch":
        field = GF(args.q)
        code = bch_code(field, extension_degree=args.m, designed_distance=args.delta)
    else:
        code = reed_muller_code(args.r, args.m)
    print(code)
    print(f"parameters: [{code.length}, {code.dimension}, {code.minimum_distance()}]")
    print("generator matrix:")
    for row in code.generator_matrix:
        print(row)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
