"""Small command-line front end."""

from __future__ import annotations

import argparse

from .constructors import hamming_code, reed_solomon_code
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

    args = parser.parse_args(argv)
    field = GF(args.q)
    if args.command == "hamming":
        code = hamming_code(field, args.r)
    else:
        code = reed_solomon_code(field, length=args.n, dimension=args.k)
    print(code)
    print(f"parameters: [{code.length}, {code.dimension}, {code.minimum_distance()}]")
    print("generator matrix:")
    for row in code.generator_matrix:
        print(row)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

