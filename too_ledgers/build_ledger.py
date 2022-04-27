#!/usr/bin/env python

from tooledger import TooLedgerMaker
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter

if __name__ == '__main__':
    p = ArgumentParser(description='ToO Ledger maker',
                       formatter_class=ArgumentDefaultsHelpFormatter)
    p.add_argument('inputfile', nargs=1,
                   help='Input TOO file.')
    p.add_argument('-o', '--output', dest='output', default='ToO-input.ecsv',
                   help='ToO ledger output filename.')
    args = p.parse_args()

    # Generate the ToO ledger.
    ledger = TooLedgerMaker()
    ledger.build_too_ledger(args.inputfile[0], args.output)

