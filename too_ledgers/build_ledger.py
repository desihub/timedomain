from tooparse import TooLedgerMaker
from argparse import ArgumentParser


if __name__ == '__main__':
    p = ArgumentParser(description='ToO Ledger maker')
    p.add_argument('inputfile', nargs=1, help='Input TOO file.')
    args = p.parse_args()

    # Generate the ToO ledger.
    ledger = TooLedgerMaker()
    ledger.build_too_ledger(args.inputfile[0], 'test.ecsv')

