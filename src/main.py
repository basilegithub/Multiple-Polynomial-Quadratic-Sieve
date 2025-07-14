# This si the main script

import argparse
import QS

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--n', required=True, help='number to factor')
    
    args = parser.parse_args()
    
    n = args.n
    
    QS.QS(n)