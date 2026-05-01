import argparse
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd


def main():
    parser = argparse.ArgumentParser(description='Create an example read-count QC plot')
    parser.add_argument('--in', dest='inp', required=True)
    parser.add_argument('--out', dest='out', required=True)
    args = parser.parse_args()

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(args.inp, sep='\t')
    ax = df.plot(x=df.columns[0], y=df.columns[1], kind='bar')
    ax.set_title('Example read counts')
    plt.tight_layout()
    plt.savefig(out_path)
    print('Saved', out_path)


if __name__ == '__main__':
    main()
