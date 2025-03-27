import argparse
import pandas as pd
import sys

def main():
    usage_examples = (
        "Usage examples:\n"
        "  python detect.py <*.out>\n"
        "  optional case:\n"
        "  python detect.py <*.out> <delta-chi-squred(default:160)> <obsgroup_number(default:0)> \n"
    )

    if len(sys.argv) == 1:
        print("No arguments provided.\n")
        print(usage_examples)
        sys.exit(1)

    parser = argparse.ArgumentParser(
        description="Filter events by chi-square threshold.",
        epilog=usage_examples,
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    parser.add_argument("outfile", help=".out file path")
    parser.add_argument("chi2_threshold", nargs="?", type=float, default=160.0)
    parser.add_argument("obsgroup", type=int, default=0, help="ObsGroup index (default: 0)")
    args = parser.parse_args()

    df = pd.read_csv(args.outfile, delim_whitespace=True)

    col = f"ObsGroup_{args.obsgroup}_chi2"
    if not {"EventID", col}.issubset(df.columns):
        raise KeyError(f"Missing required column(s): 'EventID' and/or '{col}'")

    filtered = df[df[col] > args.chi2_threshold][["EventID", col]]

    print(f"Detections with {col} > {args.chi2_threshold}:")
    for eid, chi2 in zip(filtered["EventID"], filtered[col]):
        print(f"{int(eid)} |  {chi2}")

if __name__ == "__main__":
    main()
