#!/usr/bin/env python3
import argparse
import re


TRANSCRIPT_RE = re.compile(r'transcript_id "([^"]+)"')


def load_ids(path):
    ids = set()
    with open(path) as handle:
        for line in handle:
            value = line.strip()
            if value:
                ids.add(value)
    return ids


def main():
    parser = argparse.ArgumentParser(
        description="Filter a GTF/GFF-like file to records whose transcript_id is in an ID list."
    )
    parser.add_argument("--ids", required=True, help="Text file with one transcript ID per line.")
    parser.add_argument("--gtf", required=True, help="Input candidate GTF.")
    parser.add_argument("--output", required=True, help="Filtered output GTF.")
    args = parser.parse_args()

    keep_ids = load_ids(args.ids)
    with open(args.gtf) as infile, open(args.output, "w") as outfile:
        for line in infile:
            if line.startswith("#"):
                outfile.write(line)
                continue
            match = TRANSCRIPT_RE.search(line)
            if match and match.group(1) in keep_ids:
                outfile.write(line)


if __name__ == "__main__":
    main()
