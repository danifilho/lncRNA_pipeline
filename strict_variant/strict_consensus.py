#!/usr/bin/env python3
"""Strict variant of the lncRNA consensus (EXPERIMENT, does not touch the faithful pipeline).

Two changes vs the upstream verbatim intersection:
  1. CPAT is parsed CORRECTLY -> filter on the real coding_prob (< cutoff), instead of the
     upstream read_delim column-misalignment that makes CPAT match ~everything.
  2. Isoforms are collapsed to one (longest) transcript per StringTie gene (MSTRG.x).

Reads the per-tool outputs of a finished run; writes ID lists + a counts summary.
"""
import argparse, os

def feelnc_ids(p):
    return {l.strip() for l in open(p) if l.strip()}

def cpat_noncoding(p, cutoff):
    out = set()
    with open(p) as f:
        f.readline()                      # header (5 names; rows have ID + 5 values)
        for line in f:
            c = line.rstrip("\n").split("\t")
            if len(c) < 6:
                continue
            try:
                if float(c[-1]) < cutoff:  # real coding_prob is the LAST field
                    out.add(c[0])          # ID is the unnamed first field
            except ValueError:
                pass
    return out

def lncfinder_noncoding(p):
    out = set()
    with open(p) as f:
        f.readline()
        for line in f:
            c = line.rstrip("\n").split("\t")
            if len(c) >= 2 and c[1] == "NonCoding":   # col0=ID, col1=label (row.names shift)
                out.add(c[0])
    return out

def boost_lncrna(p):
    out = set()
    with open(p) as f:
        f.readline()
        for line in f:
            c = line.rstrip("\n").split("\t")
            if len(c) >= 2 and c[1].strip() == "1":
                out.add(c[0])
    return out

def protein_hits(p, evalue, pident):
    out = set()
    with open(p) as f:
        for line in f:
            c = line.rstrip("\n").split("\t")
            if len(c) < 11:
                continue
            try:
                if float(c[2]) > pident and float(c[10]) < evalue:
                    out.add(c[0])
            except ValueError:
                pass
    return out

def fasta_lengths(p):
    L, cur = {}, None
    with open(p) as f:
        for line in f:
            if line.startswith(">"):
                cur = line[1:].split()[0]
                L[cur] = 0
            elif cur is not None:
                L[cur] += len(line.strip())
    return L

ap = argparse.ArgumentParser()
ap.add_argument("--feelnc", required=True)
ap.add_argument("--cpat", required=True)
ap.add_argument("--lncfinder", required=True)
ap.add_argument("--boost", required=True)
ap.add_argument("--diamond", required=True)
ap.add_argument("--fasta", required=True)      # candidate_transcript.fasta (for lengths)
ap.add_argument("--cpat-cutoff", type=float, default=0.46)
ap.add_argument("--diamond-evalue", type=float, default=1e-5)
ap.add_argument("--diamond-pident", type=float, default=60.0)
ap.add_argument("--outdir", required=True)
a = ap.parse_args()
os.makedirs(a.outdir, exist_ok=True)

fe = feelnc_ids(a.feelnc)
cp = cpat_noncoding(a.cpat, a.cpat_cutoff)
lf = lncfinder_noncoding(a.lncfinder)
bo = boost_lncrna(a.boost)
pr = protein_hits(a.diamond, a.diamond_evalue, a.diamond_pident)

consensus = (fe & cp & lf & bo) - pr

lengths = fasta_lengths(a.fasta)
genes = {}
for tid in consensus:
    g = tid.rsplit(".", 1)[0]            # MSTRG.x.y -> MSTRG.x
    L = lengths.get(tid, 0)
    if g not in genes or L > genes[g][1]:
        genes[g] = (tid, L)
gene_reprs = sorted(v[0] for v in genes.values())   # longest isoform per gene

with open(os.path.join(a.outdir, "strict_final_transcript_ids.txt"), "w") as fh:
    fh.write("\n".join(sorted(consensus)) + "\n")
with open(os.path.join(a.outdir, "strict_final_gene_ids.txt"), "w") as fh:
    fh.write("\n".join(gene_reprs) + "\n")

with open(os.path.join(a.outdir, "per_tool_counts.txt"), "w") as fh:
    fh.write(f"FEELnc_candidates={len(fe)}\n")
    fh.write(f"CPAT_noncoding_REAL(coding_prob<{a.cpat_cutoff})={len(cp)}\n")
    fh.write(f"LncFinder_noncoding={len(lf)}\n")
    fh.write(f"PlantLncBoost_lncRNA={len(bo)}\n")
    fh.write(f"DIAMOND_protein_hits(pident>{a.diamond_pident},e<{a.diamond_evalue})={len(pr)}\n")
    fh.write(f"STRICT_consensus_transcripts={len(consensus)}\n")
    fh.write(f"STRICT_consensus_genes(longest_isoform)={len(gene_reprs)}\n")

print(f"strict transcripts={len(consensus)}  strict genes={len(gene_reprs)}")
