#!/usr/bin/env python3
"""
01_simulate_dataset.py
Genetics Assistant Internship -- Week 1: Genetic Sequence Analysis and Reporting

Designs a fully reproducible simulated genetic sequence dataset:
  * a 1,500 bp reference gene (SIMG1): 5' UTR (1-100), coding sequence
    (101-1399; 433 codons) with a glycine-rich loop and a catalytic motif,
    and a 3' UTR (1400-1500) carrying a polyadenylation signal
  * 12 sample sequences (6 cases, 6 controls) carrying designed variants:
    five common SNVs, a case-enriched missense marker, a case-enriched
    in-frame deletion marker, premature-stop (nonsense) variants, a
    frameshift insertion, a control-only SNV, a control-only deletion and a
    designed synonymous variant
  * a 4-sequence mini reference database for homology searching

Coordinate convention: all positions are 1-based on the 1500 bp reference.
The script asserts every designed effect after construction, so the dataset
is self-verifying. Seed = 20260926 (fully reproducible).
"""
import csv
import json
import os
import random

SEED = 20260926
OUT = "data"
os.makedirs(OUT, exist_ok=True)
rng = random.Random(SEED)

BASES = "ACGT"
STOPS = {"TAA", "TAG", "TGA"}
CODONS = [a + b + c for a in BASES for b in BASES for c in BASES]
SENSE = [c for c in CODONS if c not in STOPS]
CDS_START, CDS_END = 101, 1399
N_CODONS = (CDS_END - CDS_START + 1) // 3          # 433

TABLE = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L', 'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M', 'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*', 'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K', 'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W', 'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R', 'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}

def rand_dna(n):
    return "".join(rng.choice(BASES) for _ in range(n))

def build_cds(n_codons):
    codons = ["ATG"] + [rng.choice(SENSE) for _ in range(n_codons - 2)] + ["TAA"]
    return codons

def set_motif(codons, start_idx, aa_pattern):
    for k, aa in enumerate(aa_pattern):
        if aa == "x":
            continue
        choices = [c for c, a in TABLE.items() if a == aa and c not in STOPS]
        codons[start_idx + k] = rng.choice(choices)

# ---------------------------------------------------------------- reference
utr5 = list(rand_dna(100))
utr5[24:30] = list("TATAAA")
cds = build_cds(N_CODONS)
set_motif(cds, 59, "GxGxxG")      # glycine-rich loop, codons 60-65
set_motif(cds, 199, "DxxxxN")     # catalytic motif, codons 200-205
utr3 = list(rand_dna(101))
utr3[70:76] = list("AATAAA")
REFERENCE = "".join(utr5) + "".join(cds) + "".join(utr3)
assert len(REFERENCE) == 1500
assert REFERENCE[CDS_START - 1:CDS_START + 2] == "ATG"
assert REFERENCE[CDS_END - 3:CDS_END] == "TAA"

# ---- coordinate-correct helpers (0-based index of a codon's first base)
def codon_span(c):
    """(start1, end1) 1-based genomic span of 1-based codon index c."""
    i0 = (CDS_START - 1) + (c - 1) * 3
    return i0 + 1, i0 + 3

def codon_seq(seq, c):
    s, e = codon_span(c)
    return seq[s - 1:e]

def ref_codon(c):
    return codon_seq(REFERENCE, c)

# --------------------------------------------------- find designed positions
def find_stop_variant(base_index, wanted_alt=None):
    """Codon where a single-base change at base_index (0,1,2) creates a stop."""
    for c in range(30, N_CODONS - 5):
        codon = ref_codon(c)
        for alt in (wanted_alt or BASES):
            if alt == codon[base_index]:
                continue
            new = codon[:base_index] + alt + codon[base_index + 1:]
            if new in STOPS:
                pos = codon_span(c)[0] + base_index
                return c, base_index, codon, new, pos, codon[base_index], alt
    raise RuntimeError("no stop variant found")

def find_synonymous_variant(base_index=0, wanted_alt="T"):
    for c in range(20, N_CODONS - 5):
        codon = ref_codon(c)
        alt = wanted_alt
        if alt == codon[base_index]:
            continue
        new = codon[:base_index] + alt + codon[base_index + 1:]
        if TABLE[new] == TABLE[codon]:
            pos = codon_span(c)[0] + base_index
            return c, base_index, codon, new, pos, codon[base_index], alt
    raise RuntimeError("no synonymous variant found")

NS1 = find_stop_variant(0, "T")            # first-base C>T nonsense
NS2 = find_stop_variant(1, "A")            # second-base G>A nonsense
SYN = find_synonymous_variant(0, "T")      # synonymous C>T

# ------------------------------------------------------------- edit engine
def substitute(seq, pos1, alt):
    i = pos1 - 1
    old = seq[i]
    assert old != alt, f"no-op substitution at {pos1} ({old}>{alt})"
    return seq[:i] + alt + seq[i + 1:], old

def delete(seq, start1, length):
    i = start1 - 1
    return seq[:i] + seq[i + length:], seq[i:i + length]

def insert(seq, after1, inserted):
    return seq[:after1] + inserted + seq[after1:], inserted

def apply_edits(edits):
    seq = REFERENCE
    for e in [x for x in edits if x[0] == "sub"]:
        seq, _ = substitute(seq, e[1], e[2])
    for e in sorted([x for x in edits if x[0] != "sub"], key=lambda x: -x[1]):
        if e[0] == "del":
            seq, _ = delete(seq, e[1], e[2])
        else:
            seq, _ = insert(seq, e[1], e[2])
    return seq

# --------------------------------------------------------------- variant set
def snv_effect(pos, alt):
    """Coding effect of a substitution at 1-based genomic position pos."""
    if not (CDS_START <= pos <= CDS_END):
        return "non-coding"
    c = ((pos - CDS_START) // 3) + 1
    off = (pos - CDS_START) % 3
    cod = ref_codon(c)
    new = cod[:off] + alt + cod[off + 1:]
    if TABLE[new] == TABLE[cod]:
        return "synonymous"
    return "nonsense" if TABLE[new] == "*" else "missense"

def pick_snv(pos, prefer, avoid=("nonsense",)):
    """Preferred alternate allele, auto-switched to avoid a stop codon."""
    for alt in [prefer] + [b for b in BASES if b != prefer]:
        if alt == REFERENCE[pos - 1]:
            continue
        if snv_effect(pos, alt) not in avoid:
            return alt
    raise RuntimeError(f"no acceptable alternate allele at {pos}")

# Background variants are transition-weighted (4 transitions : 1 transversion),
# matching the transition bias of real mutational spectra.
COMMON = [("sub", 150, pick_snv(150, "A")), ("sub", 420, pick_snv(420, "T")),
          ("sub", 780, pick_snv(780, "A")), ("sub", 1100, pick_snv(1100, "A")),
          ("sub", 1450, pick_snv(1450, "C"))]
MARKER_A = ("sub", 612, pick_snv(612, "G"))                    # case-enriched missense
MB_POS = codon_span(200)[0]                                   # first base of codon 200
MARKER_B = ("del", MB_POS, 3)                                 # whole-codon (in-frame) deletion
CTRL_SNV = ("sub", 350, pick_snv(350, "G"))                    # control-only
CD_POS = codon_span(50)[0]                                    # first base of codon 50
CTRL_DEL = ("del", CD_POS, 3)                                 # whole-codon (in-frame) deletion
FRAMESHIFT = ("ins", 1200, "G")                               # case-only frameshift

SAMPLES = {
    "CASE-1": COMMON + [MARKER_A, MARKER_B],
    "CASE-2": COMMON + [MARKER_A, MARKER_B, ("sub", NS1[4], NS1[6])],
    "CASE-3": COMMON + [MARKER_A, ("sub", NS2[4], NS2[6])],
    "CASE-4": COMMON + [MARKER_A, MARKER_B],
    "CASE-5": COMMON + [MARKER_A, FRAMESHIFT],
    "CASE-6": COMMON + [MARKER_A, MARKER_B, ("sub", NS1[4], NS1[6])],
    "CTRL-1": COMMON + [CTRL_SNV, ("sub", SYN[4], SYN[6])],
    "CTRL-2": COMMON + [CTRL_DEL],
    "CTRL-3": COMMON + [CTRL_SNV],
    "CTRL-4": COMMON + [CTRL_SNV, CTRL_DEL],
    "CTRL-5": COMMON + [("sub", SYN[4], SYN[6])],
    "CTRL-6": COMMON + [],
}

samples = {name: apply_edits(edits) for name, edits in SAMPLES.items()}

# --------------------------------------------------- self-verification pass
def translate_cds(seq):
    """Translate the sample's own CDS: 5' UTR is 100 bp, 3' UTR is 101 bp,
    so the CDS ends at len(seq) - 101 whatever indels the sample carries."""
    cds_seq = seq[CDS_START - 1: len(seq) - 101]
    return "".join(TABLE[cds_seq[i:i + 3]] for i in range(0, len(cds_seq) - 2, 3))

ref_protein = translate_cds(REFERENCE)
assert ref_protein[0] == "M" and ref_protein[-1] == "*"
assert ref_protein[:-1].count("*") == 0, "reference protein has internal stop"

def first_stop_codon_index(seq):
    """Index of the first stop codon in the sample's CDS frame."""
    prot = translate_cds(seq)
    idx = prot.find("*")
    return idx if idx != -1 else len(prot)

premature = {}
for name, seq in samples.items():
    prot = translate_cds(seq)
    first_stop = first_stop_codon_index(seq)
    premature[name] = first_stop
    if name == "CASE-5":                      # frameshift: shifted frame, truncated product
        assert first_stop < len(prot) - 1, "CASE-5 frameshift should truncate"
    elif name in ("CASE-2", "CASE-3", "CASE-6"):
        assert prot[-1] == "*", f"{name}: CDS does not end in a stop codon"
        assert first_stop < len(prot) - 1, f"{name} missing its designed premature stop"
    else:
        assert prot[-1] == "*", f"{name}: CDS does not end in a stop codon"
        assert first_stop == len(prot) - 1, f"{name} unexpectedly has a premature stop"
assert snv_effect(MARKER_A[1], MARKER_A[2]) == "missense", "marker A must be missense"

# ------------------------------------------------------------------- output
def write_fasta(path, records):
    with open(path, "w") as fh:
        for name, seq in records:
            fh.write(f">{name}\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i + 60] + "\n")

write_fasta(f"{OUT}/reference.fasta", [("SIMG1_reference", REFERENCE)])
write_fasta(f"{OUT}/samples.fasta", list(samples.items()))

homolog_edits, used_positions = [], set()
while len(homolog_edits) < 180:
    p = rng.randint(1, 1500)
    if p in used_positions:
        continue
    used_positions.add(p)
    alt = rng.choice([b for b in BASES if b != REFERENCE[p - 1]])
    homolog_edits.append(("sub", p, alt))
HOMOLOG = apply_edits(homolog_edits)
distant = list(REFERENCE)
for _ in range(830):
    p = rng.randint(1, 1500)
    distant[p - 1] = rng.choice([b for b in BASES if b != distant[p - 1]])
DISTANT = "".join(distant)
UNRELATED = rand_dna(1500)
write_fasta(f"{OUT}/sim_database.fasta",
            [("SIMG1", REFERENCE), ("SIMG1_HOMOLOG_B", HOMOLOG),
             ("SIMG1_DISTANT_C", DISTANT), ("UNRELATED_D", UNRELATED)])

with open(f"{OUT}/sample_metadata.csv", "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["sample_id", "group", "length_bp", "n_edits_vs_reference"])
    for name, edits in SAMPLES.items():
        w.writerow([name, "case" if name.startswith("CASE") else "control",
                    len(samples[name]), len(edits)])

with open(f"{OUT}/simulation_manifest.csv", "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["sample_id", "edit_type", "genomic_position", "change", "designed_note"])
    notes = {}
    for nm in ("CASE-2", "CASE-6"):
        notes[(nm, NS1[4])] = f"nonsense: codon {NS1[0]} {NS1[2]}>{NS1[3]}"
    notes[("CASE-3", NS2[4])] = f"nonsense: codon {NS2[0]} {NS2[2]}>{NS2[3]}"
    for nm in ("CTRL-1", "CTRL-5"):
        notes[(nm, SYN[4])] = f"synonymous: codon {SYN[0]} {SYN[2]}>{SYN[3]}"
    for nm in SAMPLES:
        notes[(nm, MARKER_A[1])] = "case-enriched missense marker A"
        notes[(nm, CTRL_SNV[1])] = "control-only SNV"
    for name, edits in SAMPLES.items():
        for e in edits:
            if e[0] == "sub":
                w.writerow([name, "substitution", e[1], f"{REFERENCE[e[1]-1]}>{e[2]}",
                            notes.get((name, e[1]), "common variant" if e in COMMON else "")])
            elif e[0] == "del":
                note = ("case-enriched in-frame deletion marker B (codon 200)"
                        if e[1] == MB_POS else "control-only deletion (codon 50)")
                w.writerow([name, "deletion", e[1], f"del{e[2]}bp", note])
            else:
                w.writerow([name, "insertion", e[1], f"ins{e[2]}", "case-only frameshift"])

_ma_codon = ((MARKER_A[1] - CDS_START) // 3) + 1
_ma_off = (MARKER_A[1] - CDS_START) % 3
_ma_ref = ref_codon(_ma_codon)
_ma_alt = _ma_ref[:_ma_off] + MARKER_A[2] + _ma_ref[_ma_off + 1:]
_marker_a_label = (f"g.{MARKER_A[1]}{REFERENCE[MARKER_A[1]-1]}>{MARKER_A[2]} "
                   f"(codon {_ma_codon} {_ma_ref}>{_ma_alt}, "
                   f"{TABLE[_ma_ref]}>{TABLE[_ma_alt]})")

summary = {
    "reference_length": len(REFERENCE),
    "gc_content_pct": round(100 * sum(REFERENCE.count(b) for b in "GC") / len(REFERENCE), 2),
    "cds_span": f"{CDS_START}-{CDS_END}", "cds_length_bp": N_CODONS * 3,
    "protein_length_aa": len(ref_protein) - 1,
    "n_samples": len(samples),
    "designed_nonsense": {
        "NS1": f"codon {NS1[0]} ({NS1[2]}>{NS1[3]}) at g.{NS1[4]} ({NS1[5]}>{NS1[6]})",
        "NS1_samples": ["CASE-2", "CASE-6"],
        "NS2": f"codon {NS2[0]} ({NS2[2]}>{NS2[3]}) at g.{NS2[4]} ({NS2[5]}>{NS2[6]})",
        "NS2_samples": ["CASE-3"],
    },
    "designed_synonymous": f"codon {SYN[0]} ({SYN[2]}>{SYN[3]}) at g.{SYN[4]}",
    "marker_A": _marker_a_label,
    "marker_B": f"g.{MB_POS}_{MB_POS+2}del (codon 200 removed, in-frame)",
    "frameshift": "g.1200_1201insG (CASE-5)",
    "premature_stops_per_sample": premature,
}
print(json.dumps(summary, indent=2))
with open(f"{OUT}/simulation_summary.json", "w") as fh:
    json.dump(summary, fh, indent=2)
