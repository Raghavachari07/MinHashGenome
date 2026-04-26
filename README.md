# MinHashGenome

A command-line tool in Java that estimates similarity between microbial genomes using **k-mer hashing and MinHash sketches** — no sequence alignment needed.

Given a set of bacterial genome files (FASTA format), it outputs a pairwise distance matrix. Genomes that are biologically close (e.g., two *E. coli* strains) will have a low distance score; unrelated species will score near 1.0.

---

## How it works

### Step 1 — Parse FASTA
Reads each `.fasta` / `.fna` file and concatenates all sequence lines into one string, skipping header lines that start with `>`.

### Step 2 — Canonical k-mer counting
Slides a window of length `k` (default: 21) across the genome sequence. For each k-mer, it also computes its **reverse complement** and keeps whichever is lexicographically smaller. This ensures that a sequence and its complement are treated as the same genomic feature.

### Step 3 — MinHash sketch
Instead of storing all k-mers (millions per genome), the tool hashes each canonical k-mer using **FNV-1a 64-bit hashing** and keeps only the `s` smallest hash values (default: s=1000). This compact set of s values is the **MinHash sketch** — a fixed-size fingerprint of the genome.

### Step 4 — Jaccard similarity
Two sketches are compared by merging their sorted hash arrays and counting how many of the combined-smallest values appear in both. This estimates the **Jaccard similarity** (shared k-mers / total k-mers). Distance = 1 − similarity.

---

## Data structures used

| Structure | Where | Why |
|---|---|---|
| `StringBuilder` | FASTA parser, reverse complement | Efficient string building without repeated concatenation |
| `PriorityQueue<Long>` (max-heap) | MinHash sketch builder | Maintains the s smallest hashes seen so far in O(log s) per insertion |
| `long[]` (sorted array) | Sketch storage | Enables O(s) linear merge for Jaccard computation |
| `List<String>`, `List<long[]>` | Main method | Stores genome names and their sketches in order |

The key algorithmic insight is the **max-heap of fixed size s**: if a new hash is smaller than the current maximum in the heap, swap it in. This runs in O(L · log s) time where L is genome length — much faster than sorting all k-mers.

---

## Usage

```bash
# Compile
javac MinHashGenome.java

# Run on genome files
java MinHashGenome ecoli.fna salmonella.fna klebsiella.fna

# Custom k and sketch size
java MinHashGenome *.fna --k 31 --sketch 2000
```

---

## Example output

```
genome        ecoli    salmonella    klebsiella
ecoli         0.0000   0.2100        0.7800
salmonella    0.2100   0.0000        0.8100
klebsiella    0.7800   0.8100        0.0000
```

Low distance = genomes are similar. *E. coli* and *Salmonella* are both Enterobacteriaceae so they cluster closer together than *Klebsiella*.

---

## Getting genome files

Download bacterial reference genomes free from NCBI RefSeq:

```bash
# Install NCBI datasets CLI, then:
./datasets download genome accession GCF_000005845.2 --filename ecoli.zip   # E. coli K-12
./datasets download genome accession GCF_000006945.2 --filename salmonella.zip
./datasets download genome accession GCF_000240185.1 --filename klebsiella.zip
```

Or browse at: https://www.ncbi.nlm.nih.gov/datasets/genome/

---

## Requirements

- Java 11 or higher
- No external dependencies — standard library only
