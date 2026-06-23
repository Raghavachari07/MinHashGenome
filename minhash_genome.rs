use std::collections::BinaryHeap;
use std::env;
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::Path;

// minhash genome distance — no alignment, just kmer sketches
// usage: minhash_genome [--k 21] [--s 1000] genome1.fasta genome2.fasta ...

fn main() {
    let args: Vec<String> = env::args().skip(1).collect();

    if args.len() < 2 {
        eprintln!("usage: minhash_genome [--k <int>] [--s <int>] file1.fasta file2.fasta ...");
        std::process::exit(1);
    }

    let mut k = 21usize;
    let mut s = 1000usize;
    let mut files = vec![];

    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "--k" => { i += 1; k = args[i].parse().expect("--k needs an integer"); }
            "--s" => { i += 1; s = args[i].parse().expect("--s needs an integer"); }
            f     => files.push(f.to_string()),
        }
        i += 1;
    }

    eprintln!("k={k} s={s}");

    let (sketches, names): (Vec<_>, Vec<_>) = files.iter().map(|path| {
        let name = Path::new(path)
            .file_name().unwrap().to_str().unwrap()
            .trim_end_matches(".fasta")
            .trim_end_matches(".fna")
            .to_string();

        let seq = read_fasta(path).unwrap_or_else(|e| panic!("couldn't read {path}: {e}"));
        eprintln!("  {name}: {} bp", seq.len());

        (sketch(&seq, k, s), name)
    }).unzip();

    // distance matrix
    print!("\t");
    for n in &names { print!("{n}\t"); }
    println!();

    for i in 0..names.len() {
        print!("{}\t", names[i]);
        for j in 0..names.len() {
            print!("{:.3}\t", 1.0 - jaccard(&sketches[i], &sketches[j], s));
        }
        println!();
    }
}

fn read_fasta(path: &str) -> std::io::Result<String> {
    let mut seq = String::new();
    for line in BufReader::new(File::open(path)?).lines() {
        let line = line?;
        let line = line.trim();
        if !line.starts_with('>') && !line.is_empty() {
            seq.push_str(&line.to_uppercase());
        }
    }
    Ok(seq)
}

fn sketch(seq: &str, k: usize, s: usize) -> Vec<u64> {
    let bytes = seq.as_bytes();
    let mut heap: BinaryHeap<u64> = BinaryHeap::new(); // max-heap; tracks s smallest

    for w in bytes.windows(k) {
        if w.contains(&b'N') { continue; }

        let h = fnv1a(&canonical(w));

        if heap.len() < s {
            heap.push(h);
        } else if h < *heap.peek().unwrap() {
            heap.pop();
            heap.push(h);
        }
    }

    let mut out: Vec<u64> = heap.into_vec();
    out.sort_unstable();
    out
}

fn canonical(kmer: &[u8]) -> Vec<u8> {
    let rc = revcomp(kmer);
    if kmer <= rc.as_slice() { kmer.to_vec() } else { rc }
}

fn revcomp(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|&b| match b {
        b'A' => b'T', b'T' => b'A',
        b'C' => b'G', b'G' => b'C',
        _    => b'N',
    }).collect()
}

fn fnv1a(s: &[u8]) -> u64 {
    s.iter().fold(0xcbf29ce484222325u64, |h, &b| {
        h.wrapping_mul(0x00000100000001B3) ^ b as u64
    })
}

// merge-walk the two sorted sketches and count shared hashes
fn jaccard(a: &[u64], b: &[u64], s: usize) -> f64 {
    let limit = s.min(a.len()).min(b.len());
    let (mut shared, mut total, mut i, mut j) = (0, 0, 0, 0);

    while total < limit && i < a.len() && j < b.len() {
        match a[i].cmp(&b[j]) {
            std::cmp::Ordering::Equal   => { shared += 1; total += 1; i += 1; j += 1; }
            std::cmp::Ordering::Less    => { total += 1; i += 1; }
            std::cmp::Ordering::Greater => { total += 1; j += 1; }
        }
    }
    total += (limit - total).min(a.len() - i);
    total += (limit - total).min(b.len() - j);

    if total == 0 { 0.0 } else { shared as f64 / total as f64 }
}
