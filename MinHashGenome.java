import java.io.*;
import java.util.*;
import java.nio.file.*;

// Tool to compare bacterial genomes using minhash sketches
// No alignment needed - just hashing kmers and comparing
// Usage: java MinHashGenome genome1.fasta genome2.fasta ...

public class MinHashGenome {

    static int K = 21;        // kmer length, 21 is standard for bacterial genomes
    static int S = 1000;      // sketch size

    public static void main(String[] args) throws IOException {

        if (args.length < 2) {
            System.out.println("Need at least 2 fasta files to compare");
            System.out.println("Usage: java MinHashGenome file1.fasta file2.fasta ...");
            return;
        }

        List<String> files = new ArrayList<>();
        for (int i = 0; i < args.length; i++) {
            if (args[i].equals("--k")) {
                K = Integer.parseInt(args[++i]);
            } else if (args[i].equals("--s")) {
                S = Integer.parseInt(args[++i]);
            } else {
                files.add(args[i]);
            }
        }

        System.out.println("k=" + K + " sketch_size=" + S);

        List<long[]> sketches = new ArrayList<>();
        List<String> names = new ArrayList<>();

        for (String f : files) {
            String name = f.replaceAll(".*/", "").replaceAll("\\.fna|\\.fasta", "");
            System.out.print("reading " + name + "... ");

            String seq = readFasta(f);
            long[] sk = minhash(seq);
            sketches.add(sk);
            names.add(name);
            System.out.println(seq.length() + " bp");
        }

        // print distance matrix
        System.out.print("\t");
        for (String n : names) System.out.print(n + "\t");
        System.out.println();

        for (int i = 0; i < names.size(); i++) {
            System.out.print(names.get(i) + "\t");
            for (int j = 0; j < names.size(); j++) {
                double dist = 1.0 - jaccard(sketches.get(i), sketches.get(j));
                System.out.printf("%.3f\t", dist);
            }
            System.out.println();
        }
    }

    // read fasta and return sequence as one big string
    static String readFasta(String path) throws IOException {
        StringBuilder seq = new StringBuilder();
        BufferedReader br = new BufferedReader(new FileReader(path));
        String line;
        while ((line = br.readLine()) != null) {
            line = line.trim();
            if (line.startsWith(">") || line.isEmpty()) continue;
            seq.append(line.toUpperCase());
        }
        br.close();
        return seq.toString();
    }

    // build minhash sketch from sequence
    // keep S smallest hash values seen across all kmers
    static long[] minhash(String seq) {
        // max-heap so we can efficiently track the S smallest values
        PriorityQueue<Long> heap = new PriorityQueue<>(Collections.reverseOrder());

        for (int i = 0; i <= seq.length() - K; i++) {
            String kmer = seq.substring(i, i + K);
            if (kmer.contains("N")) continue;

            // use canonical kmer (smaller of kmer and reverse complement)
            String canon = canonical(kmer);
            long h = hash(canon);

            if (heap.size() < S) {
                heap.add(h);
            } else if (h < heap.peek()) {
                heap.poll();
                heap.add(h);
            }
        }

        long[] sketch = new long[heap.size()];
        int i = 0;
        for (long h : heap) sketch[i++] = h;
        Arrays.sort(sketch);
        return sketch;
    }

    // return whichever is smaller: kmer or its reverse complement
    // this way ACGT and its rc map to the same thing
    static String canonical(String kmer) {
        String rc = revcomp(kmer);
        return kmer.compareTo(rc) < 0 ? kmer : rc;
    }

    static String revcomp(String s) {
        char[] rc = new char[s.length()];
        for (int i = 0; i < s.length(); i++) {
            char c = s.charAt(s.length() - 1 - i);
            if      (c == 'A') rc[i] = 'T';
            else if (c == 'T') rc[i] = 'A';
            else if (c == 'C') rc[i] = 'G';
            else if (c == 'G') rc[i] = 'C';
            else               rc[i] = 'N';
        }
        return new String(rc);
    }

    // FNV-1a hash, works well for short strings like kmers
    static long hash(String s) {
        long h = 0xcbf29ce484222325L;
        for (int i = 0; i < s.length(); i++) {
            h ^= s.charAt(i);
            h *= 0x00000100000001B3L;
        }
        return h;
    }

    // estimate jaccard similarity by merging the two sorted sketch arrays
    // count how many of the S smallest values appear in both
    static double jaccard(long[] a, long[] b) {
        int size = Math.min(S, Math.min(a.length, b.length));
        int shared = 0, total = 0;
        int i = 0, j = 0;

        while (total < size && i < a.length && j < b.length) {
            if (a[i] == b[j])      { shared++; total++; i++; j++; }
            else if (a[i] < b[j]) { total++; i++; }
            else                   { total++; j++; }
        }
        while (total < size && i < a.length) { total++; i++; }
        while (total < size && j < b.length) { total++; j++; }

        return total == 0 ? 0 : (double) shared / total;
    }
}
