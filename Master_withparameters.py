import tkinter as tk
from tkinter import filedialog, messagebox, scrolledtext
from Bio import SeqIO
import primer3
import csv
import os

def load_genome(fasta_file):
    """Load the reference genome into a SeqRecord dict."""
    return SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))

def parse_gff3(gff_file):
    """Yield gene feature dicts with locus_tag, start, end, strand, and gene name."""
    genes = []
    with open(gff_file) as handle:
        for line in handle:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.strip().split("\t")
            if len(fields) < 9 or fields[2] != "gene":
                continue
            chrom, _, _, start, end, _, strand, _, attr = fields
            locus_tag = "unknown"
            gene_name = None
            for token in attr.split(";"):
                if token.startswith("locus_tag="):
                    locus_tag = token.split("=")[1]
                elif token.startswith("gene="):
                    gene_name = token.split("=")[1]
            genes.append({
                "chrom": chrom,
                "locus_tag": locus_tag,
                "gene_name": gene_name,
                "start": int(start) - 1,    # GFF3 1-based → python 0-based
                "end": int(end),            # end is exclusive in python
                "strand": 1 if strand == "+" else -1
            })
    return genes

def calc_gc(seq):
    """Calculate GC % for a DNA sequence string."""
    seq = seq.upper()
    gc = seq.count('G') + seq.count('C')
    return round(100.0 * gc / len(seq), 2) if seq else 0.0

def extract_sequences_to_fasta(genome, gene_list, output_file, flank_size=100):
    """Extract flanks and gene region for each gene and write to FASTA with banners."""
    with open(output_file, "w") as out_f:
        for gene in gene_list:
            chrom = gene["chrom"]
            gene_name = gene["gene_name"] if gene["gene_name"] else gene["locus_tag"]
            locus = gene["locus_tag"]
            strand = gene["strand"]
            start = gene["start"]
            end = gene["end"]
            seq_obj = genome.get(chrom)
            if not seq_obj:
                continue

            # Upstream (5' flank)
            five_seq = seq_obj.seq[max(0, start - flank_size):start]
            # Downstream (3' flank)
            three_seq = seq_obj.seq[end:min(len(seq_obj.seq), end + flank_size)]
            # Full gene sequence
            gene_seq = seq_obj.seq[start:end]

            if strand == -1:
                five_seq = five_seq.reverse_complement()
                three_seq = three_seq.reverse_complement()
                gene_seq = gene_seq.reverse_complement()

            def write_record(identifier, seq, gene_name, region_type):
                out_f.write(f"; ---------- {region_type} ----------\n")
                out_f.write(f">{identifier} | gene={gene_name}\n")
                seq_str = str(seq)
                for i in range(0, len(seq_str), 60):
                    out_f.write(seq_str[i:i+60] + "\n")
                out_f.write("\n")

            write_record(f"{locus}_5prime_flank", five_seq, gene_name, "5prime_flank")
            write_record(f"{locus}_3prime_flank", three_seq, gene_name, "3prime_flank")
            write_record(f"{locus}_gene", gene_seq, gene_name, "gene")

def run_pipeline():
    genome_path = genome_entry.get()
    gff_path = gff_entry.get()
    output_csv = output_entry.get()

    if not genome_path or not gff_path or not output_csv:
        messagebox.showerror("Error", "Please select genome FASTA, GFF3, and output CSV file.")
        return

    # --- Parameters ---
    try:
        min_len = int(min_size_entry.get())
        opt_len = int(opt_size_entry.get())
        max_len = int(max_size_entry.get())
        min_tm_val = float(min_tm.get())
        opt_tm_val = float(opt_tm.get())
        max_tm_val = float(max_tm.get())
        min_gc_val = float(min_gc.get())
        max_gc_val = float(max_gc.get())
        prod_min = int(min_prod.get())
        prod_max = int(max_prod.get())
        required_flank = int(flank_size_entry.get())
    except Exception as e:
        messagebox.showerror("Error", "Please enter valid numbers for primer parameters.")
        return

    DEFAULT_SETTINGS = {
        'PRIMER_OPT_SIZE': opt_len,
        'PRIMER_MIN_SIZE': min_len,
        'PRIMER_MAX_SIZE': max_len,
        'PRIMER_OPT_TM': opt_tm_val,
        'PRIMER_MIN_TM': min_tm_val,
        'PRIMER_MAX_TM': max_tm_val,
        'PRIMER_MIN_GC': min_gc_val,
        'PRIMER_MAX_GC': max_gc_val,
        'PRIMER_NUM_RETURN': 1,
        'PRIMER_PRODUCT_SIZE_RANGE': [[prod_min, prod_max]]
    }

    # --- Load files and extract genes ---
    try:
        genome = load_genome(genome_path)
    except Exception as e:
        messagebox.showerror("Error", f"Failed to load genome FASTA: {e}")
        return

    try:
        gene_list = parse_gff3(gff_path)
    except Exception as e:
        messagebox.showerror("Error", f"Failed to parse GFF3: {e}")
        return

    # --- 1. Output extracted gene sequences as FASTA with region separation ---
    extracted_fasta = os.path.splitext(output_csv)[0] + "_extracted_sequences.fasta"
    extract_sequences_to_fasta(genome, gene_list, extracted_fasta, flank_size=required_flank)

    # --- 2. Primer design step with Tm/GC reporting ---
    n_success = 0
    with open(output_csv, "w", newline="") as outcsv:
        writer = csv.writer(outcsv)
        # Write parameter summary row
        writer.writerow([
            f"Parameters: Min Size={min_len}, Opt Size={opt_len}, Max Size={max_len}, "
            f"Min Tm={min_tm_val}, Opt Tm={opt_tm_val}, Max Tm={max_tm_val}, "
            f"Min GC={min_gc_val}, Max GC={max_gc_val}, "
            f"Min Product={prod_min}, Max Product={prod_max}, Flank Size={required_flank}"
        ])
        writer.writerow([
            "Gene Name",
            "Forward Primer (5' flank)", "Fwd Tm", "Fwd GC%",
            "Reverse Primer (3' flank)", "Rev Tm", "Rev GC%",
            "Amplified Product Length", "Status"
        ])
        for gene in gene_list:
            locus = gene["locus_tag"]
            gene_name = gene["gene_name"] if gene["gene_name"] else locus
            chrom_seq = genome[gene["chrom"]].seq
            start, end, strand = gene["start"], gene["end"], gene["strand"]

            # Flanking regions for primer design
            five_seq = chrom_seq[max(0, start - required_flank):start]
            three_seq = chrom_seq[end:min(len(chrom_seq), end + required_flank)]
            if strand == -1:
                five_seq = five_seq.reverse_complement()
                three_seq = three_seq.reverse_complement()

            # Forward primer from 5' flank
            primers_5 = primer3.bindings.design_primers(
                {'SEQUENCE_ID': gene_name + '_5prime_flank', 'SEQUENCE_TEMPLATE': str(five_seq)},
                DEFAULT_SETTINGS
            )
            fwd_primer = primers_5.get('PRIMER_LEFT_0_SEQUENCE')

            # Reverse primer from 3' flank
            primers_3 = primer3.bindings.design_primers(
                {'SEQUENCE_ID': gene_name + '_3prime_flank', 'SEQUENCE_TEMPLATE': str(three_seq)},
                DEFAULT_SETTINGS
            )
            rev_primer = primers_3.get('PRIMER_RIGHT_0_SEQUENCE')

            product_len = (end - start) + 2 * required_flank
            ok = fwd_primer and rev_primer

            # Calculate per-primer stats (Tm, GC)
            fwd_tm = round(primer3.calcTm(fwd_primer), 2) if fwd_primer else ""
            fwd_gc = calc_gc(fwd_primer) if fwd_primer else ""
            rev_tm = round(primer3.calcTm(rev_primer), 2) if rev_primer else ""
            rev_gc = calc_gc(rev_primer) if rev_primer else ""

            writer.writerow([
                gene_name,
                fwd_primer if fwd_primer else "No primer found", fwd_tm, fwd_gc,
                rev_primer if rev_primer else "No primer found", rev_tm, rev_gc,
                product_len if ok else "N/A",
                "OK" if ok else "No suitable primers"
            ])
            if ok:
                n_success += 1

    # --- GUI Feedback ---
    results_text.delete("1.0", tk.END)
    results_text.insert(
        tk.END,
        f"✅ Primer design complete for {len(gene_list)} genes\n"
        f"{n_success}/{len(gene_list)} had suitable primer pairs.\n"
        f"Output:\n• CSV: {os.path.abspath(output_csv)}\n"
        f"• Extracted genes: {os.path.abspath(extracted_fasta)}"
    )

# --- GUI setup ---
root = tk.Tk()
root.title("PCR Primer Design – with Gene Extraction FASTA")

tk.Label(root, text="Genome FASTA").grid(row=0, column=0, sticky="w")
genome_entry = tk.Entry(root, width=50); genome_entry.grid(row=0, column=1)
tk.Button(root, text="Browse",
    command=lambda: genome_entry.delete(0, tk.END) or genome_entry.insert(0, filedialog.askopenfilename())
).grid(row=0, column=2)

tk.Label(root, text="GFF3 File").grid(row=1, column=0, sticky="w")
gff_entry = tk.Entry(root, width=50); gff_entry.grid(row=1, column=1)
tk.Button(root, text="Browse",
    command=lambda: gff_entry.delete(0, tk.END) or gff_entry.insert(0, filedialog.askopenfilename())
).grid(row=1, column=2)

tk.Label(root, text="Output CSV").grid(row=2, column=0, sticky="w")
output_entry = tk.Entry(root, width=50); output_entry.insert(0, "primers.csv")
output_entry.grid(row=2, column=1)
tk.Button(root, text="Browse",
    command=lambda: output_entry.delete(0, tk.END) or output_entry.insert(0, filedialog.asksaveasfilename(defaultextension=".csv"))
).grid(row=2, column=2)

tk.Label(root, text="Min Size").grid(row=3, column=0)
min_size_entry = tk.Entry(root, width=5); min_size_entry.insert(0,"18"); min_size_entry.grid(row=3, column=1)
tk.Label(root, text="Opt Size").grid(row=3, column=2)
opt_size_entry = tk.Entry(root, width=5); opt_size_entry.insert(0,"20"); opt_size_entry.grid(row=3, column=3)
tk.Label(root, text="Max Size").grid(row=3, column=4)
max_size_entry = tk.Entry(root, width=5); max_size_entry.insert(0,"25"); max_size_entry.grid(row=3, column=5)

tk.Label(root, text="Min Tm").grid(row=4, column=0)
min_tm = tk.Entry(root, width=5); min_tm.insert(0,"58.0"); min_tm.grid(row=4, column=1)
tk.Label(root, text="Opt Tm").grid(row=4, column=2)
opt_tm = tk.Entry(root, width=5); opt_tm.insert(0,"60.0"); opt_tm.grid(row=4, column=3)
tk.Label(root, text="Max Tm").grid(row=4, column=4)
max_tm = tk.Entry(root, width=5); max_tm.insert(0,"60.0"); max_tm.grid(row=4, column=5)

tk.Label(root, text="Min GC").grid(row=5, column=0)
min_gc = tk.Entry(root, width=5); min_gc.insert(0,"20.0"); min_gc.grid(row=5, column=1)
tk.Label(root, text="Max GC").grid(row=5, column=2)
max_gc = tk.Entry(root, width=5); max_gc.insert(0,"50.0"); max_gc.grid(row=5, column=3)

tk.Label(root, text="Min Product").grid(row=5, column=4)
min_prod = tk.Entry(root, width=5); min_prod.insert(0,"20"); min_prod.grid(row=5, column=5)
tk.Label(root, text="Max Product").grid(row=6, column=0)
max_prod = tk.Entry(root, width=5); max_prod.insert(0,"5000"); max_prod.grid(row=6, column=1)
tk.Label(root, text="Flank Size (bp)").grid(row=6, column=2)
flank_size_entry = tk.Entry(root, width=5); flank_size_entry.insert(0,"100"); flank_size_entry.grid(row=6, column=3)

tk.Button(root, text="Run", command=run_pipeline, bg="lightgreen").grid(row=7, column=0, columnspan=6, pady=10)
results_text = scrolledtext.ScrolledText(root, width=80, height=15)
results_text.grid(row=8, column=0, columnspan=6)
root.mainloop()
