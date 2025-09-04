
import tkinter as tk
from tkinter import filedialog, messagebox, scrolledtext, ttk
from Bio import SeqIO
from Bio.Seq import Seq
import primer3
import csv
import os
import re

def load_genome(fasta_file):
    """Load the reference genome into a SeqRecord dict."""
    return SeqIO.to_dict(SeqIO.parse(fasta_file, "fasta"))

def load_cds_sequences(cds_file):
    """Load CDS sequences from a FASTA file."""
    cds_dict = {}
    try:
        for record in SeqIO.parse(cds_file, "fasta"):
            # Extract gene name from header (handle different formats)
            header = record.description
            gene_name = None

            # Try different patterns to extract gene name
            if 'gene=' in header:
                gene_name = header.split('gene=')[1].split()[0].strip('[]')
            elif 'locus_tag=' in header:
                gene_name = header.split('locus_tag=')[1].split()[0].strip('[]')
            else:
                # Use the record ID as fallback
                gene_name = record.id.split('|')[0] if '|' in record.id else record.id

            if gene_name:
                cds_dict[gene_name.lower()] = {
                    'sequence': str(record.seq),
                    'id': record.id,
                    'description': record.description
                }
    except Exception as e:
        print(f"Error loading CDS file: {e}")

    return cds_dict

def parse_gff3_full(gff_file):
    """Parse full GFF3 file and return all genes."""
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
                "start": int(start) - 1,
                "end": int(end),
                "strand": 1 if strand == "+" else -1
            })
    return genes

def filter_genes_by_names(all_genes, target_names):
    """Filter genes by user-specified names (case-insensitive)."""
    if not target_names.strip():
        return all_genes, [], []

    # Parse target names (comma or newline separated)
    names = [name.strip().lower() for name in re.split(r'[,\n\r]+', target_names) if name.strip()]

    filtered_genes = []
    found_names = set()

    for gene in all_genes:
        gene_name_lower = (gene["gene_name"] or "").lower()
        locus_tag_lower = (gene["locus_tag"] or "").lower()

        # Check if gene name or locus tag matches any target name
        for target_name in names:
            if target_name in [gene_name_lower, locus_tag_lower]:
                filtered_genes.append(gene)
                found_names.add(target_name)
                break

    # Report which names were found/not found
    not_found = set(names) - found_names
    return filtered_genes, list(found_names), list(not_found)

def calc_gc(seq):
    """Calculate GC % for a DNA sequence string."""
    if not seq:
        return 0.0
    seq = seq.upper()
    clean_seq = re.sub(r'[^ATGC]', '', seq)
    if not clean_seq:
        return 0.0
    gc = clean_seq.count('G') + clean_seq.count('C')
    return round(100.0 * gc / len(clean_seq), 2)

def calc_tm_safe(seq):
    """Calculate Tm using primer3 with error handling."""
    if not seq:
        return 0.0
    try:
        return round(primer3.calc_tm(seq), 2)
    except:
        try:
            return round(primer3.calcTm(seq), 2)
        except:
            return 0.0

def clean_sequence(seq):
    """Clean sequence to contain only ATGC characters."""
    if not seq:
        return ""
    return re.sub(r'[^ATGC]', '', seq.upper())

def search_primer_specificity(forward_primer, reverse_primer, genome_fasta, max_product_size=5000):
    """Robust primer specificity search."""
    try:
        fwd_clean = clean_sequence(forward_primer)
        rev_clean = clean_sequence(reverse_primer)

        if not fwd_clean or not rev_clean or len(fwd_clean) < 10 or len(rev_clean) < 10:
            return "Invalid primers"

        try:
            genome = SeqIO.to_dict(SeqIO.parse(genome_fasta, "fasta"))
        except:
            return "Genome loading failed"

        if not genome:
            return "Empty genome"

        total_amplicons = 0

        try:
            rev_rc = str(Seq(rev_clean).reverse_complement())
        except:
            return "Reverse complement failed"

        for chrom_id, record in genome.items():
            try:
                seq_str = clean_sequence(str(record.seq))
                if not seq_str:
                    continue

                # Simple exact matching
                fwd_positions = []
                start = 0
                while True:
                    pos = seq_str.find(fwd_clean, start)
                    if pos == -1:
                        break
                    fwd_positions.append(pos)
                    start = pos + 1

                rev_positions = []
                start = 0
                while True:
                    pos = seq_str.find(rev_rc, start)
                    if pos == -1:
                        break
                    rev_positions.append(pos)
                    start = pos + 1

                for fwd_pos in fwd_positions:
                    for rev_pos in rev_positions:
                        if fwd_pos < rev_pos:
                            product_size = rev_pos - fwd_pos + len(rev_clean)
                            if 50 <= product_size <= max_product_size:
                                total_amplicons += 1

            except:
                continue

        if total_amplicons == 0:
            return "No amplicons found"
        elif total_amplicons == 1:
            return "Specific (1 amplicon)"
        else:
            return f"Non-specific ({total_amplicons} amplicons)"

    except Exception as e:
        return f"Analysis failed: {str(e)[:50]}"

def extract_sequences_to_fasta(genome, gene_list, output_file, flank_size=100):
    """Extract flanks and gene region for each gene and write to FASTA."""
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

            five_seq = seq_obj.seq[max(0, start - flank_size):start]
            three_seq = seq_obj.seq[end:min(len(seq_obj.seq), end + flank_size)]
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

def preview_genes():
    """Preview which genes will be processed based on the filter."""
    input_mode = mode_var.get()

    if input_mode == "genome_gff":
        gff_path = gff_entry.get()
        if not gff_path:
            messagebox.showerror("Error", "Please select a GFF3 file first.")
            return

        try:
            all_genes = parse_gff3_full(gff_path)
            target_names = gene_filter_text.get("1.0", tk.END)

            if target_names.strip():
                filtered_genes, found_names, not_found = filter_genes_by_names(all_genes, target_names)
            else:
                filtered_genes = all_genes
                found_names = []
                not_found = []

            show_gene_preview(all_genes, filtered_genes, found_names, not_found)

        except Exception as e:
            messagebox.showerror("Error", f"Failed to preview genes: {e}")

    elif input_mode == "cds_only":
        cds_path = cds_entry.get()
        if not cds_path:
            messagebox.showerror("Error", "Please select a CDS FASTA file first.")
            return

        try:
            cds_dict = load_cds_sequences(cds_path)
            target_names = gene_filter_text.get("1.0", tk.END)

            if target_names.strip():
                names = [name.strip().lower() for name in re.split(r'[,\n\r]+', target_names) if name.strip()]
                filtered_cds = {k: v for k, v in cds_dict.items() if k in names}
                found_names = list(filtered_cds.keys())
                not_found = [name for name in names if name not in cds_dict]
            else:
                filtered_cds = cds_dict
                found_names = []
                not_found = []

            show_cds_preview(cds_dict, filtered_cds, found_names, not_found)

        except Exception as e:
            messagebox.showerror("Error", f"Failed to preview CDS: {e}")

def show_gene_preview(all_genes, filtered_genes, found_names, not_found):
    """Show gene preview window."""
    preview_window = tk.Toplevel(root)
    preview_window.title("Gene Preview")
    preview_window.geometry("600x400")

    preview_text = scrolledtext.ScrolledText(preview_window, width=70, height=20)
    preview_text.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

    preview_content = f"📊 GENE PREVIEW\n"
    preview_content += f"=" * 50 + "\n"
    preview_content += f"Total genes in GFF: {len(all_genes)}\n"
    preview_content += f"Genes to process: {len(filtered_genes)}\n\n"

    if found_names:
        preview_content += f"✅ Found genes ({len(found_names)}): {', '.join(found_names[:10])}\n"
        if len(found_names) > 10:
            preview_content += f"... and {len(found_names) - 10} more\n"
        preview_content += "\n"

    if not_found:
        preview_content += f"❌ Not found ({len(not_found)}): {', '.join(not_found[:10])}\n"
        if len(not_found) > 10:
            preview_content += f"... and {len(not_found) - 10} more\n"
        preview_content += "\n"

    preview_content += f"GENES TO PROCESS:\n"
    preview_content += f"-" * 30 + "\n"

    for i, gene in enumerate(filtered_genes[:20]):
        gene_name = gene["gene_name"] or "N/A"
        preview_content += f"{i+1:2d}. {gene['locus_tag']} ({gene_name})\n"

    if len(filtered_genes) > 20:
        preview_content += f"... and {len(filtered_genes) - 20} more genes\n"

    preview_text.insert(tk.END, preview_content)
    preview_text.config(state=tk.DISABLED)

def show_cds_preview(all_cds, filtered_cds, found_names, not_found):
    """Show CDS preview window."""
    preview_window = tk.Toplevel(root)
    preview_window.title("CDS Preview")
    preview_window.geometry("600x400")

    preview_text = scrolledtext.ScrolledText(preview_window, width=70, height=20)
    preview_text.pack(fill=tk.BOTH, expand=True, padx=10, pady=10)

    preview_content = f"📊 CDS PREVIEW\n"
    preview_content += f"=" * 50 + "\n"
    preview_content += f"Total CDS in file: {len(all_cds)}\n"
    preview_content += f"CDS to process: {len(filtered_cds)}\n\n"

    if found_names:
        preview_content += f"✅ Found CDS ({len(found_names)}): {', '.join(found_names[:10])}\n"
        if len(found_names) > 10:
            preview_content += f"... and {len(found_names) - 10} more\n"
        preview_content += "\n"

    if not_found:
        preview_content += f"❌ Not found ({len(not_found)}): {', '.join(not_found[:10])}\n"
        if len(not_found) > 10:
            preview_content += f"... and {len(not_found) - 10} more\n"
        preview_content += "\n"

    preview_content += f"CDS TO PROCESS:\n"
    preview_content += f"-" * 30 + "\n"

    for i, (gene_name, cds_info) in enumerate(list(filtered_cds.items())[:20]):
        seq_len = len(cds_info['sequence'])
        preview_content += f"{i+1:2d}. {gene_name} ({seq_len} bp)\n"

    if len(filtered_cds) > 20:
        preview_content += f"... and {len(filtered_cds) - 20} more CDS\n"

    preview_text.insert(tk.END, preview_content)
    preview_text.config(state=tk.DISABLED)

def toggle_input_mode():
    """Toggle between genome+GFF and CDS-only modes."""
    mode = mode_var.get()

    if mode == "genome_gff":
        # Show genome+GFF inputs
        genome_label.grid()
        genome_entry.grid()
        genome_browse.grid()
        gff_label.grid()
        gff_entry.grid()
        gff_browse.grid()

        # Hide CDS input
        cds_label.grid_remove()
        cds_entry.grid_remove()
        cds_browse.grid_remove()

        # Show flank size
        flank_label.grid()
        flank_size_entry.grid()

    else:  # cds_only
        # Hide genome+GFF inputs
        genome_label.grid_remove()
        genome_entry.grid_remove()
        genome_browse.grid_remove()
        gff_label.grid_remove()
        gff_entry.grid_remove()
        gff_browse.grid_remove()

        # Show CDS input
        cds_label.grid()
        cds_entry.grid()
        cds_browse.grid()

        # Hide flank size (not needed for CDS)
        flank_label.grid_remove()
        flank_size_entry.grid_remove()

def run_pipeline():
    """Main pipeline execution."""
    input_mode = mode_var.get()
    output_csv = output_entry.get()

    if not output_csv:
        messagebox.showerror("Error", "Please specify an output CSV file.")
        return

    # Get parameters
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
        required_flank = int(flank_size_entry.get()) if input_mode == "genome_gff" else 0
    except ValueError:
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
        'PRIMER_PRODUCT_SIZE_RANGE': [[prod_min, prod_max]],
        'PRIMER_MAX_POLY_X': 5,
        'PRIMER_MAX_NS_ACCEPTED': 1,
        'PRIMER_MAX_SELF_ANY': 8.0,
        'PRIMER_MAX_SELF_END': 3.0,
        'PRIMER_PAIR_MAX_COMPL_ANY': 8.0,
        'PRIMER_PAIR_MAX_COMPL_END': 3.0
    }

    # Clear results
    results_text.delete("1.0", tk.END)
    results_text.insert(tk.END, "🔄 Starting primer design...\n")
    root.update()

    try:
        if input_mode == "genome_gff":
            run_genome_gff_pipeline(DEFAULT_SETTINGS, output_csv, required_flank)
        else:
            run_cds_pipeline(DEFAULT_SETTINGS, output_csv)
    except Exception as e:
        messagebox.showerror("Error", f"Pipeline failed: {e}")

def run_genome_gff_pipeline(DEFAULT_SETTINGS, output_csv, required_flank):
    """Run pipeline using genome FASTA + GFF3."""
    genome_path = genome_entry.get()
    gff_path = gff_entry.get()
    target_names = gene_filter_text.get("1.0", tk.END)

    if not genome_path or not gff_path:
        messagebox.showerror("Error", "Please select genome FASTA and GFF3 files.")
        return

    # Load data
    try:
        genome = load_genome(genome_path)
        results_text.insert(tk.END, "✅ Genome loaded\n")
        root.update()
    except Exception as e:
        messagebox.showerror("Error", f"Failed to load genome: {e}")
        return

    try:
        all_genes = parse_gff3_full(gff_path)
        if target_names.strip():
            gene_list, found_names, not_found = filter_genes_by_names(all_genes, target_names)
            results_text.insert(tk.END, f"✅ Filtered {len(gene_list)} genes from {len(all_genes)} total\n")
            if found_names:
                results_text.insert(tk.END, f"Found: {', '.join(found_names[:5])}{'...' if len(found_names) > 5 else ''}\n")
            if not_found:
                results_text.insert(tk.END, f"❌ Not found: {', '.join(not_found[:5])}{'...' if len(not_found) > 5 else ''}\n")
        else:
            gene_list = all_genes
            results_text.insert(tk.END, f"✅ Processing all {len(gene_list)} genes\n")
        root.update()
    except Exception as e:
        messagebox.showerror("Error", f"Failed to parse GFF3: {e}")
        return

    if not gene_list:
        messagebox.showerror("Error", "No genes found to process!")
        return

    # Extract sequences
    results_text.insert(tk.END, "🔄 Extracting sequences...\n")
    root.update()

    extracted_fasta = os.path.splitext(output_csv)[0] + "_extracted_sequences.fasta"
    extract_sequences_to_fasta(genome, gene_list, extracted_fasta, flank_size=required_flank)

    # Design primers
    results_text.insert(tk.END, "🔄 Designing primers...\n")
    root.update()

    n_success = 0

    with open(output_csv, "w", newline="") as outcsv:
        writer = csv.writer(outcsv)
        writer.writerow([
            f"Parameters: Min Size={DEFAULT_SETTINGS['PRIMER_MIN_SIZE']}, "
            f"Opt Size={DEFAULT_SETTINGS['PRIMER_OPT_SIZE']}, "
            f"Max Size={DEFAULT_SETTINGS['PRIMER_MAX_SIZE']}, "
            f"Min Tm={DEFAULT_SETTINGS['PRIMER_MIN_TM']}, "
            f"Opt Tm={DEFAULT_SETTINGS['PRIMER_OPT_TM']}, "
            f"Max Tm={DEFAULT_SETTINGS['PRIMER_MAX_TM']}, "
            f"Min GC={DEFAULT_SETTINGS['PRIMER_MIN_GC']}, "
            f"Max GC={DEFAULT_SETTINGS['PRIMER_MAX_GC']}, "
            f"Min Product={DEFAULT_SETTINGS['PRIMER_PRODUCT_SIZE_RANGE'][0][0]}, "
            f"Max Product={DEFAULT_SETTINGS['PRIMER_PRODUCT_SIZE_RANGE'][0][1]}, "
            f"Flank Size={required_flank}, Genes processed: {len(gene_list)}"
        ])

        writer.writerow([
            "Gene Name",
            "Forward Primer", "Fwd Tm", "Fwd GC%",
            "Reverse Primer", "Rev Tm", "Rev GC%",
            "Product Length", "Specificity Check", "Status"
        ])

        for i, gene in enumerate(gene_list):
            if i % 10 == 0:
                results_text.insert(tk.END, f"🔄 Processing gene {i+1}/{len(gene_list)}...\n")
                results_text.see(tk.END)
                root.update()

            gene_name = gene["gene_name"] if gene["gene_name"] else gene["locus_tag"]
            chrom_seq = genome[gene["chrom"]].seq
            start, end, strand = gene["start"], gene["end"], gene["strand"]

            # Extract gene sequence for primer design
            gene_seq = chrom_seq[start:end]
            if strand == -1:
                gene_seq = gene_seq.reverse_complement()

            # Try primer design
            primers = None
            try:
                primers = primer3.bindings.design_primers(
                    {'SEQUENCE_ID': gene_name, 'SEQUENCE_TEMPLATE': str(gene_seq)},
                    DEFAULT_SETTINGS
                )
            except:
                pass

            fwd_primer = primers.get('PRIMER_LEFT_0_SEQUENCE') if primers else None
            rev_primer = primers.get('PRIMER_RIGHT_0_SEQUENCE') if primers else None
            product_len = primers.get('PRIMER_PAIR_0_PRODUCT_SIZE', 'N/A') if primers else 'N/A'

            ok = fwd_primer and rev_primer

            fwd_tm = calc_tm_safe(fwd_primer) if fwd_primer else "N/A"
            fwd_gc = calc_gc(fwd_primer) if fwd_primer else "N/A"
            rev_tm = calc_tm_safe(rev_primer) if rev_primer else "N/A"
            rev_gc = calc_gc(rev_primer) if rev_primer else "N/A"

            specificity_result = "Not tested"
            if check_specificity.get() and ok:
                specificity_result = search_primer_specificity(
                    fwd_primer, rev_primer, genome_path, 
                    max_product_size=DEFAULT_SETTINGS['PRIMER_PRODUCT_SIZE_RANGE'][0][1]
                )

            if ok:
                n_success += 1
                status = "OK"
            else:
                status = "No suitable primers"
                fwd_primer = "No primer found"
                rev_primer = "No primer found"

            writer.writerow([
                gene_name, fwd_primer, fwd_tm, fwd_gc,
                rev_primer, rev_tm, rev_gc,
                product_len, specificity_result, status
            ])

    results_text.insert(tk.END, f"✅ Primer design complete!\n")
    results_text.insert(tk.END, f"📊 {n_success}/{len(gene_list)} genes had suitable primers\n")
    results_text.insert(tk.END, f"📄 Output: {os.path.abspath(output_csv)}\n")
    results_text.insert(tk.END, f"🧬 Sequences: {os.path.abspath(extracted_fasta)}\n")

def run_cds_pipeline(DEFAULT_SETTINGS, output_csv):
    """Run pipeline using CDS FASTA only."""
    cds_path = cds_entry.get()
    target_names = gene_filter_text.get("1.0", tk.END)

    if not cds_path:
        messagebox.showerror("Error", "Please select a CDS FASTA file.")
        return

    # Load CDS sequences
    try:
        all_cds = load_cds_sequences(cds_path)
        results_text.insert(tk.END, f"✅ Loaded {len(all_cds)} CDS sequences\n")
        root.update()
    except Exception as e:
        messagebox.showerror("Error", f"Failed to load CDS file: {e}")
        return

    # Filter CDS if needed
    if target_names.strip():
        names = [name.strip().lower() for name in re.split(r'[,\n\r]+', target_names) if name.strip()]
        filtered_cds = {k: v for k, v in all_cds.items() if k in names}
        found_names = list(filtered_cds.keys())
        not_found = [name for name in names if name not in all_cds]

        results_text.insert(tk.END, f"✅ Filtered {len(filtered_cds)} CDS from {len(all_cds)} total\n")
        if found_names:
            results_text.insert(tk.END, f"Found: {', '.join(found_names[:5])}{'...' if len(found_names) > 5 else ''}\n")
        if not_found:
            results_text.insert(tk.END, f"❌ Not found: {', '.join(not_found[:5])}{'...' if len(not_found) > 5 else ''}\n")
    else:
        filtered_cds = all_cds
        results_text.insert(tk.END, f"✅ Processing all {len(filtered_cds)} CDS\n")

    root.update()

    if not filtered_cds:
        messagebox.showerror("Error", "No CDS sequences found to process!")
        return

    # Design primers
    results_text.insert(tk.END, "🔄 Designing primers from CDS...\n")
    root.update()

    n_success = 0

    with open(output_csv, "w", newline="") as outcsv:
        writer = csv.writer(outcsv)
        writer.writerow([
            f"Parameters: Min Size={DEFAULT_SETTINGS['PRIMER_MIN_SIZE']}, "
            f"Opt Size={DEFAULT_SETTINGS['PRIMER_OPT_SIZE']}, "
            f"Max Size={DEFAULT_SETTINGS['PRIMER_MAX_SIZE']}, "
            f"Min Tm={DEFAULT_SETTINGS['PRIMER_MIN_TM']}, "
            f"Opt Tm={DEFAULT_SETTINGS['PRIMER_OPT_TM']}, "
            f"Max Tm={DEFAULT_SETTINGS['PRIMER_MAX_TM']}, "
            f"Min GC={DEFAULT_SETTINGS['PRIMER_MIN_GC']}, "
            f"Max GC={DEFAULT_SETTINGS['PRIMER_MAX_GC']}, "
            f"Min Product={DEFAULT_SETTINGS['PRIMER_PRODUCT_SIZE_RANGE'][0][0]}, "
            f"Max Product={DEFAULT_SETTINGS['PRIMER_PRODUCT_SIZE_RANGE'][0][1]}, "
            f"Mode=CDS_ONLY, CDS processed: {len(filtered_cds)}"
        ])

        writer.writerow([
            "Gene Name",
            "Forward Primer", "Fwd Tm", "Fwd GC%",
            "Reverse Primer", "Rev Tm", "Rev GC%",
            "Product Length", "Specificity Check", "Status"
        ])

        for i, (gene_name, cds_info) in enumerate(filtered_cds.items()):
            if i % 10 == 0:
                results_text.insert(tk.END, f"🔄 Processing CDS {i+1}/{len(filtered_cds)}...\n")
                results_text.see(tk.END)
                root.update()

            sequence = cds_info['sequence']

            # Try primer design
            primers = None
            try:
                primers = primer3.bindings.design_primers(
                    {'SEQUENCE_ID': gene_name, 'SEQUENCE_TEMPLATE': sequence},
                    DEFAULT_SETTINGS
                )
            except:
                pass

            fwd_primer = primers.get('PRIMER_LEFT_0_SEQUENCE') if primers else None
            rev_primer = primers.get('PRIMER_RIGHT_0_SEQUENCE') if primers else None
            product_len = primers.get('PRIMER_PAIR_0_PRODUCT_SIZE', 'N/A') if primers else 'N/A'

            ok = fwd_primer and rev_primer

            fwd_tm = calc_tm_safe(fwd_primer) if fwd_primer else "N/A"
            fwd_gc = calc_gc(fwd_primer) if fwd_primer else "N/A"
            rev_tm = calc_tm_safe(rev_primer) if rev_primer else "N/A"
            rev_gc = calc_gc(rev_primer) if rev_primer else "N/A"

            # Note: Specificity testing needs genome file, not available in CDS-only mode
            specificity_result = "N/A (CDS mode)"

            if ok:
                n_success += 1
                status = "OK"
            else:
                status = "No suitable primers"
                fwd_primer = "No primer found"
                rev_primer = "No primer found"

            writer.writerow([
                gene_name, fwd_primer, fwd_tm, fwd_gc,
                rev_primer, rev_tm, rev_gc,
                product_len, specificity_result, status
            ])

    results_text.insert(tk.END, f"✅ Primer design complete!\n")
    results_text.insert(tk.END, f"📊 {n_success}/{len(filtered_cds)} CDS had suitable primers\n")
    results_text.insert(tk.END, f"📄 Output: {os.path.abspath(output_csv)}\n")

# === GUI SETUP ===
root = tk.Tk()
root.title("🧬 Ultimate PCR Primer Design Suite")
root.geometry("900x800")

# Input mode selection
mode_frame = ttk.LabelFrame(root, text="Input Mode", padding="10")
mode_frame.grid(row=0, column=0, columnspan=6, sticky="ew", padx=10, pady=5)

mode_var = tk.StringVar(value="genome_gff")

tk.Radiobutton(mode_frame, text="Genome FASTA + GFF3 (Full Pipeline)", 
               variable=mode_var, value="genome_gff", command=toggle_input_mode).grid(row=0, column=0, sticky="w")
tk.Radiobutton(mode_frame, text="CDS FASTA Only (Direct Primer Design)", 
               variable=mode_var, value="cds_only", command=toggle_input_mode).grid(row=0, column=1, sticky="w")

# File selection frame
file_frame = ttk.LabelFrame(root, text="Input Files", padding="10")
file_frame.grid(row=1, column=0, columnspan=6, sticky="ew", padx=10, pady=5)

# Genome + GFF inputs (default visible)
genome_label = tk.Label(file_frame, text="Genome FASTA")
genome_label.grid(row=0, column=0, sticky="w")
genome_entry = tk.Entry(file_frame, width=50)
genome_entry.grid(row=0, column=1)
genome_browse = tk.Button(file_frame, text="Browse", 
    command=lambda: genome_entry.delete(0, tk.END) or genome_entry.insert(0, filedialog.askopenfilename()))
genome_browse.grid(row=0, column=2)

gff_label = tk.Label(file_frame, text="GFF3 File")
gff_label.grid(row=1, column=0, sticky="w")
gff_entry = tk.Entry(file_frame, width=50)
gff_entry.grid(row=1, column=1)
gff_browse = tk.Button(file_frame, text="Browse",
    command=lambda: gff_entry.delete(0, tk.END) or gff_entry.insert(0, filedialog.askopenfilename()))
gff_browse.grid(row=1, column=2)

# CDS input (initially hidden)
cds_label = tk.Label(file_frame, text="CDS FASTA")
cds_entry = tk.Entry(file_frame, width=50)
cds_browse = tk.Button(file_frame, text="Browse",
    command=lambda: cds_entry.delete(0, tk.END) or cds_entry.insert(0, filedialog.askopenfilename()))

# Output file
tk.Label(file_frame, text="Output CSV").grid(row=2, column=0, sticky="w")
output_entry = tk.Entry(file_frame, width=50)
output_entry.insert(0, "ultimate_primers.csv")
output_entry.grid(row=2, column=1)
tk.Button(file_frame, text="Browse",
    command=lambda: output_entry.delete(0, tk.END) or output_entry.insert(0, filedialog.asksaveasfilename(defaultextension=".csv"))
).grid(row=2, column=2)

# Gene filtering frame
filter_frame = ttk.LabelFrame(root, text="Gene/CDS Selection (Case Insensitive)", padding="10")
filter_frame.grid(row=2, column=0, columnspan=6, sticky="ew", padx=10, pady=5)

tk.Label(filter_frame, text="Enter gene names (comma or line separated):").grid(row=0, column=0, sticky="w")
tk.Label(filter_frame, text="Leave empty to process ALL sequences", fg="blue").grid(row=0, column=1, sticky="w")

gene_filter_text = scrolledtext.ScrolledText(filter_frame, width=60, height=4)
gene_filter_text.grid(row=1, column=0, columnspan=2, sticky="ew", pady=5)
gene_filter_text.insert(tk.END, "# Examples:\n# sulA, opgH, galU\n# mrcB, rpoS")

tk.Button(filter_frame, text="🔍 Preview Selection", command=preview_genes, 
          bg="lightblue", font=("Arial", 10, "bold")).grid(row=2, column=0, sticky="w", pady=5)

# Primer parameters frame
param_frame = ttk.LabelFrame(root, text="Primer Parameters", padding="10")
param_frame.grid(row=3, column=0, columnspan=6, sticky="ew", padx=10, pady=5)

# Size parameters
tk.Label(param_frame, text="Min Size").grid(row=0, column=0)
min_size_entry = tk.Entry(param_frame, width=5)
min_size_entry.insert(0, "18")
min_size_entry.grid(row=0, column=1)

tk.Label(param_frame, text="Opt Size").grid(row=0, column=2)
opt_size_entry = tk.Entry(param_frame, width=5)
opt_size_entry.insert(0, "20")
opt_size_entry.grid(row=0, column=3)

tk.Label(param_frame, text="Max Size").grid(row=0, column=4)
max_size_entry = tk.Entry(param_frame, width=5)
max_size_entry.insert(0, "25")
max_size_entry.grid(row=0, column=5)

# Tm parameters
tk.Label(param_frame, text="Min Tm").grid(row=1, column=0)
min_tm = tk.Entry(param_frame, width=5)
min_tm.insert(0, "55.0")
min_tm.grid(row=1, column=1)

tk.Label(param_frame, text="Opt Tm").grid(row=1, column=2)
opt_tm = tk.Entry(param_frame, width=5)
opt_tm.insert(0, "60.0")
opt_tm.grid(row=1, column=3)

tk.Label(param_frame, text="Max Tm").grid(row=1, column=4)
max_tm = tk.Entry(param_frame, width=5)
max_tm.insert(0, "65.0")
max_tm.grid(row=1, column=5)

# GC and product parameters
tk.Label(param_frame, text="Min GC%").grid(row=2, column=0)
min_gc = tk.Entry(param_frame, width=5)
min_gc.insert(0, "20.0")
min_gc.grid(row=2, column=1)

tk.Label(param_frame, text="Max GC%").grid(row=2, column=2)
max_gc = tk.Entry(param_frame, width=5)
max_gc.insert(0, "70.0")
max_gc.grid(row=2, column=3)

tk.Label(param_frame, text="Min Product").grid(row=2, column=4)
min_prod = tk.Entry(param_frame, width=5)
min_prod.insert(0, "100")
min_prod.grid(row=2, column=5)

tk.Label(param_frame, text="Max Product").grid(row=3, column=0)
max_prod = tk.Entry(param_frame, width=5)
max_prod.insert(0, "2000")
max_prod.grid(row=3, column=1)

flank_label = tk.Label(param_frame, text="Flank Size (bp)")
flank_label.grid(row=3, column=2)
flank_size_entry = tk.Entry(param_frame, width=5)
flank_size_entry.insert(0, "200")
flank_size_entry.grid(row=3, column=3)

check_specificity = tk.BooleanVar()
check_specificity.set(True)
tk.Checkbutton(param_frame, text="Test primer specificity", 
               variable=check_specificity, font=("Arial", 10, "bold")).grid(row=3, column=4, columnspan=2)

# Run button
tk.Button(root, text="🚀 Design Primers", command=run_pipeline, 
          bg="lightgreen", font=("Arial", 14, "bold"), height=2).grid(row=4, column=0, columnspan=6, pady=15)

# Results frame
results_frame = ttk.LabelFrame(root, text="Results", padding="10")
results_frame.grid(row=5, column=0, columnspan=6, sticky="ew", padx=10, pady=5)

results_text = scrolledtext.ScrolledText(results_frame, width=90, height=15)
results_text.pack(fill=tk.BOTH, expand=True)

# Initialize GUI state
toggle_input_mode()

# Configure column weights for resizing
for i in range(6):
    root.columnconfigure(i, weight=1)

root.mainloop()
