# PCR Primer Designer (GUI)

**Automated PCR primer design for annotated genomes** using BioPython, Primer3, and Tkinter.  
Extracts gene and flanking regions from a FASTA/GFF3. Designs primers within genes, with full reporting of melting temperature, GC content, and design parameters. Results are exported to detailed CSV and FASTA files.

---

## Features

- **Graphical User Interface (GUI):**
  - Select input FASTA and GFF3 files
  - Set all key primer design parameters (size, Tm, GC%, flank size, product range)

- **Sequence Extraction:**
  - Saves 5′ flank, 3′ flank, and internal gene sequence for each locus to FASTA
  
- **Automated Primer Design:**
  - Designs primers from internal gene regions using Primer3
  - Reports Tm and GC% for each primer
  - Records product size, design summary, and primer status per gene
  
- **Parameter Documentation:**
  - All selected design parameters are exported to the CSV file

- **Open and Extensible:**
  - Easily modifiable for different regions (flanks, whole gene, central region)
  - Compatible with most bacterial genome FASTA/GFF3 formats

---

## Installation

1. **Clone this repo:**
    ```
    git clone https://github.com/YOUR_USERNAME/pcr-primer-designer.git
    ```

2. **Install dependencies:**
    ```
    pip install biopython primer3-py
    ```

3. Launch the Python GUI:
    ```
    python Master_withparameters.py
    ```

---

## Usage

1. **Select your genome FASTA and GFF3 files**
2. **Specify primer design constraints:**  
   - Primer length, melting temperature (Tm), GC %, product size, flank size, etc.
3. Click **Run**.
4. **Outputs generated:**
    - `primers.csv`: All designed primers per gene, with Tm and GC%.
    - `*_extracted_sequences.fasta`: FASTA file with gene and flanking regions.

---

## Example Outputs

- **primers.csv**
    | Gene Name | Forward Primer (internal) | Fwd Tm | Fwd GC% | Reverse Primer (internal) | Rev Tm | Rev GC% | Amplified Product Length | Status |
    |-----------|--------------------------|--------|---------|--------------------------|--------|---------|-------------------------|--------|
    | opgH     | TACAGCCGGTTTGCACTTCT     | 59.89  | 50.0    | ACGATCGCGATTCAGCTTCT     | 59.9   | 50.0    | 133 | OK  |

- **primer_extracted_sequences.fasta**
    ```
    ; ---------- gene ----------
    >BW25113_1049_gene | gene=opgH
    ATGAATAAGACAACTGAGTACATTGACGCAATGCCCATCGCCGCAAGCGAGAAAGCGGCA...
    ```

---

## Parameter Documentation

All design parameters used are saved in the first row of the output CSV for transparency and reproducibility:


