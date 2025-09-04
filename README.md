# PCR Primer Design with Specificity Testing

A Python-based graphical user interface (GUI) tool for designing PCR primers from a reference genome and gene annotations, with automated specificity testing via in-silico PCR simulation.

---

## Features

- Load reference genome in FASTA format.
- Parse gene coordinates from GFF3 annotation files.
- Extract gene sequences with customizable upstream/downstream flanking regions.
- Design primers for each gene using Primer3 with flexible parameters:
  - Primer size (min, optimal, max)
  - Melting temperature (Tm) range
  - GC content range
  - PCR product size range
- Perform primer specificity testing by searching for potential amplicons within the genome.
- Output primer information and specificity results to CSV file.
- User-friendly Tkinter GUI for loading files, setting parameters, and running the pipeline.
- Real-time progress and status display within the GUI.

---

## Installation

Ensure Python 3 is installed. Recommended to use a virtual environment.

Install required Python packages via pip:

pip install biopython primer3-py


---

## Usage

1. **Run GUI**  
   Launch the GUI application (`enhanced_primer_gui.py`).

2. **Input Files**  
   - Select the reference genome FASTA file.
   - Select the corresponding GFF3 gene annotation file.
   - Provide the output CSV filename.

3. **Set Primer Design Parameters**  
   Adjust primer size, melting temperature, GC content, product size, flank size, etc., according to your experimental needs.

4. **Enable Specificity Testing**  
   Check or uncheck the option to test primer specificity against the genome.

5. **Start Pipeline**  
   Click the "Run" button to initiate sequence extraction, primer design, and specificity testing.

6. **View Results**  
   Monitor progress in the results window.  
   After completion, output files will be saved:
   - CSV file containing primers and specificity information.
   - FASTA file of extracted gene sequences with flanks.

---

## Output

- **CSV File** containing:
  - Gene name
  - Forward and reverse primer sequences
  - Calculated melting temperatures (Tm)
  - GC content percentages
  - Expected PCR product size
  - Specificity test status (e.g., specific, non-specific, error)
  - Status indicator for primer design success

- **FASTA File** with extracted gene sequences including user-defined flanking regions.

---

## How It Works

- The program parses the genome and annotation files.
- For each gene, it extracts the sequence and flanking regions.
- Primer3 software is used to generate primers with user-specified constraints.
- Each primer pair is checked for genome-wide specificity by locating all potential binding sites and predicted amplicons via simple in-silico PCR search.
- Results are saved and displayed for user inspection.

---

## Dependencies

- [Biopython](https://biopython.org/) for sequence parsing and manipulation.
- [primer3-py](https://pypi.org/project/primer3-py/) for primer design functionalities.
- Python standard libraries: `tkinter`, `csv`, `os`, `subprocess`, `tempfile`.

---

## Example Parameter Settings

- Min/Opt/Max primer size: 18 / 20 / 25 bases
- Min/Opt/Max Tm: 58 / 60 / 60 °C
- GC content range: 20% - 50%
- PCR product size range: 35 - 5000 bp
- Flanking sequence size: 100 bp

---

## Notes and Limitations

- The in-silico specificity check uses exact sequence matching and simple reverse complement search without allowing mismatches.
- Experimental validation is recommended to confirm primer performance.
- Specificity testing may be slow for very large genomes or high numbers of primers.

---

## How to Extend

- Add mismatch tolerance in specificity search.
- Integrate with NCBI Primer-BLAST for enhanced specificity validation.
- Export additional metrics like self-dimers, hairpins, and secondary structure scores.
- Support batch processing from command line.

---

## Contact

For questions or support, please contact the developer.

---
