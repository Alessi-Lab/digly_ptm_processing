# DIA-NN PTM Processing Script

This document explains a Python script for processing post-translational modification (PTM) data from DIA-NN proteomics experiments, specifically focusing on UniMod:121 modifications.

## Code Breakdown

```python
import pandas as pd
from uniprotparser.betaparser import UniprotSequence, UniprotParser
from curtainutils.diann import lambda_function_for_diann_ptm_single_site
from curtainutils.common import read_fasta
```

### 1. Imports
- `pandas`: For data manipulation and analysis
- `UniprotSequence`, `UniprotParser`: For parsing and standardizing UniProt accession numbers
- `lambda_function_for_diann_ptm_single_site`: Function to process PTM single-site data
- `read_fasta`: Function to read protein sequences from FASTA format files

### 2. Loading Protein Sequence Database

```python
fasta_df = read_fasta(r"C:\Users\Toan Phung\Downloads\Uniprot_Human_032021.fasta")
```

- Reads a FASTA format file containing human protein sequences
- Converts it to a DataFrame for efficient lookup of protein sequences

### 3. Loading DIA-NN Data

```python
df = pd.read_csv(r"C:\Users\Toan Phung\Downloads\Supplementary Data 3 (1).R3.txt", sep="\t")
```

- Loads DIA-NN output data from a tab-delimited file
- Contains peptide sequences, modifications, and protein identifiers

### 4. Standardizing Protein Identifiers

```python
df["parse_id"] = df["Protein.Group"].apply(lambda x: str(UniprotSequence(x, parse_acc=True)) if UniprotSequence(x, parse_acc=True).accession else x)
```

- Creates a new column `parse_id` with standardized UniProt accession numbers
- Uses `UniprotSequence` to properly parse complex accession entries
- Preserves original identifiers if they can't be parsed as UniProt accessions

### 5. Processing PTM Data

```python
cf = df.apply(lambda x: lambda_function_for_diann_ptm_single_site(x, "Modified.Sequence", "parse_id", fasta_df, "UniMod:121"), axis=1)
```

- Applies the specialized function to each row in the dataframe
- Parameters:
  - `x`: Each row of the dataframe
  - `"Modified.Sequence"`: Column containing modified peptide sequences
  - `"parse_id"`: Column containing standardized protein identifiers
  - `fasta_df`: FASTA database for protein sequences
  - `"UniMod:121"`: Specific modification to process (cysteine carbamidomethylation)
- The function maps modifications to their positions in proteins, calculates sequence windows around modifications.

### 6. Filtering Relevant Modifications

```python
zf = cf[cf["Modified.Sequence"].str.contains("UniMod:121")]
```

- Filters rows to keep only those containing the UniMod:121 modification
- Ensures only peptides with the relevant modification are included in the final output

### 7. Saving Processed Data

```python
zf.to_csv(r"C:\Users\Toan Phung\Downloads\Supplementary Data 3 (1).R3.parsed.txt")
```

- Exports the processed data to a tabulated values file.
