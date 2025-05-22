# diGly PTM Data Processing Script

This document explains a Python script for processing proteomics data from MSFragger to extract post-translational modifications (PTMs), particularly GG (diglycine) modifications.

## Code Breakdown

```python
import pandas as pd
from uniprotparser.betaparser import UniprotSequence, UniprotParser
from curtainutils.diann import lambda_function_for_diann_ptm_single_site
from curtainutils.common import read_fasta
from sequal.sequence import Sequence
from sequal.modification import Modification
import re
```

### 1. Imports
- `pandas`: For data manipulation and analysis
- `UniprotSequence`, `UniprotParser`: For parsing UniProt accession numbers
- `lambda_function_for_diann_ptm_single_site`: Function to process PTM data
- `read_fasta`: Function to read FASTA protein sequence files
- `Sequence`, `Modification`: For handling protein sequences and modifications
- `re`: For regular expression pattern matching

### 2. Data Loading

```python
fasta_df = read_fasta(r"C:\Users\Toan Phung\Downloads\Uniprot_Human_032021.fasta")
df = pd.read_csv(r"C:\Users\Toan Phung\Downloads\HEK293_MG132_6H_urea_DDA_3.txt", sep="\t")
```

- Loads a FASTA file containing protein sequences from UniProt
- Loads MSFragger output data from a tab-delimited file

### 3. Parsing UniProt Accessions

```python
df["parse_id"] = df["Master Protein Accessions"].apply(lambda x: str(UniprotSequence(x, parse_acc=True)) if UniprotSequence(x, parse_acc=True).accession else x)
```

- Creates a new column `parse_id` containing standardized UniProt accession numbers
- Uses `UniprotSequence` to properly parse complex accession entries

### 4. Setting Up Regex Patterns

```python
pattern = re.compile(r"(\d+)]$")
pattern_residue = re.compile(r"\[(\w+)\]")
```

- `pattern`: Matches position numbers at the end of modification strings
- `pattern_residue`: Matches amino acid residues enclosed in square brackets

### 5. Processing Modifications

```python
for i,r in df.iterrows():
    splitted = r["Modifications"].split(";")
    for j in splitted:
        if "GG" in j:
            match = re.search(pattern, j.strip())
            seq = Sequence(r["Sequence"])
            if match:
                mod = Modification("GG", int(match.group(1))-1)
                seq.add_modifications({int(match.group(1))-1: [mod]})
                df.at[i, "Modified.Sequence"] = seq.to_proforma()
                break
            else:
                pattern_residue_match = re.search(pattern_residue, j.strip())
                if pattern_residue_match:
                    index = r["Sequence"].index(pattern_residue_match.group(1))
                    mod = Modification("GG", index)
                    seq.add_modifications({index: [mod]})
                    df.at[i, "Modified.Sequence"] = seq.to_proforma()
```

- Iterates through each row in the dataframe
- Splits the "Modifications" column by semicolons to handle multiple modifications
- Filters for modifications containing "GG" (diglycine)
- Creates a `Sequence` object from the peptide sequence
- Handles two modification formats:
  1. When position is directly provided (e.g., "GG[123]")
  2. When amino acid is provided (e.g., "GG[K]")
- Adds the modification to the sequence object
- Converts the modified sequence to ProForma notation and stores it

### 6. Final Processing

```python
df["Protein.Group"] = df["Master Protein Accessions"]
df = df[pd.notnull(df["Modified.Sequence"])]
cf = df.apply(lambda x: lambda_function_for_diann_ptm_single_site(x, "Modified.Sequence", "parse_id", fasta_df, "GG"), axis=1)
cf.to_csv(r"C:\Users\Toan Phung\Downloads\HEK293_MG132_6H_urea_DDA_3.parsed.txt", sep="\t", index=False)
```

- Adds a "Protein.Group" column (required by downstream functions)
- Filters out rows without a modified sequence
- Applies the `lambda_function_for_diann_ptm_single_site` to each row to:
  - Extract site-specific information
  - Calculate positions in protein
  - Generate sequence windows
  - Format data for CurtainPTM
- Saves the processed data to a tab-delimited file
