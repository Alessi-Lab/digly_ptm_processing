import pandas as pd
from uniprotparser.betaparser import UniprotSequence, UniprotParser
from curtainutils.diann import lambda_function_for_diann_ptm_single_site
from curtainutils.common import read_fasta
from sequal.sequence import Sequence
from sequal.modification import Modification
import re

#### DDA #######
fasta_df = read_fasta(r"C:\Users\Toan Phung\Downloads\Uniprot_Human_032021.fasta")
df = pd.read_csv(r"C:\Users\Toan Phung\Downloads\HEK293_MG132_6H_urea_DDA_3.txt", sep="\t")
df["parse_id"] = df["Master Protein Accessions"].apply(lambda x: str(UniprotSequence(x, parse_acc=True)) if UniprotSequence(x, parse_acc=True).accession else x)

pattern = re.compile(r"(\d+)]$")
pattern_residue = re.compile(r"\[(\w+)\]")
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
df["Protein.Group"] = df["Master Protein Accessions"]
df = df[pd.notnull(df["Modified.Sequence"])]
cf = df.apply(lambda x: lambda_function_for_diann_ptm_single_site(x, "Modified.Sequence", "parse_id", fasta_df,
                                                                      "GG"), axis=1)
cf.to_csv(r"C:\Users\Toan Phung\Downloads\HEK293_MG132_6H_urea_DDA_3.parsed.txt", sep="\t", index=False)


####### DIANN PTM SINGLE SITE #######

import pandas as pd
from uniprotparser.betaparser import UniprotSequence, UniprotParser
from curtainutils.diann import lambda_function_for_diann_ptm_single_site
from curtainutils.common import read_fasta

fasta_df = read_fasta(r"C:\Users\Toan Phung\Downloads\Uniprot_Human_032021.fasta")
df = pd.read_csv(r"C:\Users\Toan Phung\Downloads\Supplementary Data 3 (1).R3.txt", sep="\t")
df["parse_id"] = df["Protein.Group"].apply(lambda x: str(UniprotSequence(x, parse_acc=True)) if UniprotSequence(x, parse_acc=True).accession else x)
cf = df.apply(lambda x: lambda_function_for_diann_ptm_single_site(x, "Modified.Sequence", "parse_id", fasta_df,
                                                                      "UniMod:121"), axis=1)
zf = cf[cf["Modified.Sequence"].str.contains("UniMod:121")]
zf.to_csv(r"C:\Users\Toan Phung\Downloads\Supplementary Data 3 (1).R3.parsed.txt")
