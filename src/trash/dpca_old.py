# DETECTION OF POTENTIALLY CONTAMINATED AREAS

import pandas as pd
import numpy as np
from tqdm import tqdm
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

    
def compute_cons_nt_prop(row):
    
    if row["Clean_depth"] == 0:
        return 1
    else:
        # print(row["Cons_nt"], row["total_reads"])
        cons_nt = row["Cons_nt"]
        if cons_nt == "-":
            return row['prop_D']
    
        elif cons_nt not in ["A", "C", "G", "T"]:
            return 1
        
        prop = "prop_" + cons_nt
        return row[prop]
    
def prepare_df(qc_df):
    qc_df.iloc[:, 8:] = qc_df.iloc[:, 8:].replace(".", "0")
    
    for col in qc_df.columns:
        if qc_df[col].dtype == "object":
            try:
                qc_df[col] = qc_df[col].astype(int)
            except ValueError:
                pass
    
    qc_df["total_A"] = qc_df["A"] + qc_df["a"]
    qc_df["total_C"] = qc_df["C"] + qc_df["c"]
    qc_df["total_G"] = qc_df["G"] + qc_df["g"]
    qc_df["total_T"] = qc_df["T"] + qc_df["t"]
    qc_df["total_I"] = qc_df["I"] + qc_df["i"]
    qc_df["total_D"] = qc_df["D"] + qc_df["d"]
    
    # Calculate the proportion of all the nucleotides
    qc_df["prop_A"] = qc_df["total_A"] / qc_df["Clean_depth"]
    qc_df["prop_C"] = qc_df["total_C"] / qc_df["Clean_depth"]
    qc_df["prop_G"] = qc_df["total_G"] / qc_df["Clean_depth"]
    qc_df["prop_T"] = qc_df["total_T"] / qc_df["Clean_depth"]
    qc_df["prop_I"] = qc_df["total_I"] / qc_df["Clean_depth"]
    qc_df["prop_D"] = qc_df["total_D"] / qc_df["Clean_depth"]
    
    qc_df["c_a_prop"] = qc_df.apply(compute_cons_nt_prop, axis=1)
    
    return qc_df

def find_potential_contamination(qc_df, het_thr=0.9, depth_thr=0.2, prop_under_depth_thr=0.4, typical_depth_thr=500):
    """
    Find potential contamination in a DataFrame based on the proportion of the consensus nucleotide.

    Parameters:
        df (pd.DataFrame): Input DataFrame containing sequencing data.
        het_threshold (float): Threshold for heterozygosity. Default is 0.9.
        depth_thr (float): Threshold for depth amplicon dropout. Default is 0.5.

    Returns:
        pd.DataFrame: DataFrame with rows that have potential contamination.
    """
    
    qc_df = prepare_df(qc_df)

    typical_depth = qc_df["Clean_depth"].median()
    depth_int_thr = typical_depth * depth_thr
    prop_pos_under_depth_thr = qc_df[qc_df["Clean_depth"] < depth_int_thr].shape[0] / qc_df.shape[0]

    # Filter rows where the proportion is below the threshold
    contaminated_rows = qc_df[
        (qc_df["c_a_prop"] < het_thr) &
        (qc_df["Clean_depth"] < depth_int_thr) &
        (qc_df["Clean_depth"] > 0) # We don't want to include rows with 0 depth
    ]
    
    # Measure the proportion of positions below the depth threshold    
    # If this proportion is too high, we can consider that the sequencing quality is too poor to draw any conclusions
    
    if prop_pos_under_depth_thr > prop_under_depth_thr or typical_depth < typical_depth_thr:
        # Empty the contaminated_rows DataFrame
        contaminated_rows = pd.DataFrame(columns=qc_df.columns)

    return contaminated_rows, typical_depth, prop_pos_under_depth_thr

def generate_result_df(het_thr, depth_thr, pud_thr, med_dep_thr, sample_list):

    result_list = []
    for tsv, filename in tqdm(sample_list):
        # Get the filename without the path
        n_samples = len(sample_list)
        # Read the file
        
        try:
            df = pd.read_csv(
                tsv,       
                sep='\t',          
                compression='gzip',   
                low_memory=False,     
                on_bad_lines='warn' 
            )
        except:
            print(f"Error reading {tsv}. Skipping...")
            continue
        
            
        # Find potential contamination
        contaminated_rows, typical_depth, prop_pos_under_depth_thr = find_potential_contamination(df,
            het_thr=het_thr,
            depth_thr=depth_thr,
            prop_under_depth_thr=pud_thr,
            typical_depth_thr=med_dep_thr
        )
        
        
        if contaminated_rows.empty:
            continue  # No need to do anything else

        # Store the number of conflicting sites, their position in the genome and the name of the read
        mut_count = len(contaminated_rows)
        positions = contaminated_rows["Ref_pos"].tolist()
        
        # We also want to build the minor allele sequence 
        # Look at the contaminated rows and get the second most frequent nucleotide in columns total_A, total_C, total_G, total_T
        counts = contaminated_rows[['total_A', 'total_C', 'total_T', 'total_G']].copy()
        cons_nuc_idx = contaminated_rows['Cons_nt'].apply(lambda x: f'total_{x}')
        # Set the most frequent nucleotide to -1
        for idx, col in cons_nuc_idx.items():
            counts.at[idx, col] = -1
        m_a_variants = counts.idxmax(axis=1).str.replace('total_', '')
        
        # Build the sequence by taking the consensus sequence (df['Cons_nt']) and replacing the positions of the minor allele
        cons_sequence = df['Cons_nt'].copy()
        cons_sequence_str = "".join(cons_sequence.tolist())
        
        # We have to make -1 to all the positions because they are 1 - 29903 but the sequence is 0 - 29902
        python_positions = [pos - 1 for pos in positions]
        cons_sequence.loc[python_positions] = m_a_variants.values
        # m_a_sequence_str = "".join(cons_sequence)
        
        # Also output the sequence with N at the positions of the minor allele
        masked_sequence = df['Cons_nt'].copy()
        masked_sequence.loc[python_positions] = "N"
        masked_sequence_str = "".join(masked_sequence.tolist())
        
        # Build the list of the potentially contaminated sites in contaminated rows by concatenating cons_nt, position and minor allele
        pot_cont_sites = []
        for i, row in contaminated_rows.iterrows():
            cons_nt = row['Cons_nt']
            pos = row['Ref_pos']
            m_a = m_a_variants[i]
            pot_cont_sites.append(f"{cons_nt}{pos}{m_a}")
        
        # Save these as a tuple in a list
        result_list.append((filename, mut_count, pot_cont_sites, int(typical_depth), float(prop_pos_under_depth_thr), masked_sequence_str, cons_sequence_str))

    result_list_df = pd.DataFrame(result_list, columns=["read_name", "het_count", "het_mutations", "typical_depth", "prop_pos_under_thr", "masked_sequence_str", "cons_sequence_str"])

    filename = "result_samples_" + str(n_samples) + "_" + str(het_thr) + "_" + str(depth_thr) + "_" + str(pud_thr) + "_" + str(med_dep_thr) + ".csv"
    
    result_list_df.to_csv("/nfs/research/goldman/anoufa/data/" + filename, index=False)
    print(f"Results saved to {filename}")
    
    return result_list_df
    
def save_sequences_to_fasta(result_list_df):
    records = []
    for filename, _, _, _, _, m_a_sequence_str, cons_sequence_str in result_list_df.itertuples(index=False):
        # Create a SeqRecord for the minor allele sequence
        minor_seq_record = SeqRecord(Seq(m_a_sequence_str), id=f"{filename}_masked", description="Masked sequence")
        records.append(minor_seq_record)
        
        # Create a SeqRecord for the consensus sequence
        cons_seq_record = SeqRecord(Seq(cons_sequence_str), id=f"{filename}_consensus", description="Consensus sequence")
        records.append(cons_seq_record)
    
    # Save all records to a fasta file
    SeqIO.write(records, "/nfs/research/goldman/anoufa/data/masked_consensus_sequences.fasta", "fasta")
    

if __name__ == "__main__":
    # Example usage
    het_thr = 1.1 # if het_thr > 1.0, it is not used
    depth_thr = 0.1 # % of the median depth under which dropout is considered
    prop_under_depth_thr = 0.2 # % of the positions under the depth threshold (sample removed if too high)
    typical_depth_thr = 900 # median depth threshold (sample removed if too low)
    
    # Generate the result DataFrame
    result_df = generate_result_df(het_thr, depth_thr, prop_under_depth_thr, typical_depth_thr, n_samples=100)
    # Save the sequences to a fasta file
    save_sequences_to_fasta(result_df)