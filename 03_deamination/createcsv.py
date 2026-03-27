#!/usr/bin/env python3

import pandas as pd
from pathlib import Path

col_substitutions = ['G>A','C>T', 'A>G', 'T>C', 'A>C', 'A>T', 'C>G', 'C>A', 'T>G', 'T>A', 'G>C', 'G>T', 'A>-', 'T>-', 'C>-', 'G>-', '->A', '->T', '->C', '->G', 'S']

def analyze():
    data_dir = Path("01_data")
    results = []
    
    sample_dirs = [d for d in data_dir.iterdir() if d.is_dir()]
    
    for sample_dir in sample_dirs:
        sample_name = sample_dir.name
        misincorp_file = sample_dir / "misincorporation.txt"
        
        if not misincorp_file.exists():
            continue
            
        try:
            df = pd.read_csv(misincorp_file, sep='\t')
            
            mask = (df['Pos'] <= 25) & (df['End'].isin(['3p', '5p']))
            df_filtered = df[mask].copy()
            
            if df_filtered.empty:
                continue
            
            for end in ['3p', '5p']:
                end_data = df_filtered[df_filtered['End'] == end]
                
                if end_data.empty:
                    continue
                
                for pos in range(1, 26):
                    pos_data = end_data[end_data['Pos'] == pos]
                    
                    if len(pos_data) != 2:
                        continue
                    
                    pos_strand = pos_data[pos_data['Std'] == '+'].iloc[0] if len(pos_data[pos_data['Std'] == '+']) > 0 else None
                    neg_strand = pos_data[pos_data['Std'] == '-'].iloc[0] if len(pos_data[pos_data['Std'] == '-']) > 0 else None
                    
                    if pos_strand is None or neg_strand is None:
                        continue
                    
                    for sub_type in col_substitutions:
                        pos_count       = pos_strand[sub_type]
                        neg_count       = neg_strand[sub_type]
                        total_pos_count = pos_strand['Total']
                        total_neg_count = neg_strand['Total']
                        combined_count  = pos_count + neg_count
                        combined_total  = total_pos_count + total_neg_count
                        
                        if combined_total == 0:
                            continue
                        
                        frequency = combined_count / combined_total
                        
                        results.append({
                            'Sample'         : sample_name,
                            'End'            : end,
                            'Position'       : pos,
                            'Substitution'   : sub_type,
                            'Pos_Count'      : pos_count,
                            'Neg_Count'      : neg_count,
                            'Combined_Count' : combined_count,
                            'Pos_Total'      : total_pos_count,
                            'Neg_Total'      : total_neg_count,
                            'Combined_Total' : combined_total,
                            'Frequency'      : frequency
                        })
                    
        except Exception as e:
            continue
    
    return results

def save_results(results):
    if not results:
        return
    
    csv_df = pd.DataFrame(results)
    csv_filename = "03_results/frequencies.csv"
    csv_df.to_csv(csv_filename, index=False)
    
    max_results = []
    for sub_type in col_substitutions:
        sub_data = [r for r in results if r['Substitution'] == sub_type]
        if sub_data:
            maxes = max(sub_data, key=lambda x: x['Frequency'])
            max_results.append({
                'Substitution'  : maxes['Substitution'],
                'Max_Frequency' : maxes['Frequency'],
                'Sample'        : maxes['Sample'],
                'End'           : maxes['End'],
                'Position'      : maxes['Position'],
                'Combined_Count': maxes['Combined_Count'],
                'Combined_Total': maxes['Combined_Total'],
                'Pos_Count'     : maxes['Pos_Count'],
                'Neg_Count'     : maxes['Neg_Count'],
                'Pos_Total'     : maxes['Pos_Total'],
                'Neg_Total'     : maxes['Neg_Total']
            })
    
    max_df = pd.DataFrame(max_results)
    max_csv_filename = "03_results/max_frequencies.csv"
    max_df.to_csv(max_csv_filename, index=False)

def main():
    results = analyze()
    save_results(results)

if __name__ == "__main__":
    main() 