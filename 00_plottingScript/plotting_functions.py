import pandas as pd
import os
import matplotlib.pyplot as plt
import numpy as np
import re
import multiprocessing as mp
from plotting_config import *

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
ROOT_DIR = os.path.dirname(os.path.dirname(SCRIPT_DIR))

def get_sample_type(sample):
    return SAMPLE_TYPES.get(sample, 'unknown')

def get_real_sample_name(sample):
    return SAMPLE_NAMES.get(sample, sample)

def extract_busco_data(file_path):
    try:
        with open(file_path, 'r') as f:
            content = f.read()
        
        results_line = None
        for line in content.split('\n'):
            if line.strip().startswith('C:') and '%' in line:
                results_line = line.strip()
                break
        
        if not results_line:
            return None
        
        c_match = re.search(r'C:([\d.]+)%', results_line)
        s_match = re.search(r'S:([\d.]+)%', results_line)
        d_match = re.search(r'D:([\d.]+)%', results_line)
        f_match = re.search(r'F:([\d.]+)%', results_line)
        m_match = re.search(r'M:([\d.]+)%', results_line)
        
        if not all([c_match, s_match, d_match, f_match, m_match]):
            return None
        
        lines = content.split('\n')
        s_count = d_count = f_count = m_count = 0
        
        for line in lines:
            if 'Complete and single-copy BUSCOs (S)' in line:
                s_count = int(line.split()[0])
            elif 'Complete and duplicated BUSCOs (D)' in line:
                d_count = int(line.split()[0])
            elif 'Fragmented BUSCOs (F)' in line:
                f_count = int(line.split()[0])
            elif 'Missing BUSCOs (M)' in line:
                m_count = int(line.split()[0])
        
        return {
            'S': float(s_match.group(1)), 'D': float(d_match.group(1)),
            'F': float(f_match.group(1)), 'M': float(m_match.group(1)),
            'S_count': s_count, 'D_count': d_count,
            'F_count': f_count, 'M_count': m_count
        }
    except Exception as e:
        print(f"Error extracting BUSCO data from {file_path}: {e}")
        return None

def get_deamination_plot_data(sample_file):    
    def get_data(end, std, col):
        try:
            df = pd.read_csv(sample_file, sep='\t')
            filtered_df = df[(df['End'] == end) & (df['Std'] == std)]
            return filtered_df[col].head(25).tolist()
        except Exception as e:
            print(f"Error reading deamination data: {e}")
            return [0] * 25
    
    pos_3p_CtoT = get_data('3p', '+', 'C>T')
    pos_3p_AtoG = get_data('3p', '+', 'A>G')
    pos_3p_Total = get_data('3p', '+', 'Total')
    neg_3p_CtoT = get_data('3p', '-', 'C>T')
    neg_3p_AtoG = get_data('3p', '-', 'A>G')
    neg_3p_Total = get_data('3p', '-', 'Total')

    total_3p_CtoT_freq = [(x + y)/(z + w) if (z + w) > 0 else 0 for x,y,z,w in zip(pos_3p_CtoT, neg_3p_CtoT, pos_3p_Total, neg_3p_Total)]
    total_3p_AtoG_freq = [(x + y)/(z + w) if (z + w) > 0 else 0 for x,y,z,w in zip(pos_3p_AtoG, neg_3p_AtoG, pos_3p_Total, neg_3p_Total)]

    total_3p_CtoT_freq.reverse()
    total_3p_AtoG_freq.reverse()

    pos_5p_CtoT = get_data('5p', '+', 'C>T')
    pos_5p_AtoG = get_data('5p', '+', 'A>G')
    pos_5p_Total = get_data('5p', '+', 'Total')
    neg_5p_CtoT = get_data('5p', '-', 'C>T')
    neg_5p_AtoG = get_data('5p', '-', 'A>G')
    neg_5p_Total = get_data('5p', '-', 'Total')

    total_5p_CtoT_freq = [(x + y)/(z + w) if (z + w) > 0 else 0 for x,y,z,w in zip(pos_5p_CtoT, neg_5p_CtoT, pos_5p_Total, neg_5p_Total)]
    total_5p_AtoG_freq = [(x + y)/(z + w) if (z + w) > 0 else 0 for x,y,z,w in zip(pos_5p_AtoG, neg_5p_AtoG, pos_5p_Total, neg_5p_Total)]

    return {
        '5p_CtoT': total_5p_CtoT_freq, '5p_AtoG': total_5p_AtoG_freq,
        '3p_CtoT': total_3p_CtoT_freq, '3p_AtoG': total_3p_AtoG_freq
    }

def plot_busco_representative(output_dir):
    print("Creating representative BUSCO plots...")
    
    busco_dir = os.path.join(output_dir, 'busco')
    os.makedirs(busco_dir, exist_ok=True)
    
    fig, axes = plt.subplots(1, 3, figsize=(24, 8))
    
    for i, sample in enumerate(REPRESENTATIVE_SAMPLES):
        ax = axes[i]
        
        if sample.startswith('WA'):
            file_path = f'00_busco/01_data/{sample}_paired_busco/short_summary.specific.viridiplantae_odb10.{sample}_paired_busco.txt'
        else:
            file_path = f'00_busco/01_data/{sample}_busco/short_summary.specific.viridiplantae_odb10.{sample}_busco.txt'
        
        try:
            busco_data = extract_busco_data(file_path)
            if busco_data is None:
                continue
                
            categories = ['S', 'D', 'F', 'M']
            percentages = [busco_data['S'], busco_data['D'], busco_data['F'], busco_data['M']]
            colors = [COLORS['busco'][cat] for cat in categories]
            
            bars = ax.bar(categories, percentages, color=colors, alpha=0.8)
            
            for bar, percentage in zip(bars, percentages):
                height = bar.get_height()
                ax.text(bar.get_x() + bar.get_width()/2., height + 0.5,
                        f'{percentage:.1f}%', ha='center', va='bottom', fontsize=18, fontweight='bold')
            
            category_labels = ['Complete\nSingle-copy', 'Complete\nDuplicated', 'Fragmented', 'Missing']
            for bar, label in zip(bars, category_labels):
                ax.text(bar.get_x() + bar.get_width()/2., -5, label, 
                        ha='center', va='top', fontsize=22, rotation=90, color='black')
            
            ax.set_title(f'{get_real_sample_name(sample)}', fontsize=18, fontweight='bold')
            ax.set_ylim(0, 100)
            ax.tick_params(axis='x', labelsize=22)
            ax.tick_params(axis='y', labelsize=22)
            ax.set_ylabel('Percentage (%)', fontsize=22, labelpad=10)
            ax.tick_params(axis='x', which='both', bottom=True, top=False, labelbottom=False)
            ax.grid(True, alpha=0.3, axis='y')
            ax.set_yticks([0, 20, 40, 60, 80, 100])
            
        except Exception as e:
            print(f"Error processing {file_path}: {e}")
            ax.text(0.5, 0.5, f'Error loading data\n{get_real_sample_name(sample)}', 
                    ha='center', va='center', transform=ax.transAxes, fontsize=12)
    
    plt.tight_layout(pad=1.0)
    output_file = os.path.join(busco_dir, 'busco_representative.png')
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Representative BUSCO plot saved to: {output_file}")

def plot_deamination_representative(output_dir):
    print("Creating representative deamination plots...")
    
    deamination_dir = os.path.join(output_dir, 'deamination')
    os.makedirs(deamination_dir, exist_ok=True)
    
    fig, axes = plt.subplots(1, 3, figsize=(24, 8))
    
    for i, sample in enumerate(REPRESENTATIVE_SAMPLES):
        ax = axes[i]
        file_path = os.path.join(ROOT_DIR, 'FINALdata', 'deamination_data', sample, 'misincorporation.txt')
        
        try:
            plot_data = get_deamination_plot_data(file_path)
            
            x_5p = list(range(1, 26))
            ax.plot(x_5p, plot_data['5p_AtoG'], color='#56B4E9', label='A>I', linewidth=2, alpha=0.8)
            ax.plot(x_5p, plot_data['5p_CtoT'], color='#F04442', label='C>U', linewidth=2, alpha=0.8)
            
            x_3p = list(range(-25, 0))
            ax.plot(x_3p, plot_data['3p_AtoG'], color='#56B4E9', linewidth=2, alpha=0.8)
            ax.plot(x_3p, plot_data['3p_CtoT'], color='#F04442', linewidth=2, alpha=0.8)
            
            ax.set_title(f'{get_real_sample_name(sample)}', fontsize=22, fontweight='bold')
            ax.set_ylim(0.00, 0.002)
            ax.set_yticks([0.0000, 0.0004, 0.0008, 0.0012, 0.0016, 0.0020])
            ax.set_xlim(-25.5, 25.5)
            
            tick_positions = [-25, -20, -15, -10, -5, 0, 5, 10, 15, 20, 25]
            tick_labels = ['25', '20', '15', '10', '5', '0', '-5', '-10', '-15', '-20', '-25']
            ax.set_xticks(tick_positions)
            ax.set_xticklabels(tick_labels)
            ax.tick_params(axis='x', labelsize=22)
            ax.tick_params(axis='y', labelsize=22)
            ax.set_xlabel('Position', fontsize=22, labelpad=10)
            if i == 0:
                ax.set_ylabel('Misincorporation Frequency', fontsize=22, labelpad=10)
                ax.legend(fontsize=22, loc='upper right')
            
            ax.grid(True, alpha=0.3)
            ax.axvline(x=0, color='gray', linestyle='--', alpha=0.5)
            
        except Exception as e:
            print(f"Error processing {file_path}: {e}")
            ax.text(0.5, 0.5, f'Error loading data\n{sample}', 
                    ha='center', va='center', transform=ax.transAxes, fontsize=18)
    
    plt.tight_layout()
    output_file = os.path.join(deamination_dir, 'deamination_representative.png')
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Representative deamination plot saved to: {output_file}")

def plot_inserts_representative(output_dir):
    print("Creating representative insert plots...")
    
    inserts_dir = os.path.join(output_dir, 'inserts')
    os.makedirs(inserts_dir, exist_ok=True)
    
    fig, axes = plt.subplots(1, 3, figsize=(24, 8))
    
    for i, sample in enumerate(REPRESENTATIVE_SAMPLES):
        ax = axes[i]
        
        if sample.startswith('WA'):
            file_path = os.path.join(ROOT_DIR, 'FINALdata', 'inserts_data', f'{sample}_paired.csv')
        else:
            file_path = os.path.join(ROOT_DIR, 'FINALdata', 'inserts_data', f'{sample}.csv')
        
        try:
            df = pd.read_csv(file_path)
            df_filtered = df[df['insert_length'] >= 1].copy()
            
            insert_length_counts = df_filtered['insert_length'].value_counts().sort_index()
            x = insert_length_counts.index.values
            y = insert_length_counts.values
            
            mask_pos = (x >= 1) & (x <= 800)
            
            if np.any(mask_pos):
                ax.plot(x[mask_pos], y[mask_pos], color='#000000', linewidth=3, label='insert')
            
            sample_name = f'{get_real_sample_name(sample)} ({get_sample_type(sample).title()})'
            ax.set_title(sample_name, fontsize=20, fontweight='bold')
            ax.set_xlim(1, 800)
            ax.set_ylim(-1, np.max(y[mask_pos]) + np.max(y[mask_pos]) * 0.1)
            
            ax.set_xticks([1] + list(np.arange(100, 800 + 1, 100)))
            ax.tick_params(axis='x', rotation=45, labelsize=22)
            ax.tick_params(axis='y', labelsize=22)
            
            ax.set_xlabel('Insert Length (bp)', fontsize=22, labelpad=10)
            if i == 0:
                ax.set_ylabel('Count', fontsize=22, labelpad=10)
            
            ax.grid(True, alpha=0.3)
            
        except Exception as e:
            print(f"Error processing {file_path}: {e}")
            ax.text(0.5, 0.5, f'Error loading data\n{get_real_sample_name(sample)}', 
                    ha='center', va='center', transform=ax.transAxes, fontsize=18)
    
    plt.tight_layout()
    output_file = os.path.join(inserts_dir, 'inserts_representative.png')
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Representative insert plot saved to: {output_file}")

def plot_transrate_representative(output_dir):
    print("Creating representative transrate plots...")
    
    transrate_dir = os.path.join(output_dir, 'transrate')
    os.makedirs(transrate_dir, exist_ok=True)
    
    fig, axes = plt.subplots(1, 3, figsize=(24, 8))
    
    for i, sample in enumerate(REPRESENTATIVE_SAMPLES):
        ax = axes[i]
        file_path = os.path.join(ROOT_DIR, 'FINALdata', 'transrate_data', 'tr2_assembly', f'{sample}.csv')
        
        try:
            df = pd.read_csv(file_path)
            
            score_columns = ['sCnuc_Harmonic', 'sCcov_Harmonic', 'sCord_Harmonic', 'sCseg_Harmonic']
            score_labels = ['sCnuc', 'sCcov', 'sCord', 'sCseg']
            
            scores = []
            for col in score_columns:
                if col in df.columns:
                    scores.append(df[col].iloc[0])
                else:
                    scores.append(0)
            
            x_pos = np.arange(len(score_labels))
            bars = ax.bar(x_pos, scores, color=[COLORS['transrate'][col] for col in score_columns], alpha=0.8)
            
            for bar, score in zip(bars, scores):
                height = bar.get_height()
                ax.text(bar.get_x() + bar.get_width()/2., height + 0.01,
                        f'{score:.3f}', ha='center', va='bottom', fontsize=18, fontweight='bold')
            
            sample_name = f'{get_real_sample_name(sample)} ({get_sample_type(sample).title()})'
            ax.set_title(sample_name, fontsize=20, fontweight='bold')
            ax.set_ylim(0, 1.1)
            
            ax.set_xticks(x_pos)
            ax.set_xticklabels(score_labels, fontsize=22)
            ax.tick_params(axis='y', labelsize=22)
            
            ax.set_ylabel('Score', fontsize=22, labelpad=10)
            ax.grid(True, alpha=0.3, axis='y')
            
        except Exception as e:
            print(f"Error processing {file_path}: {e}")
            ax.text(0.5, 0.5, f'Error loading data\n{get_real_sample_name(sample)}', 
                    ha='center', va='center', transform=ax.transAxes, fontsize=18)
    
    plt.tight_layout()
    output_file = os.path.join(transrate_dir, 'transrate_representative.png')
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Representative transrate plot saved to: {output_file}")

def plot_busco_individual_worker(args):
    sample, output_dir = args
    
    busco_dir = os.path.join(output_dir, 'busco')
    os.makedirs(busco_dir, exist_ok=True)
    
    if sample.startswith('WA'):
        file_path = f'00_busco/01_data/{sample}_paired_busco/short_summary.specific.viridiplantae_odb10.{sample}_paired_busco.txt'
    else:
        file_path = f'00_busco/01_data/{sample}_busco/short_summary.specific.viridiplantae_odb10.{sample}_busco.txt'
    
    if not os.path.exists(file_path):
        return f"Warning: {file_path} not found, skipping {sample}"
    
    try:
        busco_data = extract_busco_data(file_path)
        if busco_data is None:
            return f"Could not extract BUSCO data for {sample}"
        
        fig, ax = plt.subplots(1, 1, figsize=(8, 6))
        
        categories = ['Complete\nSingle-copy', 'Complete\nDuplicated', 'Fragmented', 'Missing']
        percentages = [busco_data['S'], busco_data['D'], busco_data['F'], busco_data['M']]
        colors = [COLORS['busco'][cat] for cat in ['S', 'D', 'F', 'M']]
        
        bars = ax.bar(categories, percentages, color=colors, alpha=0.8)
        
        for bar, percentage in zip(bars, percentages):
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height + 0.5,
                    f'{percentage:.1f}%', ha='center', va='bottom', fontsize=14, fontweight='bold')
        
        ax.set_title(f'{get_real_sample_name(sample)} BUSCO Scores', fontsize=16, fontweight='bold')
        ax.set_ylim(0, 100)
        ax.tick_params(axis='x', labelsize=12)
        ax.tick_params(axis='y', labelsize=12)
        ax.set_ylabel('Percentage (%)', fontsize=12, labelpad=10)
        ax.tick_params(axis='x', which='both', bottom=True, top=False, labelbottom=False)
        ax.grid(True, alpha=0.3, axis='y')
        ax.set_yticks([0, 20, 40, 60, 80, 100])
        
        plt.tight_layout()
        output_file = os.path.join(busco_dir, f'busco_{sample}.png')
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        return f"Successfully processed {get_real_sample_name(sample)} ({sample})"
        
    except Exception as e:
        return f"Error processing {sample}: {e}"

def plot_deamination_individual_worker(args):
    sample, output_dir = args
    
    deamination_dir = os.path.join(output_dir, 'deamination')
    os.makedirs(deamination_dir, exist_ok=True)
    
    file_path = os.path.join(ROOT_DIR, 'FINALdata', 'deamination_data', sample, 'misincorporation.txt')
    
    if not os.path.exists(file_path):
        return f"Warning: {file_path} not found, skipping {sample}"
    
    try:
        plot_data = get_deamination_plot_data(file_path)
        
        fig, ax = plt.subplots(1, 1, figsize=(10, 6))
        
        x_5p = list(range(1, 26))
        ax.plot(x_5p, plot_data['5p_AtoG'], color='#56B4E9', label='A>I', linewidth=2, alpha=0.8)
        ax.plot(x_5p, plot_data['5p_CtoT'], color='#F04442', label='C>U', linewidth=2, alpha=0.8)
        
        x_3p = list(range(-25, 0))
        ax.plot(x_3p, plot_data['3p_AtoG'], color='#56B4E9', linewidth=2, alpha=0.8)
        ax.plot(x_3p, plot_data['3p_CtoT'], color='#F04442', linewidth=2, alpha=0.8)
        
        ax.set_title(f'{get_real_sample_name(sample)} Deamination Pattern ({get_sample_type(sample).title()})', fontsize=14, fontweight='bold')
        ax.set_ylim(0.00, 0.002)
        ax.set_yticks([0.0000, 0.0004, 0.0008, 0.0012, 0.0016, 0.0020])
        ax.set_xlim(-25.5, 25.5)
        
        tick_positions = [-25, -20, -15, -10, -5, 0, 5, 10, 15, 20, 25]
        tick_labels = ['25', '20', '15', '10', '5', '0', '-5', '-10', '-15', '-20', '-25']
        ax.set_xticks(tick_positions)
        ax.set_xticklabels(tick_labels)
        ax.tick_params(axis='x', labelsize=10)
        ax.tick_params(axis='y', labelsize=10)
        ax.set_xlabel('Position', fontsize=12, labelpad=10)
        ax.set_ylabel('Misincorporation Frequency', fontsize=12, labelpad=10)
        ax.legend(fontsize=10, loc='upper right')
        
        ax.grid(True, alpha=0.3)
        ax.axvline(x=0, color='gray', linestyle='--', alpha=0.5)
        
        plt.tight_layout()
        output_file = os.path.join(deamination_dir, f'deamination_{sample}.png')
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        return f"Successfully processed {get_real_sample_name(sample)} ({sample})"
        
    except Exception as e:
        return f"Error processing {sample}: {e}"

def plot_inserts_individual_worker(args):
    sample, output_dir = args
    
    inserts_dir = os.path.join(output_dir, 'inserts')
    os.makedirs(inserts_dir, exist_ok=True)
    
    if sample.startswith('WA'):
        file_path = os.path.join(ROOT_DIR, 'FINALdata', 'inserts_data', f'{sample}_paired.csv')
    else:
        file_path = os.path.join(ROOT_DIR, 'FINALdata', 'inserts_data', f'{sample}.csv')
    
    if not os.path.exists(file_path):
        return f"Warning: {file_path} not found, skipping {sample}"
    
    try:
        df = pd.read_csv(file_path)
        df_filtered = df[df['insert_length'] >= 1].copy()
        
        fig, ax = plt.subplots(1, 1, figsize=(10, 6))
        
        insert_length_counts = df_filtered['insert_length'].value_counts().sort_index()
        x = insert_length_counts.index.values
        y = insert_length_counts.values
        
        mask_pos = (x >= 1) & (x <= 800)
        
        if np.any(mask_pos):
            ax.plot(x[mask_pos], y[mask_pos], color='#000000', linewidth=2, label='insert')
        
        ax.set_title(f'{get_real_sample_name(sample)} insert Length Distribution ({get_sample_type(sample).title()})', fontsize=14, fontweight='bold')
        ax.set_xlim(1, 800)
        ax.set_ylim(-1, np.max(y[mask_pos]) + np.max(y[mask_pos]) * 0.1)
        
        ax.set_xticks([1] + list(np.arange(100, 800 + 1, 100)))
        ax.tick_params(axis='x', rotation=45, labelsize=10)
        ax.tick_params(axis='y', labelsize=10)
        
        ax.set_xlabel('Insert Length (bp)', fontsize=12, labelpad=10)
        ax.set_ylabel('Count', fontsize=12, labelpad=10)
        
        ax.grid(True, alpha=0.3)
        
        plt.tight_layout()
        output_file = os.path.join(inserts_dir, f'inserts_{sample}.png')
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        return f"Successfully processed {get_real_sample_name(sample)} ({sample})"
        
    except Exception as e:
        return f"Error processing {sample}: {e}"

def plot_transrate_individual_worker(args):
    sample, output_dir = args
    
    transrate_dir = os.path.join(output_dir, 'transrate')
    os.makedirs(transrate_dir, exist_ok=True)
    
    file_path = os.path.join(ROOT_DIR, 'FINALdata', 'transrate_data', 'tr2_assembly', f'{sample}.csv')
    
    if not os.path.exists(file_path):
        return f"Warning: {file_path} not found, skipping {sample}"
    
    try:
        df = pd.read_csv(file_path)
        
        fig, ax = plt.subplots(1, 1, figsize=(8, 6))
        
        score_columns = ['sCnuc_Harmonic', 'sCcov_Harmonic', 'sCord_Harmonic', 'sCseg_Harmonic']
        score_labels = ['sCnuc', 'sCcov', 'sCord', 'sCseg']
        
        scores = []
        for col in score_columns:
            if col in df.columns:
                scores.append(df[col].iloc[0])
            else:
                scores.append(0)
        
        x_pos = np.arange(len(score_labels))
        bars = ax.bar(x_pos, scores, color=[COLORS['transrate'][col] for col in score_columns], alpha=0.8)
        
        for bar, score in zip(bars, scores):
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height + 0.01,
                    f'{score:.3f}', ha='center', va='bottom', fontsize=12, fontweight='bold')
        
        ax.set_title(f'{get_real_sample_name(sample)} Transrate2 Scores ({get_sample_type(sample).title()})', fontsize=14, fontweight='bold')
        ax.set_ylim(0, 1.1)
        
        ax.set_xticks(x_pos)
        ax.set_xticklabels(score_labels, fontsize=12)
        ax.tick_params(axis='y', labelsize=12)
        
        ax.set_ylabel('Score', fontsize=12, labelpad=10)
        ax.grid(True, alpha=0.3, axis='y')
        
        plt.tight_layout()
        output_file = os.path.join(transrate_dir, f'transrate_{sample}.png')
        plt.savefig(output_file, dpi=300, bbox_inches='tight')
        plt.close()
        
        return f"Successfully processed {get_real_sample_name(sample)} ({sample})"
        
    except Exception as e:
        return f"Error processing {sample}: {e}"

def plot_busco_individual(output_dir, num_processes=40):
    print("Creating individual BUSCO plots for all samples...")
    
    busco_dir = os.path.join(output_dir, 'busco')
    os.makedirs(busco_dir, exist_ok=True)
    
    args_list = [(sample, output_dir) for sample in SAMPLES.keys()]
    
    with mp.Pool(processes=num_processes) as pool:
        results = pool.map(plot_busco_individual_worker, args_list)
    
    for result in results:
        print(f"  {result}")

def plot_deamination_individual(output_dir, num_processes=40):
    print("Creating individual deamination plots for all samples...")
    
    deamination_dir = os.path.join(output_dir, 'deamination')
    os.makedirs(deamination_dir, exist_ok=True)
    
    args_list = [(sample, output_dir) for sample in SAMPLES.keys()]
    
    with mp.Pool(processes=num_processes) as pool:
        results = pool.map(plot_deamination_individual_worker, args_list)
    
    for result in results:
        print(f"  {result}")

def plot_inserts_individual(output_dir, num_processes=40):
    print("Creating individual insert plots for all samples...")
    
    inserts_dir = os.path.join(output_dir, 'inserts')
    os.makedirs(inserts_dir, exist_ok=True)
    
    args_list = [(sample, output_dir) for sample in SAMPLES.keys()]
    
    with mp.Pool(processes=num_processes) as pool:
        results = pool.map(plot_inserts_individual_worker, args_list)
    
    for result in results:
        print(f"  {result}")

def plot_transrate_individual(output_dir, num_processes=40):
    print("Creating individual transrate plots for all samples...")
    
    transrate_dir = os.path.join(output_dir, 'transrate')
    os.makedirs(transrate_dir, exist_ok=True)
    
    args_list = [(sample, output_dir) for sample in SAMPLES.keys()]
    
    with mp.Pool(processes=num_processes) as pool:
        results = pool.map(plot_transrate_individual_worker, args_list)
    
    for result in results:
        print(f"  {result}")

def create_concatenated_representative_plots(output_dir, plot_types):
    print("Creating concatenated representative plots image...")
    
    available_plots = []
    plot_files = {
        'deamination': os.path.join(output_dir, 'deamination', 'deamination_representative.png'),
        'inserts': os.path.join(output_dir, 'inserts', 'inserts_representative.png'),
        'transrate': os.path.join(output_dir, 'transrate', 'transrate_representative.png'),
        'busco': os.path.join(output_dir, 'busco', 'busco_representative.png')
    }
    
    for plot_type in plot_types:
        if plot_type in plot_files and os.path.exists(plot_files[plot_type]):
            available_plots.append(plot_type)
    
    if len(available_plots) < 2:
        print("  Not enough representative plots to concatenate")
        return
    
    # Reorder available_plots to ensure correct sequence: deamination, inserts, transrate, busco
    desired_order = ['deamination', 'inserts', 'transrate', 'busco']
    ordered_plots = []
    for plot_type in desired_order:
        if plot_type in available_plots:
            ordered_plots.append(plot_type)
    
    # Add any remaining plot types that weren't in the desired order
    for plot_type in available_plots:
        if plot_type not in ordered_plots:
            ordered_plots.append(plot_type)
    
    available_plots = ordered_plots
    
    fig, axes = plt.subplots(4, 3, figsize=(26, 36))
    
    # plot_labels = {
    #     'busco': 'BUSCO Scores',
    #     'deamination': 'Deamination Patterns', 
    #     'inserts': 'Insert Length Distributions',
    #     'transrate': 'Transrate2 Assembly Scores'
    # }111

    representative_samples = {
        'deamination': ['DAL192', 'WA13', 'WA22'],
        'inserts': ['DAL192', 'WA13', 'WA22'],
        'transrate': ['DAL192', 'WA13', 'WA22'],
        'busco': ['DAL192', 'WA13', 'WA22']
    }
    
    for row_idx, plot_type in enumerate(available_plots):
        if plot_type not in plot_files or not os.path.exists(plot_files[plot_type]):
            continue
        
        for col_idx, sample in enumerate(representative_samples[plot_type]):
            ax = axes[row_idx, col_idx]
            
            try:
                if plot_type == 'busco':
                    create_busco_individual_plot(ax, sample, output_dir)
                elif plot_type == 'deamination':
                    create_deamination_individual_plot(ax, sample, output_dir)
                elif plot_type == 'inserts':
                    create_inserts_individual_plot(ax, sample, output_dir)
                elif plot_type == 'transrate':
                    create_transrate_individual_plot(ax, sample, output_dir)
            except Exception as e:
                print(f"    Error creating {plot_type} plot for {sample}: {e}")
                ax.text(0.5, 0.5, f'Error loading\n{sample}', 
                        ha='center', va='center', transform=ax.transAxes, fontsize=16)
    
    for row_idx in range(len(available_plots), 4):
        for col_idx in range(3):
            axes[row_idx, col_idx].axis('off')
    
    plt.tight_layout(pad=12.0, h_pad=1.0, w_pad=1.0, rect=[0.05, 0.05, 0.95, 0.95])
    
    for row_idx, plot_type in enumerate(available_plots):
        if plot_type not in plot_files or not os.path.exists(plot_files[plot_type]):
            continue
        
        pos = axes[row_idx, 0].get_position()
        row_top_y = pos.y1 - 0.01
        
        fig.text(0.08, row_top_y, chr(65 + row_idx), 
                fontsize=48, fontweight='bold', ha='center', va='top', rotation=0)
        
        row_center_y = pos.y0 + (pos.y1 - pos.y0) / 2 
        row_titles = {
            'deamination': 'Deamination',
            'inserts': 'Insert Length', 
            'transrate': 'TransRate2',
            'busco': 'BUSCO'
        }
        
        title_x = pos.x0 - 0.08
        
        fig.text(title_x, row_center_y, row_titles[plot_type], 
                fontsize=32, fontweight='bold', ha='center', va='center', rotation=90)
    
    output_file = os.path.join(output_dir, 'representative_plots_concatenated.png')
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Concatenated representative plots saved to: {output_file}")

def create_busco_individual_plot(ax, sample, output_dir):
    if sample.startswith('WA'):
        file_path = f'00_busco/01_data/{sample}_paired_busco/short_summary.specific.viridiplantae_odb10.{sample}_paired_busco.txt'
    else:
        file_path = f'00_busco/01_data/{sample}_busco/short_summary.specific.viridiplantae_odb10.{sample}_busco.txt'
    
    if not os.path.exists(file_path):
        ax.text(0.5, 0.5, f'File not found\n{sample}', 
                ha='center', va='center', transform=ax.transAxes, fontsize=16)
        return
    
    try:
        busco_data = extract_busco_data(file_path)
        if busco_data is None:
            ax.text(0.5, 0.5, f'No BUSCO data\n{sample}', 
                    ha='center', va='center', transform=ax.transAxes, fontsize=24)
            return
        
        categories = ['S', 'D', 'F', 'M']
        category_labels = ['Complete\nSingle-copy', 'Complete\nDuplicated', 'Fragmented', 'Missing']
        percentages = [busco_data['S'], busco_data['D'], busco_data['F'], busco_data['M']]
        colors = [COLORS['busco'][cat] for cat in categories]
        
        bars = ax.bar(categories, percentages, color=colors, alpha=0.8)
        
        ax.set_xticks(range(len(categories)))
        ax.set_xticklabels(category_labels, rotation=90, ha='center')
        
        for bar, percentage in zip(bars, percentages):
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height + 0.5,
                    f'{percentage:.1f}%', ha='center', va='bottom', fontsize=20, fontweight='bold')
        
        ax.set_title(f'{get_real_sample_name(sample)}', fontsize=24, fontweight='bold')
        ax.set_ylim(0, 100)
        ax.tick_params(axis='x', labelsize=18)
        ax.tick_params(axis='y', labelsize=20)
        ax.set_ylabel('Percentage (%)', fontsize=22, labelpad=10)
        ax.grid(True, alpha=0.3, axis='y')
        ax.set_yticks([0, 20, 40, 60, 80, 100])
        
    except Exception as e:
        ax.text(0.5, 0.5, f'Error processing\n{sample}', 
                ha='center', va='center', transform=ax.transAxes, fontsize=24)

def create_deamination_individual_plot(ax, sample, output_dir):
    file_path = os.path.join(ROOT_DIR, 'FINALdata', 'deamination_data', sample, 'misincorporation.txt')
    
    if not os.path.exists(file_path):
        ax.text(0.5, 0.5, f'File not found\n{sample}', 
                ha='center', va='center', transform=ax.transAxes, fontsize=16)
        return
    
    try:
        plot_data = get_deamination_plot_data(file_path)
        
        color_blue = '#56B4E9'
        color_red = '#F04442'
        
        x_5p = list(range(1, 26))
        ax.plot(x_5p, plot_data['5p_AtoG'], color=color_blue, label='A>I', linewidth=2, alpha=0.8)
        ax.plot(x_5p, plot_data['5p_CtoT'], color=color_red, label='C>U', linewidth=2, alpha=0.8)
        
        x_3p = list(range(-25, 0))
        ax.plot(x_3p, plot_data['3p_AtoG'], color=color_blue, linewidth=2, alpha=0.8)
        ax.plot(x_3p, plot_data['3p_CtoT'], color=color_red, linewidth=2, alpha=0.8)
        
        ax.set_title(f'{get_real_sample_name(sample)}', fontsize=24, fontweight='bold')
        ax.set_ylim(0.00, 0.002)
        ax.set_yticks([0.0000, 0.0004, 0.0008, 0.0012, 0.0016, 0.0020])
        ax.set_xlim(-25.5, 25.5)
        
        tick_positions = [-25, -20, -15, -10, -5, 0, 5, 10, 15, 20, 25]
        tick_labels = ['25', '20', '15', '10', '5', '0', '-5', '-10', '-15', '-20', '-25']
        ax.set_xticks(tick_positions)
        ax.set_xticklabels(tick_labels)
        ax.tick_params(axis='x', labelsize=20)
        ax.tick_params(axis='y', labelsize=20)
        ax.set_xlabel('Position', fontsize=22, labelpad=10)
        ax.set_ylabel('Misincorporation Frequency', fontsize=22, labelpad=10)
        ax.legend(fontsize=18, loc='upper right')
        
        ax.grid(True, alpha=0.3)
        ax.axvline(x=0, color='gray', linestyle='--', alpha=0.5)
        
    except Exception as e:
        ax.text(0.5, 0.5, f'Error processing\n{sample}', 
                ha='center', va='center', transform=ax.transAxes, fontsize=24)

def create_inserts_individual_plot(ax, sample, output_dir):
    if sample.startswith('WA'):
        file_path = os.path.join(ROOT_DIR, 'FINALdata', 'inserts_data', f'{sample}_paired.csv')
    else:
        file_path = os.path.join(ROOT_DIR, 'FINALdata', 'inserts_data', f'{sample}.csv')
    
    if not os.path.exists(file_path):
        ax.text(0.5, 0.5, f'File not found\n{sample}', 
                ha='center', va='center', transform=ax.transAxes, fontsize=16)
        return
    
    try:
        df = pd.read_csv(file_path)
        df_filtered = df[df['insert_length'] >= 1].copy()
        
        lower_bound = 1
        upper_bound = 800
        
        insert_length_counts = df_filtered['insert_length'].value_counts().sort_index()
        x = insert_length_counts.index.values
        y = insert_length_counts.values
        
        mask_pos = (x >= lower_bound) & (x <= upper_bound)
        
        if np.any(mask_pos):
            ax.plot(x[mask_pos], y[mask_pos], color='#000000', linewidth=2, label='insert')
        
        ax.set_title(f'{get_real_sample_name(sample)}', fontsize=24, fontweight='bold')
        ax.set_xlim(lower_bound, upper_bound)
        ax.set_ylim(-1, np.max(y[mask_pos]) + np.max(y[mask_pos]) * 0.1)
        
        ax.set_xticks([1] + list(np.arange(100, upper_bound + 1, 100)))
        ax.tick_params(axis='x', rotation=45, labelsize=20)
        ax.tick_params(axis='y', labelsize=20)
        
        ax.set_xlabel('Insert Length (bp)', fontsize=22, labelpad=10)
        ax.set_ylabel('Count', fontsize=22, labelpad=10)
        
        ax.grid(True, alpha=0.3)
        
    except Exception as e:
        ax.text(0.5, 0.5, f'Error processing\n{sample}', 
                ha='center', va='center', transform=ax.transAxes, fontsize=24)

def create_transrate_individual_plot(ax, sample, output_dir):
    file_path = os.path.join(ROOT_DIR, 'FINALdata', 'transrate_data', 'tr2_assembly', f'{sample}.csv')
    
    if not os.path.exists(file_path):
        ax.text(0.5, 0.5, f'File not found\n{sample}', 
                ha='center', va='center', transform=ax.transAxes, fontsize=16)
        return
    
    try:
        df = pd.read_csv(file_path)
        
        score_columns = ['sCnuc_Harmonic', 'sCcov_Harmonic', 'sCord_Harmonic', 'sCseg_Harmonic']
        score_labels = ['sCnuc', 'sCcov', 'sCord', 'sCseg']
        
        scores = []
        for col in score_columns:
            if col in df.columns:
                scores.append(df[col].iloc[0])
            else:
                scores.append(0)
        
        x_pos = np.arange(len(score_labels))
        bars = ax.bar(x_pos, scores, color=[COLORS['transrate'][col] for col in score_columns], alpha=0.8)
        
        for bar, score in zip(bars, scores):
            height = bar.get_height()
            ax.text(bar.get_x() + bar.get_width()/2., height + 0.01,
                    f'{score:.3f}', ha='center', va='bottom', fontsize=20, fontweight='bold')
        
        ax.set_title(f'{get_real_sample_name(sample)}', fontsize=24, fontweight='bold')
        ax.set_ylim(0, 1.1)
        
        ax.set_xticks(x_pos)
        ax.set_xticklabels(score_labels, fontsize=20)
        ax.tick_params(axis='y', labelsize=20)
        
        ax.set_ylabel('Score', fontsize=22, labelpad=10)
        ax.grid(True, alpha=0.3, axis='y')
        
    except Exception as e:
        ax.text(0.5, 0.5, f'Error processing\n{sample}', 
                ha='center', va='center', transform=ax.transAxes, fontsize=24)
