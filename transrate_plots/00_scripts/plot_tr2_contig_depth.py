import pandas as pd
import matplotlib.pyplot as plt
from glob import glob
import os
import numpy as np

SAMPLE_CONFIGS = {
    'fresh': {
        'samples': ['DAL192', 'DAL193', 'DAL195', 'DAL212', 'DAL218', 'DAL224', 'DAL227'],
        'color': '#C5CF9D'
    },
    'herbaria': {
        'samples': ['WA22', 'WA25', 'WA26'],
        'color': '#BA6F6D'
    },
    'silica': {
        'samples': ['WA13', 'WA18', 'WA02', 'WA08', 'WA03', 'WA07', 'WA12', 'WA09', 'WA11', 'WA14'],
        'color': '#80B0B8'
    }
}

def get_sample_category(sample_name):
    for category, config in SAMPLE_CONFIGS.items():
        if sample_name in config['samples']:
            return category
    return None

def create_depth_boxplot(all_data):
    
    sample_type_data = {
        'fresh': [],
        'herbaria': [],
        'silica': []
    }
    
    for data in all_data:
        sample_name = data['csv_basename']
        depth_values = data['depth_count']
        
        category = get_sample_category(sample_name)
        if category:
            sample_type_data[category].append(depth_values)
    
    fig, ax = plt.subplots(figsize=(10, 8))
    
    plot_data = []
    plot_labels = []
    plot_colors = []
    
    for sample_type in ['fresh', 'herbaria', 'silica']:
        if sample_type_data[sample_type]:
            plot_data.append(sample_type_data[sample_type])
            plot_labels.append(sample_type.capitalize())
            plot_colors.append(SAMPLE_CONFIGS[sample_type]['color'])
    
    bp = ax.boxplot(plot_data, tick_labels=plot_labels, patch_artist=True)
    
    for patch, color in zip(bp['boxes'], plot_colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)
    
    ax.set_title('Contig Count by Sample Type', fontsize=16, fontweight='bold')
    ax.set_ylabel('Number of Contigs', fontsize=12)
    ax.set_xlabel('Sample Type', fontsize=12)
    
    
    ax.grid(axis='x', alpha=0.0)
    
    plt.tight_layout()
    output_path = '02_plots/tr2_contigs/contig_count.png'
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"{'Box plot':<34}{output_path}")

def create_contig_count_barchart(all_data):
    
    sorted_data = sorted(all_data, key=lambda x: x['csv_basename'])
    
    sample_names = [data['csv_basename'] for data in sorted_data]
    contig_counts = [data['depth_count'] for data in sorted_data]
    
    sample_colors = []
    for data in sorted_data:
        category = get_sample_category(data['csv_basename'])
        if category:
            sample_colors.append(SAMPLE_CONFIGS[category]['color'])
        else:
            sample_colors.append('#CCCCCC')
    
    fig, ax = plt.subplots(figsize=(14, 8))
    
    bars = ax.bar(range(len(sample_names)), contig_counts, color=sample_colors, alpha=0.8)
    
    ax.set_title('Contig Count by Sample', fontsize=16, fontweight='bold')
    ax.set_ylabel('Number of Contigs', fontsize=12)
    ax.set_xlabel('Sample', fontsize=12)
    
    ax.set_xticks(range(len(sample_names)))
    ax.set_xticklabels(sample_names, rotation=45, ha='right')
    
    
    ax.grid(axis='x', alpha=0.0)
        
    for i, (bar, count) in enumerate(zip(bars, contig_counts)):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + max(contig_counts)*0.01, 
                str(count), ha='center', va='bottom', fontsize=10)
    
    legend_elements = []
    for sample_type, config in SAMPLE_CONFIGS.items():
        legend_elements.append(plt.Rectangle((0,0),1,1, facecolor=config['color'], 
                                          alpha=0.8, label=sample_type.capitalize()))
    
    ax.legend(handles=legend_elements, loc='upper right')
    
    plt.tight_layout()
    output_path = '02_plots/tr2_contigs/contig_count_barchart.png'
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"{'Bar chart':<34}{output_path}")

def create_fragments_per_contig_barchart(all_data):    
    sorted_data = sorted(all_data, key=lambda x: x['csv_basename'])
    
    sample_names = [data['csv_basename'] for data in sorted_data]
    avg_fragments = [data['avg_fragments'] for data in sorted_data]
    avg_both_mapped = [data['avg_both_mapped'] for data in sorted_data]
    
    sample_colors = []
    for data in sorted_data:
        category = get_sample_category(data['csv_basename'])
        if category:
            sample_colors.append(SAMPLE_CONFIGS[category]['color'])
        else:
            sample_colors.append('#CCCCCC')
    
    fig, ax = plt.subplots(figsize=(16, 8))
    
    x = np.arange(len(sample_names))
    width = 0.35
    
    bars1 = ax.bar(x - width/2, avg_fragments, width, label='Fragments', color=sample_colors, alpha=0.8)
    
    bars2 = ax.bar(x + width/2, avg_both_mapped, width, label='Both Mapped', color=sample_colors, alpha=0.8, hatch='///')
    
    ax.set_title('Average Fragments and Both Mapped per Contig by Sample', fontsize=16, fontweight='bold')
    ax.set_ylabel('Average Count per Contig', fontsize=12)
    ax.set_xlabel('Sample', fontsize=12)
    
    ax.set_xticks(x)
    ax.set_xticklabels(sample_names, rotation=45, ha='right')
    
    
    ax.grid(axis='x', alpha=0.0)
        
    for i, (bar1, bar2, avg_frag, avg_both) in enumerate(zip(bars1, bars2, avg_fragments, avg_both_mapped)):
        ax.text(bar1.get_x() + bar1.get_width()/2, bar1.get_height() + max(max(avg_fragments), max(avg_both_mapped))*0.01, 
                f'{avg_frag:.1f}', ha='center', va='bottom', fontsize=9)
        ax.text(bar2.get_x() + bar2.get_width()/2, bar2.get_height() + max(max(avg_fragments), max(avg_both_mapped))*0.01, 
                f'{avg_both:.1f}', ha='center', va='bottom', fontsize=9)
    
    legend_elements = []
    for sample_type, config in SAMPLE_CONFIGS.items():
        legend_elements.append(plt.Rectangle((0,0),1,1, facecolor=config['color'], 
                                          alpha=0.8, label=sample_type.capitalize()))
    
    ax.legend(loc='upper right')
    
    plt.tight_layout()
    output_path = '02_plots/tr2_contigs/contig_avg_fragments_per_contig.png'
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"{'Fragments per contig bar chart':<34}{output_path}")

def create_fragments_ratio_barchart(all_data, csv_assemblies):    
    assembly_data = {}
    for csv_file in csv_assemblies:
        sample_name = os.path.basename(csv_file).split('.')[0]
        try:
            df = pd.read_csv(csv_file)
            if 'p_fragments_mapped' in df.columns:
                p_fragments_mapped = df['p_fragments_mapped'].iloc[0] * 100
                assembly_data[sample_name] = p_fragments_mapped
            else:
                print(f"Warning: 'p_fragments_mapped' column not found in {csv_file}")
                assembly_data[sample_name] = 0
        except Exception as e:
            print(f"Error reading {csv_file}: {e}")
            assembly_data[sample_name] = 0
    
    sorted_data = sorted(all_data, key=lambda x: x['csv_basename'])
    
    sample_names = []
    p_fragments_values = []
    
    for data in sorted_data:
        sample_name = data['csv_basename']
        if sample_name in assembly_data:
            sample_names.append(sample_name)
            p_fragments_values.append(assembly_data[sample_name])
        else:
            print(f"Warning: No assembly data found for sample {sample_name}")
    
    if not sample_names:
        print("Error: No valid assembly data found for any samples")
        return
    
    sample_colors = []
    for sample_name in sample_names:
        category = get_sample_category(sample_name)
        if category:
            sample_colors.append(SAMPLE_CONFIGS[category]['color'])
        else:
            sample_colors.append('#CCCCCC')
    
    fig, ax = plt.subplots(figsize=(14, 8))
    
    bars = ax.bar(range(len(sample_names)), p_fragments_values, color=sample_colors, alpha=0.8)
    
    ax.set_title('Percentage of Fragments Mapped by Sample', fontsize=16, fontweight='bold')
    ax.set_ylabel('Percentage of Fragments Mapped (%)', fontsize=12)
    ax.set_xlabel('Sample', fontsize=12)
    
    ax.set_xticks(range(len(sample_names)))
    ax.set_xticklabels(sample_names, rotation=45, ha='right')
    
    
    ax.grid(axis='x', alpha=0.0)
        
    for i, (bar, value) in enumerate(zip(bars, p_fragments_values)):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + max(p_fragments_values)*0.01, 
                f'{value:.1f}%', ha='center', va='bottom', fontsize=10)
    
    legend_elements = []
    for sample_type, config in SAMPLE_CONFIGS.items():
        legend_elements.append(plt.Rectangle((0,0),1,1, facecolor=config['color'], 
                                          alpha=0.8, label=sample_type.capitalize()))
    
    ax.legend(handles=legend_elements, loc='upper right')
    
    plt.tight_layout()
    output_path = '02_plots/tr2_contigs/contig_p_fragments_mapped.png'
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"{'P_fragments_mapped bar chart':<34}{output_path}")

def create_fragments_ratio_barchart_original(all_data):    
    sorted_data = sorted(all_data, key=lambda x: x['csv_basename'])
    
    sample_names = [data['csv_basename'] for data in sorted_data]
    
    ratios = []
    for data in sorted_data:
        if data['avg_both_mapped'] > 0:
            ratio = data['avg_both_mapped'] / data['avg_fragments'] 
        else:
            ratio = 0
        ratios.append(ratio)

    sample_colors = []
    for data in sorted_data:
        category = get_sample_category(data['csv_basename'])
        if category:
            sample_colors.append(SAMPLE_CONFIGS[category]['color'])
        else:
            sample_colors.append('#CCCCCC')
    
    fig, ax = plt.subplots(figsize=(14, 8))
    
    bars = ax.bar(range(len(sample_names)), ratios, color=sample_colors, alpha=0.8)
    
    ax.set_title('Fragments Mapped to Both Mapped Ratio by Sample', fontsize=16, fontweight='bold')
    ax.set_ylabel('Ratio (Fragments / Both Mapped)', fontsize=12)
    ax.set_xlabel('Sample', fontsize=12)
    
    ax.set_xticks(range(len(sample_names)))
    ax.set_xticklabels(sample_names, rotation=45, ha='right')
    
    
    ax.grid(axis='x', alpha=0.0)
        
    for i, (bar, ratio) in enumerate(zip(bars, ratios)):
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + max(ratios)*0.01, 
                f'{ratio:.2f}', ha='center', va='bottom', fontsize=10)
    
    legend_elements = []
    for sample_type, config in SAMPLE_CONFIGS.items():
        legend_elements.append(plt.Rectangle((0,0),1,1, facecolor=config['color'], 
                                          alpha=0.8, label=sample_type.capitalize()))
    
    ax.legend(handles=legend_elements, loc='upper right')
    
    plt.tight_layout()
    output_path = '02_plots/tr2_contigs/contig_fragments_both_mapped_ratio.png'
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"{'Fragments ratio bar chart':<34}{output_path}")

if __name__ == '__main__':
    csv_files = glob('01_data/tr2_contigs/*.csv')
    csv_assemblies = glob('01_data/tr2_assembly/*.csv')
    
    all_data = []
    
    for csv_file in csv_files:
        sample_name = os.path.basename(csv_file).split('.')[0]

        df = pd.read_csv(csv_file)
        row_count = len(df) - 1
        
        if 'fragments' in df.columns:
            avg_fragments = df['fragments'].mean()
        else:
            avg_fragments = 0
        
        if 'both_mapped' in df.columns:
            avg_both_mapped = df['both_mapped'].mean()
        else:
            avg_both_mapped = 0
        
        all_data.append({
            'csv_basename': sample_name,
            'depth_count': row_count,
            'avg_fragments': avg_fragments,
            'avg_both_mapped': avg_both_mapped
        })
    
    create_depth_boxplot(all_data)
    create_contig_count_barchart(all_data)
    create_fragments_per_contig_barchart(all_data)
    create_fragments_ratio_barchart(all_data, csv_assemblies)
    create_fragments_ratio_barchart_original(all_data)
