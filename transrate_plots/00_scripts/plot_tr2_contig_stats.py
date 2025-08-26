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
    """Map sample name to its category (fresh, herbaria, or silica)"""
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
    
    ax.set_title('Read Depth by Sample Type', fontsize=16, fontweight='bold')
    ax.set_ylabel('Number of Reads', fontsize=12)
    ax.set_xlabel('Sample Type', fontsize=12)
    
    plt.tight_layout()
    output_path = '02_plots/tr2_contigs/read_depth.png'
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    plt.close()
    
    print(f"Box plot saved to: {output_path}")

if __name__ == '__main__':
    csv_files = glob('01_data/tr2_contigs/*.csv')
    
    all_data = []
    
    column_names = ['name', 'length', '']

    for csv_file in csv_files:
        sample_name = os.path.basename(csv_file).split('.')[0]

        df = pd.read_csv(csv_file)
        row_count = len(df) - 1
        
        all_data.append({
            'csv_basename': sample_name,
            'depth_count': row_count
        })
        print(f"Processed {csv_file} - {row_count} data rows (excluding header)")

    
    create_depth_boxplot(all_data)
