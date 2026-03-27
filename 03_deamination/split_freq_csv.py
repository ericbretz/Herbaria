import pandas as pd
import os

df = pd.read_csv('/home/eric/scratch/Herbaria/00_deamination/03_results/frequencies.csv')

output_dir = '/home/eric/scratch/Herbaria/00_deamination/03_results/freq_categorized'
os.makedirs(output_dir, exist_ok=True)

substitution_categories = df['Substitution'].unique()

for category in substitution_categories:
    if pd.isna(category) or category == 'Substitution':
        continue
    
    category_df = df[df['Substitution'] == category]
    
    safe_category = str(category).replace('>', '_to_').replace('-', 'del').replace('->', 'ins_')
    filename = f"{safe_category}.csv"
    filepath = os.path.join(output_dir, filename)
    
    category_df.to_csv(filepath, index=False)
