import pandas as pd

all_csv = '03_results/frequencies.csv'
max_csv = '03_results/max_frequencies.csv'
ct_csv = '03_results/CT_statistics.csv'
ag_csv = '03_results/AG_statistics.csv'

all_df = pd.read_csv(all_csv)
max_df = pd.read_csv(max_csv)

CT_mean = all_df[all_df['Substitution'] == 'C>T']['Frequency'].mean()
CT_max = all_df[all_df['Substitution'] == 'C>T']['Frequency'].max()
CT_min = all_df[all_df['Substitution'] == 'C>T']['Frequency'].min()
CT_std = all_df[all_df['Substitution'] == 'C>T']['Frequency'].std()
CT_median = all_df[all_df['Substitution'] == 'C>T']['Frequency'].median()

AG_mean = all_df[all_df['Substitution'] == 'A>G']['Frequency'].mean()
AG_max = all_df[all_df['Substitution'] == 'A>G']['Frequency'].max()
AG_min = all_df[all_df['Substitution'] == 'A>G']['Frequency'].min()
AG_std = all_df[all_df['Substitution'] == 'A>G']['Frequency'].std()
AG_median = all_df[all_df['Substitution'] == 'A>G']['Frequency'].median()

ct_stats = []
ag_stats = []

ct_stats.append({
    'Metric': 'Average',
    'Value': CT_mean,
    'Description': 'Mean frequency of C>T substitutions'
})
ct_stats.append({
    'Metric': 'Median',
    'Value': CT_median,
    'Description': 'Median frequency of C>T substitutions'
})
ct_stats.append({
    'Metric': 'Max',
    'Value': CT_max,
    'Description': 'Maximum frequency of C>T substitutions'
})
ct_stats.append({
    'Metric': 'Min',
    'Value': CT_min,
    'Description': 'Minimum frequency of C>T substitutions'
})
ct_stats.append({
    'Metric': 'Std',
    'Value': CT_std,
    'Description': 'Standard deviation of C>T substitution frequencies'
})

ag_stats.append({
    'Metric': 'Average',
    'Value': AG_mean,
    'Description': 'Mean frequency of A>G substitutions'
})
ag_stats.append({
    'Metric': 'Median',
    'Value': AG_median,
    'Description': 'Median frequency of A>G substitutions'
})
ag_stats.append({
    'Metric': 'Max',
    'Value': AG_max,
    'Description': 'Maximum frequency of A>G substitutions'
})
ag_stats.append({
    'Metric': 'Min',
    'Value': AG_min,
    'Description': 'Minimum frequency of A>G substitutions'
})
ag_stats.append({
    'Metric': 'Std',
    'Value': AG_std,
    'Description': 'Standard deviation of A>G substitution frequencies'
})

sample_count = all_df['Sample'].nunique()

ct_df = pd.DataFrame(ct_stats)
ct_df.to_csv(ct_csv, index=False)

ag_df = pd.DataFrame(ag_stats)
ag_df.to_csv(ag_csv, index=False)

