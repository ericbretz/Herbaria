import shutil
from glob import glob

samples = ['DAL192','DAL193', 'DAL195', 'DAL212', 'DAL218',
           'DAL224', 'DAL227', 'WA02', 'WA03', 
           'WA07', 'WA08', 'WA09', 'WA11', 
           'WA12', 'WA13', 'WA14', 'WA18', 
           'WA22', 'WA25', 'WA26']

csv_assemblies = glob('../*/transrate2/nuclear/assembly.csv')
csv_assemblies.sort()
sample_names   = [x.split('/')[1] for x in csv_assemblies]
sample_names.sort()
zipped         = zip(samples, csv_assemblies)
data_path      = '01_data/tr2_assembly'

print(len(csv_assemblies))

# sys.exit()
for sample, csv_assembly in zipped:
    shutil.copy(csv_assembly, f'{data_path}/{sample}.csv')