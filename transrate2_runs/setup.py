import os

dal   = ['DAL192', 'DAL193', 'DAL195', 'DAL212', 'DAL218', 'DAL224', 'DAL227']
wa    = ['WA02', 'WA03', 'WA07', 'WA08', 'WA09', 'WA11', 'WA12', 'WA13', 'WA14', 'WA18', 'WA22', 'WA25', 'WA26']

for d in dal:
    os.makedirs(d, exist_ok=True)

for w in wa:
    os.makedirs(w, exist_ok=True)