import os
import sys
import json
import argparse
import subprocess
import shutil

class Herbaria:
    def __init__(self):
        self.sampleCurrent :int = 0
        self.sampleTotal   :int = 0
        self.sampleName    :str = ''
        self.sampleOrganism:str = ''
        self.pathDict      :dict = {}
        self.pathIter      :bool = False

        self.argsSample    :list = []
        self.argsOrganism  :list = []
        self.progress_data :dict = {}

    def parseArgs(self):
        parser = argparse.ArgumentParser(description='Herbaria')
        parser.add_argument('-s', '--sample',   type=str, help='Sample Type (all, dal, wa, none)', default='all')
        parser.add_argument('-o', '--organism', type=str, help='Organelle Type (all, chloro, mito, nuclear, none)', default='all')
        parser.add_argument('-p', '--path',     action='store_true', help='Reiterate File Paths', default=False)
        args = parser.parse_args()

        if args.path:
            self.pathIter = True

        if args.sample not in ['all', 'dal', 'wa', 'none']:
            print('Error: Invalid Sample Type (all, dal, wa, none)')
            sys.exit()

        if any(org not in ['all', 'chloro', 'mito', 'nuclear', 'none'] for org in args.organism.split(',')):
            print('Error: Invalid Organism Type (all, chloro, mito, nuclear, none)')
            sys.exit()

        if args.sample == 'none' or args.organism == 'none':
            self.argsSample = []
            self.argsOrganism = []
        elif args.sample and args.sample != 'all':
            self.argsSample = [s for s in args.sample.split(',') if s in ['wa', 'dal']]
        else:
            self.argsSample = ['wa', 'dal']

        if args.organism and args.organism != 'all':
            self.argsOrganism = [o for o in args.organism.split(',') if o in ['chloro', 'mito', 'nuclear']]
        else:
            self.argsOrganism = ['chloro', 'mito', 'nuclear']


    def iterateDir(self):
        if os.path.exists('pathDict.json') and not self.pathIter:
            self.pathDict = json.load(open('pathDict.json'))
            return
        for root, dirs, files in os.walk('.'):
            for d in dirs:
                if d.startswith('DAL') or d.startswith('WA'):
                    self.pathDict[d] = {'assembly': '', 'chloro': [], 'mito': [], 'nuclear': []}
                    for file in os.listdir(os.path.join(root, d)):
                        if file.endswith('.cds.fa'):
                            self.pathDict[d]['assembly'] = file
                        elif file.startswith('chloro'):
                            self.pathDict[d]['chloro'].append(file)
                        elif file.startswith('mito'):
                            self.pathDict[d]['mito'].append(file)
                        elif file.startswith('nuclear'):
                            self.pathDict[d]['nuclear'].append(file)
        for v in self.pathDict.values():
            for key in ['chloro', 'mito', 'nuclear']:
                v[key].sort()
        with open('pathDict.json', 'w') as json_file:
            json.dump(self.pathDict, json_file, indent=4)

    def getTotal(self):
        for k,v in self.pathDict.items():
            if any(k.startswith(s.upper()) for s in self.argsSample):
                for key,var in v.items():
                    if k in self.progress_data and key in self.progress_data[k]:
                        continue
                    if key in self.argsOrganism and var:
                        self.sampleTotal += 1
                    
    def iterateRuns(self):
        for key, value in self.pathDict.items():
            for sample in self.argsSample:
                if key.startswith(sample.upper()) and value['assembly'] != '':
                    for organism in self.argsOrganism:
                        if key in self.progress_data and organism in self.progress_data[key]:
                            continue
                        self.sampleCurrent += 1
                        self.sampleName = key
                        self.sampleOrganism = organism
                        output = os.path.join(key, 'transrate2', organism)
                        assemblyPath = os.path.join(key, value['assembly'])
                        try:
                            leftPath = os.path.join(key, value[organism][0])
                            rightPath = os.path.join(key, value[organism][1])
                        except:
                            continue
                        keys = [key, organism]
                        color_codes = {'mito': '91', 'chloro': '92', 'nuclear': '94'}
                        color_code = color_codes.get(organism, '0')
                        print(f'{key:<10}\033[{color_code}m{organism:<10}\033[0m', end='\r')
                        self.transrateRun(assemblyPath, leftPath, rightPath, output, keys)
                        print(f'{key:<10}\033[{color_code}m{organism:<10}\033[0m{self.sampleCurrent}/{self.sampleTotal}')

    def transrateRun(self, assembly, left, right, output, keys):
        if not os.path.exists(output):
            os.makedirs(output)
        cmd = ['transrate2', '-a', assembly, '-l', left, '-r', right, '-o', output, '-t', str(24), '-s']
        trRun = subprocess.Popen(cmd, shell=False, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout, stderr = trRun.communicate()
        if trRun.returncode != 0:
            print(f'{keys[0]:<10}{keys[1]:<10}FAILED')
        else:
            self.saveProgress(keys[0], keys[1])
        self.transrateCleanup(output)

    def transrateCleanup(self, output):
        for root, dirs, files in os.walk(output):
            for file in files:
                if file.endswith(('.sam', '.bam', '.bai', '.csv', '.fa')):
                    os.rename(os.path.join(root, file), os.path.join(output, file))
        for root, dirs, files in os.walk(output):
            for d in dirs:
                if d != 'logs':
                    shutil.rmtree(os.path.join(root, d))
            for file in files:
                if not file.endswith(('.sam', '.bam', '.bai', '.csv', '.fa')) and 'logs' not in root:
                    os.remove(os.path.join(root, file))

    def saveProgress(self, key, organism):
        progress_data = {}
        if os.path.exists('progress.json'):
            with open('progress.json', 'r') as progress_file:
                progress_data = json.load(progress_file)
        
        if key in progress_data:
            if organism not in progress_data[key]:
                progress_data[key].append(organism)
        else:
            progress_data[key] = [organism]
        
        with open('progress.json', 'w') as progress_file:
            json.dump(progress_data, progress_file, indent=4)

    def loadProgress(self):
        progress_data = {}
        if os.path.exists('progress.json'):
            with open('progress.json', 'r') as progress_file:
                progress_data = json.load(progress_file)
        self.progress_data = progress_data

            
        


main = Herbaria()
main.parseArgs()    #iterate arguments
main.loadProgress() #load progress data
main.iterateDir()   #create a dictionary of file paths (-p flag to reiterate)
main.getTotal()     #print total number of samples
main.iterateRuns()  #iterate through sample(WA,DAL), and organism(chloro,mito,nuclear) types for transrate2 run