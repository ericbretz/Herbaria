import os
import sys
import glob
import json
import subprocess
from concurrent.futures import ThreadPoolExecutor

class Damage:
    def __init__(self):
        self.files_bam   = glob.glob(os.path.join('*', 'transrate2', 'nuclear', '*.postSample.sorted.bam'))
        self.files_fasta = glob.glob(os.path.join('*', '*.fa'))

        self.dir_output = 'deamination'
        self.done_count = 0
        os.makedirs(self.dir_output, exist_ok=True)

        self.sample_files = {}

        for file_bam in self.files_bam:
            sample_name = file_bam.split('/')[-1].split('.')[0]
            self.sample_files[sample_name] = {'bam': file_bam}

        for file_fasta in self.files_fasta:
            sample_name = file_fasta.split('/')[-1].split('.')[0]
            if sample_name in self.sample_files:
                self.sample_files[sample_name]['fasta'] = file_fasta

    def run_mapdamage(self, sample, files):
        output = os.path.join(self.dir_output, sample)
        os.makedirs(output, exist_ok=True)
        bam_file = files['bam']
        fasta_file = files['fasta']

        command = f"mapDamage -i {bam_file} -r {fasta_file} -d {output} --merge-libraries"
        results = subprocess.run(command, shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
        if results.returncode != 0:
            print(f"Error running mapDamage2.0 for {sample}")
            sys.exit()
        self.done_count += 1
        print(f"{self.done_count} / {len(self.sample_files)} -- {sample}")

    def mapDamage_threading(self):
        threads = len(self.sample_files)
        with ThreadPoolExecutor(max_workers=threads) as executor:
            futures = [executor.submit(self.run_mapdamage, sample, files) for sample, files in self.sample_files.items()]
            for future in futures:
                future.result()

if __name__ == '__main__':
    damage = Damage()
    damage.mapDamage_threading()