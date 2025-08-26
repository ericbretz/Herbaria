'''
This script analyzes reads by measuring the distances between them.

Key points:
- Processes BAM files by breaking them into chunks
- Uses multiprocessing to run multiple chunks simultaneously
- Each process handles its own chunk of references independently
- For each reference, reads are processed in batches of 1000 to manage memory
- For each read:
  - Finds its mate
  - Calculates soft-clipped positions
  - Determines insert and overlap lengths
- Results are combined and saved to CSV
- Progress bar shows completion status

The parallel processing makes it much faster than single thread processing.
'''

'''
ref_name       - Reference name
read_name      - Read name
read_start     - Read start position
read_end       - Read end position
mate_start     - Mate start position
mate_end       - Mate end position
read_length    - Read length
mate_length    - Mate length
insert_length  - Insert length
overlap_length - Overlap length


reference:              -----GATCGTCGTAAACTGACTAGTACTTGACGATGATGACGTAGCAGTTTGTCACAGGTTT
read:                   atgctGATCGTCGTAAAC
mate:                   |    |           |            GACGATGATGACGTAGCA
overlap: 0              |    |           |            |                |
insert: 12              |    |           |............|                |
read start: 5           |    ^           |            |                |
read end: 17            |                ^            |                |
read adjusted start: 0  ^                             |                |
mate start: 30                                        ^                |
mate end: 47                                                           ^
read total length: 18   |________________|
mate total length: 18                                 |________________|


"adjusted" just means accounting for soft-clipped bases.

'''

import os
import sys
import pandas as pd
import numpy as np
import pysam
import glob
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing as mp
from tqdm import tqdm
import psutil

class ReadDistance:
    def __init__(self):
        # Glob is pretty neat, its going to find all the bam files for us so we dont have to manually iterate every directory.
        self.files_bam      = glob.glob("**/*sorted.bam", recursive=True)
        self.index_files    = glob.glob("**/*sorted.bam.bai", recursive=True)
        # csv columns
        self.pandas_columns = ['ref_name', 'read_name', 'read_start', 'read_end', 'mate_start', 'mate_end', 
                               'read_length', 'mate_length', 'insert_length', 'overlap_length']
        # save directory
        self.save_dir       =  'read_distance'
        
        # Calculate number of workers based the total memory we have
        total_memory = psutil.virtual_memory().total
        self.max_workers    = min(mp.cpu_count() - 1, max(1, int(total_memory / (1024 * 1024 * 1024))))  # Leave one core free
        
        # Calculate batch size based on how much of that total memory is actually available
        available_memory    = psutil.virtual_memory().available
        self.batch_size     = min(10000, max(1000, int(available_memory / (1024 * 1024 * 10))))  # 10MB per batch
        
        # Calculate references per process, 100MB per reference, 10 references per process
        # Each reference is going to be around 50-100mb
        self.refs_per_process = max(1, min(10, int(available_memory / (1024 * 1024 * 100))))     # 100MB per reference, 10 references per process
        
        os.makedirs(self.save_dir, exist_ok=True)
        print(f"Using {self.max_workers} workers")
        print(f"References per process: {self.refs_per_process}")
        print(f"Batch size: {self.batch_size} reads per batch")

    def process_read_batch(self, reads, bam):
        '''Process a batch of reads'''
        # Build a dictionary of reads by query_name
        reads_by_name = {}
        for read in reads:
            if not read.is_proper_pair:
                continue
            reads_by_name.setdefault(read.query_name, []).append(read)

        n_pairs = sum(1 for v in reads_by_name.values() if len(v) == 2)
        data = {
            'ref_name'      : np.empty(n_pairs, dtype=object),
            'read_name'     : np.empty(n_pairs, dtype=object),
            'read_start'    : np.zeros(n_pairs, dtype=np.int32),
            'read_end'      : np.zeros(n_pairs, dtype=np.int32),
            'mate_start'    : np.zeros(n_pairs, dtype=np.int32),
            'mate_end'      : np.zeros(n_pairs, dtype=np.int32),
            'read_length'   : np.zeros(n_pairs, dtype=np.int32),
            'mate_length'   : np.zeros(n_pairs, dtype=np.int32),
            'insert_length' : np.zeros(n_pairs, dtype=np.int32),
            'overlap_length': np.zeros(n_pairs, dtype=np.int32)
        }

        valid_count = 0
        for read_name, pair in reads_by_name.items():
            # Find exactly one read1 and one read2
            read1s = [r for r in pair if r.is_read1]
            read2s = [r for r in pair if r.is_read2]
            if len(read1s) != 1 or len(read2s) != 1:
                continue  # skip multimapped or ambiguous pairs

            read = read1s[0]
            mate = read2s[0]

            # The rest of the calculations (as before)
            read_soft_start     = read.query_alignment_start
            read_soft_end       = read.query_length - read.query_alignment_end
            mate_soft_start     = mate.query_alignment_start
            mate_soft_end       = mate.query_length - mate.query_alignment_end

            read_total_length   = read.query_alignment_length + read_soft_start + read_soft_end
            mate_total_length   = mate.query_alignment_length + mate_soft_start + mate_soft_end

            read_adjusted_start = read.reference_start - read_soft_start
            read_adjusted_end   = read.reference_end + read_soft_end
            mate_adjusted_start = mate.reference_start - mate_soft_start
            mate_adjusted_end   = mate.reference_end + mate_soft_end

            read_start = min(read_adjusted_start, read_adjusted_end)
            read_end   = max(read_adjusted_start, read_adjusted_end)
            mate_start = min(mate_adjusted_start, mate_adjusted_end)
            mate_end   = max(mate_adjusted_start, mate_adjusted_end)

            if read_end < mate_start:
                insert_length = mate_start - read_end
                overlap_length = 0
            elif mate_end < read_start:
                insert_length = read_start - mate_end
                overlap_length = 0
            else:
                insert_length = 0
                overlap_length = min(read_end, mate_end) - max(read_start, mate_start)

            if (read.reference_start > read.reference_end) != (mate.reference_start > mate.reference_end):
                insert_length = -insert_length

            if 'N' in read.cigarstring or 'N' in mate.cigarstring:
                continue

            if abs(read.reference_start - mate.reference_start) > 1000:
                continue

            data['ref_name'][valid_count]       = read.reference_name
            data['read_name'][valid_count]      = read.query_name
            data['read_start'][valid_count]     = read_adjusted_start
            data['read_end'][valid_count]       = read_adjusted_end
            data['mate_start'][valid_count]     = mate_adjusted_start
            data['mate_end'][valid_count]       = mate_adjusted_end
            data['read_length'][valid_count]    = read_total_length
            data['mate_length'][valid_count]    = mate_total_length
            data['insert_length'][valid_count]  = insert_length
            data['overlap_length'][valid_count] = overlap_length

            valid_count += 1

        if valid_count > 0:
            return pd.DataFrame({k: v[:valid_count] for k, v in data.items()})
        return pd.DataFrame(columns=self.pandas_columns)

    def process_reference_batch(self, file_bam, references):
        '''Process multiple references in a single process'''
        try:
            # Open BAM file
            bam            = pysam.AlignmentFile(file_bam, "rb")
            all_batches    = []
            
            for reference in references:
                try:
                    # Get all reads for the reference
                    reads       = list(bam.fetch(reference=reference))
                    total_reads = len(reads)
                    
                    # Process reads in batches
                    for i in range(0, total_reads, self.batch_size):
                        batch    = reads[i:i + self.batch_size]
                        batch_df = self.process_read_batch(batch, bam)
                        if not batch_df.empty:
                            all_batches.append(batch_df)
                except Exception as e:
                    print(f"\nError processing reference {reference}: {str(e)}")
                    continue
            
            bam.close()
            
            # Combine all of the batches
            if all_batches:
                return pd.concat(all_batches, ignore_index=True)
            return pd.DataFrame(columns=self.pandas_columns)
            
        except Exception as e:
            print(f"\nError in process processing references {references}: {str(e)}")
            if 'bam' in locals():
                bam.close()
            return pd.DataFrame(columns=self.pandas_columns)

    def get_read_distance(self):
        for file_bam in self.files_bam:
            # Check if output file already exists
            bam_name = file_bam.split('/')[-1].split('.')[0]
            output_file = f'{self.save_dir}/{bam_name}.csv'
            
            if os.path.exists(output_file):
                print(f'\nSkipping {file_bam} - output file already exists at {output_file}')
                continue
                
            print(f'\nProcessing file: {file_bam}')
            
            # Open BAM file just to get references
            with pysam.AlignmentFile(file_bam, "rb") as bam:
                references = bam.references
            
            # Create batches of references for each process
            ref_batches = [references[i:i + self.refs_per_process] 
                         for i in range(0, len(references), self.refs_per_process)]
            
            # Create a list to store all results
            all_results = []
            
            # Process reference batches in parallel
            with ProcessPoolExecutor(max_workers=self.max_workers) as executor:
                # Submit all reference batches
                future_to_batch = {
                    executor.submit(self.process_reference_batch, file_bam, batch): batch 
                    for batch in ref_batches
                }
                
                # Just a little progress bar to keep track of things
                with tqdm(total=len(ref_batches), desc="Processing reference batches", unit="batch") as pbar:
                    for future in as_completed(future_to_batch):
                        batch = future_to_batch[future]
                        try:
                            result_df = future.result()
                            if not result_df.empty:
                                all_results.append(result_df)
                        except Exception as e:
                            print(f"\nError processing batch {batch}: {str(e)}")
                        pbar.update(1)
            
            # Combine all the results
            if all_results:
                final_df = pd.concat(all_results, ignore_index=True)
                final_df.to_csv(output_file, index=False)
                print(f"\nSaved {len(final_df)} processed reads to {output_file}")
            else:
                print(f"\nNo valid reads found in {file_bam}")

if __name__ == "__main__":
    # Set the start method for multiprocessing
    # You want to use spawn for stuff like this because fork is going to cause problems
    mp.set_start_method('spawn')
    read_distance = ReadDistance()
    read_distance.get_read_distance()