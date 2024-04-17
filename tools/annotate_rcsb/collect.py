import os
import pandas as pd
from glob import glob
from tqdm import tqdm

path_to_pdbs = '/home/olebedenko/PDB/PDB_2023/structures/all/pdb/'
paths_to_pdbs = glob(os.path.join(path_to_pdbs, 'pdb*.ent.gz'))
original_pdb_ids = [os.path.basename(path_to_pdb).split('.')[0][3:] for path_to_pdb in paths_to_pdbs]

annotated_csv_files = glob('annotations/all_data_*csv')
annotation_df = pd.concat((pd.read_csv(csv_file) for csv_file in tqdm(annotated_csv_files)), ignore_index=True)
annotation_df.to_csv(os.path.join('annotations', 'all.csv'), index=False)

problematic_csv_files = glob('annotations/all_problem_structures_*csv')
problems_df = pd.concat((pd.read_csv(csv_file) for csv_file in tqdm(problematic_csv_files)), ignore_index=True)
problems_df.to_csv(os.path.join('annotations', 'all_problems.csv'), index=False)

resulted_pdb_ids = sorted(set(annotation_df['pdb_id'].to_list() + problems_df['pdb_id'].to_list()))

if original_pdb_ids == resulted_pdb_ids:
    print(f'All input structures were proccessed. {len(problems_df)} errors were wrote to annotations/all_problems.csv')
else:
    print(
    f'{len(original_pdb_ids) - len(resulted_pdb_ids)} structures were not proccessed. {len(problems_df)} errors were wrote to annotations/all_problems.csv')
