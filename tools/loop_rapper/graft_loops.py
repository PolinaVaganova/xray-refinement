import gemmi
from glob import glob
import os
import pandas as pd


def graft_loop(loop_pdb, path_to_original_pdb, chain_id):
    # take every loop and open pdb files
    loop_st = gemmi.read_structure(loop_pdb)[0]
    protein_st = gemmi.read_structure(path_to_original_pdb)[0]

    # create structure for grafting
    grafted_structure = gemmi.Structure()
    grafted_model = gemmi.Model('1')
    grafted_chain = gemmi.Chain(chain_id)

    for residue in protein_st.find_chain(chain_id):
        if residue.seqid.num < int(gap_start_res_id) and residue.name != 'HOH':
            grafted_chain.add_residue(residue)

    for residue in loop_st.find_chain(chain_id):
        grafted_chain.add_residue(residue)

    for residue in protein_st.find_chain(chain_id):
        if residue.seqid.num > int(gap_end_res_id) and residue.name != 'HOH':
            grafted_chain.add_residue(residue)

    grafted_model.add_chain(grafted_chain)
    grafted_structure.add_model(grafted_model)

    return grafted_structure


def merge_grafted_pdb(path_to_grafted_pdb, annotation_df, pdb_id):
    chain_ids = annotation_df.loc[annotation_df['pdb_id'] == pdb_id]['chain_ids'].values[0].split('_')
    pdb_paths_chain1 = glob(os.path.join(path_to_grafted_pdb, f'*{chain_ids[0]}.pdb'))
    pdb_paths_chain2 = glob(os.path.join(path_to_grafted_pdb, f'*{chain_ids[1]}.pdb'))

    for idx1, pdb_chain1 in enumerate(pdb_paths_chain1):
        for idx2, pdb_chain2 in enumerate(pdb_paths_chain2):
            if idx1 == idx2:
                with open(pdb_chain1, 'r') as chain1:
                    chain1_contents = chain1.read()

                with open(pdb_chain2, 'r') as chain1:
                    chain2_contents = chain1.read()

                merged_pdbs = chain1_contents + '\n' + chain2_contents

                with open(os.path.join(path_to_grafted_pdb, f'{pdb_id}_rapper_{idx1 + 1:03d}_merged.pdb'),
                          'w') as output_file:
                    output_file.write(merged_pdbs)


if __name__ == "__main__":
    # read inputs
    path_to_out_dir = 'grafted_loop'
    pdb_id = '1nko'
    path_to_annotation_df = '/home/polina/xray_refinment/1_filter_pdb/all_filtered_supercell.csv'
    path_to_original_pdb = f'/home/polina/xray_refinment/2_refinment/rcsb_data/pdb/{pdb_id}.pdb'

    # open df with data
    annotation_df = pd.read_csv(path_to_annotation_df, sep=',')

    # iterate over gaps in data
    gap_num = annotation_df[annotation_df['pdb_id'] == pdb_id]['gap_num'].astype(int).values[0]
    for gap_idx in range(gap_num):
        gap_start_end_res_ids = \
            annotation_df[annotation_df['pdb_id'] == pdb_id]['gap_start_end_ids'].values[0].split('_')[
                gap_idx]
        gap_start_res_id = gap_start_end_res_ids.split('-')[0]
        gap_end_res_id = gap_start_end_res_ids.split('-')[1]

        # get chain id
        chain_id = annotation_df.loc[annotation_df['pdb_id'] == pdb_id]['chain_ids'].values[0].split('_')[
            gap_idx]

        paths_to_loops = glob(os.path.join(f'{pdb_id}_20_{chain_id}', 'looptest-*.pdb'))

        for idx, loop_pdb in enumerate(paths_to_loops):
            grafted_structure = graft_loop(loop_pdb, path_to_original_pdb, chain_id)
            os.makedirs(path_to_out_dir, exist_ok=True)
            grafted_structure.write_pdb(
                os.path.join(path_to_out_dir, f'{pdb_id}_rapper_{idx + 1:03d}_{chain_id}.pdb'),
                gemmi.PdbWriteOptions(cryst1_record=False))
            
    if gap_num > 1:
        merge_grafted_pdb(path_to_out_dir, annotation_df, pdb_id)
