import os
import gemmi
from typing import Union
from hybrid_36 import hy36encode
from glob import glob
import pandas as pd
from tqdm import tqdm
import sys

AnyPath = Union[str, bytes, os.PathLike]

def renumber_residues(st: gemmi.Structure) -> gemmi.Structure:
    chain_renum = 1
    first_res_num = None

    result = gemmi.Structure()
    for model in st:
        new_model = gemmi.Model(model.name)
        for chain in model:
            new_chain = gemmi.Chain(hy36encode(2, chain_renum))
            chain_renum += 1
            for residue in chain:
                if first_res_num is None:
                    first_res_num = residue.seqid.num
                residue.seqid.num = residue.seqid.num - first_res_num + 1
                new_chain.add_residue(residue)

            new_model.add_chain(new_chain)
        result.add_model(new_model)
    result.cell = st.cell
    result.spacegroup_hm = st.spacegroup_hm
    return result

def read_pdb(path: AnyPath) -> gemmi.Structure:
    return gemmi.read_pdb(str(path), split_chain_on_ter=True)


def copy_b_factor(structure,
                  original_b_factor_dict):
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    if original_b_factor_dict.get(residue.seqid.num):
                        if original_b_factor_dict[residue.seqid.num].get(atom.name):
                            atom.b_iso = original_b_factor_dict[residue.seqid.num][atom.name]
                        else:
                            atom.b_iso = -1
                    else:
                        atom.b_iso = -1
    return structure


if __name__ == '__main__':
    # set correct paths
    option = 'modeller'
    # path_to_csv = '/home/polina/xray_refinment/1_filter_pdb/rapper_1u36.csv'
    path_to_original_pdbs = '/home/polina/xray_refinment/2_refinment/rcsb_data/pdb'
    data_input_dir = '/home/polina/xray_refinment/2_refinment/xray-refinement/data/input/'
    path_to_modelled_loops = os.path.join('/home/polina/xray_refinment/2_refinment/loop_modelling/', option)
    path_to_original_mtz = '/home/polina/xray_refinment/2_refinment/rcsb_data/mtz'
    pdb_ids = ['1K33_1|Chain.B99990009', '1U3Y_1|Chain.B99990028', '1NXX_1|Chain.B99990074']

    # pdb_df = pd.read_csv(path_to_csv, sep=',')
    # pdb_ids = pdb_df['pdb_id']

    for pdb_id in pdb_ids:
        # set correct paths
        # path_to_out_dir = os.path.join(data_input_dir, f'{pdb_id}_{option}')
        path_to_out_dir = os.path.join(data_input_dir, f'{pdb_id}')

        if not os.path.exists(path_to_out_dir):
            os.makedirs(path_to_out_dir)


        # if option == 'rapper':
        #     path_to_modelled_loops = os.path.join(path_to_modelled_loops, pdb_id, f'{pdb_id}_20')
        #     path_to_modelled_pdbs = glob(os.path.join(path_to_modelled_loops, f'{pdb_id}_rapper_*.pdb'))

        # elif option == 'modeller':
        #     # fix this path later
        #     path_to_modelled_pdbs = glob(os.path.join(path_to_modelled_loops, '*Chain*.pdb'))

        # print(len(path_to_modelled_pdbs))
        # path_to_modelled_pdbs.sort()

        pdb_code = pdb_id[:4].lower()
        path_to_modelled_pdbs = (os.path.join(path_to_modelled_loops, pdb_code, f'{pdb_id}.pdb'))
        path_to_original_pdb = os.path.join(path_to_original_pdbs, f'{pdb_code}.pdb')

        st_original = read_pdb(path_to_original_pdb)
        st_original_renumber = renumber_residues(st_original)

        st_original_entities = st_original_renumber.entities

        # add headers
        original_header = st_original_renumber.make_pdb_headers()

        original_b_factor_dict = {}

        # for modeller option structure renumbering is needed
        for model in st_original_renumber:
            for chain in model:
                for residue in chain:
                    original_b_factor_dict[residue.seqid.num] = {}
                    for atom in residue:
                        original_b_factor_dict[residue.seqid.num].update({atom.name: atom.b_iso})


        # for idx, path_to_modelled_pdb in tqdm(enumerate(path_to_modelled_pdbs)):
        path_to_fout = os.path.join(path_to_out_dir, f'{pdb_id}.pdb')
        # print(path_to_modelled_pdbs)
        # # exit()
        st_modelled = read_pdb(path_to_modelled_pdbs)
        st_modelled_b_factors = copy_b_factor(st_modelled, original_b_factor_dict)
        # st_modelled_b_factors.add_entity_types(overwrite=True)
        st_modelled.entities = st_original_renumber.entities
        st_modelled_b_factors_str = st_modelled_b_factors.make_pdb_string(gemmi.PdbWriteOptions(cryst1_record=False, end_record=True, ter_records=True))
        with open(path_to_fout, "w") as fout:
            fout.writelines(original_header)
            fout.writelines(st_modelled_b_factors_str)

        if not os.path.islink(os.path.join(path_to_out_dir, f'{pdb_code}.mtz')):
            os.symlink(os.path.join(path_to_original_mtz, f'{pdb_code}.mtz'), os.path.join(path_to_out_dir, f'{pdb_code}.mtz'))