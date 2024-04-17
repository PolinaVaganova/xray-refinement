import os
from annotation_utils import Monomer
import pandas as pd
from multiprocessing import Process, Queue
import glob

path_to_pdbs = '/home/olebedenko/PDB/PDB_2023/structures/all/pdb/'
# path_to_monomer_csv = '/home/olebedenko/PDB/PDB_2023/all/annotation/data/protein_monomer/prot.csv'
path_to_crystal_csv = 'crystal.idx'


def get_annotation(paths_to_pdbs, result_queue):
    annotation_df = pd.DataFrame(
        columns=['pdb_id', 'chain_num', 'chain_ids', 'gap_num', 'gap_len', 'gap_start_end_ids', 'is_gap_in_ss',
                 'is_gap_on_terminal_end', 'hetatm', 'unk', 'a', 'b', 'c', 'alpha', 'beta',
                 'gamma', 'sp.gp.', 'Z'])
    problem_structures = pd.DataFrame(columns=['pdb_id', 'exception'])

    # open index data
    colnames = ['pdb_id', 'CRYST1', 'a', 'b', 'c', 'alpha', 'beta', 'gamma', 'sp.gp.', 'Z']
    crystal_df = pd.read_csv(path_to_crystal_csv, skiprows=4, sep=' \s{1,}|\t', header=None, names=colnames)
    crystal_df['pdb_id'] = crystal_df['pdb_id'].str.lower()
    crystal_df = crystal_df.drop(columns=['CRYST1'])

    # annotate monomers
    for idx, path_to_pdb in enumerate(paths_to_pdbs):
        try:
            st = Monomer(path_to_pdb)
            pdb_id = os.path.basename(path_to_pdb).split('.')[0][3:]
            hetatm = st.count_hetatm()
            unk = st.count_unk()
            st = st.clean()
            # get information about gaps
            gaps_start_end_ids = st.get_gaps_start_end()
            gap_num = len(gaps_start_end_ids)

            # chains number
            chain_num = st.count_chains()

            # chaind ids
            chain_ids = st.get_chain_ids()


            # find is gap rids in ss structure and write gap length
            is_gap_in_ss = []
            ss_rids = st.get_ss_rids()
            gap_lengths = []
            is_gap_on_terminal_end = []
            for gap in gaps_start_end_ids:
                gap_rids = [rid for rid in range(gap[0], gap[1] + 1)]
                if set(gap_rids) & set(ss_rids):
                    is_gap_in_ss.append(True)
                else:
                    is_gap_in_ss.append(False)
                gap_lengths.append(gap[1] - gap[0] + 1)

                # find out is gap between ss or not
                min_ss_rid = min(ss_rids)
                max_ss_rid = max(ss_rids)

                if gap_rids[0] < min_ss_rid or gap_rids[0] > max_ss_rid:
                    is_gap_on_terminal_end.append(True)
                else:
                    is_gap_on_terminal_end.append(False)

            # get info from crystal_df
            pdb_id_crystal = crystal_df[crystal_df['pdb_id'] == pdb_id]
            a, b, c, alpha, beta, gamma, sp_gp, Z = pdb_id_crystal['a'].values[0], pdb_id_crystal['b'].values[0], \
                pdb_id_crystal['c'].values[0], pdb_id_crystal[
                'alpha'].values[0], pdb_id_crystal['beta'].values[0], pdb_id_crystal[
                'gamma'].values[0], pdb_id_crystal['sp.gp.'].values[0], \
                pdb_id_crystal['Z'].values[0]


            # write results to dataframe
            gap_ids_str = "_".join("-".join(map(str, ids)) for ids in gaps_start_end_ids)
            is_gap_in_ss_str = "_".join(map(str, is_gap_in_ss))
            is_gap_on_terminal_end = "_".join(map(str, is_gap_on_terminal_end))
            gap_lengths_str = "_".join(map(str, gap_lengths))
            chain_ids_str = "_".join(map(str, chain_ids))

            annotation_df.loc[idx] = [pdb_id, chain_num, chain_ids_str, gap_num, gap_lengths_str,
                                      gap_ids_str, is_gap_in_ss_str, is_gap_on_terminal_end, hetatm, unk, a, b, c,
                                      alpha, beta, gamma, sp_gp, Z]

        except Exception as e:
            problem_structures.loc[idx] = [pdb_id, e]

    result_queue.put((annotation_df, problem_structures))


if __name__ == "__main__":
    paths_to_pdbs = glob.glob(os.path.join(path_to_pdbs, 'pdb*.ent.gz'))

    num_chunks = 20
    chunk_size = (len(paths_to_pdbs) // num_chunks) + 1

    chunks = [paths_to_pdbs[i:(i + chunk_size)] for i in range(0, len(paths_to_pdbs), chunk_size)]

    result_queue = Queue()

    processes = []

    for chunk in chunks:
        p = Process(target=get_annotation, args=(chunk, result_queue))
        processes.append(p)
        p.start()

    for p in processes:
        annotation_df, problem_structures = result_queue.get()
        annotation_df.to_csv(f'annotations/all_data_{p.pid}.csv', index=False)
        problem_structures.to_csv(f'annotations/all_problem_structures_{p.pid}.csv', index=False)

    for p in processes:
        p.join()
