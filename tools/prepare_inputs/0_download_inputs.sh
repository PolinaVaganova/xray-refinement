#!/bin/bash

csv_file=/home/polina/xray_refinment/1_filter_pdb/run_16.04_1nko.csv
run_name=1un2
error_filename="$run_name.failed.list"
touch $error_filename


# prepare cif files
mkdir -p cif
cd cif
tail -n +2 $csv_file | cut -d ',' -f 1 | while IFS= read -r pdb_id; do
    if ! [ -f "$pdb_id-sf.cif" ] && ! [ -f "$pdb_id.cif" ]; then
        echo "Download CIF file.."
        if wget "https://files.rcsb.org/download/$pdb_id-sf.cif"; then
            echo "CIF file downloaded"
        elif wget "https://files.rcsb.org/download/$pdb_id.cif"; then
            echo "CIF file downloaded"
        else
            echo "Failed download cif $pdb_id" >> "../$error_filename"
        fi
    fi    
done 


# prepare mtz files
mkdir -p ../mtz
cd ../mtz
tail -n +2 $csv_file | cut -d ',' -f 1 | while IFS= read -r pdb_id; do
    if ! [ -f "${pdb_id}.mtz" ]; then
        if [ -f "../cif/${pdb_id}-sf.cif" ]; then
            if phenix.cif_as_mtz "../cif/${pdb_id}-sf.cif" --ignore_bad_sigmas --merge --output_file_name="${pdb_id}-sf.mtz"; then
                echo "Turm to sf-mtz"
            elif phenix.cif_as_mtz "../cif/${pdb_id}.cif" --ignore_bad_sigmas --merge --output_file_name="${pdb_id}-sf.mtz"; then
                echo "Turm to sf-mtz"
            else
                echo "Failed turn to sf-mtz $pdb_id" >> "../$error_filename"
            fi
        fi
        if [ -f "${pdb_id}-sf.mtz" ]; then
            if phenix.reflection_file_converter --expand_to_p1 "${pdb_id}-sf.mtz" --write_mtz_amplitudes --mtz_root_label="FOBS" --label="FOBS" --generate_r_free_flags --non_anomalous --mtz "${pdb_id}.mtz"; then
                echo "Turn to mtz"
            elif phenix.reflection_file_converter --expand_to_p1 "${pdb_id}-sf.mtz" --write_mtz_amplitudes --mtz_root_label="IOBS" --label="IOBS" --generate_r_free_flags --non_anomalous --mtz "${pdb_id}.mtz"; then
                echo "Turn to mtz"
            else
                echo "Failed turn to mtz $pdb_id" >> "../$error_filename"
            fi
        fi
    fi
done


# prepare pdb files
mkdir -p ../pdb
cd ../pdb
tail -n +2 $csv_file | cut -d ',' -f 1 | while IFS= read -r pdb_id; do
    if ! [ -f "${pdb_id}.pdb" ]; then
        if gunzip -c "/home/olebedenko/PDB/PDB_2023/structures/all/pdb/pdb${pdb_id}.ent.gz" > "${pdb_id}.pdb"; then
            echo "pdb unpacked"
        else
            echo "Failed unpack pdb $pdb_id" >> "../$error_filename"
        fi
    fi
done


# prepare fasta files
mkdir -p ../fasta
cd ../fasta
tail -n +2 $csv_file | cut -d ',' -f 1 | while IFS= read -r pdb_id; do
    if ! [ -f "$pdb_id" ]; then
        if wget "https://www.rcsb.org/fasta/entry/$pdb_id"; then
            echo "fasta file downloaded"
        else
            echo "Failed download fasta $pdb_id" >> "../$error_filename"
        fi
    fi
done
