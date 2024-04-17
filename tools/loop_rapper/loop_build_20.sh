#!/bin/bash

rapper ../params.xml model-loops \
    --pdb "1nko.pdb" \
    --map "1nko_2mFo-DFc_map.ccp4" \
    --chain-id "A" \
    --models "20" \
    --cryst-d-high "2.19" \
    --seq "GNDIS" \
    --start "69" \
    --stop "73" \
    --use-CCP4i-file-name "false" \
    --edm-fit "true" \
    --sidechain-mode "smart" \
    --runs-dir "1nko_20_A" \
    --sidechain-radius-reduction "0.75" \
    --WRITE-INDIVIDUAL-MODELS "true" \
    --enforce-strict-anchor-geometry "true" \
    --sidechain-library /home/olebedenko/tools/CCP4/ccp4-8.0/share/rapper/data/richardson.lib