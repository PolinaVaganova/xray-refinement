#!/bin/bash


source venv/bin/activate
source /opt/amber/amber.sh
export PYTHONPATH=$PYTHONPATH:$(pwd)

python init.py
