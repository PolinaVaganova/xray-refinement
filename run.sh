#!/bin/bash


source venv/bin/activate

export AMBERHOME="/opt/amber/"
export PATH="$AMBERHOME/bin:$PATH"

python run_locally.py
