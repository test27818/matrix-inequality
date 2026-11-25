#!/bin/sh
set -eu
cd "$(dirname "$0")"
export PYTHONDONTWRITEBYTECODE=1
python3 deep_continuum_bounds.py
python3 deep_continuum_layer.py
python3 deep_continuum_frequency.py
