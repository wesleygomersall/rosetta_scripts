#!/bin/bash

source ~/miniforge3/etc/profile.d/conda.sh # for conda to work
conda activate pdbfixer

pdbfixer \
    /projects/parisahlab/wesg/af3_data/20250920_AF3_output/comdn230-c16-13/comdn230-c16-13_model.pdb \
    --output=/projects/parisahlab/wesg/af3_data/20250920_AF3_output/comdn230-c16-13/comdn230-c16-13_model_fixed.pdb

pdbfixer \
    /projects/parisahlab/wesg/af3_data/20250920_AF3_output/comdn230-c16-14/comdn230-c16-14_model.pdb \
    --output=/projects/parisahlab/wesg/af3_data/20250920_AF3_output/comdn230-c16-14/comdn230-c16-14_model_fixed.pdb

pdbfixer \
    /projects/parisahlab/wesg/af3_data/20250920_AF3_output/comdn230-c16-15/comdn230-c16-15_model.pdb \
    --output=/projects/parisahlab/wesg/af3_data/20250920_AF3_output/comdn230-c16-15/comdn230-c16-15_model_fixed.pdb

pdbfixer \
    /projects/parisahlab/wesg/af3_data/20250920_AF3_output/comdn230-c16-16/comdn230-c16-16_model.pdb \
    --output=/projects/parisahlab/wesg/af3_data/20250920_AF3_output/comdn230-c16-16/comdn230-c16-16_model_fixed.pdb

pdbfixer \
    /projects/parisahlab/wesg/af3_data/20250920_AF3_output/comdn230-c16-17/comdn230-c16-17_model.pdb \
    --output=/projects/parisahlab/wesg/af3_data/20250920_AF3_output/comdn230-c16-17/comdn230-c16-17_model_fixed.pdb

pdbfixer \
    /projects/parisahlab/wesg/af3_data/20250920_AF3_output/comdn230-c16-18/comdn230-c16-18_model.pdb \
    --output=/projects/parisahlab/wesg/af3_data/20250920_AF3_output/comdn230-c16-18/comdn230-c16-18_model_fixed.pdb

pdbfixer \
    /projects/parisahlab/wesg/af3_data/20250920_AF3_output/comdn230-c16/comdn230-c16_model.pdb \
    --output=/projects/parisahlab/wesg/af3_data/20250920_AF3_output/comdn230-c16/comdn230-c16_model_fixed.pdb

