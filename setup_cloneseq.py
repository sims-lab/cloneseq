#!/usr/bin/env python
"""
Automate the setup of cloneseq analysis.
Creates the necessary directory structure,
and links in the appropriate files into data subdirectory.
"""

import pathlib
import shutil
import textwrap

if __name__ == '__main__':

    # plate_name = 'human-Plate7'
    # plate_name = 'human-Plate8'
    # plate_name = 'human-Plate9'
    # plate_name = 'human-Plate10'

    # plate_name = 'mouse-Plate5'
    # plate_name = 'mouse-Plate6'
    # plate_name = 'mouse-Plate11'
    # plate_name = 'mouse-Plate12'
    # plate_name = 'mouse-Plate13'
    # plate_name = 'mouse-Plate14'

    project_dir = pathlib.Path('/ifs/projects/proj093/')
    plate_dir = project_dir.joinpath(f'analysis/cloneseq/{plate_name}')
    plate_data_dir = plate_dir.joinpath('data')
    plate_data_dir.mkdir(parents=True, exist_ok=True)

    # gather and link data sets
    data_dir = project_dir.joinpath('data')
    datasets = data_dir.glob(f'{plate_name}*.gz')

    for dataset in datasets:
        dst = plate_data_dir.joinpath(dataset.name)
        shutil.copy(dataset, dst, follow_symlinks=False)

    # create pipeline.yml file
    pipeline_yml_text = f"""
    kmer_filtering:
        kmer_size: 10
        min_kmer_matches: 20

    filename_vector_fasta: /ifs/projects/proj093/backup/vector/vector.fa

    input_fastq_glob: /ifs/projects/proj093/analysis/cloneseq/{plate_name}/data/*.gz
    """
    pipeline_yml_text = textwrap.dedent(pipeline_yml_text)
    pipeline_yml = plate_dir.joinpath('pipeline.yml')
    pipeline_yml.write_text(pipeline_yml_text)
