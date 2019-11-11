#!/usr/bin/env python
"""
Automate the setup of cloneseq analysis.
Creates the necessary directory structure,
and links in the appropriate files into the data subdirectory.
"""

import re
import shutil
import pathlib
import textwrap
import pandas as pd


def setup_folders(project, pipeline, plate):
    """Creates the following folder structure for RNASeq data analysis:

    project
        /analysis
            /pipeline
                /plate
                    /data

    and returns a dictionary that maps each element to a pathlib Path object.

    Arguments:
    ----------
    project : str
         path/to/project
    pipeline : str
         pipeline name
    plate : str
        plate name

    Returns:
    --------
    paths : dict str : pathlib.Path
        Mapping of folder names to paths.
    """
    # construct all paths
    project_dir = pathlib.Path(project)
    pipeline_dir = project_dir.joinpath(f'analysis/{pipeline}')
    plate_dir = pipeline_dir.joinpath(plate)
    data_dir = plate_dir.joinpath('data')

    # make all directories
    data_dir.mkdir(parents=True, exist_ok=True)

    return dict(project=project_dir,
                pipeline=pipeline_dir,
                plate=plate_dir,
                data=data_dir
    )


def link_in_data(src_data_dir, dst_data_dir, mapping_info):
    src_name_to_dst_name = remap_file_names(mapping_info)
    create_symlinks(src_data_dir, dst_data_dir, src_name_to_dst_name)


def remap_file_names(path_to_email:str, lane:str='lane1') -> dict:
    """Map sample IDs to human readable file names.

    The file names assigned by the sequencing facility are not human
    readable. Based on the data provided in the email from the facility,
    we create a mapping from the sample ID to a human readable name with
    the following format:

    <sample ID>_<strand>.fastq.gz : <plate>-<well>-<lane>.fastq.<strand>.gz

    Arguments:
    ----------
    path_to_email : str
        Path to the email sent by the sequencing facility.
        For a sample text, see parse_email.
    lane : str (default lane1)
        The sequencing lane. This information is not contained in the email,
        but should be consistent for all files within one email as each email
        maps to a single batch.

    Returns:
    --------
    mapping : dict
        Mapping from sample IDs to human readable file names.

    TODO:
    -----
    - Handle non-stranded sequencing data.
    """

    dataframe = parse_email(path_to_email)

    mapping = dict()
    problems = []
    for ii, row in dataframe.iterrows():
        if (row['plate'] != '') and (row['well'] != ''):
            for strand in [1, 2]:
                src = f"{row['sample_id']}_{strand}.fastq.gz"
                dst = f"{row['plate']}-{row['well']}-{lane}.fastq.{strand}.gz"
                mapping[src] = dst
        else:
            problems.append(f"{row['sample_id']} ({row['sample_name']})")

    if problems:
        msg = 'Could not map the following samples:\n'
        msg += '\n'.join(problems)
        import warnings
        warnings.warn(msg)

    return mapping


def parse_email(file_path):
    """
    Given an email with the following format:

    WTCHG_693615_70015361 Sample:Plate15_A10 Count:12577357
    WTCHG_693615_70025362 Sample:Plate15_B10 Count:4835494
    WTCHG_693615_70105370 Sample:Plate15_B11 Count:7296271
    WTCHG_693615_70115371 Sample:Plate15_C11_lysis_ctl Count:732017
    WTCHG_693615_70125372 Sample:Plate15_D11 Count:1453627
    WTCHG_693615_70195379 Sample:Plate15_C12 Count:1842041
    WTCHG_693615_70205380 Sample:Plate15_D12_ceph10pgRNA_ctl Count:756343
    WTCHG_693615_70215381 Sample:Plate15_E12 Count:65603
    WTCHG_693615_71345014 Sample:single_GM12878_ctl_14 Count:776734

    return a pandas dataframe with the following format:

    sample_id,             sample_name,                 plate,   well, count
    WTCHG_693615_70015361, Plate15_A10,                 Plate15, A10,  12577357
    WTCHG_693615_70025362, Plate15_B10,                 Plate15, B10,  4835494
    WTCHG_693615_70105370, Plate15_B11,                 Plate15, B11,  7296271
    WTCHG_693615_70115371, Plate15_C11_lysis_ctl,       Plate15, C11,  732017
    WTCHG_693615_70125372, Plate15_D11,                 Plate15, D11,  1453627
    WTCHG_693615_70195379, Plate15_C12,                 Plate15, C12,  1842041
    WTCHG_693615_70205380, Plate15_D12_ceph10pgRNA_ctl, Plate15, D12,  756343
    WTCHG_693615_70215381, Plate15_E12,                 Plate15, E12,  65603
    WTCHG_693615_71345014, single_GM12878_ctl_14,              ,    ,  776734
    """

    with open(file_path) as f:
        lines = f.readlines()

    # remove all lines that do not appear to contain a sample ID;
    # this should remove
    # a) empty lines
    # b) the header
    lines = [line for line in lines if 'WTCHG' in line]

    # parse remaining lines using regex
    data = []
    problems = []
    for line in lines:
        sample_id    = re.search('^WTCHG_[0-9]+_[0-9]+', line).group(0)
        sample_name  = re.search('Sample:([A-z0-9_]+)', line).group(1)
        count_as_str = re.search('Count:([0-9]+)', line).group(1)

        # Not all sample names follow the <plate>_<well><rest> format,
        # as, for example, shown in the last line of the sample text.
        # re.search may thus return None, which we need to handle.
        match = re.search('(Plate[0-9]{1,2})_([A-Z][0-9]{1,2})', line)
        if match:
            plate, well = match.groups()
        else:
            plate, well = '',''
            problems.append(line)

        data.append((sample_id, sample_name, plate, well, int(count_as_str)))

    if problems:
        msg = "Could not determine plate and/or well assignement in the following lines:\n"
        msg += ''.join(problems)
        import warnings
        warnings.warn(msg)

    return pd.DataFrame(data,
                        columns=['sample_id',
                                 'sample_name',
                                 'plate',
                                 'well',
                                 'count'])


def create_symlinks(src_data_dir, dst_data_dir, src_name_to_dst_name):
    """Create symlinks in the project working directory to the raw sequencing data in the backup folder."""
    assert src_data_dir.exists(), "Source directory does not exist."
    assert dst_data_dir.exists(), "Target directory does not exist."
    for src_name, dst_name in src_name_to_dst_name.items():
        src = src_data_dir.joinpath(src_name)
        assert src.exists(), f"File {src} does not exist."
        dst = dst_data_dir.joinpath(dst_name)
        # shutil.copy(src, dst, follow_symlinks=False)
        dst.symlink_to(src)


def create_cloneseq_yml(plate_dir):
    pipeline_yml_text = f"""
    kmer_filtering:
        kmer_size: 10
        min_kmer_matches: 20

    filename_vector_fasta: /ifs/projects/proj093/backup/vector/vector.fa

    input_fastq_glob: {plate_dir}/data/*.gz
    """
    pipeline_yml_text = textwrap.dedent(pipeline_yml_text)
    pipeline_yml = plate_dir.joinpath('pipeline.yml')
    pipeline_yml.write_text(pipeline_yml_text)


if __name__ == '__main__':

    # plate_name = 'mouse-Plate15'
    # path_to_email = '/ifs/projects/proj093/backup/2019-07-15-plate15/email_plate_15.txt'
    # path_to_data = '/ifs/projects/proj093/backup/2019-07-15-plate15/190626_K00181_0151_AH7MTJBBXY/'

    # plate_name = 'mouse-Plate16'
    # path_to_email = '/ifs/projects/proj093/backup/2019_09_19_plates_16_and_17/email_plate_16.txt'
    # path_to_data = '/ifs/projects/proj093/backup/2019_09_19_plates_16_and_17/190913_K00181_0175_BHF27LBBXY/'

    # plate_name = 'mouse-Plate17'
    # path_to_email = '/ifs/projects/proj093/backup/2019_09_19_plates_16_and_17/email_plate_17.txt'
    # path_to_data = '/ifs/projects/proj093/backup/2019_09_19_plates_16_and_17/190913_K00181_0175_BHF27LBBXY/'

    # plate_name = 'mouse-Plate18'
    # path_to_email = '/ifs/projects/proj093/backup/2019_09_19_plate_18/email_plate_18.txt'
    # path_to_data = '/ifs/projects/proj093/backup/2019_09_19_plate_18/190913_K00181_0175_BHF27LBBXY/'

    data_info = pd.read_csv('/ifs/projects/proj093/analysis/data_info.csv')
    for ii, dataset in data_info.iterrows():
        if dataset['processed'] is not True:

            plate_name = dataset['plate_name']
            path_to_data = dataset['path_to_data']
            path_to_email = dataset['path_to_email']

            project_paths = setup_folders(
                project='/ifs/projects/proj093/',
                pipeline='cloneseq',
                plate=plate_name,
            )

            # gather and link data sets
            link_in_data(src_data_dir=pathlib.Path(path_to_data),
                         dst_data_dir=project_paths['data'],
                         mapping_info=path_to_email)

            # create pipeline.yml file
            create_cloneseq_yml(project_paths['plate'])
