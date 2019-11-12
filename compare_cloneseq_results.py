#!/usr/bin/env python
"""Compare the output of the cloneseq pipeline for different run (and
presumably versions of the code).
"""

import pathlib
import hashlib
import pandas as pd

from sys import stdout
from functools import partial

def force_pathlib_path(func):
    def wrapped_func(path, *args, **kwargs):
        pathlib_path = _force_pathlib_path(path)
        return func(pathlib_path, *args, **kwargs)
    return wrapped_func


def _force_pathlib_path(path):
    if isinstance(path, pathlib.Path):
        return path
    else:
        return pathlib.Path(path)


@force_pathlib_path
def _get_paths(topdir, include, exclude):
    paths = topdir.glob(include)
    if exclude:
        paths = set(paths) - set(topdir.glob(exclude))
    return paths


def _get_relative_path(parent, child):
    relative_path = child.relative_to(parent)
    if isinstance(relative_path, pathlib.PosixPath):
        return relative_path
    else:
        return relative_path.as_posix_path()


def compare_files(file1, file2, func):
    return func(file1) == func(file2)


@force_pathlib_path
def get_size(file_path):
    return file_path.stat().st_size


@force_pathlib_path
def get_total_lines(file_path):
    return len(file_path.open().readlines())


@force_pathlib_path
def get_hash(file_path):
    return hashlib.md5(file_path.read_bytes()).hexdigest()



def compare_directories(dir1, dir2,
                        include='*', exclude=None,
                        filter_paths_with=None,
                        compare_files_with=get_hash,
                        verbose=True):
    """
    TODO
    """
    dir1 = _force_pathlib_path(dir1)
    dir2 = _force_pathlib_path(dir2)

    paths1 = _get_paths(dir1, include, exclude)
    paths2 = _get_paths(dir2, include, exclude)

    if filter_paths_with:
        paths1 = filter_paths_with(paths1)
        paths2 = filter_paths_with(paths2)

    relative_paths1 = [_get_relative_path(dir1, path) for path in paths1]
    relative_paths2 = [_get_relative_path(dir2, path) for path in paths2]

    matches, mismatches = [], []
    potential_matches = set(relative_paths1) & set(relative_paths2) # intersection
    for ii, candidate in enumerate(potential_matches):
        if verbose:
            stdout.write(f"\r    {ii+1} / of {len(potential_matches)} potential matches processed.")
            stdout.flush()
        if compare_files_with(dir1.joinpath(candidate)) == \
           compare_files_with(dir2.joinpath(candidate)):
            matches.append(candidate)
        else:
            mismatches.append(candidate)
    if verbose:
        stdout.write("\n")
        stdout.flush()

    errors = set(relative_paths1) ^ set(relative_paths2) # symmetric difference
    errors1 = [dir1.joinpath(error) for error in errors if error in relative_paths1]
    errors2 = [dir2.joinpath(error) for error in errors if error in relative_paths2]
    errors = errors1 + errors2

    return matches, mismatches, errors


def compare_dataframes(filepath1, filepath2, sep='\t', *args, **kwargs):
    df1 = pd.read_csv(filepath1, sep=sep, *args, **kwargs)
    df2 = pd.read_csv(filepath2, sep=sep, *args, **kwargs)

    if df1.equals(df2):
        return df1, None, None

    else:
        lines1 = pathlib.Path(filepath1).open().readlines()
        lines2 = pathlib.Path(filepath2).open().readlines()
        # exclude headers
        lines1 = lines1[1:]
        lines2 = lines2[1:]
        _, unique_to_file1, unique_to_file2 = \
            _compare_lines(lines1, lines2)
        indices1 = [ii for ii, line in enumerate(lines1) if line in unique_to_file1]
        indices2 = [ii for ii, line in enumerate(lines2) if line in unique_to_file2]

        unique_to_df1 = df1.iloc[indices1]
        unique_to_df2 = df2.iloc[indices2]

        return None, unique_to_df1, unique_to_df2


def _compare_lines(lines1, lines2):
    set1 = set(lines1)
    set2 = set(lines2)

    common_lines = set1 & set2 # intersection
    unique_to_file1 = set1 - set2
    unique_to_file2 = set2 - set1
    return common_lines, unique_to_file1, unique_to_file2


def exclude_files_with_extension(filepaths, extensions):
    return [filepath for filepath in filepaths \
            if pathlib.Path(filepath).suffix not in extensions]


if __name__ == '__main__':

    analysis_dir = pathlib.Path("/ifs/projects/proj093/analysis")
    top_dir1 = analysis_dir.joinpath('cloneseq.backup_20191111/mouse-Plate17')
    top_dir2 = analysis_dir.joinpath('cloneseq/mouse-Plate17')

    # check if output summary files are identical;
    # this is much cheaper than comparing all directories,
    # so we do it before we start an in-depth analysis of all files in all subdirectories
    summary1 = top_dir1.joinpath('summary.tsv')
    summary2 = top_dir2.joinpath('summary.tsv')
    common_rows, unique_to_df1, unique_to_df2 = compare_dataframes(summary1, summary2)

    if (len(unique_to_df1) == 0) and (len(unique_to_df2) == 0):
        print('The summary.tsv tables are identical.')
        print('Skipping in-depth analysis of all other files created by the pipeline upstream.')

    else:
        print(f"The summary.tsv files differ between the two folders!")
        print("Rows unique to 1st table:")
        print(unique_to_df1)
        print()
        print("Rows unique to 2nd table:")
        print(unique_to_df2)
        print()

        # iterate over subfolders and assert equality of all files in the folder;
        # stop when a file is found that differs between the folders, or
        # if that file does not exist in one of the folders
        print("Comparing subdirectories...")
        subdirs = [
            'vector.dir',
            'fastq.dir',
            'mapped.dir',
            'merged_bam.dir',
            'deduplicated.dir',
            'deduplicatedclone_codes.dir',
        ]

        exclude_log_and_err_files = partial(exclude_files_with_extension,
                                            extensions=['.log', '.err'])

        for subdir in subdirs:
            print(f"  {subdir}:")
            matches, mismatches, errors = compare_directories(
                top_dir1.joinpath(subdir),
                top_dir2.joinpath(subdir),
                filter_paths_with = exclude_log_and_err_files,
                compare_files_with = get_size,
            )
            if mismatches:
                print(f"There are mismatches in the {subdir} subdirectories:")
                for mismatch in mismatches:
                    print(mismatch)
                print("Aborting further analysis.")
                break

            if errors:
                print(f"There are errors in the {subdir} subdirectories:")
                for err in errors:
                    print(err)
                print("Aborting further analysis.")
                break
