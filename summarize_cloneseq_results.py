#!/usr/bin/env python
"""
Summarize pipeline output across data sets.

# TODO:
- handle seperator in _convert_filepaths_to_dataframes
"""

import re
import pathlib
import numpy as np
import pandas as pd
import StyleFrame

from functools import partial


def convert_filepaths_to_dataframes(func):
    """Given a potentially mixed list of file paths or dataframes, replace
    filepaths by their corresponing dataframes before passing the list
    to the function.
    """
    def wrapped_func(dataframes_or_filenames, *args, **kwargs):
        dataframes = [_convert_filepath_to_dataframe(obj) for obj in dataframes_or_filenames]
        return func(dataframes, *args, **kwargs)
    return wrapped_func


def convert_filepath_to_dataframe(func):
    """Given an object that is either a file path or dataframe, replace
    the filepath by its corresponing dataframe before passing the object
    to the function.
    """
    def wrapped_func(obj, *args, **kwargs):
        dataframe = _convert_filepath_to_dataframe(obj)
        return func(dataframe, *args, **kwargs)
    return wrapped_func


def _convert_filepath_to_dataframe(obj, *args, **kwargs):
    if isinstance(obj, pd.DataFrame):
        return obj
    elif isinstance(obj, (str, pathlib.Path)):
        return pd.read_csv(obj, *args, **kwargs)

_convert_filepath_to_dataframe = partial(_convert_filepath_to_dataframe, sep='\t')


@convert_filepaths_to_dataframes
def concatenate_tables(dataframes):
    return pd.concat(dataframes, axis=0, ignore_index=True, sort=False)


@convert_filepaths_to_dataframes
def write_to_excel_workbook(dataframes,
                            output_filepath,
                            sheet_names=None,
                            *args, **kwargs):
    """Save multiple dataframes as individual sheets in a single excel workbook.
    """
    if not sheet_names:
        sheet_names = [f'Sheet{ii}' for ii, _ in enumerate(dataframes)]

    style = StyleFrame.Styler(wrap_text=False,
                              horizontal_alignment='general',
                              font='Courier',
    )
    with StyleFrame.StyleFrame.ExcelWriter(output_filepath) as writer:
        for df, sheet_name in zip(dataframes, sheet_names):
            columns = tuple(df.columns.values)
            sf = StyleFrame.StyleFrame(df)
            sf.apply_column_style(columns, style)
            sf.to_excel(writer, sheet_name=sheet_name, best_fit=columns)


def sort_strings_ending_in_numbers(strings):
    """Sort a list strings ending in numbers.

    ['bcd3', 'abc1', 'abc11', abc2'] -> ['abc1', 'abc2', 'abc11', 'bcd3']

    Adapted from:
    -------------
    # stackoverflow.com/a/4318182/2912349
    """
    def key(mixed):
        string, num = re.search(r"^(\D+)(\d+)", mixed).groups()
        return string, int(num)
    return sorted(strings, key=key)


def argsort_strings_ending_in_numbers(strings):
    """Return the indices that would sort a given list of strings ending in numbers.

    unsorted = ['bcd3', 'abc1', 'abc11', abc2']
    sorted = [unsorted[idx] for idx in argsort(unsorted)]
    print(sorted) # ['abc1', 'abc2', 'abc11', 'bcd3']

    Adapted from:
    -------------
    # stackoverflow.com/a/4318182/2912349
    # stackoverflow.com/a/6979121/2912349
    """

    def key(mixed):
        string, num = re.search(r"^(\D+)(\d+)", mixed).groups()
        return string, int(num)
    return [x for x, y in sorted(enumerate(strings), key = lambda z: key(z[1]))]


@convert_filepath_to_dataframe
def get_barcode_recovery_statistics(dataframe):
    recovered = np.array([is_barcode(value, 12) for value in dataframe['barcode']])
    total_recovered = np.sum(recovered == 1)
    total_partials = np.sum(recovered == 0.5)
    total_missing = np.sum(recovered == 0.)
    total_samples = len(dataframe)
    recovery_rate = (total_recovered + total_partials) / total_samples
    return total_recovered, total_partials, total_missing, total_samples, recovery_rate


def is_barcode(candidate, expected_length):
    """
    Classifies candidate barcode sequences into

    - barcodes (1.0): sequences containing only ATGC bases and having the right length
    - partials (0.5): sequences containing some ATGC bases, and not necessarily being of the right length
    - nonsense (0.0): sequences containing only N bases; length is irrelevant
    """
    if np.all([char in 'ATGC' for char in candidate]) and len(candidate) == expected_length:
        return 1.
    elif np.all([char == 'N' for char in candidate]):
        return 0.
    else:
        import warnings
        warnings.warn(f'Barcode is only a partial match: {candidate}.')
        return 0.5


if __name__ == '__main__':

    # --------------------------------------------------------------------------------
    # concatenate all cloneseq summary.tsv tables into one excel workbook

    project_dir = pathlib.Path('/ifs/projects/proj093/analysis/cloneseq/')
    filepaths = list(project_dir.glob('*/summary.tsv'))
    plate_names = [p.parent.name for p in filepaths]

    # sort plates by plate name and number;
    # plate names are assumed to follow the format <organism>-Plate<number>, e.g. mouse-Plate15
    order = argsort_strings_ending_in_numbers(plate_names)
    filepaths = [filepaths[idx] for idx in order]
    plate_names = [plate_names[idx] for idx in order]

    # add a table that is a concatenation of all other tables
    summary = concatenate_tables(filepaths)
    filepaths.append(summary)
    plate_names.append('summary')

    write_to_excel_workbook(filepaths,
                            project_dir.joinpath('summary.xlsx'),
                            sheet_names=plate_names)

    # --------------------------------------------------------------------------------
    # compute barcode recovery rate per table

    total_recovered, total_partials, total_missing, total_samples, recovery_rate = \
        zip(*[get_barcode_recovery_statistics(filepath) for filepath in filepaths])

    df = pd.DataFrame(dict(
        plate           = plate_names,
        total_recovered = total_recovered,
        total_partials  = total_partials,
        total_missing   = total_missing,
        total_samples   = total_samples,
        recovery_rate   = recovery_rate,
    ))

    df.to_csv(project_dir.joinpath('barcode_recovery_statistics.xlsx'))
    print(df)
