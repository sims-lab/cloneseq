#!/usr/bin/env python
"""
Summarize pipeline output across data sets.

# TODO:
- handle seperator in _convert_filepaths_to_dataframes
"""

import re
import pathlib
import pandas as pd
import StyleFrame

from functools import partial


def convert_filepaths_to_dataframes(func):
    def wrapped_func(dataframes_or_filenames, *args, **kwargs):
        dataframes = _convert_filepaths_to_dataframes(dataframes_or_filenames)
        return func(dataframes, *args, **kwargs)
    return wrapped_func


def _convert_filepaths_to_dataframes(dataframes_or_filepaths, *args, **kwargs):
    """Given a potentially mixed list of file paths or dataframes,
    read in data where necessary and return a list of dataframes."""
    output = []
    for obj in dataframes_or_filepaths:
        if isinstance(obj, pd.DataFrame):
            output.append(obj)
        elif isinstance(obj, (str, pathlib.Path)):
            df = pd.read_csv(obj, *args, **kwargs)
            output.append(df)
    return output

_convert_filepaths_to_dataframes = partial(_convert_filepaths_to_dataframes, sep='\t')


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


if __name__ == '__main__':

    # --------------------------------------------------------------------------------
    # concatenate all cloneseq summary.tsv tables into one excel workbook

    project_dir = pathlib.Path('/ifs/projects/proj093/analysis/cloneseq/')
    filepaths = list(project_dir.glob('*/summary.tsv'))
    sheet_names = [p.parent.name for p in filepaths]

    # sort sheets by plate name and number;
    # sheet names are assumed to follow the format <organism>-Plate<number>, e.g. mouse-Plate15
    order = argsort_strings_ending_in_numbers(sheet_names)
    filepaths = [filepaths[idx] for idx in order]
    sheet_names = [sheet_names[idx] for idx in order]

    # add a table that is a concatenation of all other tables
    summary = concatenate_tables(filepaths)
    filepaths.append(summary)
    sheet_names.append('summary')

    write_to_excel_workbook(filepaths,
                            project_dir.joinpath('summary.xlsx'),
                            sheet_names=sheet_names)
