"""
This module:
- Imports the Sample and Probe classes from prep_progr_as4_classes.py.
- Gathers the values of command line arguments passed.
- Unpacks data from relevant files.
- Has functions defined for getting intersection and difference
in probes of given samples.
- Can handle any number of structure acronyms,
and will calculate the results of above mentioned functions
for all of the samples that were matched by the listed acronyms.
- Creates the sample objects as elements of a list,
where each stores the same number of probe objects
as the number of experiments,
where a probe has one expression value,
which is the value at the intersection of
the probe's row index and the sample's column index
from file 'MicroarrayExpression.csv'.
- Modifies the probe objects stored
in each sample object's 'probes' attribute
to contain only probes that have their expression
above the value defined by the respective command line argument.
- Prints all Sample objects created, along with their attributes.
- Prints results of functions calculating intersection and difference,
as sets of probe IDs.
"""

import sys
import os
import getopt
import math
from functools import reduce
from typing import Optional, Union, Any

import pandas as pd

from prep_progr_as4_classes import Probe, Sample

NEWL = '\n'
external_parameters = getopt.getopt(sys.argv[1:], "")

cutoff_value = int(external_parameters[1][0])
above_background: Optional[bool] = None
sample_acronyms: list[str]
try:
    above_background = bool(int(external_parameters[1][-1]))
    sample_acronyms = external_parameters[1][1:-1]
except ValueError:
    sample_acronyms = external_parameters[1][1:]

def cast_cell_value(value_: Any) -> Union[int, str, None]:
    """
    For 'chromosome' column, where values may be of multiple types,
    this function converts the cell value to
    None if it is an empty value,
    integer if it can be cast as that,
    and leaves it unchanged if it is an arbitrary string.
    """

    if isinstance(value_, float) and math.isnan(value_):
        return None
    try:
        return int(value_)
    except ValueError:
        return value_

def get_intersection_in_probes(*samples: Sample) -> set[int]:
    """
    From any number of Sample objects, returns the IDs of probes
    that were found in all Sample objects' 'probes' attribute.
    If not having performed the filtering of Probe objects per Sample
    based on minimum expression
    and whether the probe is above background,
    this function will most likely return all of the IDs.
    """

    print("Probes shared between regions:")
    return reduce(
        lambda accumulator_, next_: accumulator_.intersection(next_),
        map(set, [[
            probe_.probe_id for probe_ in sample_.probes
        ] for sample_ in samples])
    )

def get_difference_in_probes(*samples: Sample) -> dict[str, set[int]]:
    """
    From any number of Sample objects, returns the IDs of probes
    that are unique to each Sample object.
    The result is cast as a dictionary,
    where the keys contain information about the Sample object,
    and their values contain the set of unique probe IDs.
    If not having performed the filtering of Probe objects per Sample
    based on minimum expression
    and whether the probe is above background,
    this function will most likely return none of the IDs
    for all Samples.
    """

    sys.stdout = open(os.devnull, 'w')
    # This time the function will not print text
    t_shared_expressions: set[float] = get_intersection_in_probes(*samples)
    sys.stdout = sys.__stdout__
    return {
        f"Unique probe IDs for REGION: {sample_.structure_name}"
        f"; ID: {sample_.structure_id}":
        set(map(
            lambda probe_: probe_.probe_id,
            sample_.probes
        )) - t_shared_expressions
        for sample_ in samples
    }

probes_data: pd.DataFrame = pd.read_csv(
    "C:\\Users\\andri\\Documents\\Hanze\\Programming 1\\Probes.csv"
)
samples_data: pd.DataFrame = pd.read_csv(
    "C:\\Users\\andri\\Documents\\Hanze\\Programming 1"
    "\\SampleAnnot.csv"
)
probe_samples_data: pd.DataFrame = pd.read_csv(
    "C:\\Users\\andri\\Documents\\Hanze\\Programming 1"
    "\\MicroarrayExpression.csv",
    header=None
).iloc[:, 1:]
above_background_map: pd.DataFrame = pd.read_csv(
    "C:\\Users\\andri\\Documents\\Hanze\\Programming 1\\PACall.csv",
    header=None
).iloc[:, 1:]

selected_samples: list[Sample] = [Sample(
    int(sample_row_[0]),
    sample_acronym_,
    sample_row_[5],
    int(sample_row_[6]),
    [Probe(
        int(probes_data.loc[probe_row_index_, 'probe_id']),
        int(probes_data.loc[probe_row_index_, 'gene_id']),
        probes_data.loc[probe_row_index_, 'gene_name'],
        cast_cell_value(probes_data.loc[probe_row_index_, 'chromosome']),
        float(probe_samples_data.iat[probe_row_index_, sample_index_]),
        bool(int(above_background_map.iat[probe_row_index_, sample_index_]))
    ) for probe_row_index_ in range(probes_data.shape[0])]
) for sample_acronym_ in sample_acronyms
for sample_row_, sample_index_ in zip(
    *(lambda row_: (row_.values, row_.index))(samples_data[
        samples_data['structure_acronym'] == sample_acronym_
    ])
)]

for sample_object_ in selected_samples:
    if above_background is not None:
        sample_object_.probes = \
            sample_object_.get_probes_with_expression_greater_than(
                cutoff_value, above_background
            )
    else:
        sample_object_.probes = \
            sample_object_.get_probes_with_expression_greater_than(
                cutoff_value
            )

for sample_object_ in selected_samples:
    print(sample_object_)
print(get_intersection_in_probes(*selected_samples))
print(get_difference_in_probes(*selected_samples))
