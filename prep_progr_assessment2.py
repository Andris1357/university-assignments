# Runs on Python 3.10

import pandas as pd
from typing import Union
from functools import reduce

new_line = "\n"

# TASK 1 ###

# Groups probe IDs by gene ID
probes_data: pd.DataFrame = pd.read_csv("C:\\Users\\andri\\Documents\\Hollandia\\Programming\\Probes.csv")
gene_ids: set[int] = set(probes_data.loc[:, "gene_id"])
probes_grouped_by_gene: dict[int, list[int]] = dict((
    key_, list(probes_data[probes_data["gene_id"] == key_]["probe_id"].values)
) for key_ in gene_ids)
del probes_data

# Retrieves average expression for each probe and casts that into a dictionary where keys are the probe IDs
probes_expression_samples: pd.DataFrame = pd.read_csv(
    "C:\\Users\\andri\\Documents\\Hollandia\\Programming\\MicroarrayExpression.csv",
    header= None
)
probes_expression_averages: Union[pd.DataFrame, dict[int, float]] = pd.DataFrame(
    {"probe_id": probes_expression_samples.iloc[:, 0]}
).join(
    pd.DataFrame({"avg": probes_expression_samples.iloc[:, 1:].mean(axis=1)})
)
probes_expression_averages = dict(
    (key_, val_) for key_, val_ in zip(probes_expression_averages.iloc[:, 0], probes_expression_averages.iloc[:, 1])
)

# Retrieves and outputs the ID of probe with the highest expression average for each gene ID
temp__highest_average_value: float
genes_highest_average_expression: dict[int, dict[str, Union[float, list[int]]]] = dict((
    key_, {
        "highest_avg": (temp__highest_average_value := max(
            map(probes_expression_averages.get, probes_grouped_by_gene[key_])
        )),
        "respective_probe_id(s)": [x for x in probes_expression_averages if probes_expression_averages[x] == temp__highest_average_value]
    } # Respective id(s) is of type list, because there may be an edge case where multiple equivalent values exist that are the highest
) for key_ in probes_grouped_by_gene)
print(f"ID of probe with the highest average expression for each gene{new_line}{genes_highest_average_expression}")

# TASK 2 ###

"""
Retrieves those expression value columns from MicroarrayExpression.csv, 
    where the column index matches the row index of the rows from Probes.csv, 
    where the column "structure_acronym" has a value of "LHM" or "PHA", respectively.
Retrieves the index of all rows from the above selected expression value columns,
    where the expression value is higher than 15,
    for each sample column,
    for each region.
Deallocates variable housing MicroarrayExpression.csv.
"""
samples: pd.DataFrame = pd.read_csv("C:\\Users\\andri\\Documents\\Hollandia\\Programming\\SampleAnnot.csv")
LHM_sample_columns: pd.DataFrame = pd.DataFrame({
    f"expression_values_column_{index_}": probes_expression_samples.iloc[:, index_] for index_ in samples[samples["structure_acronym"] == "LHM"].index
})
LHM_sample_filtered_row_numbers: dict[str, list[int]] = {
    column_header_: list(LHM_sample_columns[LHM_sample_columns[column_header_] > 15].index) for column_header_ in LHM_sample_columns.columns
}
PHA_sample_columns: pd.DataFrame = pd.DataFrame({
    f"expression_values_column_{index_}": probes_expression_samples.iloc[:, index_] for index_ in samples[samples["structure_acronym"] == "PHA"].index
})
PHA_sample_filtered_row_numbers: dict[str, list[int]] = {
    column_header_: list(PHA_sample_columns[PHA_sample_columns[column_header_] > 15].index) for column_header_ in PHA_sample_columns.columns
}
del probes_expression_samples

"""
Retrieves the IDs of those probes that are present in all samples from both brain regions.
Retrieves probe IDs for each sample that do not exist in other samples of the same group.
"""
probes_shared_across_LHM_samples: set[int] = reduce(
    lambda accumulator_, next_: set(accumulator_).intersection(set(next_)),
    map(lambda items_: items_[1], LHM_sample_filtered_row_numbers.items()) # .items() returns a sequence of tuples, where in each tuple key is at index 0, values at index 1.
)
probes_shared_across_PHA_samples: set[int] = reduce(
    lambda accumulator_, next_: set(accumulator_).intersection(set(next_)),
    map(lambda items_: items_[1], PHA_sample_filtered_row_numbers.items())
)
probes_unique_between_LHM_samples: set[int] = reduce(
    lambda accumulator_, next_: set(accumulator_) ^ set(next_), # "^" when applied to sets, returns the items that are not found in both sets
    map(lambda items_: items_[1], LHM_sample_filtered_row_numbers.items())
)
probes_unique_between_PHA_samples: set[int] = reduce(
    lambda accumulator_, next_: set(accumulator_) ^ set(next_),
    map(lambda items_: items_[1], PHA_sample_filtered_row_numbers.items())
)

"""
Retrieves the union of the IDs of probes present across the samples taken from brain regions LHM or PHA.
Retrieves and outputs the IDs of those probes that are present in both brain regions.
Retrieves and outputs the IDs of those probes that are unique for each brain region.
"""
all_probes_across_LHM_samples: set[int] = reduce(
    lambda accumulator_, next_: set(accumulator_).union(set(next_)),
    map(lambda items_: items_[1], LHM_sample_filtered_row_numbers.items())
)
all_probes_across_PHA_samples: set[int] = reduce(
    lambda accumulator_, next_: set(accumulator_).union(set(next_)),
    map(lambda items_: items_[1], PHA_sample_filtered_row_numbers.items())
)
probes_shared_between_regions = all_probes_across_LHM_samples.intersection(all_probes_across_PHA_samples)
probes_unique_between_regions: dict[str, set[int]] = dict((
    key_, vals_ ^ probes_shared_between_regions
) for key_, vals_ in zip(
    ["LHM_unique", "PHA_unique"],
    [all_probes_across_LHM_samples, all_probes_across_PHA_samples]
))
print(
    f"IDs of probes shared between LHM and PHA regions:{new_line}{probes_shared_between_regions}",
    f"IDs of probes unique for LHM and PHA regions:{new_line}{probes_unique_between_regions}"
)