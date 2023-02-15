"""
This module:
- Unpacks the data from DaMN.xlsx, and retrieves
    - Unique treatments
    - Unique inoculation treatments
    - Rows of individuals that were not killed by scientist handling
    errors, and the location of these within the source table
- Imports definition of Flea class.
- Creates objects from all rows in the file DaMN.xlsx,
where an object's parameters are passed into the constructor
in the order of the columns from left to right.
- Creates a list of dictionaries, where a dictionary is returned
from calling function "get_average_metrics_per_treatment" with
a certain column header that was included among the questions.
- Prints attributes of one of the Flea objects as an example,
besides all dictionaries from the above step
"""

import pandas as pd
import numpy as np

from prep_progr_as5_class import Flea

NEW_LINE = '\n'

flea_data: pd.DataFrame = pd.read_excel(
    "C:\\Users\\andri\\Documents\\Hanze\\Programming 1\\DaMN.xlsx",
    sheet_name="Sheet1"
)
treatments: set[str] = set(flea_data['Treatment'])
inoculation_treatments: set[str] = set(flea_data['In_treatment'])
successful_specimen_indices = flea_data['H_errors'] != 1
successful_specimen = flea_data[successful_specimen_indices]

def get_average_metric_per_treatment(
    flea_objects_: list[Flea], metric_column_: str, inoculation_: bool
) -> dict[str, float]:
    """
    Can be applied to those four columns, for which average values have
    to be calculated based on the task description.
    For grouping based on inoculation treatment, only the column
    can to be used that stores whether an individual got infected.
    Maps the possible inputs (headers of above mentioned 4 columns)
    to the respective output descriptions, and names of Flea object
    attributes.
    Returns a dictionary that will have
    - elements for all groupings (number of distinct treatments)
    - keys describing
        - The basis of the grouping
        - Which column the averages were calculated from
    - a floating point value as the average of the column that was
    chosen for the calculation
    Only those individuals are taken into account that did not die
    due to scientist handling errors.
    """

    assert metric_column_ in flea_data.columns.values[[4, 5, 7, 8]], \
        f"The function only handles data from columns" \
        f" {flea_data.columns.values[[4, 5, 7, 8]]}"
    if inoculation_:
        assert metric_column_ == "Inf_mets", \
            "Only column 'Inf_mets' should be calculated" \
            "when grouping based on inoculation treatment"

    t_metric_description_map_: dict[str, str] = dict((
        column_, description_
    ) for column_, description_ in zip(
        flea_data.columns.values[[4, 5, 7, 8]],
        ["mortality (in days lived)", "offspring count",
         "infection_rate", "spore yield"]
    ))
    t_metric_attribute_map_: dict[str, str] = dict((
        column_, attribute_
    ) for column_, attribute_ in zip(
        flea_data.columns.values[[4, 5, 7, 8]],
        ["life_duration", "offspring_count",
         "was_infected", "spore_yield"]

    ))

    return dict((
        f"Average {t_metric_description_map_[metric_column_]} for"
        f"{' inoculation' if inoculation_ else ''} treatment {key_}",
        sum(
            getattr(
                flea_, t_metric_attribute_map_[metric_column_]
            ) for flea_ in filter(
                lambda object_: getattr(
                    object_,
                    'treatment' if not inoculation_ else 'inoculation_treatment'
                ) == key_ and not object_.was_handling_error_victim
                , flea_objects_
            )
        ) / successful_specimen[successful_specimen[
            'Treatment' if not inoculation_ else 'In_treatment'
        ] == key_].shape[0]
    ) for key_ in (
        treatments if not inoculation_ else inoculation_treatments
    ))


flea_objects: list[Flea] = [Flea(*map(
    lambda cell_: bool(int(cell_[1])) if cell_[0] in [7,10] else cell_[1],
    zip(range(11), flea_data.iloc[row_, :])
)) for row_ in range(flea_data.shape[0])]

average_metrics_per_treatment: list[dict[str, float]] = [
    get_average_metric_per_treatment(
        np.asarray(flea_objects)[[successful_specimen_indices]].tolist(),
        column_,
        column_ == 'Inf_mets'
    ) for column_ in [
        'Age_death', 'Total_juv', 'Spore_yield', 'Inf_mets'
    ]
]

print(f"Example individual:{NEW_LINE * 2}{flea_objects[5]}{NEW_LINE}")
print("Average metrics grouped by treatments:\n")
for dict_ in average_metrics_per_treatment:
    print(dict_)
