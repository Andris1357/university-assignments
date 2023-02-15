import sys
import getopt
from typing import Generator, Union, Any, Callable
import csv
import functools as fn
import pandas

from calculate_timeseries_metrics import \
    calculate_standard_deviation, calculate_historical_volatility,\
    calculate_moving_average, calculate_moving_average_crossover_divergence

BACKSLASH = '\\'
external_parameters = getopt.getopt(sys.argv[1:], "")

input_file_path: str = external_parameters[1][0]
separator: str = external_parameters[1][1]
no_header: bool = bool(int(external_parameters[1][2]))
short_moving_average_length: int = int(external_parameters[1][3])
calculation_sequence_length: int = int(external_parameters[1][4])
is_moving_average_exponential: bool = bool(int(external_parameters[1][5]))
output_file_name: str = external_parameters[1][6]
# Values can be passed in any order: vol, ma, macd, std
calculation_references: list[str] = external_parameters[1][7:]

assert all(
    reference_ in [
        'vol', 'ma', 'macd', 'std'
    ] for reference_ in calculation_references
), "Only calculations ['vol', 'ma', 'macd', 'std'] can be specified."

if 'macd' in calculation_references:
    assert short_moving_average_length < calculation_sequence_length,\
        "The fourth argument should refer to a calculation sequence " \
        "length that is longer than that of the third argument."

def flatten_list(list_of_lists: list[list[Any]]) -> list[Any]:
    """Arranges all items of lists from within a list
    to be placed in one list"""
    return fn.reduce(
        lambda accumulator_, next_: accumulator_ + next_,
        list_of_lists
    )

def generate_dataframe_rows(
    *timeseries_columns_: pandas.Series(float),
    **metrics_: Union[Callable[..., float], fn.partial]
) -> Generator[list[Union[float, None]], None, None]:
    """Generates output file rows one by one, using the functions
    that were selected for calculating timeseries metrics.
    Using a generator is useful,
    because the method csv.writer.writerow
    inserts those one by one anyways"""

    for window_start_, window_end_ in enumerate(range(
        calculation_sequence_length - 1,
        max(map(len, timeseries_columns_))
    )):
        yield flatten_list(
            [[timeseries_columns_[column_group_].iloc[window_end_]] + [
                 metrics_[argument_key_](
                    timeseries_columns_[column_group_].iloc[window_start_: window_end_]
                 ) for argument_key_ in metrics_
            ] for column_group_ in range(len(timeseries_columns_))]
        )

# Maps functions to their abbreviations
map_calculation_functions: dict[
    str, Union[Callable[..., float], fn.partial]
] = {
    'vol': calculate_historical_volatility,
    'ma': calculate_moving_average if not is_moving_average_exponential
    else fn.partial(calculate_moving_average, exponential_=True),
    'macd': fn.partial(
        calculate_moving_average_crossover_divergence,
        short_length_=short_moving_average_length
    ),
    'std': calculate_standard_deviation
}

input_file: pandas.DataFrame = pandas.read_csv(
    input_file_path, sep=separator, header=(None if no_header else 0)
)
yieldRow: Generator[list[float], None, None] = generate_dataframe_rows(
    *[input_file.iloc[:, column_index_] for column_index_ in range(input_file.shape[1])],
    **dict((key_, map_calculation_functions[key_]) for key_ in calculation_references)
)

with open(
    f"{input_file_path[:input_file_path.rfind(BACKSLASH)]}{BACKSLASH}"
    f"{output_file_name}.csv",
    'w', newline='', encoding='UTF8'
) as f:
    writer = csv.writer(f)

    # Creates headers
    writer.writerow(flatten_list(
        [[
            t_current_column := f"Column{column_index_}" if no_header
            else input_file.columns.values[column_index_],
            *[f"{t_current_column}_{metric_}" for metric_ in calculation_references]
        ] for column_index_ in range(input_file.shape[1])]
    ))

    """Inserts rows that will not have metrics associated
    due to residing at a row index less than the calculation length"""
    for row_index_ in range(calculation_sequence_length - 1):
        writer.writerow(flatten_list(
            [[
                float(input_file.iat[row_index_, column_index_]),
                *[None for _ in calculation_references]
            ] for column_index_ in range(input_file.shape[1])]
        ))

    for timestep_ in range(input_file.shape[0] - calculation_sequence_length + 1):
        writer.writerow(next(yieldRow))

print(f"Generated file {output_file_name}")
