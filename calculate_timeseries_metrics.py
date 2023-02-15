import statistics as stat
import pandas as pd

def calculate_moving_average(
    sequence_: pd.Series, exponential_: bool = False
) -> float:
    """Calculates the average value for a sequence of numbers.
    If the option exponential_ is set to True,
    it will calculate the average in a way that the values
    will be weighted (multiplied) by an ascending weight sequence.
    The weights in the sequence will add up to the same sum,
    as if we multiplied the length of the sequence by one
    (which represents the weights used for a simple average,
    where there is no difference for the weights).
    Having the weights in an ascending order renders
    the most recent timeseries value of the highest importance,
    and the earliest one of the lowest."""

    t_values_count = len(sequence_)
    if not exponential_:
        return sum(map(float, sequence_)) / t_values_count

    return sum(
        value_ * (index_ + 1) * (
            1 / (t_values_count / 2 + 0.5)
        ) for index_, value_ in enumerate(map(float, sequence_))
    ) / t_values_count

def calculate_standard_deviation(sequence_: pd.Series) -> float:
    """Calculates the standard deviation for a sequence of numbers."""
    return stat.stdev(map(float, sequence_))

def calculate_moving_average_crossover_divergence(
    sequence_: pd.Series, short_length_: int
) -> float:
    """Calculates the difference between the values of two
    moving averages with differing length, for a sequence of numbers.
    It yields a negative (positive) value
    when the value of the short moving average is below (above)
    that of the longer one,
    indicating that the timeseries
    is trending downwards (upwards) in the short term.
    It yields a near-zero value
    when the moving averages are crossing each other,
    indicating a reversal of the previous trend"""
    return calculate_moving_average(sequence_) - calculate_moving_average(
        sequence_.iloc[len(sequence_) - short_length_:]
    )

def get_deviation(current_value_: float, earlier_value_: float) -> float:
    """Calculates the percentage change between two timeseries values,
    defaulting to 100% if the earlier value equals zero."""
    if earlier_value_ == 0.:
        return 1.
    return current_value_ / earlier_value_ - 1

def calculate_historical_volatility(sequence_: pd.Series) -> float:
    """Calculates the historical volatility of a sequence of numbers,
    meaning what is the average percentage change between two values
    measured at consecutive timestamps.
    This metric operates on a length that is 1 shorter
    than the length of other metrics,
    because for the first value of the sequence,the change between it
    and the previous value cannot be calculated."""

    t_values_count: int = len(sequence_)
    t_average_deviation: float = sum(
        get_deviation(
            float(sequence_.iloc[x]), float(sequence_.iloc[x - 1])
        ) for x in range(1, t_values_count)
    ) / (t_values_count - 1)

    return (sum((
        get_deviation(
            float(sequence_.iloc[x]), float(sequence_.iloc[x - 1])
        ) - t_average_deviation
    ) ** 2 for x in range(
        1, t_values_count
    )) / (t_values_count - 1)) ** 0.5
