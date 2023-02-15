"""
This module:
- Contains class definitions of Sample and Probe objects
- Contains class definition of custom exception
"""

from types import NoneType
from typing import Union

NEWL = '\n'

class ParametersUnfilled(Exception):
    """
    Exception that is raised when
    either Probe or Sample constructors
    receive parameters that are not of the expected type.
    It does not check for unprovided arguments,
    because Python will not let an instance be constructed anyways,
    if any of the constructor's arguments are missing.
    """

    def __init__(self, **incorrect_parameter_types_: type) -> None:
        self.message = NEWL.join(
            "Parameter" + index_
            + " has been provided a value of incorrect type" + repr(
                incorrect_parameter_types_[index_]
            ) for index_ in incorrect_parameter_types_
        )
        super().__init__(self.message)

class Probe(object):
    """
    Describes a Probe object, where
    - Parameters 'probe_id_', 'gene_id_', 'gene_name_', 'chromosome_'
    will be provided from file 'Probes.csv'.
    - Parameter 'expression_'
    will be provided from file 'MicroarrayExpression.csv'.
    - Parameter 'above_background_'
    will be provided from file 'PACall.csv'.

    If there are any arguments, the type of which
    did not meet their predescribed type,
    the exception 'ParametersUnfilled' will be raised.
    When called by built-ins 'print' and 'repr',
    prints attributes of Probe instance,
    based on which it can be identified.
    """

    created_probe_count: int = 0

    def __init__(
        self,
        probe_id_: int,
        gene_id_: int,
        gene_name_: str,
        chromosome_: Union[str, int, None],
        expression_: float,
        above_background_: bool
    ) -> None:
        t_argument_types: list[list[type]] = [
            [int], [int], [str], [str, int, NoneType], [float], [bool]
        ]
        t_incorrect_parameters: filter = filter(
            lambda param_: param_[1] not in t_argument_types[param_[0]],
            enumerate(map(
                type,
                [probe_id_, gene_id_, gene_name_,
                 chromosome_, expression_, above_background_]
            ))
        )
        if len(list(t_incorrect_parameters)) > 0:
            raise ParametersUnfilled(
                **{param_[0]: param_[1] for param_ in t_incorrect_parameters}
            )

        self.probe_id = probe_id_
        self.gene_id = gene_id_
        self.gene_name = gene_name_
        self.chromosome = chromosome_
        self.expression = expression_
        self.above_background = above_background_
        Probe.created_probe_count += 1

    def __str__(self) -> str:
        return f"Probe instance{NEWL}{'-' * 10}{NEWL}" \
               f"Probe id: {self.probe_id}{NEWL}" \
               f"Corresponding gene id: {self.gene_id}{NEWL}" \
               f"Corresponding gene name: {self.gene_name}{NEWL}" \
               f"Corresponding chromosome id: {self.chromosome}"


class Sample(object):
    """
    Describes a Sample object, where
    - Parameters 'structure_id_', 'structure_acronym_',
    'structure_name_', 'polygon_id_'
    will be provided from file 'SampleAnnot.csv'.
    - ID of corresponding polygon is provided because there are cases
    when the first 3 attributes match for a sample.
    - Parameter 'probes_'
    will be assigned lists of Probe objects
    that were created with the expression values found in the column
    that belongs to this sample,
    based on its row index in 'SampleAnnot.csv'.

    If there are any arguments, the type of which
    did not meet their predescribed type,
    the exception 'ParametersUnfilled' will be raised.
    The class has a method defined to get the list of those Probes,
    that have an expression greater than a target value,
    and are optionally also above background.
    When called by built-ins 'print' and 'repr',
    prints attributes of Sample instance,
    based on which it can be identified.
    """

    created_sample_count: int = 0

    def __init__(
        self,
        structure_id_: int,
        structure_acronym_: str,
        structure_name_: str,
        polygon_id_: int,
        probes_: list[Probe]
    ) -> None:
        t_incorrect_parameters: filter = filter(
            lambda param_: param_[1] != [
                int, str, str, int, list
            ][param_[0]] and (True if param_[0] != 4 else not all(
                isinstance(probe_, Probe) for probe_ in probes_
            )),
            enumerate(map(
                type,
                [structure_id_, structure_acronym_,
                 structure_name_, polygon_id_, probes_]
            ))
        )
        if len(list(t_incorrect_parameters)) > 0:
            raise ParametersUnfilled(
                **{param_[0]: param_[1] for param_ in t_incorrect_parameters}
            )

        self.structure_id = structure_id_
        self.structure_acronym = structure_acronym_
        self.structure_name = structure_name_
        self.polygon_id = polygon_id_
        self.probes = probes_
        Sample.created_sample_count += 1

    def __str__(self) -> str:
        return f"Sample instance{NEWL}{'-' * 10}{NEWL}" \
               f"Structure id: {self.structure_id}{NEWL}" \
               f"Structure acronym: {self.structure_acronym}{NEWL}" \
               f"Structure name: {self.structure_name}{NEWL}" \
               f"Polygon id: {self.polygon_id}"

    def get_probes_with_expression_greater_than(
        self,
        cutoff_: int,
        above_background_: bool = False
    ) -> list[Probe]:
        """
        From a list of Probe objects, returns those
        that have their expression above the value
        specified by 'cutoff_'.
        If default 'above_background_' is overridden as True,
        it excludes Probes that were above background,
        based on PACall.csv.
        """

        return list(filter(
            lambda probe_: probe_.expression > cutoff_ and (
                True if not above_background_ else probe_.above_background
            ),
            self.probes
        ))
