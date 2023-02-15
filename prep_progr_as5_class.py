"""
Houses definition of the Flea class.
"""

from typing import Union

NEW_LINE = '\n'

class Flea:
    """
    Defines an individual from the experiment, where one
    can be instantiated by spreading all values from a certain row.
    When an object created from this class has "repr" or "print"
    called upon, it lists the names and values of its attributes.
    """

    def __init__(self, *parameters_) -> None:
        assert len(parameters_) == 11, "Incorrect amount of columns passed"

        self.specimen_id: int = parameters_[0]
        self.treatment: str = parameters_[1]
        self.nanoplastic_treatment: str = parameters_[2]
        self.inoculation_treatment: str = parameters_[3]
        self.life_duration: int = parameters_[4]
        self.offspring_count: int = parameters_[5]
        self.birth_count: int = parameters_[6]
        self.was_infected: bool = parameters_[7]
        self.spore_yield: int = parameters_[8]
        self.first_birth_at: int = parameters_[9]
        self.was_handling_error_victim: bool = parameters_[10]

    def __str__(self) -> str:
        t_instance_attributes: list[str] = list(filter(
            lambda attribute_: '__' not in attribute_,
            dir(self)
        ))
        t_attribute_value_map: dict[str, Union[int, str, bool]] = dict((
            attribute_, getattr(self, attribute_)
        ) for attribute_ in t_instance_attributes[1:])

        return f"Attributes of individual with ID " \
               f"{getattr(self, t_instance_attributes[0])}:{NEW_LINE}"\
               f"{t_attribute_value_map}"
