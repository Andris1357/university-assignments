"""
The module converts a .gb file to a format which would be the result of using the GenBank Extractor

It defines a GenBankParser class, which handles the conversion of the input file
The name of the output file will be the name of the input file
with a .txt extension instead of ".gb" and postfixed with "_features".
"""

import sys
from typing import Union, Final, Optional
import re
import getopt
from functools import reduce
import numpy as np

class Feature:
    """
    Defines a feature object that will be used to construct the elements of the output file.
    In the output file, they will be printed in the format ">{type_} /{attribute_}\n{sequence_}",
        where the each 60 characters in the sequence will be separated by a newline.
    """
    
    def __init__(
        self,
        sequence_: str,
        attribute_: dict[str, str],
        type_: str
    ) -> None:
        self.sequence = sequence_
        self.attribute = attribute_
        self.type = type_

class GenBankParser:
    """
    Defines an object that can convert a .gb file into the desired format.
    It has a class attribute that is used to translate the characters to their counterpart,
        when parsing a sequence that is prefixed with "complement".
    A GenBankParser object is initiated by passing in the path of the input file. It extracts:
        - The whole text from the file
        - The definition from the text
        - The "ORIGIN" section from the text that stores the source sequence
        - All specifications from the "FEATURES" section,
            where for each of them it is identified, whether
    """
    
    complement_map: Final[dict[str, str]] = (
        lambda lowercase_map_: lowercase_map_ | {
            key_.upper(): lowercase_map_[key_].upper() for key_ in lowercase_map_
        }
    )({'a': 't', 't': 'a', 'g': 'c', 'c': 'g'})

    def __init__(
        self,
        file_path_: str
    ) -> None:
        with open(file_path_, 'r') as file_:
            self.source: str = file_.read()
        self.definition: str = self.extractDefinition()
        self.specifications: list[dict[str, Union[str, bool, list[list[int]]]]] = []
        self.features: list[Feature] = []
        self.origin: str = self.extractOrigin()

    def extractDefinition(self) -> str:
        """
        Returns the "DEFINITION" section from the text of the input file.
        """
        return re.search(
            r'(?<=DEFINITION).*(?=\n)',
            self.source
        ).group().strip()

    def extractOrigin(self) -> str:
        """
        Returns the "ORIGIN" section from the input file in a way
            where only the characters are extracted and re-concatenated.
        """
        return ''.join(re.findall(
            r'[agtc]+',
            self.source[self.source.index("ORIGIN"):]
        ))

    def getSlices(self, specification_: str) -> list[list[int]]:
        """
        Extracts the index range from a specification text segment if it has one listed,
            and all ranges if it features a "join" clause;
            it only extracts one index if it is not a range.
        """
        if "join(" not in specification_:
            return [*map(int, re.search(
                r'\d+\.{0,2}\>?\d*',
                specification_
            ).group().replace('>', '').split('..'))]
        else:
            temp_truncated_spec = specification_[
                specification_.index("join(") + 5 :
            ].replace(' ', '').replace('>', '')
            return [
                list(map(int, segment_.strip('\n').split('..')))
                for segment_ in temp_truncated_spec[:temp_truncated_spec.index(')')].split(',')
            ]
    
    def extractSpecifications(self) -> None:
        """
        When the GenBankParser object's "source" attribute already stores a string,
            this method can be used to fill the object's "specifications" attribute
            with a list of dictionaries that contain the definition
            of how to create the Feature objects.
        """
        temp_specifications: str = self.source[
            re.search(
                r'FEATURES\s+Location\/Qualifiers',
                self.source
            ).span()[1] : self.source.index("ORIGIN")
        ]
        
        while len(temp_specifications) > 1:
            temp_current_match: Optional[re.Match] = re.search(
                r'"\n\s+[^\s\/]+',
                temp_specifications
            )
            
            # If there is no match, that is because we reached the last specification
            if temp_current_match is None:
                temp_current_match = re.search(r'"\n$', temp_specifications)
            
            temp_specification: str = temp_specifications[
                : temp_current_match.span()[0] + 1
            ].strip()
            temp_is_complement: bool = "complement" in temp_specification
            
            self.specifications.append({
                'type': re.search(r"\b\S*\b", temp_specification).group(),
                'attribute': re.search(
                    r'\/.*=.*"\n',
                    temp_specification
                ).group().lstrip('/').rstrip('\n'),
                'span': self.getSlices(
                    temp_specification if not temp_is_complement else temp_specification.lstrip(
                        "complement("
                    ).rstrip(')')
                ),
                'is_complement': temp_is_complement
            })
            
            """This way the quotation mark will be escaped and the pattern will
            find the beginning of the next specification section"""
            temp_specifications = temp_specifications[temp_current_match.span()[0] + 1 :]

    def sliceOrigin(self, slices_: Union[list[list[int]], list[int]], format_option_: str) -> str:
        """
        It takes the index specifiers of a specification, and extracts those characters
            from the sample sequence that match the specified index ranges.
        If the format option is uppercased,
            it casts the extracted characters mentioned above to uppercase
            and appends the characters before the first indexes of the ranges as lowercase
        """
        if type(slices_[0]) == int: # It is not a nested list, it contains only one span
            temp_sequence: str = self.origin[
                slices_[0] - 1
                : (slices_[1] if len(slices_) > 1 else slices_[0])
            ]
            if format_option_ == "uppercased":
                temp_sequence = self.origin[: slices_[0]] + temp_sequence.upper()
            return temp_sequence
        
        else: # If if contains multiple slices, we know it was prefixed with "join"
            temp_selected_indexes: list[int] = reduce(
                lambda accumulator_, next_: accumulator_ + next_,
                [[*range(
                    range_[0] - 1,
                    range_[1] if len(range_) > 1 else range_[0]
                )] for range_ in slices_]
            )
            if format_option_ == "uppercased":
                return ''.join(*np.asarray([*self.origin])[[temp_selected_indexes]])
            else:
                return ''.join(
                    character_ if index_ not in temp_selected_indexes else character_.upper()
                    for index_, character_ in enumerate(
                        self.origin[: temp_selected_indexes[-1] + 1]
                    )
                )
    
    def extractFeatures(self, format_option_: str) -> None:
        """
        Creates Feature objects from the specifications the GenBankParser object currently stores.
        If "complement" is True for the given specification,
            the sequence will be translated and reversed.
        """
        assert format_option_ in ["uppercased", "separated"], \
            "Option may only be specified as 'separated' or 'uppercase'"
        
        for specification_ in self.specifications:
            temp_feature: Feature = Feature(
                self.sliceOrigin(specification_['span'], format_option_),
                (lambda key_value_: {
                    key_value_[0]: key_value_[1].strip('"')
                })(specification_['attribute'].split('=')),
                specification_['type']
            )
            
            if specification_['is_complement']:
                temp_complement: list[str] = [*map(
                    lambda character_: GenBankParser.complement_map[character_],
                    temp_feature.sequence
                )]
                temp_complement.reverse()
                temp_feature.sequence = ''.join(temp_complement)

            self.features.append(temp_feature)

    def writeFile(self, input_path_: str, format_option_: str) -> None:
        """
        Creates file according to the specified format
        """
        with open(
            input_path_.rstrip('.gb') + "_features"
            + "uppercased" if format_option_ == "uppercased" else '' + '.txt',
            'w',
            newline='',
            encoding='UTF-8'
        ) as output_:
            output_.write("GenBank Feature Extractor results" + '\n' + self.definition + '\n' * 2)
            for feature_ in self.features:
                output_.write(
                    f'>{feature_.type} /{list(feature_.attribute.keys())[0]}='
                    f'"{list(feature_.attribute.values())[0]}"' + '\n'
                )
                output_.write(''.join(map(
                    lambda character_index_: feature_.sequence[character_index_ : (
                        character_index_ + 60 if character_index_ + 60 <= len(
                            feature_.sequence
                        ) else -1
                    )] + '\n',
                    range(0, len(feature_.sequence), 60)
                )) + '\n')

external_parameters: tuple[list[tuple[str, str]], list[str]] = getopt.getopt(sys.argv[1:], "")
input_file_path: Final[str] = external_parameters[1][0]
format_option: Final[str] = external_parameters[1][1]

genbank_parser: GenBankParser = GenBankParser(input_file_path)
genbank_parser.extractSpecifications()
genbank_parser.extractFeatures(format_option)
genbank_parser.writeFile(input_file_path, format_option)

print(f"Created file at location '{input_file_path.rstrip('.gb') + '_features'}"
      f"{'uppercased' if format_option == 'uppercased' else ''}.txt'")
