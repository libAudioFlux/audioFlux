from enum import Enum

__all__ = [
    'ReassignType'
]


class ReassignType(Enum):
    """
    Reassign Type
    """
    ALL = 0
    FRE = 1
    TIME = 2

    NONE = 3
