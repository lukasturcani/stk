import typing
from collections import abc


_T = typing.TypeVar('_T')
OneOrMany = typing.Union[_T, abc.Iterable[_T]]
