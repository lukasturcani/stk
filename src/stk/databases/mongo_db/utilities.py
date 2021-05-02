"""
MongoDB Utilities
=================

"""


class HashableDict(dict):
    def __hash__(self):
        return hash((
            frozenset(self),
            frozenset(_to_hashable(self.values())),
        ))

    def __eq__(self, other):
        return super().__eq__(other)


def _to_hashable(iterable):
    for item in iterable:
        if isinstance(item, list):
            yield tuple(_to_hashable(subitem) for subitem in item)
        else:
            yield item
