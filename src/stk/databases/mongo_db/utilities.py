"""
MongoDB Utilities
=================

"""


class HashableDict(dict):
    def __hash__(self):
        try:
            return hash((
                frozenset(self),
                frozenset(_to_hashable(self.values())),
            ))
        except Exception:
            print((
                (self, ),
                _to_hashable(self.values()),
            ))
            raise

    def __eq__(self, other):
        return super().__eq__(other)


def _to_hashable(item):
    if isinstance(item, list):
        return tuple(_to_hashable(subitem) for subitem in item)
    else:
        return item
