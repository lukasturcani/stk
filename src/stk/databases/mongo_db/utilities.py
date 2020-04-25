class HashableDict(dict):
    def __hash__(self):
        return hash((frozenset(self), frozenset(self.values())))

    def __eq__(self, other):
        return super().__eq__(other)
