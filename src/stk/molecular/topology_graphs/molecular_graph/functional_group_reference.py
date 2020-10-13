class FunctionalGroupReference:
    def __init__(self, node_id, functional_group_id):
        self._node_id = node_id
        self._functional_group_id = functional_group_id

    def get_node_id(self):
        return self._node_id

    def get_functional_group_id(self):
        return self._functional_group_id
