class MockMongoClient:
    def __getattr__(self, name):
        return MockMongoDatabase()

    def __getitem__(self, name):
        return MockMongoDatabase()


class MockMongoDatabase:
    def __getattr__(self, name):
        return MockMongoCollection()

    def __getitem__(self, name):
        return MockMongoCollection()


class MockMongoCollection:
    def __init__(self):
        self._documents = []

    def insert_one(self, document):
        self._documents.append(document)

    def find_one(self, key):
        for document in self._documents:
            if self._is_match(document, key):
                return document

    @staticmethod
    def _is_match(document, key):
        for key_ in key:
            if key_ not in document:
                return False
            if document[key_] != key[key_]:
                return False

        return True
