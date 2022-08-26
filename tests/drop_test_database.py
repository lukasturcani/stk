"""
A utility for dropping the databases created by running the tests.

This script is not part of the test suite. It is meant to be run
manually by developers when they feel the need to.

"""


import pymongo


def main():
    databases_to_drop = (
        "_stk_pytest_database",
        "_stk_test_database_for_testing",
    )
    client = pymongo.MongoClient()
    for database in databases_to_drop:
        client.drop_database(database)


if __name__ == "__main__":
    main()
