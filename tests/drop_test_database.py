"""
A utility for dropping the database running the tests creates.

This script is not part of the test suite. It is meant to the run
manually by developers.

"""


import pymongo


def main():
    pymongo.MongoClient().drop_database('_stk_pytest_database')


if __name__ == '__main__':
    main()
