import pytest
import pymongo


def pytest_addoption(parser):
    parser.addoption(
        '--mongodb_uri',
        action='store',
        default='mongodb://localhost:27017/',
    )


@pytest.fixture
def mongo_client(pytestconfig):
    return pymongo.MongoClient(pytestconfig.getoption('mongodb_uri'))
