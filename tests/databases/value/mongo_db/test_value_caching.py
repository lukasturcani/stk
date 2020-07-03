import stk
import pymongo


def test_put_caching():
    collection = '_test_put_caching'
    database_name = '_test_put_caching'
    client = pymongo.MongoClient()
    client.drop_database(database_name)

    database = stk.ValueMongoDb(
        mongo_client=client,
        collection=collection,
        database=database_name,
    )
    molecule = stk.BuildingBlock('CCC')
    database.put(molecule, 43)
    database.put(molecule, 43)

    cache_info = database._put.cache_info()
    assert cache_info.hits == 1
    assert cache_info.misses == 1

    database.put(molecule, 40)
    cache_info = database._put.cache_info()
    assert cache_info.hits == 1
    assert cache_info.misses == 2


def test_get_caching():
    collection = '_test_get_caching'
    database_name = '_test_get_caching'
    client = pymongo.MongoClient()
    client.drop_database(database_name)

    database = stk.ValueMongoDb(
        mongo_client=client,
        collection=collection,
        database=database_name,
    )
    molecule = stk.BuildingBlock('CCC')
    database.put(molecule, 43)
    database.get(molecule)
    database.get(molecule)

    cache_info = database._get.cache_info()
    assert cache_info.hits == 1
    assert cache_info.misses == 1
