def pytest_add_option(parser):
    parser.addoption(
        '--mongodb_uri',
        action='store',
        default='mongodb://localhost:27017/',
    )


def pytest_generate_tests(metafunc):
    mongodb_uri = metafunc.config.option.mongodb_uri
    if 'mongodb_uri' in metafunc.fixturenames:
        metafunc.parametrize('mongodb_uri', [mongodb_uri])
