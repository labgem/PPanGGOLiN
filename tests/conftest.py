import pytest

def pytest_addoption(parser):
    parser.addoption(
        "--cpu",
        action="store",
        default="2",
        help="Number of CPUs to use in functional tests"
    )
    parser.addoption(
        "--full",
        action="store_true",
        default=False,
        help="Run full test suite (unit + functional)"
    )
    parser.addoption(
        "--update-golden",
        action="store_true",
        default=False,
        help="Update golden hashes JSON instead of just testing.",
    )


@pytest.fixture(scope="module")
def num_cpus(request):
    return request.config.getoption("--cpu")


@pytest.fixture
def update_golden(request):
    return request.config.getoption("--update-golden")


def pytest_collection_modifyitems(config, items):
    if config.getoption("--full"):
        # Run all tests (unit + functional)
        return

    # By default, skip functional tests
    skip_functional = pytest.mark.skip(reason="Skipping functional tests by default")
    for item in items:
        if "functional_tests" in str(item.fspath):
            item.add_marker(skip_functional)
