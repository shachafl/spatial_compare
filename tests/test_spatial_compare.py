from spatial_compare import SpatialCompare, get_column_ordering
import pandas as pd
import pathlib

SC_DIR = pathlib.Path().absolute()
print(SC_DIR)

XENIUM_STEM = "xenium.h5ad"
TEST_DIR = SC_DIR.joinpath("tests").joinpath("data")
XENIUM_DIR = TEST_DIR.joinpath(XENIUM_STEM)
TEST_DF_RECORDS = [dict(a=1, b=0.5, c=0.1), dict(a=.5, b=.8, c=.9), dict(a=0.1, b=1., c=0.1)]
TEST_DF = pd.DataFrame.from_records(TEST_DF_RECORDS)
TEST_DF.index  = ["a", "b", "c"]


def test_get_column_ordering():
    ordered_columns = get_column_ordering(TEST_DF, ordered_rows=["a", "b", "c"])
    print(ordered_columns)
    assert ordered_columns == ["a", "c", "b"]