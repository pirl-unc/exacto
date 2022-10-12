import os


def get_data_path(name: str) -> str:
    """
    Return the absolute path to a file in the test/data directory.

    Parameters
    ----------
    name    :   Name of file.

    Returns
    -------
    Absolute path to a file in the test/data directory.
    """
    return os.path.join(os.path.dirname(__file__), "data", name)
