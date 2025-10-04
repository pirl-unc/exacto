from importlib.metadata import version, PackageNotFoundError

try:
    __version__ = version("exacto")
except PackageNotFoundError:
    __version__ = "unknown"