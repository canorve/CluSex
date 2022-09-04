import sys



from clusex.lib.join import joinsexcat 
from clusex.lib.make import CatArSort
from clusex.lib.make import MakeMask
from clusex.lib.make import MakeSatBox 

from clusex.lib.ds9 import ds9kron
from clusex.lib.sky import SkyCal 

if sys.version_info[:2] >= (3, 8):
    # TODO: Import directly (no need for conditional) when `python_requires = >= 3.8`
    from importlib.metadata import PackageNotFoundError, version  # pragma: no cover
else:
    from importlib_metadata import PackageNotFoundError, version  # pragma: no cover

try:
    # Change here if project is renamed and does not equal the package name
    dist_name = "CluSex"
    __version__ = version(dist_name)
except PackageNotFoundError:  # pragma: no cover
    __version__ = "unknown"
finally:
    del version, PackageNotFoundError
