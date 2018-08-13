# ztf_tiling

## Data:
data sources (in `data/`):
 * `ztf_fields.txt`: list of ZTF field centers for primary and secondary grids.  From XXX?
 * `ZTF_quadrantcenters.dat`: centers of individual quadrants relative to the pointing center in the TAN projection **OLD**
 * `ZTF_quadrantcenters_20180813.dat`: centers of individual quadrants, updated 2018-08-13 (by @dlakaplan)

## Installation:
Standard `python ./setup.py install`

## Requirements:
 * astropy
 * healpy (not actually for most functionality)

## Usage:
Basic check to see if a position is within a tile:
```python
from ztf_tiling import ztf_tiling
from astropy import units as u
from astropy.coordinates import SkyCoord

# test a position
c=SkyCoord(33.30833*u.deg, 32.73167*u.deg)
# get the tile that corresponds roughly to this position
Z=ztf_tiling.get_tile(651)
print(Z.inside(c))
```
should return 21.  To figure out which CCD/quadrant:
```
print('CCD=%d' % ((Z.inside(c)-1)//4+1))
print('Quadrant=%d' % (Z.inside(c)-4*(Z.inside(c)-1)//4))
```
should return 6 and 1

