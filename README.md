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

### Basic check to see if a position is within a tile:
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

### Look at a all of the pixels in a healpix map and see what the probability for a tile is:
```python
from ztf_tiling import ztf_tiling
import healpy as hp

nside=256
npix=hp.nside2npix(nside)

# coordinates for each pixel in the map
RA,Dec=hp.pix2ang(nside, np.arange(npix),lonlat=True)
Z=ztf_tiling.get_tile(651)
# this can be slow for a large value of npix
on=Z.inside(RA,Dec)
# now plot the result, rotating to the center of the tile:
hp.mollview(on,rot=(Z.RA.value,Z.Dec.value,0),xsize=4000)
```
Note that this will appear a bit blotchy toward the edges becase of the coarse pixelization.


