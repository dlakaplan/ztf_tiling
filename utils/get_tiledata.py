
"""
Download a bunch of ZTF metadata for different tiles
uses the ztf_tile list of tile centers (just the first set) to identify the Dec ranges
and then finds the most recent pointing that has a tile number in the appropriate range
writes the metadata to a separate csv file for each tile

Can use this to e.g., update the tile quadrant centers or other pointing information

DLK 2018-08-13
"""


from __future__ import print_function
from ztfquery import query
zquery = query.ZTFQuery()

from ztf_tiling import ztf_tiling
from astropy.table import Table
from astropy import units as u, constants as c
from astropy.coordinates import SkyCoord
import os
import numpy as np


tiles=Table.read(os.path.join(ztf_tiling.data._datadir,'ztf_fields.dat'),format='ascii.commented_header')
# only set 1
tiles=tiles[:np.argmin(np.diff(tiles['Dec']))+1]

Decs=np.unique(tiles['Dec'])
Decs=Decs[Decs>-30]
for dec in Decs:
    fields=tiles[tiles['Dec']==dec]['ID']
    zquery.load_metadata(sql_query="fid=2 and ccdid=1 and qid=1 and field between %d and %d order by obsjd DESC" % (fields.min(),fields.max()))
    if len(zquery.metatable > 0):
        print(dec,zquery.metatable['field'][0])
        field=zquery.metatable['field'][0]
        zquery.load_metadata(sql_query="fid=2 and field=%d and obsjd>%f" % (zquery.metatable['field'][0],zquery.metatable['obsjd'][0]-35./86400))
        t=Table.from_pandas(zquery.metatable)
        t.write('ztf_%04d.csv' % field, format='ascii.csv')
    else:
        print(dec,'None')
        

        
                         
