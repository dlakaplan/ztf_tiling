"""
Look at the per-tile csv files downloaded previously (by get_tiledata.py)
examine each and determine the per-quadrant pointing offsets
these are then averaged (along with updated pixel scale information)
and written to a file of the form:

ZTF_quadrantcenters_YYMMDD.dat

which can be used by ztf_tile package
"""

from ztf_tiling import ztf_tiling
from astropy.table import Table,Column
from astropy import units as u, constants as c
from astropy.coordinates import SkyCoord
import os,glob,datetime
import numpy as np
from astropy.wcs import WCS
from astropy.io import fits
import astropy.wcs

tiles=Table.read(os.path.join(ztf_tiling.data._datadir,'ztf_fields.dat'),format='ascii.commented_header')
# only set 1
tiles=tiles[:np.argmin(np.diff(tiles['Dec']))+1]

files=sorted(glob.glob('ztf_????.csv'))
# tile center
RA=np.zeros(len(files))
Dec=np.zeros(len(files))
scale=np.zeros((len(files),64))*u.arcsec
rot=np.zeros((len(files),64))*u.deg
# quadrant centers
ra=np.zeros((len(files),64))*u.deg
dec=np.zeros((len(files),64))*u.deg
x=np.zeros((len(files),64))*u.deg
y=np.zeros((len(files),64))*u.deg


for j,file in enumerate(files):
    tile=int(file.split('_')[1].split('.')[0])
    #if not tile==722:
    #    continue
    RA[j]=tiles[tile-1]['RA']
    Dec[j]=tiles[tile-1]['Dec']
    Z=ztf_tiling.ZTFtile(RA[j]*u.deg,Dec[j]*u.deg)
    data=Table.read(file,format='ascii.csv')
    w=[]
    
    for ccd in range(16):
        for q in range(4):
            try:
                k=np.argwhere((data['ccdid']==(ccd+1)) & (data['qid']==(q+1)))[0][0]
            except IndexError:
                continue
            w.append(WCS(naxis=2))
            i=ccd*4+q
            w[-1].wcs.crpix=(data[k]['crpix1'],data[k]['crpix2'])
            w[-1].wcs.ctype=['RA---TAN','DEC--TAN']
            w[-1].wcs.crval=(data[k]['crval1'],data[k]['crval2'])
            w[-1].wcs.cd=[[data[k]['cd11'],data[k]['cd12']],
                          [data[k]['cd21'],data[k]['cd22']]]
            #ff=w[-1].calc_footprint(axes=(3072,3080))
            scale[j,i]=astropy.wcs.utils.proj_plane_pixel_scales(w[-1])[0]*3600*u.arcsec
            m=w[-1].pixel_scale_matrix
            rot[j,i]=-((np.arctan2(m[0,1],m[0,0])*u.rad).to(u.deg)+180*u.deg)
            ra[j,i]=data[k]['crval1']*u.deg
            dec[j,i]=data[k]['crval2']*u.deg

    x[j],y[j]=Z.celestial2projplane(ra[j],dec[j])

x[ra==0]=np.nan
y[ra==0]=np.nan
scale[ra==0]=np.nan
rot[ra==0]=np.nan
rot[rot<-100*u.deg]+=360*u.deg
qnum=np.arange(64)
ccd=qnum/4+1
q=qnum-(ccd-1)*4+1

t=Table([Column(ccd,name='CCD'),Column(q,name='Quadrant'),
         Column(np.nanmedian(x,axis=0),name='DX'),
         Column(np.nanmedian(y,axis=0),name='DY'),
         Column(np.nanmedian(scale,axis=0),name='scale')])
t.write('ZTF_quadrantcenters_%s.dat' % datetime.datetime.now().strftime('%Y%m%d'),format='ascii.commented_header')
