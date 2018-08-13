# Copyright (C) 2018 David Kaplan
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation; either version 2 of the License, or (at your
# option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
# Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.

from __future__ import print_function
import numpy as np
from astropy.table import Table
from astropy import constants as c, units as u
from astropy.coordinates import SkyCoord
from astropy.wcs import WCS
import astropy
import os
from . import data

##################################################
# utility functions
##################################################

##################################################
# figure out points inside a polygon
def inside(X,Y,x0,y0):
    """
    result=inside(X,Y,x0,y0)
    see if the points (x0,y0) is inside the convex polygon defined by (X,Y)
    X,Y and x0,y0 should be 1D arrays
    (X,Y) should be in order, but can be clockwise or counter-clockwise
    based on:
    http://demonstrations.wolfram.com/AnEfficientTestForAPointToBeInAConvexPolygon/
    """
    X=X-np.array([x0]).T
    # account for wraps in RA
    if (X>180).any():
        X[X>=180]-=360
    if (X<-180).any():
        X[X<=-180]+=360
    # we won't be able to deal with wraps in Dec
    # so we don't want to do this for very close to the poles
    Y=Y-np.array([y0]).T
    n=X.shape[1]
    aplus=np.ones((X.shape[0]),dtype=np.bool)
    aminus=np.ones((X.shape[0]),dtype=np.bool)
    for i1 in range(n):
        i2=(i1+1) % n
        a=X[:,i2]*Y[:,i1]-X[:,i1]*Y[:,i2]
        aplus=aplus & (a>0)
        aminus=aminus & (a<0)
    return aplus | aminus


##################################################
# order arrays into a convex polygon
def order(X,Y,x0,y0):
    """
    result=order(X,Y,x0,y0)

    X,Y and x0,y0 should be 1D arrays
    based on:
    http://demonstrations.wolfram.com/AnEfficientTestForAPointToBeInAConvexPolygon/
    """
    X=X-x0
    # account for wraps in RA
    if (X>180).any():
        X[X>=180]-=360
    if (X<-180).any():
        X[X<=-180]+=360
    # we won't be able to deal with wraps in Dec
    # so we don't want to do this for very close to the poles
    Y=Y-y0
    angles=np.arctan2(Y,X)
    return np.argsort(angles)


##################################################
# return the ZTFtile object associated with a given tile number
def get_tile(number, missing=None):
    """
    Z=get_tile(number, missing=None)

    return a ZTFtile object corresponding to the number in the tile grid
    """
    try:
        tiles=Table.read(os.path.join(data._datadir,'ztf_fields.dat'),format='ascii.commented_header')
    except IOError:
        raise IOError('Unable to read ZTF field file %s' % os.path.join(data._datadir,'ztf_fields.dat'))
    try:
        return ZTFtile(tiles[number]['RA']*u.deg,tiles[number]['Dec']*u.deg,number=number,missing=missing)
    except IndexError:
        print('Tile %s is not available' % number)
        return None

def metadata2wcs(metadata):
    """
    wcs=metadata2wcs(metadata)
    converts ZTF metadata records
    like those at https://github.com/MickaelRigault/ztfquery
    to wcs objects (one for each row)
    """
    w=[]
    for i in range(len(metadata)):
        w.append(naxis=2)
        w[-1].wcs.ctype=['RA---TAN','DEC--TAN']
        w[-1].wcs.crval=[metadata['crval1'][i],metadata['crval2'][i]]
        w[-1].wcs.crpix=[metadata['crpix1'][i],metadata['crpix2'][i]]
        w[-1].wcs.cd=[[metadata['cd11'][i],metadata['cd12'][i]],
                      [metadata['cd21'][i],metadata['cd22'][i]]]
    if len(metadata)==1:
        return w[0]
    return w
##################################################
# class ZTFtile
##################################################
class ZTFtile:

    _quadrantcenters='ZTF_quadrantcenters_20180813.dat'
    # make sure it can find the table of quadrant centers
    try:
        t=Table.read(os.path.join(data._datadir,_quadrantcenters),format='ascii.commented_header')
    except IOError:
        raise IOError('Unable to read ZTF quadrant center file %s' % os.path.join(data._datadir,_quadrantcenters))
    # Cartesian coordinates in the plane of projection
    quadrant_x=t['DX']*u.deg
    quadrant_y=t['DY']*u.deg
    quadrant_size=np.array([3072,3080])
    try:
        quadrant_scale=t['scale']*u.arcsec
    except KeyError:
        # updated from 1"/pixel to 1.0126 based on data
        # this is the mean across all Dec ranges and quadrants
        # DLK, 2018-08-13
        quadrant_scale=1.0126*u.arcsec*np.ones(64)

    def __init__(self, RA, Dec, number=None, missing=None, rotation='model'):
        """
        z=ZTFtile(RA, Dec, number, missing, rotation='model')
        creates a ZTFtile object centered at the given coordinates

        coordinates can be Quantities (with units) or floats (assumed degrees)

        number is just used for informational purposes (bookkeeping)

        missing lists the missing quadrants (potentially None or [])
        if missing is not None, should be a list of quadrants (0--63) or CCD,quadrant pairs (1-16,1-4)
        they will be removed from the quadrant-by-quadrant WCS

        rotation is one of 'model' or 'none'
        If 'none', then all quadrants have 0 rotation
        if 'model', assumes rotation is X*tan(Dec)
        where X is X offset in projection plane of quadrant center and Dec is quadrant center Dec
        

        """
        if isinstance(RA, astropy.units.quantity.Quantity):            
            self.RA=RA
        else:
            self.RA=RA*u.deg
        if isinstance(Dec, astropy.units.quantity.Quantity):            
            self.Dec=Dec
        else:
            self.Dec=Dec*u.deg
        assert rotation.lower() in ['model','none',None]
        self.rotation=rotation.lower()

        self.number=number

        self.quadrant_RA,self.quadrant_Dec=self.quadrant_centers()
        if self.rotation == 'model':
            self.quadrant_rotation=self.quadrant_x*np.tan(self.quadrant_Dec)
        else:
            self.quadrant_rotation=0*u.deg*np.ones(len(self.quadrant_x))
        
        self.missing_quadrants=[]
        if missing is not None and len(missing)>0:
            for m in missing:
                if isinstance(m,int) or len(m)==1:
                    # assume this is 0 origin
                    self.missing_quadrants.append(m)
                elif len(m)==2:
                    # this assumes CCD and quadrant numbers are 1-origin
                    self.missing_quadrants.append(4*(m[0]-1)+(m[1]-1))
        self.good_quadrants=np.setdiff1d(np.arange(len(self.quadrant_RA)),
                                         np.array(self.missing_quadrants))
        
        self._wcs=[]
        for i in range(len(self.quadrant_RA)):
            if i in self.good_quadrants:
                self._wcs.append(self.quadrant_WCS(self.quadrant_RA[i],self.quadrant_Dec[i],self.quadrant_scale[i],self.quadrant_rotation[i]))
            else:
                self._wcs.append(None)

    def quadrant_WCS(self, quadrant_RA, quadrant_Dec, scale, rotation):
        """
        w=quadrant_WCS(quadrant_RA, quadrant_Dec, scale, rotation)
        returns a WCS object that is specific to the quadrant RA and Dec specified
        scale, rotation  also specified

        CD = scale * [[-cos(rotation),sin(rotation)],
                      [-sin(rotation),-cos(rotation)]]

        size are determined by class variables
        """
        
        w=WCS(naxis=2)
        w.wcs.crpix=(self.quadrant_size+1.0)/2
        w.wcs.ctype=['RA---TAN','DEC--TAN']
        w.wcs.crval=[quadrant_RA.value,
                     max(min(quadrant_Dec.value,90),-90)]
        w.wcs.cd=[[-scale.to(u.deg).value*np.cos(rotation),-scale.to(u.deg).value*np.sin(rotation)],
                  [scale.to(u.deg).value*np.sin(rotation),-scale.to(u.deg).value*np.cos(rotation)]]
        return w

    def projplane2celestial(self, x, y):
        """
        alpha,delta=projplane2celestial(x,y)
        convert from projection plane (TAN) to celestial RA,Dec
        """

        # convert to Native longitude (phi) and latitude (theta)
        # need intermediate step (Rtheta) to deal with TAN projection
        # Calabretta & Greisen (2002), Eqn. 14,15
        phi=np.arctan2(x,-y)
        Rtheta=np.sqrt(x**2+y**2)
        # now to theta using the TAN projection
        # Calabrett & Greisen (2002), Eqn. 55
        theta=np.arctan2(1,Rtheta.to(u.rad).value)*u.rad
        
        # central point of projection (native pole)
        alphap=self.RA
        deltap=self.Dec
        # Native longitude/latitude of fiducial point
        # appropriate for zenithal projections including TAN
        phi0=0*u.deg
        theta0=90*u.deg
        # Native longitude/latitue of celestial pole
        # for delta0<theta0 then phip should be 180 deg
        # and theta0=90 deg
        phip=180*u.deg
        thetap=self.Dec
        
        # Celestial longitude/latitude
        # Calabretta & Greisen (2002), Eqn. 2
        alpha=alphap+np.arctan2(-np.cos(theta)*np.sin(phi-phip),
                                np.sin(theta)*np.cos(deltap)-np.cos(theta)*np.sin(deltap)*np.cos(phi-phip))
        
        delta=(np.arcsin(np.sin(theta)*np.sin(deltap)+np.cos(theta)*np.cos(deltap)*np.cos(phi-phip))).to(u.deg)
        # fix RA wraps
        alpha[alpha<0*u.deg]+=360*u.deg
        return alpha,delta

    def celestial2projplane(self, alpha, delta):
        # central point of projection (native pole)
        alphap=self.RA
        deltap=self.Dec
        # Native longitude/latitude of fiducial point
        # appropriate for zenithal projections including TAN
        phi0=0*u.deg
        theta0=90*u.deg
        # Native longitude/latitue of celestial pole
        # for delta0<theta0 then phip should be 180 deg
        # and theta0=90 deg
        phip=180*u.deg
        thetap=self.Dec

        # convert to native long, lat
        # Calabretta & Greisen, 2002
        # Eqn. 5
        phi = phip + np.arctan2(-np.cos(delta)*np.sin(alpha-alphap),
                                np.sin(delta)*np.cos(deltap) - np.cos(delta)*np.sin(deltap)*np.cos(alpha-alphap))
        theta = np.arcsin(np.sin(delta)*np.sin(deltap)+np.cos(delta)*np.cos(deltap)*np.cos(alpha-alphap))

        # now to projection plane
        # C&G 2002, Eqn. 54
        Rtheta=(180*u.deg/np.pi)/np.tan(theta)
        # and usine Eqn. 12,13
        x=(Rtheta*np.sin(phi)).to(u.deg)
        y=-(Rtheta*np.cos(phi)).to(u.deg)

        return x,y
        

    def quadrant_centers(self):
        """
        alpha,delta=quadrant_centers()
        return celestial coordinates for the centers of each quadrant
        given the pointing center RA,Dec
        """
        return self.projplane2celestial(self.quadrant_x,self.quadrant_y)

    def _inside(self, RA, Dec):
        """
        _inside(RA,Dec)
        returns whether a given RA,Dec (floats in degrees) is inside the FOV

        returns an integer for each coordinate.
        This will be 0 if it is not inside.  Otherwise it will map to the quadrant:
        n=4*(CCD-1)+Quadrant
        with CCD=[1,..16]
        and Quadrant=[1,..4]
        """
        on=np.zeros_like(RA,dtype=np.int64)
        for i in range(len(self._wcs)):
            if self._wcs[i] is not None:
                footprint=self._wcs[i].calc_footprint(axes=self.quadrant_size)
            on[inside(footprint[:,0],footprint[:,1],RA,Dec)]=i+1
            
        return on

    def _inside_nogaps(self, RA, Dec):
        """
        _inside_nogaps(RA,Dec)
        returns whether a given RA,Dec (floats in degrees) is inside the rough FOV
        ignores chipgaps
        """
        # get all of the corners of the quadrants
        footprint=np.zeros((4*len(self._wcs),2))
        for i in range(len(self._wcs)):
            footprint[(4*i):(4*(i+1)),:]=self._wcs[i].calc_footprint(axes=self.quadrant_size)
        corners=SkyCoord(footprint[:,0]*u.deg,footprint[:,1]*u.deg)
        # distances from these corners to the central pointing position
        corner_distances=corners.separation(SkyCoord(self.RA,self.Dec))
        # find the 4 furthest points - these will define the outer boundaries
        outer_corners=footprint[np.argsort(corner_distances)[-4:],:]
        # sort them by angle to make sure they are convex
        corner_order=order(outer_corners[:,0],outer_corners[:,1],
                           self.RA.value,self.Dec.value)[::-1]
        outer_corners=outer_corners[corner_order,:]
        on=inside(outer_corners[:,0],outer_corners[:,1],RA,Dec)
        return (on>0)

    def inside(self, *args):
        """
        inside(*args)
        returns whether a given coordinate is inside the FOV
        coordinate can be SkyCoord
        or separate RA,Dec
        RA,Dec can be Quantity (with units) or floats (in degrees)

        returns an integer for each coordinate.
        This will be 0 if it is not inside.  Otherwise it will map to the quadrant:
        n=4*(CCD-1)+Quadrant
        with CCD=[1,..16]
        and Quadrant=[1,..4]
        """

        if len(args)==1:
            if isinstance(args[0],astropy.coordinates.sky_coordinate.SkyCoord):
                return self._inside([args[0].ra.value],[args[0].dec.value])
        elif len(args)==2:
            if isinstance(args[0],astropy.units.quantity.Quantity):
                return self._inside([args[0].value],[args[1].value])
            else:
                return self._inside(args[0],args[1])
            
    def inside_nogaps(self, *args):
        """
        inside_nogaps(*args)
        returns whether a given coordinate is inside the rough FOV
        ignoring the chip gaps
        coordinate can be SkyCoord
        or separate RA,Dec
        RA,Dec can be Quantity (with units) or floats (in degrees)
        """

        if len(args)==1:
            if isinstance(args[0],astropy.coordinates.sky_coordinate.SkyCoord):
                return self._inside([args[0].ra.value],[args[0].dec.value])
        elif len(args)==2:
            if isinstance(args[0],astropy.units.quantity.Quantity):
                return self._inside_nogaps([args[0].value],[args[1].value])
            else:
                return self._inside_nogaps(args[0],args[1])

    def region(self, filename, color='green'):
        s="icrs;point(%s,%s) # point=diamond text={Center} color={%s}\n" % (self.RA.value,self.Dec.value,color)        
        
        for i in range(len(self._wcs)):
            if self._wcs[i] is not None:
                f=self._wcs[i].calc_footprint(axes=self.quadrant_size)
                x=np.append(f[:,0],f[0,0])
                y=np.append(f[:,1],f[0,1])
                for j in range(len(x)-1):
                    s+="icrs;line(%s,%s,%s,%s) # line=0 0 color={%s}\n" % (x[j],y[j],x[j+1],y[j+1],color)
        if filename is not None and len(filename)>0:
            try:
                fo=open(filename,'w')
            except IOError:
                raise IOError("Could not open region file %s for writing" % filename)
            fo.write(s)
        return s

class HP_coverage:

    def __init__(self, hpfile, nside=256):
        """
        HP_coverage(hpfile, nside=256)
        healpix coverage map
        input healpix file and desired nside for degradation (if needed)
        """

        import healpy

        self.nside=nside
        self.npix=healpy.nside2npix(self.nside)

        # read the healpix map
        try:
            b,h=healpy.read_map(hpfile,h=True,nest=False)
        except IOError:
            raise IOError('Unable to read healpix file %s' % hpfile)
        
        # convert the header to a dictionary
        h={x[0]: x[1] for x in h}
        self.nside_in=h['NSIDE']
        self.npix_in=healpy.nside2npix(self.nside_in)

        # degrade if necessary
        if self.nside < self.nside_in:
            self.hpmap=healpy.ud_grade(b, self.nside, power=-2)
        else:
            self.hpmap=b
            self.nside=self.nside_in
            self.npix=healpy.nside2npix(self.nside)

        # coordinates for each pixel in the map
        RA,Dec=healpy.pix2ang(self.nside, np.arange(self.npix),lonlat=True)
        self.RA=RA*u.deg
        self.Dec=Dec*u.deg

    def tile_coverage(self, ztftiles):
        """
        covered=tile_coverage(ztftiles)
        returns an array the same size as the healpix map
        composed of booleans to show which pixels are covered
        ztftiles is list of ZTFtile objects()
        """

        try:
            l=len(ztftiles)
        except:
            ztftiles=[ztftiles]

        covered=np.zeros(len(self.hpmap),dtype=np.bool)
        for itile in range(len(ztftiles)):
            covered[ztftiles[itile].inside(self.RA,self.Dec)]=True
        return covered
    
    def tile_probability(self, ztftiles):
        """
        probability=tile_probability(tile_RA,tile_Dec)
        ztftiles is list of ZTFtile objects()
        returns total covered probability
        """

        covered=self.tile_coverage(ztftiles)
        return (self.hpmap*covered).sum()

    def get_tile_values(self, hpmap, ztftiles):
        """
        values=get_tile_values(hpmap, ztftiles)
        gets the values used for sorting the tile contributions
        this is implemented as the integral of the probability for the given map over each tile
        ztftiles is list of ZTFtile objects()
        """
        
        values=np.zeros(len(ztftiles))
        for itile in range(len(ztftiles)):
            values[itile]=hpmap[ztftiles[itile].inside(self.RA,self.Dec)].sum()
        return values
    
    def find_tiles(self, ztftiles, probability_target=0.9, verbose=False):
        """
        tiles,probability=find_tiles(ztftiles, probability_target=0.9, verbose=False)

        give ztftiles (list of ZTFtile objects())

        will find the combination to achieve the total probability target
        """

        # starting values and order
        tile_values=self.get_tile_values(self.hpmap, ztftiles)
        tile_order=np.argsort(tile_values)[::-1]

        tiles=[]
        # this is a copy of the map that we can modify
        hpmapc=np.copy(self.hpmap)
        covered=np.zeros(len(self.hpmap),dtype=np.bool)
        summed_probability=0
        while summed_probability < probability_target and len(tiles)<len(ztftiles):
            # add it to the list
            tiles.append(tile_order[0])

            individual_prob=(hpmapc[ztftiles[tile_order[0]].inside(self.RA,self.Dec)].sum())
            if verbose:
                print('Adding tile %d (%f,%f): individual probability = %.3f, total probability = %.3f' % (tile_order[0],
                                                                                                           ztftiles[tile_order[0]].RA.value,
                                                                                                           ztftiles[tile_order[0]].Dec.value,
                                                                                                           individual_prob,
                                                                                                           summed_probability+individual_prob))
            covered[ztftiles[tile_order[0]].inside(self.RA,self.Dec)]=True
            summed_probability+=individual_prob
            hpmapc[ztftiles[tile_order[0]].inside(self.RA,self.Dec)]=0
                
            # redo the priorities to account for the new 
            # probability values
            tile_values=self.get_tile_values(hpmapc, 
                                             ztftiles)
            # priorities for each of those
            tile_order=np.argsort(tile_values)[::-1]

        return tiles,summed_probability

    def percentile(self, level):
        """
        value=percentile(level)
        determines the value of the map above which there is level of total probabiilty
        """
        m=np.sort(self.hpmap)[::-1]
        msum=np.cumsum(m)
        include=msum<level*m.sum()
        # add one more point just to be sure
        if include.sum()<len(m):
            include[include.sum()]=True
        return m[include][-1]

    def percentile_area(self, level):
        """
        area=percentile_area(level)
        determines the area of the map above which there is level of total probabiilty        
        """
        
        return (self.hpmap>=self.percentile(level)).sum()*healpy.nside2pixarea(self.nside,degrees=True)*u.deg**2
