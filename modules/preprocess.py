# -*- coding: utf-8 -*-
'''
Author: Harry Addison
Created: 18/01/23

Calculations for applying corrections to spectra
before any analysis is carried out.
'''


import numpy as np
from astropy.coordinates import SkyCoord
from dustmaps.sfd import SFDQuery
from extinction import fitzpatrick99, remove
# from gpr_model import make_model, model_values


def norm_flux(spec, error=True):
    '''
    Normalise the given spectrum's flux and flux error values according
    to the maximum flux value.

    Input:
        > "spec" = spectrum data to be normalise, should contain at
                   least 2 columns named "FLUX" and "ERR_FLUX".

    Output:
        > "spec" = Normalised spectrum data in the same format as
                   the input "spec".
    '''

    # Find maximum flux.
    max_flux = spec["flux"].max()

    # Normalise the flux and flux errors.
    spec["flux"] /= max_flux

    if error == True:
        spec["flux_err"] /= max_flux

    return spec


def gal_extinct_remove(spec, ra, dec):
    '''
    Obtain the line of sight Milky Way extinction at the given
    coordinates using a dust map. Correct the given spectrum for the
    obtained Galactic line of sight extinction.

    Inputs:
        > "spec" = Spectrum to apply extinction correction to.
        > "coords" = coordinates used to obtain the line of sight
                     extinction.

    Output:
        > "spec" = Extinction corrected spectrum.
    '''

    # Store the coordinates in astropy SkyCoord object.
    coord = SkyCoord(ra, dec, frame='icrs')

    # Query the Schlegel, Finkbeiner & Davis dust map at the given
    # coords for extinction (ebv).
    ebv = float(SFDQuery()(coord))

    # Find the extinction at the wavelengths of spectrum.
    r_v = 3.1
    wave = np.array(spec["wave"], dtype=float)
    wl_extinct = fitzpatrick99(wave, a_v=(r_v * ebv),
                               r_v=r_v, unit="aa")

    # Apply correction for ebv to spectrum.
    spec["flux"] = remove(wl_extinct, spec["flux"])

    return spec
