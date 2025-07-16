# -*- coding: utf-8 -*-
'''
Author: Harry Addison
Created: 02/12/22

Input and output functions:

    > Load a 4MOST L1 spectrum from a fits file and
      make sure its in a more useable format i.e
      dataframe.
    > Load data for a given spectral region. Data
      includes the name, rest wavelength,
      lower bound wavelength range, upper bound
      wavelength range.
'''

import argparse
import textwrap
import yaml
import astropy.units as u
from astropy.table import QTable
from astropy.table import Table


def parser():
    '''
    Function that allows the user to parse a initfile path via the
    command line.
    '''

    help_msg = "A help message and description of the code"  #TODO

    parser = argparse.ArgumentParser(prog="TiDES Spectral Analysis",
                                     formatter_class=argparse.RawDescriptionHelpFormatter,
                                     description=textwrap.dedent(help_msg))

    parser.add_argument("-i",
                        "--initfile",
                        default=None,
                        help=("yaml file containing input frames "
                              "and other procedure arguments"))
    return parser


def load_config(initfile=None):
    '''
    Obtain the config parameters from the given yaml file. If no file path provided
    try to see if a yaml file was parsed via the command line.
    '''

    if initfile == None:
        arg_parser = parser()
        args = arg_parser.parse_args()
        initfile = args.initfile
        print(f"\nInput parameter file: {initfile}\n\n")

        # If no parameter file was provided print help message
        if initfile is None:
            arg_parser.print_help()
            exit()

    # Open the initfile and read in the parameters
    with open(str(initfile), "r") as file:
        config = yaml.safe_load(file)

    return config


def load_spec(spec_fn):
    '''
    Load in the spectrum dataframe from the given file.
    Inputs:
    > "spec_fn" = spectrum data file path.
    Outputs:
    > "spec" = table/dataframe of the spectrum data. Has column names:
               "WAVE" (wavelength), "FLUX", "ERR_FLUX" (flux error).
    '''

    spec = Table.read(spec_fn)

    # Each column contains one row with an array of that column's data
    # Change format of data so each data point is a row.
    # Also remove unessasary columns.
    wl = spec["WAVE"][0] * spec["WAVE"].unit
    flux = spec["FLUX"][0] * spec["FLUX"].unit
    flux_err = spec["ERR_FLUX"][0] * spec["ERR_FLUX"].unit

    spec = Table(data=[wl, flux, flux_err],
                 names=["WAVE", "FLUX", "ERR_FLUX"])

    return spec


def load_template(fn):
    '''
    Load in the spectrum dataframe from the given file.
    Inputs:
    > "fn" = spectrum data file path.
    Outputs:
    > "template" = table/dataframe of the spectrum data. Has column names:
                   "WAVE" (wavelength), "FLUX", "ERR_FLUX" (flux error).
    '''

    template = Table.read(fn)

    return template


def load_sn_info(file_path):
    '''
    Obtain the information of the supernovae from the file at the
    given file path.

    Inputs:
        > "file_path" = path corresponding to the file containing the
                        information of the given supernova.

    Outputs:
        > "info" = table containing the information of the supernovae.
    '''

    data = QTable.read(file_path)

    return data


def save_spec(spec, save_path):
    '''
    Save the processed spectrum to a fits file.

    Inputs:
    > "spec" = Table of the spectrum.
    > "save_path" = path to save the fits file to.
    '''

    spec.write(save_path, format="fits", overwrite=True)

    return


def load_telluric(file_path, z=None):
    '''
    Load in the telluric absorption regions. Apply a de-redshift to
    convert the wavelength regions into "rest" frame (optional).

    Input:
    > "file_path" = path to the directory where the telluric absorption
                    region wavelength ranges are saved.
    > "z" = Redshift to use to convert wavelength regions to "rest"
            frame. If "z" = None then the wavelengths are not
            de-redshifted into "rest" frame.

    Output:
    > "regions" = table containing the telluric absorption feature
                  blue-ward and red-ward wavelength bounds.
    '''

    regions = QTable.read(file_path, format="csv", comment="#", header_start=0)

    # De-redshift wavelengths to "rest" frame.
    if z is not None:
        regions["blue_wl"] /= (z + 1)
        regions["red_wl"] /= (z + 1)

    return regions
