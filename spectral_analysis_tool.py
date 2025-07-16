'''
Author: Harry Addison

Code to analyse SNe spectra. 
'''

import os
import sys
import matplotlib.pyplot as plt
from matplotlib.backend_bases import MouseButton
import numpy as np
from astropy.table import QTable
import modules.input_output as io
from modules.functions import vel, calc_pew, interpolate
from modules.preprocess import norm_flux, gal_extinct_remove


# Load in input parameters
params = io.load_config()

c = 299792.458 # km/s

# Load the info of the produced SNe spectra
sne_info = io.load_sn_info(params["sn_info_path"])

# Load in the features to be measured
if params["features"] == "Ia":
    fn = f"{os.path.dirname(__file__)}/data/SN_Ia_features.ecsv"
else:
    fn = params["features"]
features = QTable.read(fn, format="ascii.ecsv")

# Load in measurements table
# First see if there is a file for previously analysed spectra
# if not then load the output from the l1 spectrum code and add columns.
try:
    sne_measure = io.load_sn_info(params["sn_measure_path"])
except FileNotFoundError:
    names = ["id"]
    names += list(np.array(list(([f"{feat}_vel", f"{feat}_pew"] for feat in features["feature"]))).flatten())
    sne_measure = QTable(names=names)

# Print the instructions on how to analyse the spectrum.
print('''
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      Program user instructions:
      ~~~~~~~~~~~~~~~~~~~~~~~~~~

      This program is for semi-automated spectral analysis of supernovae
      spectra.
      The user will input a suitable pseudo continuum for a certain
      feature and then a Gaussian will be fitted to the absorption
      feature. 
      The velocity will then be obtained from the fitted Gaussian
      through the assumption that the Doppler shift and broadening are
      are the dominant contibutors of the absorption feature shift and
      broadening.

    
      For each feature follow these instructions:
        
        1. Select on the spectrum plot the two end points for the 
           pseudo continuum.
           Please select the lower wavelength point followed by the 
           upper wavelength point.
           If the absorption feature does not exist, place the points 
           at negative fluxes.
                > spacebar: place point at cursor location
                > backspace: remove the last point placed

        2. Press middle mouse button once happy with the point placement.

        3. Wait for a Gaussian profile to be fitted to the feature.

        4. Check that the Gaussian fit is suitable. 
                > type in terminal "y" if fit is good.
                > type in terminal "n" if fit is bad.

        5. > If bad fit: repeat steps 1-4 until the fit is good.
           > If good fit: continue to next absorption feature/spectrum.
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
      ''')

# Run through each of the spectra produced
for i, id in enumerate(sne_info["id"]):

    # If the spectrum has previously been analysed then skip it.
    if id not in sne_measure["id"]:

        print("Analysing SNe ID: ", id)
        row = [id]

        # Obtain the template spectrum.
        template = io.load_template("%s/template_%s.fits"%(params["sn_spec_dir"],
                                                           id))

        # Remove Galactic line of sight extinction.
        template = gal_extinct_remove(template, sne_info["ra"][i], sne_info["dec"][i])

        # Convert to rest wavelength
        template["wave"] /= (sne_info["z"][i] + 1)

        # Normalise the spectrum
        template = norm_flux(template, error=False)

        # For each of the features provided in the input parameter file,
        # measure their properties.
        for feat in features:

            # convert the wavelength space of the feature into velocity
            template["vel"] = -vel(template["wave"], feat["rest_wl"])

            # Plot the spectrum in velocity space
            plt.figure()
            plt.plot(template["vel"], template["flux"])

            plt.xlabel(r"$\rm{Velocity~relative~to~H}\alpha~\rm{rest}$")
            plt.ylabel("Normalised Flux")

            plt.title("SNe ID: %s\nFeature to measure: %s"
                    %(id, feat["feature"]))
            plt.xlim(-20000, 20000)

            mask = (template["vel"]>-20000) & (template["vel"]<20000)
            plt.ylim(min(template[mask]["flux"]), max(template[mask]["flux"]))
            plt.get_current_fig_manager().full_screen_toggle()

            # Get the user to select 2 points in which the maxima and
            # minima of the feature's P-cygni profile exists.
            selected_points = False

            while selected_points is False:
                print("\n\n Please select only 2 end points in which"
                      "the minima of the %s absorption feature exists"
                      "\n\n"%(feat["feature"]))

                region = plt.ginput(n=100, timeout=0,
                                    show_clicks=True,
                                    mouse_add=None,
                                    mouse_pop=None,
                                    mouse_stop=MouseButton.MIDDLE)

                if len(region) != 2:
                    continue

                if len(region) == 2:
                    selected_points = True

                    plt.close()

            # reformat the points data into velocity and flux columns.
            region = np.array(region).transpose()

            # Check that the feature exists (i.e +ve continuum fluxes).
            # if doesn't exist then add nan values to master table.
            if region[1][0] < 0 or region[1][1] < 0:
                row.append(np.nan)  # velocity
                row.append(np.nan)  # pew

            # if exists then find the velocity from the minima.
            else:
                # trim data spectrum to within given wavelengths
                mask = ((template["vel"] > region[0][0]) &
                        (template["vel"] < region[0][1]))

                feat_spec = template[mask]

                # Measure velocity from minimum.
                min_flux_ind = np.argmin(feat_spec["flux"])
                velocity = feat_spec["vel"][min_flux_ind]

                # Measure pew
                continuum = interpolate(feat_spec["wave"], feat_spec,
                            col_names=["wave", "flux"])

                pew = calc_pew(feat_spec, continuum, keys=["wave", "flux"])

                # Add the best fit velocity and pew to the row.
                # Invert the velocity. Originally defined it so that a redshifted
                # line was == +vel while blueshifted line == -vel. 
                # Want the definition of vel to be relative to the centre of
                # the SN and so +vel == moving outwards from the centre.
                # As blueshifted line == moving outwards from the centre of SN
                # (i.e towards us), make this the +vel. Therefore need to invert
                # the original vel.
                row.append(-velocity)
                row.append(pew)
                print(row)

        # Add the row to measurements table.
        sne_measure.add_row(row)
        sne_measure.pprint(max_lines=-1, max_width=-1)

        # Save the updated table
        print("\n\nUpdating and saving the velocity file\n\n")
        sne_measure.write(params["sn_measure_path"], overwrite=True)
