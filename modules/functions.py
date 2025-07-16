from astropy.table import Table
from scipy.integrate import simpson



def vel(min_wl, rest_wl):
    '''
    Relativistic doppler shift equation rearranged to get the velocity.
    '''

    c = 299792.458 # km/s

    wl_ratio_2 = (min_wl / rest_wl)**2

    vel = -c * (wl_ratio_2 - 1) / (wl_ratio_2 + 1)

    return vel


def interpolate(new_x, data, col_names=["x", "y"]):
    '''
    Linear interpolation between two points with error.

    Inputs:
    > "new_x" = array of x values to obtain y and y_err data for.
    > "data" = table/dataframe that is sorted in accending x is
               stuctured:
                    > data[0] = x values
                    > data[1] = y values
    > "col_names" = names to call the columns of the output table.

    Outputs:
    > "new_data" = table/dataframe of the interpolated data. Has column
                   names given by the input "col_names".
    '''

    # Min and max x data points
    p1 = data[:][0]
    p2 = data[:][-1]

    # Interpolate using y = mx + c and solving simultaneous equations.
    # (Has been simplified to one line.)
    new_y = (((p1[1] - p2[1]) / (p1[0] - p2[0])) * (new_x - p1[0])) + p1[1]
    
        # Creating the output table of the interpolated data.
    new_data = Table(data=[new_x, new_y],
                     names=col_names)

    return new_data


def calc_pew(spec, continuum, keys=["x", "y"]):
    integrand = 1 - (spec[keys[1]].value / continuum[keys[1]].value)
    pew = simpson(integrand, spec[keys[0]].value)
    return pew
