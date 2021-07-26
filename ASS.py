""" This Amazing Scaling Script (ASS) allows the user to linearly scale a
spectrum to match synthetic photometry obtained from s3/SMS
(https://github.com/cinserra/S3). """

# Need to enter a python2.7/iraf environment
# Loading in necessary libraries
import numpy as np
import math
import matplotlib.pyplot as plt
import spectres
import os
import shutil

""" User inputs section """
spectrum_name = "some-spectrum-filename-123.dat"

# These are the filters to base the linear scaling function on [can only be 2!]
filter_list = ["PS1/rs", "PS1/is"]
# These are the additional filters to plot in the final comparison plot
# - can be as many as desired, but do not need to repeat those in list above
full_filter_list = ["PS1/gs", "PS1/zs"]

# name for the file used by SMS to read in list of spectra
spec_filename = "spectra.txt"
# name for the file used by SMS to read in list of filters
filters_filename = "filters.txt"

# specify the path to the metadata folder with the different filters
metadata_path = "../S3-master/src/s3/metadata"
""" end of user inputs section """


def _mkdir(newdir):
    """ this function is for creating the required directory for the output
    files """
    """works the way a good mkdir should :)
        - already exists, silently complete
        - regular file in the way, raise an exception
        - parent directory(ies) does not exist, make them as well
    """
    if os.path.isdir(newdir):
        pass
    elif os.path.isfile(newdir):
        raise OSError(
            "a file with the same name as the desired "
            "dir, '%s', already exists." % newdir
        )
    else:
        head, tail = os.path.split(newdir)
        if head and not os.path.isdir(head):
            _mkdir(head)
        # print "_mkdir %s" % repr(newdir)
        if tail:
            os.mkdir(newdir)
    return newdir


def spectres_fn(spectrum, binning_factor):
    """ This function is used to rebin the spectrum. """
    binned_steps = spectrum[1::binning_factor, 0]

    #  This rebins the spectrum by the requested amount.
    spectrum_binned = spectres.spectres(
        binned_steps, spectrum[:, 0], spectrum[:, 1]
    )

    # This splices the separate 1D arrays together
    binned_spectrum = np.c_[binned_steps, spectrum_binned]
    # This drops any nan-related rows from the array
    binned_spectrum = binned_spectrum[~np.isnan(binned_spectrum).any(axis=1)]

    return binned_spectrum


def read_phot_obs(filter_list):
    """ This function is used to read in the OBSERVED photometry values from
    the keyboard; to be input by the user.
    For clarity:
        The photometric values are the absolute values that have been observed.
        The synthetic values are those produced by SMS.py from the uncalibrated
        spectra. """

    # Loading in the photometric filter magnitudes for the spectra
    observed_photometry = np.zeros(len(filter_list))
    observed_photometry_errors = np.zeros(len(filter_list))
    for aa in range(len(filter_list)):
        observed_photometry[aa] = float(
            input(
                "What is the "
                + filter_list[aa]
                + " PHOTOMETRIC magnitude for this epoch?\t\t"
            )
        )
        observed_photometry_errors[aa] = float(
            input("And its error?\t\t\t\t\t\t\t")
        )

    return observed_photometry, observed_photometry_errors


def read_phot_syn(filter_list):
    """ This function is used to read in the SYNTHETIC photometry values from
    the keyboard; to be input by the user.
    For clarity:
        The photometric values are the absolute values that have been observed.
        The synthetic values are those produced by SMS.py from the uncalibrated
        spectra. """

    # Loading in the synthetic filter magnitudes for the spectra
    synthetic_photometry = np.zeros(len(filter_list))
    for aa in range(len(filter_list)):
        synthetic_photometry[aa] = float(
            input(
                "What is the "
                + filter_list[aa]
                + " SYNTHETIC magnitude for this epoch?\t\t\t"
            )
        )

    return synthetic_photometry


def read_filter_files(filter_list, metadata_path):
    """ This function reads the filter file from the metadata and extracts the
    effective peak, and the range of wavelength coverage of the filter.
    Note: the wavelength range is not currentyl used for anything. """

    filter_properties = np.zeros((len(filter_list), 3))
    for aa in range(len(filter_list)):

        file = open(metadata_path + "/" + filter_list[aa] + ".txt", "r")
        lines = file.readlines()

        filter_peak = lines[2].strip("\n")
        filter_min_lam = lines[4].split(" ")[0]
        filter_max_lam = lines[-1].split(" ")[0]

        filter_properties[aa] = filter_peak, filter_min_lam, filter_max_lam

        file.close()

    return filter_properties


def write_spec_list_file(filter_list, file_string, spec_filename):
    """ This function writes the name of the intermediate spectrum filename
    to the spec_filename file for use with SMS. """

    file = open(spec_filename, "w")

    for aa in range(len(filter_list)):
        file.write(file_string)
        file.write("\n")

    file.close()

    return


def write_filter_list_file(filter_list, filters_filename):
    """ This function writes the name of the filters under investigation
    to the filters_filename file for use with SMS. """

    file = open(filters_filename, "w")

    for aa in range(len(filter_list)):
        file.write(filter_list[aa])
        file.write("\n")

    file.close()

    return


def plot_filter_bandpasses(filter, metadata_path):
    """ This function reads the filter file from the metadata and extracts the
    full bandpass information of the filter (for plotting purposes). """

    file = open(metadata_path + "/" + filter + ".txt", "r")
    lines = file.readlines()

    array = lines[4:]

    waves = []
    flux = []
    for aa in range(len(array)):
        waves.append(float(array[aa].split(" ")[0]))
        flux.append(float(array[aa].split(" ")[1]))

    file.close()

    return np.c_[waves, flux]


""" Start of the main body of code """
# We begin with an observed spectrum which is not flux-calibrated.
# We also have observed photometry that covers the spectral range.

# Throughout these code comments, I will for simplicity consider the g and r
# band filters, but the user can choose any two filters to use.

#  Splitting the spectrum filename and the extension. Necessary for naming the
# intermediate calibration files later.
specname = spectrum_name[:-4]
extension = spectrum_name[-4:]

# Loading in the spectrum to be flux-calibrated.
spec = np.loadtxt(spectrum_name)

# Some code to re-bin the spectrum - not really needed at this point.
"""
# loading in the spectra to be calibrated
spec_raw = np.loadtxt(spectrum_name)

# this bins the spectrum for better resolution
spec = spectres_fn(spec_raw, 10)

print len(spec_raw), len(spec)

plt.plot(spec_raw[:,0], spec_raw[:,1])
plt.plot(spec[:,0], spec[:,1], alpha = 0.5)

plt.show()
"""

# Creating files for use with SMS - these file contain lists of spectra and
# filters for computing synthetic magnitudes.
write_spec_list_file(filter_list, spectrum_name, spec_filename)
write_filter_list_file(filter_list, filters_filename)

# Running the functions to obtain user inputs for observed and synthetic
# photometry. Code first reads in the observed photometry (g_obs and r_obs), and
# the synthetic photometry obtained from SMS for the spectrum (g_syn and r_syn).
observed_photometry, observed_photometry_errors = read_phot_obs(filter_list)
synthetic_photometry = read_phot_syn(filter_list)

# Then we calculate the differences in observed bands
# (delta_obs = g_obs - r_obs),
# and the differences in synthetic bands (delta_syn = g_syn and r_syn).
delta_phot_obs = observed_photometry[0] - observed_photometry[1]
delta_phot_syn = synthetic_photometry[0] - synthetic_photometry[1]

# Then we compute the differences between these.
delta_phot_obs_phot_syn = delta_phot_obs - delta_phot_syn

# Here we read in some filter band info from the metadata files.
filter_band_info = read_filter_files(filter_list, metadata_path)

# Calculating the log of the linear scaling functions.
log_delta_phot_obs_phot_syn = pow(10, (delta_phot_obs_phot_syn / -2.5))

# From this, we calculate the gradient (m) of the scaling function.
#  Note: this gradient is not the m in mx + c, but actually m/c!
gradient = (log_delta_phot_obs_phot_syn - 1) / (
    filter_band_info[0, 0]
    - (filter_band_info[1, 0] * log_delta_phot_obs_phot_syn)
)

# Then we apply this scaling slope to the original observed spectrum.
scaled_spec = spec[:, 1] * ((gradient * spec[:, 0]) + 1)

# Creating a temporary directroy to store any intermiediate files creates
# - keeps the working directory clean.
temp_folder = _mkdir("temp_spectra")
file_string_one = temp_folder + "/" + specname + "_corr_1.dat"
# The code then combines the flux_scaled array with waves and saves for running
# through SMS.
# Note: this spectrum has been linearly scaled, but is still offset from the
# true values; i.e. it now has the right shape/slope, but wrong 'height'.
np.savetxt(file_string_one, np.c_[spec[:, 0], scaled_spec[:]])

#  Over-writing the files to use with SMS - they now contain information
# referring to the intermediate scaled spectrum created above.
write_spec_list_file(filter_list, file_string_one, spec_filename)
write_filter_list_file(filter_list, filters_filename)

# Prompt the user to re-run SMS and then input the new synthetic mags.
print "\nRe-run SMS\n"
# The code now reads in the new/updated synthetic photometry for this
# semi-scaled spectrum (g_syn and r_syn values overwritten).
# Note: gere we can add in more photometry points for comparison if needed. If
# added, they will be used to assist with determining the offset needed to
# complete the calibration.
synthetic_photometry = read_phot_syn(filter_list)

# The code computes the difference between the observed and synthetic
# photometry (delta_g = g_obs - g_syn, delta_r = r_obs - r_syn, others can be
# computed if additional photometry has been supplied).
delta_obs_syn = np.zeros(len(filter_list))
for bb in range(len(delta_obs_syn)):
    delta_obs_syn[bb] = observed_photometry[bb] - synthetic_photometry[bb]

# Code calculates the average of these delta values (delta_av).
av_delta_obs_syn = np.average(delta_obs_syn)

# From this, it computes the constant (c) of the scaling function.
scaling_constant = pow(10, (av_delta_obs_syn / -2.5))

# Using these values for m and c, the code then determines the final corrected
# flux values. It does this by applying them to the ORIGINAL, UNSCALED spectrum
# (NOT the intermediate, semi-scaled spectrum!).
spec_final = (scaling_constant * ((gradient * spec[:, 0]) + 1)) * spec[:, 1]

# Saving the fully flux-calibrated spectrum.
file_string_two = temp_folder + "/" + specname + "_corr_2.dat"
np.savetxt(file_string_two, np.c_[spec[:, 0], spec_final[:]])

# Prompt the user to re-run SMS a final time to obtain errors for the
# flux-calibration.
print "\nRe-run SMS (with full filter list)\n"
# Running the functions to obatin user inputs for observed and synthetic
# photometry again. Here we can read in additional photometry for including in
# our error calculations and plot.
extended_photometry = read_phot_obs(full_filter_list)
observed_photometry = np.append(observed_photometry, extended_photometry[0])
observed_photometry_errors = np.append(
    observed_photometry_errors, extended_photometry[1]
)

full_filter_list = np.append(filter_list, full_filter_list)

#  Over-writing the files to use with SMS - they now contain information
# referring to the final scaled spectrum created above.
write_spec_list_file(full_filter_list, file_string_two, spec_filename)
write_filter_list_file(full_filter_list, filters_filename)

#  Reading in this final set of synthetic photometry.
synthetic_photometry = read_phot_syn(full_filter_list)

# Here the code computes the errors in the flux-calibrated spectrum.
for cc in range(len(full_filter_list)):
    error = observed_photometry[cc] - synthetic_photometry[cc]
    print "\nError in ", full_filter_list[cc], ": ", error

# Print errors to the screen for user vetting.
print "\nAre these errors acceptable?"
verdict = float(input("0=NO, 1=YES:\t"))

# If the errors are acceptable, then the flux-calibrated spectrum is saved in
# the curretn working directory; if not, the spectrum is not saved.
# In this instance, the user may choose to change e.g. the 2 filter bands used
# for thelinear scaling process.
if verdict == 1:
    file_string = specname + "_flux_corrected.dat"
    np.savetxt(file_string, np.c_[spec[:, 0], spec_final[:]])
    # Print the linear scaling equation used.
    print "\nEquation used for scaling (mx + c):"
    print "({:.3E} * x) + {:.3E}".format(
        gradient * scaling_constant, scaling_constant
    )
    # Print the observed and synthetic photometry in table format for copying
    # across to a text file/README for posterity.
    print "Here are the observed and final flux-clibrated photometry points:"
    print "Obs.\tSyn."
    for cc in range(len(full_filter_list)):
        print observed_photometry[cc], "\t", synthetic_photometry[cc]

elif verdict == 0:
    print "Well what is wrong with them?"
else:
    print "Not a valid option!"


# Do some housekeeping here and remove all intermediate files.
shutil.rmtree(temp_folder)
os.remove(spec_filename)
os.remove(filters_filename)

# Plotting up the old and new spectra for visual comparison, along with
# photometry points
lightspeed = 3e18  # Angstroms

filter_band_info = read_filter_files(full_filter_list, metadata_path)

# Here we are converting our photometry points (and errors) into flux space for
# plotting alongside the spectra.
observed_photometry_flux = np.zeros(len(observed_photometry))
observed_photometry_flux_errors = np.zeros((2, len(observed_photometry)))
synthetic_photometry_flux = np.zeros(len(synthetic_photometry))
for dd in range(len(observed_photometry)):
    observed_photometry_flux[dd] = (
        lightspeed / (filter_band_info[dd, 0] ** 2)
    ) * (10 ** -((observed_photometry[dd] + 48.60) / 2.5))
    observed_photometry_flux_errors[0, dd] = (
        lightspeed / (filter_band_info[dd, 0] ** 2)
    ) * (
        10
        ** -(
            (observed_photometry[dd] + observed_photometry_errors[dd] + 48.60)
            / 2.5
        )
    )
    observed_photometry_flux_errors[1, dd] = (
        lightspeed / (filter_band_info[dd, 0] ** 2)
    ) * (
        10
        ** -(
            (observed_photometry[dd] - observed_photometry_errors[dd] + 48.60)
            / 2.5
        )
    )

    synthetic_photometry_flux[dd] = (
        lightspeed / (filter_band_info[dd, 0] ** 2)
    ) * (10 ** -((synthetic_photometry[dd] + 48.60) / 2.5))

# Calculating the observed photometry errors (negative and positive errors)
observed_photometry_flux_errors[0, :] = (
    observed_photometry_flux[:] - observed_photometry_flux_errors[0, :]
)
observed_photometry_flux_errors[1, :] = (
    observed_photometry_flux_errors[1, :] - observed_photometry_flux[:]
)

# Plotting everything up.
plt.plot(
    spec[:, 0], spec[:, 1], color="black", alpha=0.5, label="Original",
)
plt.plot(
    spec[:, 0], spec_final[:], color="blue", alpha=0.5, label="Flux-Calibrated"
)

plt.plot(
    filter_band_info[:, 0],
    observed_photometry_flux[:],
    marker="*",
    linestyle="",
    color="red",
    markersize=15,
    label="Observed Photometry",
)
plt.errorbar(
    filter_band_info[:, 0],
    observed_photometry_flux[:],
    observed_photometry_flux_errors,
    linestyle="",
    ecolor="red",
)

plt.plot(
    filter_band_info[:, 0],
    synthetic_photometry_flux[:],
    marker="x",
    linestyle="",
    color="purple",
    markersize=15,
    label="Synthetic Photometry",
)

#  Plotting the bandpasses of the relevant filters under investigation.
# The relative strengths have been scaled to the observed photometry values.
for aa in range(len(full_filter_list)):
    filter_bandpass = plot_filter_bandpasses(
        full_filter_list[aa], metadata_path
    )
    plt.plot(
        filter_bandpass[:, 0],
        filter_bandpass[:, 1]
        / np.max(filter_bandpass[:, 1])
        * observed_photometry_flux[aa],
        alpha=0.4,
        linewidth=2.0,
        color="black",
    )


plt.legend(loc="best")
plt.xlabel("Wavelength [$\mathrm{\AA}$]")
plt.ylabel("Flux [erg s$^{-1}$ cm$^{-2}$ $\mathrm{\AA}$]")
plt.show()
