# SMS-LinearScalingFunction

James Gillanders (07/2021): Here I present a simple piece of code that can be used alongside the the s3/SMS code (https://github.com/cinserra/S3) to linearly scale an observed spectrum to match observed photometry.

ASS.py
------

This program can be used to linearly scale spectra. Below I detail the theory behind how the code works, in a general sense. For more information on how specific things in the code work, see the comments contined within the ASS.py file.

Code theory - why does it work?

We have an observed spectrum, with some flux (F_obs).

Note: for simplicity in the discussion below, I will present the calibration as if we are using the g and r band filters, but in reality, any two user-specified filters can be used for the calibration.

We want to flux-calibrate these values to match observed photometry. To do this, we will apply a linear scaling function (m*x + c) to F_obs:

    F_corr = ((m * lam) + b) * F_obs   [1]

How do we convert from flux space to magnitude space?
We use the following equation:

    Flux = (c / lam_eff**2) * 10**(-(mag + 48.6) / 2.5)   [2]
   
given the spectral flux densities are expressed per unit wavelength, and where c is the speed of light, and lam_eff is the effective peak wavelength of the filter band.

To go the other way, it is simply:

    mag = -2.5 * log10(Flux * (lam_eff**2 / c)) - 48.6   [3]

So our equation [1] above can be rewritten for mag space instead of flux space:

    mag = -2.5 * log10((m * lam_eff) + b) + mag_obs   [4]

(there is a lot of term cancellation here...)

Rearranging to separate out our terms into two different 'measurables':

    mag_phot = -2.5 * log10(b[(m * lam_eff / b) +  1]) + mag_spec   [5]
    
where mag_phot is the photometric magnitude, and mag_spec is the synthetic photometric point obtained from the spectrum, using SMS.

    mag_phot = -2.5 * log10(b) -2.5 * log10([m * lam_eff / b] + 1) + mag_spec   [6]

    mag_phot = -2.5 * log10(b) -2.5 * log10([d * lam_eff] + 1) + mag_spec   [7]

where d = m / b

Now we calculate delta_phot (g_phot - r_phot) and delta_spec (g_spec - r_spec) and use these to deduce the current gradient (delta_spec) of the spectrum, and
also what the gradient should be (delta_phot):

    g_phot = -2.5 * log10(b) -2.5 * log10([d * lam_eff_g] + 1) + g_spec   [8]
    
    r_phot = -2.5 * log10(b) -2.5 * log10([d * lam_eff_r] + 1) + r_spec   [9]

Combining these, we get:

    delta_phot = delta_spec -2.5 * log10{([d * lam_eff_g] + 1) / ([d * lam_eff_r] + 1)}   [10]

Rearranging and trying to isolate d:

    ([d * lam_eff_g] + 1) / ([d * lam_eff_r] + 1) = 10**([delta_phot - delta_spec] / -2.5)   [11]

    d = [10**([delta_phot - delta_spec] / -2.5) - 1] / [lam_eff_g - (lam_eff_r * [10**([delta_phot - delta_spec] / -2.5)])]   [12]

(there is a lot of term cancellation here...)

This d value is now used to scale the spectrum. It contains the necessary information to correct the gradient of the spectrum. Note: it is not the 'true' value of the gradient of our correction; recall that d = m / b

But we do not have a value for our scaling constant, b. How can we scale the spectrum without our offset?

Correct - we do not have a value for b, but what we can do is apply the gradient correction to the spectrum, re-run through SMS to get updated synthetic photometry points, and then use these (and the observed photometry points) to determine the offset from the true values. Then, we can take a simple average of these offsets to determine our best value for b.

So, what does this look like?

Well, using our equation [1] for correcting flux:

    F_corr = ((m * lam) + b) * F_obs

Or:

    F_corr = b[(m * lam / b) + 1] * F_obs
    
    F_corr = b[(d * lam) + 1] * F_obs

So what we do is assume that b = 1. Then we have:

    F_corr = [(d * lam) + 1] * F_obs   [13]

So we apply this to the unscaled spectrum. Here, F_obs is the array/list of unscaled flux values, lam is the array/list of wavelength values of the spectrum, and F_corr is the array/list of the corrected flux values.

This scaled spectrum is then output to a file, for the user to re-run through SMS to obtain updated values for g_spec and r_spec.

With these updated values, determining the value for b is a trivial task. We compute delta_g (g_phot - g_spec) and delta_r (r_phot - r_spec), and, taking a
simple average, we can calculate b. Looking at this equation again:

    mag_phot = -2.5 * log10(b) -2.5 * log10([d * lam_eff] + 1) + mag_spec,   [14]
    
where -2.5 * log10([d * lam_eff] + 1) + mag_spec is equivalent to the new photometry values read in (i.e. the new values of g_spec and r_spec).

This means our equation then becomes:

    mag_phot = -2.5 * log10(b) + mag_spec,   [15]
    
where we stress that this mag_spec is the NEW/updated synthetic photometric value.

From this:

    g_phot - g_spec = delta_g = -2.5 * log10(b)   [16]

And:

    r_phot - r_spec = delta_r = -2.5 * log10(b)   [17]

These will give slightly different values for b in most cases. So, instead of computing different values for b, we average our offsets, and use that:

    delta_av = (delta_g + delta_r) / 2   [18]

So:

    delta_av = -2.5 * log10(b)   [19]

Finally:

    b = 10**(delta_av / -2.5)   [20]

Now we have values for d and b, and we can produce our final flux-calibrated spectrum (F_corr = b[(d * lam) +  1] * F_obs). We then save this, and run it through SMS again to determine the final errors in our flux calibration.

The code then prints the errors in the calculation for the user to determine if the calibration was successful or not.

Finally, the code plots the unscaled and scaled version of the spectra, and also over-plots the true photometric photometry, and the final synthetic photometry.

Additionally, the m*x + c scaling equation is output, as well as a neat table of observed versus synthetic photometry values for easy exporting to a file.

