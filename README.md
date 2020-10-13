This procedure file contains the multiple functions necessary to analyze 
fluorescence microscopy images and extract molecular number and brightness information from
time series (temporal brightness) or still snapshots (spatial brightness).

To be run, the procedure requires IgorPro (wavemetrics) 7 or higher version. To run the procedure file, just open it from Igor and compile it.

The functions actually dealing with brightness calculation are
1. temporal_brighntess() which calulates brightness values from each pixel of a time series. The user has the option to select
a polygonal ROI.This function can accomodate series of images acquired both with photon counting detectors (recommended) or
analog detectors (upon providing Koffset, S and sigma0). All data treatment operations consist in the calculation of mean and variances.
Upon user selection (doboxcar=1) the function can call
2. boxcar2(), a function implementing a boxcar filter (as described in Trullo et al. 2013, DOI 10.1002/jemt.22277) to remove the contribution of slow (relative to the boxcar window) 
fluctiations to the variance of the time series, and hence the brightness. 
The function contains also the code to provide a Gaussianity score for the selected ROI. This score is determined as the ratio
of the actual score over the threshold (As defined upon performing the Kolmogorov-Smirnov (KS) goodness-of-fit test for two continuous distributions
(that of actual pixels in the ROI and that of the pixel of an 'ideal' ROI, Gauss distributed with the same mean and variance as the first). The higher the score 
the more Gaussian the distribution of the pixels in the ROI.
3. SpIDA_photoncounting() Calculates one brightness value from a single image (or ROI thereof) based on the relationship between the variance
and the mean of the pixel intensities in the region. The code works for images acquired using photon counting detectors (recommended). 
4. Spida_analog()
This function deals with the analog case using the simplest approach. It requires preliminary determination of detector's offset, S-factor and sigma0. 
Once these parameters are known, each pixel's intensity value is corrected for the offset and divided by the S-factor. Then the variance of a ROI is calculated.
To this variance, the variance due to detector's dark counts (sigma0^2) is subtracted. The spatial brightness arises from dividing this subtracted value by the corrected average intensity within the ROI.
The function contains also the code to provide a Gaussianity score for the selected ROI. This score is determined as the ratio
of the actual score over the threshold (As defined upon performing the Kolmogorov-Smirnov (KS) goodness-of-fit test for two continuous distributions
(that of actual pixels in the ROI and that of the pixel of an 'ideal' ROI, Gauss distributed with the same mean and variance as the first). The higher the score 
the more Gaussian the distribution of the pixels in the ROI.
For Spatial Brightness calculations with analog detectors we refer here also to Godin et al. PNAS 2011, doi/10.1073/pnas.1018658108 and to the related code available under https://neurophotonics.ca/software
	
The concepts contained in the above three functions can be exported to any other language/code, and provide the simple foundation of Molecular Brightness calculation.

The other functions, under Auxiliary Functions simply deal with the polygonal selection of the ROIs and with the graphical display of the results of the analysis, and are specific to IgorPro programming language. 

Paolo Annibale, 2020
