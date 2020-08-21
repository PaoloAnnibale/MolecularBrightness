{\rtf1\ansi\ansicpg1252\cocoartf1561\cocoasubrtf610
{\fonttbl\f0\fmodern\fcharset0 Courier;}
{\colortbl;\red255\green255\blue255;\red0\green0\blue0;}
{\*\expandedcolortbl;;\cssrgb\c0\c0\c0;}
\paperw11900\paperh16840\margl1440\margr1440\vieww10800\viewh8400\viewkind0
\deftab720
\pard\pardeftab720\sl280\partightenfactor0

\f0\fs24 \cf2 \expnd0\expndtw0\kerning0
\outl0\strokewidth0 \strokec2 This procedure file contains the multiple functions necessary to analyze \
fluorescence microscopy images and extract molecular number and brightness information from\
time series (temporal brightness) or still snapshots (spatial brightness)\
\
The functions actually dealing with brightness calculation are\
1. temporal_brighntess() which calulates brightness values from each pixel of a time series. The user has the option to select\
a polygonal ROI.This function can accomodate series of images acquired both with photon counting detectors (recommended) or\
analog detectors (upon providing Koffset, S and sigma0). All data treatment operations consist in the calculation of mean and variances.\
Upon user selection (doboxcar=1) the function can call\
2. boxcar2(), a function implementing a boxcar filter (as described in Trullo et al. 2013, DOI 10.1002/jemt.22277) to remove the contribution of slow (relative to the boxcar window) \
fluctiations to the variance of the time series, and hence the brightness. \
The function contains also the code to provide a Gaussianity score for the selected ROI. This score is determined as the ratio\
of the actual score over the threshold (As defined upon performing the Kolmogorov-Smirnov (KS) goodness-of-fit test for two continuous distributions\
(that of actual pixels in the ROI and that of the pixel of an 'ideal' ROI, Gauss distributed with the same mean and variance as the first). The higher the score \
the more Gaussian the distribution of the pixels in the ROI.\
3. SpIDA_photoncounting() Calculates one brightness value from a single image (or ROI thereof) based on the relationship between the variance\
and the mean of the pixel intensities in the region. The code works for images acquired using photon counting detectors (recommended). For the analog\
case we refer here to Godin et al. PNAS 2011, doi/10.1073/pnas.1018658108 and to the related code available under https://neurophotonics.ca/software\
The function contains also the code to provide a Gaussianity score for the selected ROI. This score is determined as the ratio\
of the actual score over the threshold (As defined upon performing the Kolmogorov-Smirnov (KS) goodness-of-fit test for two continuous distributions\
(that of actual pixels in the ROI and that of the pixel of an 'ideal' ROI, Gauss distributed with the same mean and variance as the first). The higher the score \
the more Gaussian the distribution of the pixels in the ROI.\
	\
The concepts contained in the above three functions can be exported to any other language/code, and provide the simple foundation of Molecular Brightness calculation.\
\
The other functions, under Auxiliary Functions simply deal with the polygonal selection of the ROIs and with the graphical display of the results of the analysis, and are specific to IgorPro programming language. \
\
Paolo Annibale, 2020}