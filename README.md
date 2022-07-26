# SeisPol
Dependences: MATLAB Signal Processing Toolbox
SeisPol takes in two to three data streams (fileEast, fileNorth, fileVert) formatted in MATLAB
structures following the format output by rdmseed (CIT). Inputs may be either miniSEED files or .mat
variable files containing only the miniSEED structure variable. It is suggested to provide .mat files for
input to improve processing speed. It is also suggested to properly include metadata within the miniSEED
files provided to the SeisPol or rdmseed functions, as SeisPol will pull important metadata (such as
sampling rate and station information) from the structure variable. Time-constraint inputs are divided into:
● Period: MATLAB duration variable describing the length of the tapered-cosine window
● Window: MATLAB timerange variable describing the total analysis length.
These MATLAB variable formats are utilized for their verbosity and flexibility.
Input data streams are demeaned, detrended, then filtered using a butterworth filter in the range of
2-8 Hz, which is the target frequency window for our purposes. This bandwidth may be changed as
needed. Butterworth filtering was selected as it presents the best trade-off between attenuation within the
passband and rejection outside the passband. The default filter runs at order n = 4 to improve the rejection
qualities. If components are different lengths, components smaller than the largest will be padded with 0
values to match the longest time-series. The function does not currently support the identification of gaps
in data at this time.
Each period is extracted from the time-series using a Tukey tapered-cosine window. This window
shape was selected on account of its relatively even frequency response. A taper of 50\% was selected to
emphasize the center of the analytic period, but allow for adjacent windows to reconstruct the full signal
when combined, as each period overlaps by 50% with the adjacent periods. This has the advantage of
building additional correlation between analysis periods, providing additional support for results.
The central component of SeisPol is the polarity detection algorithm. First, a cross-correlation
matrix is derived from the target period by the equation
where M is the matrix of signal data with columns [East, North, Vertical]. From C, an eigenvalue problem
is solved to derive the principal axes of particle motion. Should the components be correlated, the
eigenvalue solutions will present the axes of an ellipse or ellipsoid that contains the particle motion of the
seismic signal. The eigenvalues are ordered , with corresponding eigenvectors . The
eigenvalues are used to calculate the Rectilinear Polarization Ratio Dlp & Planar Polarization Ratio Dpp:
These values give an overall sense of rectilinear or planar polarization at the point of signal collection.
From the direction cosines 𝜷, 𝜸, and 𝜹 of the eigenvectors, we compute:
Representing the incident seismic angle 𝜙 and azimuthal seismic angle 𝛩. The function calculates one
data point for each Period within the Window.
SeisPol returns a variable containing each data point within the Window and the divisions for each
time segment, allowing for later plotting. SeisPol also returns a plot of the East-component within the
Window, along with each data point from each Period, allowing for quick analysis after processing. This
enables easy tuning of the Period and Window hyperparameters.
