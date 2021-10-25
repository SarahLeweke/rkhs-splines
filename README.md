# rkhs-splines
Vector-valued and scalar-valued approximation splines in reproducing kernel Hilbert spaces for the reconstruction of neuronal current from MEG and EEG measurements.

## Summary
`spline_MEG` is a function reconstructing the neuronal current inside the brain from given magnetoencephalography (MEG) or electroencephalography (EEG) measurements.
The underlying model of the problem is given by:

 * The multiple-shell model with a ball modelling the human brain and three shells modelling the brain fluid, the scull, and the skalp, [[10]](#10).
 * The relation between the sought neuronal current and the measurements is derived from quasi-static Maxwell's equation, [[11]](#11).
 * The measurement positions are chosen accordingly to measurements from a clinical trial, [[5]](#5).

`spline_MEG` can solve three different inverse problems, that is, the reconstruction of the original vector-valued neuronal current from either given MEG or EEG measurements as well as the 
reconstruction of a scalar-valued component of the neuronal current from given MEG measurements. 
In all cases, the corresponding functional inverse problem is formulated as a Fredholm integral equation of the first kind and stated in [[4]](#4),[[6]](#6),[[7]](#7) and additionally in [[3]](#3) for the scalar-valued problem.
Besides, the vector-valued neuronal current can be calculated from the reconstructed scalar-valued component.

In order to solve the inverse problems, splines based on reproducing kernel Hilbert spaces are constructed, see [[1]](#1)-[[3]](#3),[[9]](#9) (scalar case) and [[4]](#4),[[7]](#7) (vector case).
The splines can be used for the inversion of real data as well as for the inversion of generated synthetic data. 
An inversion without regularization, the minimization of the regularized Tikhonov-functional via splines [[3]](#3),[[4]](#4), or a regularized penalized basis pursuit can be chosen [[12]](#12).  
Five parameter choice methods are implemented for finding the best regularization parameter, namely an automatic as well as a manual L-curve method, the discrepancy principle, the quasi-optimality criterion, and generalized cross validation. 
If an exact solution is known, the method resulting in the smallest normalized root mean square error on the plotting grid is chosen. 

For the visualization of the results, MATLAB plots can be produced and data can be exported for TikZ plots. 
For the sake of speed, all plots can be disabled. 

`spline_MEG` can generate its own synthetic data based on the forward operator applied to orthonormal basis functions or splines. 
In order to avoid the inverse crime, different parameters for constructing the data as well as for the inversion can be specified. 
Besides, there are two methods for calculating the spline matrix available.
On the one hand, the faster method is based on the truncated singular value decomposition.
On the other hand, a spline-free method based on quasi-Monte-Carlo integration and numerical differentiation is implemented.
In all possible settings, the data can be corrupted with additive white Gaussian noise.

## Technical Requirements

 * The code is produced and tested only with MATLAB R2021a, see [[8]](#8).
 * For the visualization of the L-curve in MATLAB, the add-on **labelpoints** version 4.1.1. is used.
   Adam Danz (2021). [labelpoints](https://www.mathworks.com/matlabcentral/fileexchange/46891-labelpoints), MATLAB Central File Exchange. Retrieved October 20, 2021.
 * Parallel computing requires the Parallel Computing Toolbox and is, hence, disabled by default. 
    Only generating the synthetic test data via the quasi-Monte-Carlo integration can be performed in parallel.

## Installation 

 * Install [labelpoints](https://www.mathworks.com/matlabcentral/fileexchange/46891-labelpoints) toolbox via MATLAB → Environment → Add-ons → Get Add-ons. Alternatively, the usage of this toolbox can be commented out in `solve_inverseSplineProblem.m` line 46.
 * Install Parallel Computing Toolbox, if the usage is desired. 
 * Save the folders `code`, `solution`, `data`, and `pics` in one path.
 * Open MATLAB and go into the folder `code`.
 * Start `spline_MEG`.

## Workflow

The main function `spline_MEEG` runs without additional input arguments, since the default values are determined in the function argument validation.
There, all possible values for the input arguments and their structures are listed. Explanations for the particular input arguments are listed in the next table.

| argument | description |
| ---|---|
| `data_generatingCase`  | type of exact solution for generating synthetic data (splines, onb, vector splines) or loading real data|
| `plot_program` | visualization of plots via MATLAB, LaTeX, or none|
| `h` | nodal width parameter for approximation spline|
| `spline_sequenceCase` | sequence used for building the approximation spline, particular sequences are defined mathematically in line 49-67 |
| `spline_truncationIndex`  | truncation degree of approximation spline series|
| `verify_closedSeriesRepresentations`  | if true, particular closed series representation are checked against summation of truncated series|
| `verify_basisMathFunctions` | if true, math functions such as Clenshaw algorithms or implementation of spherical harmonics / Legendre polynomials are verified|
| `reg_method` | regularization method for inversion such as Tikhonov minimization, regularized penalized basis pursuit, or none  |
| `rel_regPar`  | vector of relative regularization parameters, which will be multiplied by the absolute value of the spline matrix in line 131 |
| `prob_case`  | reconstructed inverse problem |
| `recalc` | if true, recalculation of the used splines instead of loading existing ones is forced|
| `noiseLevel`  | level for additive white Gaussian noise in percent points|
| `visualize_splineMatrix` | if true, a visualization of the used spline matrix is plotted|


The model parameters, such as the radii of the shells and their conductivity as well as the positions of the measurement points, are collected in `solution/deviceFile.mat` separated by the MEG and EEG application. 
Therein, the particular quantities as well as  their units are stored.
In the case of the EEG problem, the used coefficients β<sub>n</sub><sup>(L)</sup> can be found in `device.EEG`.
Note that these coefficients depend on the conductivites and the radii of the model. If these parameters change, the coefficients  β<sub>n</sub><sup>(L)</sup>  must be recalculated according to [Sec. 4 [7]](#7).

## Parallel Computing

For running the quasi-Monte-Carlo integration in parallel, choose `data_generatingCase = 2splines_avoidIC` and set the number of integration points in  `num_integration_point` in line 47 in `generate_splineData.m` according to the desired accuracy, then uncomment line 217-218.
Additionally, set `parallel_computing = true` in the function argument validation of `integragte_quasiMonteCarlo_ball_scalarMEG.m`. 
For acceptable accuracy within our synthetic test cases, we required 4.5e6 integration points. 
The precalculated spline matrices can be found in `/solution` depending on the nodal width *h*.


## Reproducibility

### Approximation Splines
The used sequences for generating the splines are coded in line 50-69 of `spline_MEEG`. Even though we have tested several nodal widths *h* within our numerical tests, the best results were achieved with *h*=0.85.
This is also the default value of the nodal width in the code. 
The used regularization parameter for the MEG synthetic spline test cases are stated in Table 1 and for the EEG synthetic spline test case in Table 2 of [[4]](#4).
The numerical results of [[4]](#4) are entirely produced with this code.

### Synthetic Data
All used parameters are stated in the paper [[4]](#4) as well as coded in `generate_splineData.m` for reproducing the non-noisy data. 
Since no fixed seed was used for noising the data, the used data sets are stored in the folder `data` depending on the test case and the noise level.
If a particular data set should be loaded, uncomment the lines 123-124 in `spline_MEEG`. Insert the precise path in line 123. 
For the calculation of the corresponding exact solution, it is necessary to state `data_generatingCase` and `prob_case` according to the used data set.

### Real Data
The used real data sets can be found in the folder `data`. Note that `data_real_scalarMEG` and `data_real_vectorMEG` consists of exactly the same data. 
There are two data files required due to the structure of the implementation. 

## References

<a id="1">[1]</a> 
Amirbekyan, A. (2007). 
The Application of Reproducing Kernel Based Spline Approximation to Seismic Surface
and Body Wave Tomography: Theoretical Aspects and Numerical Results. PhD thesis. 
University of Kaiserlsautern, Department of Mathematics, Geomathematics Group.

<a id="2">[2]</a> 
Amirbekyan, A., Michel, V. (2008). 
Splines on the three-dimensional ball and their application to seismic
body wave tomography. Inverse Probl. 24, 015022.

<a id="3">[3]</a> 
Fokas, A. S., Hauk, O., V. Michel (2012). 
Electro-magneto-encephalography for the three-shell model:
numerical implementation via splines for distributed current in spherical geometry. Inverse Probl.
28, 035009.

<a id="4">[4]</a> 
Hauk, O. Leweke, S., Michel, V. (upcomming). Vector-valued Spline Method for the Spherical Multiple-shell
Electro-magnetoencephalography Problem.

<a id="5">[5]</a> 
Hauk, O. (2013) MEG and EEG measurement data. Personal Communication. MRC Cognition and Brain
Sciences Unit, Cambridge, UK.

<a id="6">[6]</a> 
Leweke, S. (2018) The Inverse Magneto-electroencephalography Problem for the Spherical Multiple-shell Model–
Theoretical Investigations and Numerical Aspects. PhD thesis. University of Siegen, Department of
Mathematics, Geomathematics Group. 

<a id="7">[7]</a> 
Leweke, S. Michel, V. and Fokas, A. S. (2020).
Electro-magnetoencephalography for a spherical multiple-
shell model: novel integral operators with singular-value decompositions. Inverse Probl. 36, 035003.

<a id="8">[8]</a> 
MATLAB (2021). version 9.10.0 (R2021a). Natick, Massachusetts: The MathWorks Inc.

<a id="9">[9]</a> 
Michel, V. (2013).
Lectures on Constructive Approximation. Fourier, Spline, and Wavelet Methods on the Real
Line, the Sphere, and the Ball. New York: Birkhäuser.

<a id="10">[10]</a> 
de Munck, J. C. (1988). 
The potential distribution in a layered anisotropic spheroidal volume conductor. 
J. Appl. Phys. 64, 464–470.

<a id="11">[11]</a> 
Plonsey, R. (1969). Biomagnetic Phenomena. New York: McGraw-Hill.

<a id="12">[12]</a> 
Simeoni, M. M. J.-A. (2021). Functional penalised basis pursuit on spheres. In: Applied and Computational
Harmonic Analysis 53, 1–53.
