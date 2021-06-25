clear P % Parameters struct

% This example starts in the same way as example 1. To demonstrate that in
% media with uniform refractive index n_0, either the FFTBPM or FDBPM
% solver can be used for propagation. The output of example 1 is propagated
% 2e-4 m through a medium of uniform refractive index 1.45 with both
% solvers. You can see that the results are identical by comparing figures
% 2 and 3.

% The example also demonstrates the use of non-standard color maps for the
% plots.

%% Part 1 run with FDBPM
%% General and solver-related settings
P.name = mfilename;
P.useAllCPUs = true; % If false, BPM-Matlab will leave one processor unused. Useful for doing other work on the PC while simulations are running.
P.useGPU = false; % (Default: false) Use CUDA acceleration for NVIDIA GPUs

%% Visualization parameters
P.figNum = 1;
P.updates = 30;            % Number of times to update plot. Must be at least 1, showing the final state.
P.downsampleImages = false; % Due to a weird MATLAB bug, MATLAB may crash when having created imagesc (or image) plots with dimensions larger than roughly 2500x2500 and then calling mex functions repeatedly. This flag will enable downsampling to 500x500 of all data before plotting, hopefully avoiding the issue.
P.displayScaling = 1;  % Zooms in on figures. Set to 1 for no zooming.  

%The colormap options for the different subplots are 
%1: GBP, 2: HSV, 3:parula, 4: gray, 5: cividis
P.Intensity_colormap = 5; 
P.Phase_colormap = 2; 
P.n_colormap = 4; 

%% Resolution-related parameters (check for convergence)
P.Lx_main = 20e-6;        % [m] x side length of main area
P.Ly_main = 20e-6;        % [m] y side length of main area
P.Nx_main = 200;          % x resolution of main area
P.Ny_main = 200;          % y resolution of main area
P.padfactor = 1.5;  % How much absorbing padding to add on the sides of the main area (1 means no padding, 2 means the absorbing padding on both sides is of thickness Lx_main/2)
P.dz_target = 1e-6; % [m] z step size to aim for
P.alpha = 3e14;             % [1/m^3] "Absorption coefficient" per squared unit length distance out from edge of main area

%% Problem definition
P.lambda = 1000e-9; % [m] Wavelength
P.n_background = 1.45; % [] (may be complex) Background refractive index, (in this case, the cladding)
P.n_0 = 1.46; % [] reference refractive index
P.Lz = 2e-3; % [m] z propagation distances for this segment

P.shapes = [ 0 0 5e-6  2  1.46]; % See the readme file for details

P.E = @calcInitialE; % Defined at the end of this file. See the readme file for details

% Run solver
P = FD_BPM(P);

E = P.E; % Store this for use in the FFTBPM

%% Part 2 run with FDBPM
P.figNum = 2;
P.figTitle = 'FD BPM';

P.n_0 = 1.45;
P.shapes = []; % Remove the shape

P.Lx_main = 100e-6;        % [m] x side length of main area
P.Ly_main = 100e-6;        % [m] y side length of main area
P.Nx_main = 500;          % x resolution of main area
P.Ny_main = 500;          % y resolution of main area
P.alpha = 8e13;             % [1/m^3] "Absorption coefficient" per squared unit length distance out from edge of main area

P.Lz = 2e-4; % [m] z propagation distances for this segment

% Run solver
P = FD_BPM(P);

%% Part 2 run with FFTBPM for comparison
P.E = E; % Use E field output from part 1
P.figNum = 3;
P.figTitle = 'FFT BPM';

% Run solver
[E_out] = FFT_BPM(P);

%% USER DEFINED E-FIELD INITIALIZATION FUNCTION
function E = calcInitialE(X,Y,Eparameters) % Function to determine the initial E field. Eparameters is a cell array of additional parameters such as beam size
w_0 = 2.5e-6;
offset = 2.5e-6;
amplitude = exp(-((X-offset).^2+Y.^2)/w_0^2);
phase = zeros(size(X));
E = amplitude.*exp(1i*phase); % Electric field
end