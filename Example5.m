clear P % Parameters struct

% This example shows a basic case of launching a Gaussian beam
% into the cylindrical GRINTECH lens (one dimensional -y- version of GT-CFRL-100-025-20-CC (810) ) 
% realized using FDBPM and the subsequent propagation of the E field from GRIN lens's output in
% air using FFTBPM

%% Part 1 run with FDBPM
%% General and solver-related settings
P.name = mfilename;
P.useAllCPUs = true;
P.useGPU = false;

%% Visualization parameters
P.figNum = 1;
P.figTitle = 'In GRIN Lens';
P.saveVideo = false; % To save the field intensity and phase profiles at different transverse planes
P.updates = 100;            % Number of times to update plot. Must be at least 1, showing the final state.
P.downsampleImages = false; % Due to a weird MATLAB bug, MATLAB may crash when having created imagesc (or image) plots with dimensions larger than roughly 2500x2500 and then calling mex functions repeatedly. This flag will enable downsampling to 500x500 of all data before plotting, hopefully avoiding the issue.
P.displayScaling = 1;  % Zooms in on figures 1 & 3a,b. Set to 1 for no zooming.  

%% Resolution-related parameters (check for convergence)
P.Lx_main = 1.0e-3;        % [m] x side length of main area
P.Ly_main = 0.7e-3;        % [m] y side length of main area
P.Nx_main = 1400;          % x resolution of main area
P.Ny_main = 1000;          % y resolution of main area
P.padfactor = 1;  % How much absorbing padding to add on the sides of the main area (1 means no padding, 2 means the absorbing padding on both sides is of thickness Lx_main/2)
P.dz_target = 20e-6; % [m] z step size to aim for
P.alpha = 3e14;             % [1/m^3] "Absorption coefficient" per squared unit length distance out from edge of main area

%% Problem definition
P.lambda = 810e-9; % [m] Wavelength
P.n_cladding = 1.0; % [] Cladding refractive index
P.n_0 = 1.521;
P.Lz = 6.05e-3; % [m] z propagation distances for this segment

% In the shapes 2D array, each row is a shape such as a core in a fiber.
% Column 1 are the x coordinates, column 2 are the y coordinates, column 3
% are radii, column 4 are the types of the shapes, column 5 are the peak
% refractive indices and column 6 is the g parameter, only needed if any of
% the shapes are GRIN lenses.

% Shape types are 1: Circular step-index disk, 2: Antialiased circular
% step-index disk, 3: Parabolic graded index disk, 4: GRIN lens focusing in
% both x and y, 5: GRIN lens focusing only in y.
P.shapes = [ 0 0 0.5e-3  5  1.521 260];

% P.E can be either a function that takes X, Y and Eparameters as inputs
% and provides the complex E field as output, or it can be a struct with 3
% fields: a 'field' field which is the complex E-field matrix, and 'Lx' and
% 'Ly' fields that describe the side lengths of the provided E matrix. In
% the case of a struct, the provided E field will be adapted to the new
% grid using the interp2 function.
P.E = @calcInitialE; % Defined at the end of this file

% Run solver
P = FD_BPM(P);

%% Output E from FDBPM propagating in air with FFTBPM
P.figNum = 2;
P.figTitle = 'In air';
P.saveVideo = false; % To save the field intensity and phase profiles at different transverse planes
P.updates = 30;            % Number of times to update plot. Must be at least 1, showing the final state.

P.Lx_main = 2.5e-3;        % [m] x side length of main area
P.Ly_main = 0.7e-3;        % [m] y side length of main area
P.Nx_main = 1500;          % x resolution of main area
P.Ny_main = 400;          % y resolution of main area
P.Lz = 1e-2; % [m] z propagation distances for this segment

P.E = P.Efinal;

% Run solver
FFT_BPM(P);

%% USER DEFINED E-FIELD INITIALIZATION FUNCTION
function E = calcInitialE(X,Y,Eparameters) % Function to determine the initial E field. Eparameters is a cell array of additional parameters such as beam size
w_0 = 5e-6;
offset = 0;
amplitude = exp(-((X-offset).^2+Y.^2)/w_0^2);
phase = zeros(size(X));
E = amplitude.*exp(1i*phase); % Electric field
end