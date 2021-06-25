clear P % Parameters struct

% This example shows how to define a refractive index profile using the P.n
% input instead of the P.shapes input. P.n is a handle to a user-defined
% function at the end of the model file (this file), in which the user can
% specify the refractive index profile using a series of analytical
% expressions. In this example, we define the refractive index profile to
% be an ellipsoidal core with part of one side cut off.

%% General and solver-related settings
P.name = mfilename;
P.useAllCPUs = false; % If false, BPM-Matlab will leave one processor unused. Useful for doing other work on the PC while simulations are running.
P.useGPU = false; % (Default: false) Use CUDA acceleration for NVIDIA GPUs

%% Visualization parameters
P.updates = 100;            % Number of times to update plot. Must be at least 1, showing the final state.
P.plotEmax = 1; % Max of color scale in the intensity plot, relative to the peak of initial intensity

%% Resolution-related parameters (check for convergence)
P.Lx_main = 20e-6;        % [m] x side length of main area
P.Ly_main = 20e-6;        % [m] y side length of main area
P.Nx_main = 200;          % x resolution of main area
P.Ny_main = 200;          % y resolution of main area
P.padfactor = 2;  % How much absorbing padding to add on the sides of the main area (1 means no padding, 2 means the absorbing padding on both sides is of thickness Lx_main/2)
P.dz_target = 1e-6; % [m] z step size to aim for
P.alpha = 2e14;             % [1/m^3] "Absorption coefficient" per squared unit length distance out from edge of main area

%% Problem definition
P.lambda = 1000e-9; % [m] Wavelength
P.n_background = 1.45; % [] (may be complex) Background refractive index, (in this case, the cladding)
P.n_0 = 1.47; % [] reference refractive index
P.Lz = 2e-3; % [m] z propagation distances for this segment

P.n.func = @calcRI; % Defined at the end of this file. See the readme file for details

P.E = @calcInitialE; % Defined at the end of this file. See the readme file for details

% Run solver
P = FD_BPM(P);

%% USER DEFINED E-FIELD INITIALIZATION FUNCTION
function E = calcInitialE(X,Y,Eparameters) % Function to determine the initial E field. Eparameters is a cell array of additional parameters such as beam size
w_0 = 2.5e-6;
offset = -2.5e-6;
amplitude = exp(-((X-offset).^2+Y.^2)/w_0^2);
phase = zeros(size(X));
E = amplitude.*exp(1i*phase); % Electric field
end

%% USER DEFINED RI INITIALIZATION FUNCTION
function n = calcRI(X,Y,n_background,nParameters) % This user-customizable function is the one that will be used to define the refractive index profile. You can either use analytical expressions or, e.g., load data from a file and use interpn to adapt the data to the simulation grid.
% n may be complex
n = n_background*ones(size(X)); % Start by setting all pixels to n_background
n(X.^2/2 + Y.^2 < (2e-6)^2 & X > -2e-6) = n_background + 0.02; % Equation for ellipse with part of the side cut off
end