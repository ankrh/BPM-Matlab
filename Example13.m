P = BPMmatlab.model;

% This example, like example 12, shows how to define a completely arbitrary
% (complex) refractive index profile using the P.n input, but in this
% example we specify the properties of P.n manually. In this example, we
% load a built-in MATLAB image of the moon, scale it to values that are
% realistic refractive indices and store it as the 'n' property of the P.n
% input. We choose some arbitrary values for the 'Lx' and 'Ly' side length
% properties that correspond to this 2D array. BPM-matlab will automatically
% convert the data to the simulation grid by linear interpolation.

% The mode solver also works with such arbitrary refractive index profiles,
% which we demonstrate here by finding and launching the 5th mode into this
% simulation.

% In the first segment, the "fiber" is non-twisted, so nothing happens to
% the intensity profile as the beam propagates (because the injected field
% is a mode), but the fiber is twisted in the second segment, which makes
% the intensity profile change slightly.

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
P.padfactor = 1.5;  % How much absorbing padding to add on the sides of the main area (1 means no padding, 2 means the absorbing padding on both sides is of thickness Lx_main/2)
P.dz_target = 1e-6; % [m] z step size to aim for
P.alpha = 3e14;             % [1/m^3] "Absorption coefficient" per squared unit length distance out from edge of main area

%% Problem definition
P.lambda = 1000e-9; % [m] Wavelength
P.n_background = 1.45; % [] (may be complex) Background refractive index, (in this case, the cladding)
P.n_0 = 1.46; % [] reference refractive index
P.Lz = 2e-3; % [m] z propagation distances for this segment

P.n.n = flipud(single(imread('./ExampleData/moon.tif'))).'/255*0.02 + P.n_background;
P.n.Lx = 10e-6;
P.n.Ly = 10e-6*(537/358);

P = findModes(P,5);
P.E = P.modes(5);

%% First segment, non-twisted fiber
% Run solver
P = FD_BPM(P);

%% Second segment, twisted fiber
P.twistRate = 1000;
% Run solver
P = FD_BPM(P);
