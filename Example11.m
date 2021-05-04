clear P % Parameters struct

% This example shows that refractive indices can be complex, leading to
% losses during propagation. The imaginary part of the refractive index is
% called the extinction coefficient (kappa) and is related to the
% absorption coefficient (alpha) by alpha = 4*pi*kappa/lambda. Both
% n_cladding and the refractive indices of the individual shapes can be
% complex.

% Here, light is launched into two cores. The cladding and the bottom core
% absorb strongly. As a result, the intensity can be seen to decrease
% rapidly in the bottom core, but also gradually in the top core. This is
% because the evanescent field of the top core mode extends into the
% cladding.

% The example also illustrates that the mode solver can be set to solve for
% modes of each core individually using the singleCoreModes flag. This
% makes no difference in the case where one core absorbs differently than
% the other, but if the cores' absorptions are the same then the mode
% solver will otherwise find symmetric and antisymmetric modes that extend
% into both cores.

%% General and solver-related settings
P.name = mfilename;
P.useAllCPUs = false; % If false, BPM-Matlab will leave one processor unused. Useful for doing other work on the PC while simulations are running.
P.useGPU = false; % (Default: false) Use CUDA acceleration for NVIDIA GPUs

%% Visualization parameters
P.updates = 100;            % Number of times to update plot. Must be at least 1, showing the final state.
P.plotEmax = 1; % Max of color scale in the intensity plot, relative to the peak of initial intensity

%% Resolution-related parameters (check for convergence)
P.Lx_main = 50e-6;        % [m] x side length of main area
P.Ly_main = 50e-6;        % [m] y side length of main area
P.Nx_main = 200;          % x resolution of main area
P.Ny_main = 200;          % y resolution of main area
P.padfactor = 1.5;  % How much absorbing padding to add on the sides of the main area (1 means no padding, 2 means the absorbing padding on both sides is of thickness Lx_main/2)
P.dz_target = 1e-6; % [m] z step size to aim for
P.alpha = 3e14;             % [1/m^3] "Absorption coefficient" per squared unit length distance out from edge of main area

%% Problem definition
P.lambda = 1000e-9; % [m] Wavelength

P.n_cladding = 1.45 + 1e-3i; % [] Cladding refractive index (may be complex)
P.n_0 = 1.46; % [] reference refractive index
P.Lz = 2e-3; % [m] z propagation distances for this segment

% In the shapes 2D array, each row is a shape such as a core in a fiber.
% Column 1 are the x coordinates, column 2 are the y coordinates, column 3
% are radii, column 4 are the types of the shapes, column 5 are the peak
% refractive indices (may be complex) and column 6 is the g parameter, only
% needed if any of the shapes are GRIN lenses.

% Shape types are 1: Circular step-index disk, 2: Antialiased circular
% step-index disk, 3: Parabolic graded index disk, 4: GRIN lens focusing in
% both x and y, 5: GRIN lens focusing only in y.
P.shapes = [ 0 -10e-6 5e-6  2  1.46 + 1e-3i;
             0  10e-6 5e-6  2  1.46 + 0i];

nModes = 2; % For mode finding
plotModes = true; % If true, will plot the found modes
sortByLoss = false; % If true, sorts the list of found modes in order of ascending loss. If false, sorts in order of ascending imaginary part of eigenvalue (descending propagation constant)
singleCoreModes = true; % If true, finds modes for each core/shape individually. Note that the resulting "modes" will only be true modes of the entire structure if the core-to-core coupling is negligible.
P = findModes(P,nModes,singleCoreModes,sortByLoss,plotModes);

% P.E can be either a function that takes X, Y and Eparameters as inputs
% and provides the complex E field as output, or it can be a struct with 3
% fields: a 'field' field which is the complex E-field matrix, and 'Lx' and
% 'Ly' fields that describe the side lengths of the provided E matrix. In
% the case of a struct, the provided E field will be adapted to the new
% grid using the interp2 function.
% P.E = @calcInitialE; % Defined at the end of this file
P.E = modeSuperposition(P,[1 2]);

% Run solver
P = FD_BPM(P);

%% USER DEFINED E-FIELD INITIALIZATION FUNCTION
function E = calcInitialE(X,Y,Eparameters) % Function to determine the initial E field. Eparameters is a cell array of additional parameters such as beam size
w_0 = 2.5e-6;
offset = 2.5e-6;
amplitude = exp(-((X-offset).^2+Y.^2)/w_0^2);
phase = zeros(size(X));
E = amplitude.*exp(1i*phase); % Electric field
end