clear P % Parameters struct

% This example shows how to define a completely arbitrary refractive index
% profile using the P.n input instead of the P.shapes input. P.n is a
% handle to a user-defined function at the end of the model file, in which
% the user can either specify the refractive index profile using a series
% of analytical expressions or by loading data from a file and adapting it
% to the simulation grid using interpn. In this example, we use a built-in
% MATLAB image of the moon as our refractive index profile. The mode solver
% also works with these arbitrary refractive index profiles, demonstrated
% here by launching the 5th mode into this simulation.

% In the simulation, we can see that only the phase of the field changes
% during propagation but not the intensity distribution, proving that the
% input field is indeed a mode of the moon-shaped refractive index profile.

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
P.n_background = 1.45; % [] Background refractive index
P.n_0 = 1.46; % [] reference refractive index
P.Lz = 2e-3; % [m] z propagation distances for this segment

% To define the refractive index profile, you can use either P.shapes or
% P.n. Here we use P.n, which can be either a function that takes X, Y,
% n_background and nParameters as inputs and provides the complex (or real)
% refractive index profile as output, or it can be a struct with 3 fields:
% an 'n' field which is the complex refractive index matrix, and 'Lx' and
% 'Ly' fields that describe the side lengths of the provided n matrix. In
% the case of a struct, the provided n will be adapted to the new grid
% using the interpn function.
P.n = @calcInitialn; % Defined at the end of the file

% P.E can be either a function that takes X, Y and Eparameters as inputs
% and provides the complex E field as output, or it can be a struct with 3
% fields: a 'field' field which is the complex E-field matrix, and 'Lx' and
% 'Ly' fields that describe the side lengths of the provided E matrix. In
% the case of a struct, the provided E field will be adapted to the new
% grid using the interpn function.
% P.E = @calcInitialE; % Defined at the end of this file
nModes = 5; % For mode finding
plotModes = true; % If true, will plot the found modes
sortByLoss = false; % If true, sorts the list of found modes in order of ascending loss. If false, sorts in order of ascending imaginary part of eigenvalue (descending propagation constant)
singleCoreModes = false; % If true and if P.shapes is defined, finds modes for each core/shape individually. Note that the resulting "modes" will only be true modes of the entire structure if the core-to-core coupling is negligible.
P = findModes(P,nModes,singleCoreModes,sortByLoss,plotModes);
P.E = P.modes(5);

% Run solver
P = FD_BPM(P);

%% USER DEFINED E-FIELD INITIALIZATION FUNCTION
% function E = calcInitialE(X,Y,Eparameters) % Function to determine the initial E field. Eparameters is a cell array of additional parameters such as beam size
% w_0 = 2.5e-6;
% offset = -2.5e-6;
% amplitude = exp(-((X-offset).^2+Y.^2)/w_0^2);
% phase = zeros(size(X));
% E = amplitude.*exp(1i*phase); % Electric field
% end

%% USER DEFINED RI INITIALIZATION FUNCTION
function n = calcInitialn(X,Y,n_background,nParameters) % This user-customizable function is the one that will be used to define the refractive index profile. You can either use analytical expressions or, e.g., load data from a file and use interpn to adapt the data to the simulation grid.
A = flipud(imread('moon.tif')).'; % A built-in MATLAB image of the moon
Lx = 10e-6; % Arbitrarily chosen x width of the image
[Nx,Ny] = size(A);
dx = Lx/Nx;
dy = dx;
x = dx*(-(Nx-1)/2:(Nx-1)/2);
y = dy*(-(Ny-1)/2:(Ny-1)/2);
n = n_background + 0.02*interpn(x,y,single(A),X,Y,'linear',n_background)/255;
end