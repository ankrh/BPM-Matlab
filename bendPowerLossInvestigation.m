clear P % Parameters struct

% Comparison with R. T. Schermer and J. H. Cole, "Improved Bend Loss Formula Verified for Optical Fiber by Simulation and Experiment," in IEEE Journal of Quantum Electronics, vol. 43, no. 10, pp. 899-909, Oct. 2007, doi: 10.1109/JQE.2007.903364.

%% General and solver-related settings
P.name = mfilename;

%% Visualization parameters
P.updates = 5;            % Number of times to update plot. Must be at least 1, showing the final state.

%% Resolution-related parameters (check for convergence)
P.Lx_main = 50e-6;        % [m] x side length of main area
P.Ly_main = 50e-6;        % [m] y side length of main area
P.Nx_main = 150;          % x resolution of main area
P.Ny_main = 150;          % y resolution of main area
P.padfactor = 1.5;  % How much absorbing padding to add on the sides of the main area (1 means no padding, 2 means the absorbing padding on both sides is of thickness Lx_main/2)
P.dz_target = 0.2e-6; % [m] z step size to aim for
P.alpha = 3e14;             % [1/m^3] "Absorption coefficient" per squared unit length distance out from edge of main area

%% Reproducing Fig. 4
P.lambda = 1064e-9; % [m] Wavelength
P.n_cladding = 1.52; % [] Cladding refractive index
NA = 0.1;
n_core = sqrt(NA^2 + P.n_cladding^2);
P.n_0 = n_core;
P.Lz = 1e-2; % [m] z propagation distances for this segment

P.bendingRoC = 1.24e-2;
P.bendDirection = 0;
% P.rho_e = 0;

nModes = 30; % For mode finding
plotModes = true; % If true, will plot the found modes
sortByLoss = true; % If true, sorts the list of found modes in order of ascending loss. If false, sorts in order of ascending imaginary part of eigenvalue (descending propagation constant)

% In the shapes 2D array, each row is a shape such as a core in a fiber.
% Column 1 are the x coordinates, column 2 are the y coordinates, column 3
% are radii, column 4 are the types of the shapes, column 5 are the peak
% refractive indices and column 6 is the g parameter, only needed if any of
% the shapes are GRIN lenses.

% Shape types are 1: Circular step-index disk, 2: Antialiased circular
% step-index disk, 3: Parabolic graded index disk, 4: GRIN lens focusing in
% both x and y, 5: GRIN lens focusing only in y.
P.shapes = [ 0 0 12.5e-6    2  n_core];

% P.E can be either a function that takes X, Y and Eparameters as inputs
% and provides the complex E field as output, or it can be a struct with 3
% fields: a 'field' field which is the complex E-field matrix, and 'Lx' and
% 'Ly' fields that describe the side lengths of the provided E matrix. In
% the case of a struct, the provided E field will be adapted to the new
% grid using the interp2 function.
P = findModes(P,nModes,sortByLoss,plotModes);

%% Reproducing Fig 7, blue curve:
P.Lz = 1e-2; % [m] z propagation distances for this segment

P.lambda = 1320e-9; % [m] Wavelength
P.n_cladding = 1.447; % [] Cladding refractive index
NA = 0.117;
n_core = sqrt(NA^2 + P.n_cladding^2);

P.n_0 = n_core;

P.bendingRoC = 0.5e-2;
P.bendDirection = 0;
% P.rho_e = 0;

P.shapes = [ 0 0 4.1e-6    2  n_core];

plotModes = false; % If true, will plot the found modes
P = findModes(P,nModes,sortByLoss,plotModes);
P.E = P.modes(1);

% Run solver
P = FD_BPM(P);

bendloss = -log(P.powers(end))/P.Lz*(-10*log10(exp(-1))) % [dB/m]

% Results for wavelength 1320 nm: (Nx=150, dz=0.2e-6)
% RoC [m]   bendloss [dB/m]
% 1.5e-2     0.0055
% 1.25e-2   0.0289
% 1e-2      0.6303
% 0.75e-2   15.5092
% 0.5e-2    272.3184

%% Reproducing Fig 7, red curve:
P = clearData(P);

P.Lz = 1e-2; % [m] z propagation distances for this segment

P.lambda = 1550e-9; % [m] Wavelength
P.n_cladding = 1.440; % [] Cladding refractive index
NA = 0.117;
n_core = sqrt(NA^2 + P.n_cladding^2);
P.n_0 = n_core;

P.bendingRoC = 1.5e-2;
P.bendDirection = 0;
% P.rho_e = 0;

P.shapes = [ 0 0 4.1e-6    2  n_core];

plotModes = false; % If true, will plot the found modes
P = findModes(P,nModes,sortByLoss,plotModes);
P.E = P.modes(1);

% Run solver
P = FD_BPM(P);

bendloss = -log(P.powers(end))/P.Lz*(-10*log10(exp(-1))) % [dB/m]

% Results for wavelength 1550 nm:
% RoC [m]   bendloss [dB/m]
% 1.5e-2    0.5450
% 1.25e-2   3.1963
% 1.0e-2    54.3441
% 0.75e-2   180.1932
% 0.5e-2    1270.3429

S = load('train');
sound(S.y.*0.3,S.Fs);

%% USER DEFINED E-FIELD INITIALIZATION FUNCTION
function E = calcInitialE(X,Y,Eparameters) % Function to determine the initial E field. Eparameters is a cell array of additional parameters such as beam size
w_0 = 2.5e-6;
offset = 2.5e-6;
amplitude = exp(-((X-offset).^2+Y.^2)/w_0^2);
phase = zeros(size(X));
E = amplitude.*exp(1i*phase); % Electric field
end