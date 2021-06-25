clear P % Parameters struct

% This example is a reproduction of the RP Fiber Power demo file "Long
% period Bragg grating.fpw", in which a fiber refractive index is
% defined with a super-Gaussian profile and a sinusoidal variation along z.
% This fiber grating functions as a mode converter between LP01 and LP03.

% The refractive index (RI) is here defined by setting P.n.func to
% a function, defined at the end of the model file, that takes X, Y, Z,
% n_background and nParameters as inputs and returns the 3D RI. The RI
% could be complex, but is in this case real. The grating has period Lambda
% and dbeta = 2*pi/Lambda. dbeta is in this case passed to the RI function
% as the first cell in the nParameters cell array.

% For 3D RI functions, the P.n.Nx, P.n.Ny, P.n.Nz, P.n.Lx and P.n.Ly fields
% are necessary but used only for calculation of the RI. They should be set
% to ensure that the RI is calculated in a large enough window and with
% high enough resolution. High values of P.n.Nx*P.n.Ny*P.n.Nz will require
% large amounts of memory (CPU or GPU).

% First, the grating length Lambda is set to infinity (dbeta = 0) and the
% modes of the unmodulated fiber are found. By manual inspection of the
% mode profiles we find that LP01 is mode 1, LP02 is mode 6 and LP03 is
% mode 15. From the theory of long period gratings we know that Lambda =
% lambda/(neff_1 - neff_15) in order to achieve phase matching, where
% lambda is the wavelength and neff_1 and neff_15 are the effective
% refractive indices of the LP01 and LP03 modes.

% After calculating the optimal Lambda and setting the corresponding dbeta
% in nParameters, FD_BPM is run. We can observe how the field's overlap with
% the LP03 mode reaches nearly 1 around z = 12.9 cm.

%% General and solver-related settings
P.name = mfilename;

%% Visualization parameters
P.updates = 100;            % Number of times to update plot. Must be at least 1, showing the final state.
P.plotEmax = 3; % Max of color scale in the intensity plot, relative to the peak of initial intensity

%% Resolution-related parameters (check for convergence)
P.Lx_main = 100e-6;        % [m] x side length of main area
P.Ly_main = 100e-6;        % [m] y side length of main area
P.Nx_main = 100;          % x resolution of main area
P.Ny_main = 100;          % y resolution of main area
P.padfactor = 1.25;  % How much absorbing padding to add on the sides of the main area (1 means no padding, 2 means the absorbing padding on both sides is of thickness Lx_main/2)
P.dz_target = 5e-6; % [m] z step size to aim for
P.alpha = 3e14;             % [1/m^3] "Absorption coefficient" per squared unit length distance out from edge of main area

%% Problem definition
P.lambda = 1060e-9; % [m] Wavelength
P.n_background = 1.44; % [] (may be complex) Background refractive index, (in this case, the cladding)
P.n_0 = 1.442; % [] reference refractive index
P.Lz = 140e-3; % [m] z propagation distances for this segment

P.n.func = @calcRI; % See the readme file for details
P.n.Nx = 200;
P.n.Ny = 200;
P.n.Nz = 1500;
P.n.Lx = 100e-6;
P.n.Ly = 100e-6;

P.nParameters = {0}; % No grating for initial mode finding

nModes = 20; % For mode finding
plotModes = true; % If true, will plot the found modes
sortByLoss = false; % If true, sorts the list of found modes in order of ascending loss. If false, sorts in order of ascending imaginary part of eigenvalue (descending propagation constant)
singleCoreModes = false; % If true, finds modes for each core/shape individually. Note that the resulting "modes" will only be true modes of the entire structure if the core-to-core coupling is negligible.
P = findModes(P,nModes,singleCoreModes,sortByLoss,plotModes);
P.modes(1).label = 'LP01';
P.modes(6).label = 'LP02';
P.modes(15).label = 'LP03';
P.calcModeOverlaps = true;

Lambda = P.lambda/(P.modes(1).neff - P.modes(15).neff)
P.nParameters = {2*pi/Lambda}; % Grating with correct period for phase matching

P.E = P.modes(1); % See the readme file for details

% Run solver
tic
P = FD_BPM(P);
toc

%% USER DEFINED RI INITIALIZATION FUNCTION
function n = calcRI(X,Y,Z,n_background,nParameters) % This user-customizable function is the one that will be used to define the refractive index. You can either use analytical expressions or, e.g., load data from a file and use interpn to adapt the data to the simulation grid.
% n may be complex
n = n_background*ones(size(X)); % Start by setting all pixels to n_background
r_co = 30e-6;
r_cl = 40e-6;
n_cl = n_background;
n_co = n_background + 0.002;
dbeta = nParameters{1};
R = sqrt(X.^2 + Y.^2);
n(R <= r_cl) = n_cl + exp(-log(2)*(R(R <= r_cl)/r_co).^8).*(n_co - n_cl).*(1 + 0.1*sin(dbeta*Z(R <= r_cl))); % Equation for super-Gaussian modulating sinusoidally in z
end
