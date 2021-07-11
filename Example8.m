clear P % Parameters struct

% This example shows that modes can be found and stored in the parameters
% struct using the findModes function, then used as initial E fields. Here,
% we use the LP31o-like mode of a non-trivial refractive index profile.

%% General and solver-related settings
P.name = mfilename;
P.useAllCPUs = true; % If false, BPM-Matlab will leave one processor unused. Useful for doing other work on the PC while simulations are running.
P.useGPU = false; % (Default: false) Use CUDA acceleration for NVIDIA GPUs

%% Visualization parameters
P.updates = 100;            % Number of times to update plot. Must be at least 1, showing the final state.

%% Resolution-related parameters (check for convergence)
P.Lx_main = 30e-6;        % [m] x side length of main area
P.Ly_main = 25e-6;        % [m] y side length of main area
P.Nx_main = 300;          % x resolution of main area
P.Ny_main = 250;          % y resolution of main area
P.padfactor = 1.5;  % How much absorbing padding to add on the sides of the main area (1 means no padding, 2 means the absorbing padding on both sides is of thickness Lx_main/2)
P.dz_target = 5e-7; % [m] z step size to aim for
P.alpha = 3e14;             % [1/m^3] "Absorption coefficient" per squared unit length distance out from edge of main area

%% Problem definition
P.lambda = 1000e-9; % [m] Wavelength
P.n_background = 1.45; % [] (may be complex) Background refractive index, (in this case, the cladding)
P.n_0 = 1.46; % [] reference refractive index
P.Lz = 2e-3; % [m] z propagation distances for this segment

P.n.func = @calcRI;

% We search for the 10 modes with effective refractive index closest to
% n_0:
P = findModes(P,10);

P.E = P.modes(9); % The 9th mode is an LP31o-like mode

% Run solver
FD_BPM(P);

%% USER DEFINED RI FUNCTIONS
function n = calcRI(X,Y,n_background,nParameters)
% n may be complex
n = n_background*ones(size(X)); % Start by setting all pixels to n_background
% Core 1 is step index:
n((X + 2.5e-6).^2 + Y.^2 < 5e-6^2) = 1.46;

% Core 2 is graded index:
corepos = [5e-6 0];
r = 2.5e-6;
R = sqrt((X - corepos(1)).^2 + (Y - corepos(2)).^2);
n(R < r) = n_background + (1.47 - n_background)*(1 - (R(R < r)/r).^2); % Equation for parabolic graded index core
end
