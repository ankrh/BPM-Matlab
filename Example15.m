clear P % Parameters struct

% This example shows that the mode finder can correctly find the LPlm
% labels even of modes of nonstandard radially symmetric reractive index
% distributions.
% 
% The refractive index is chosen to approximate that of Figure 3 of
% https://www.rp-photonics.com/lp_modes.html. In the plot on the RP
% photonics site, the black curves correspond to the LP01 and LP02 modes,
% the red to the LP11 and LP12 modes, the green to the LP21 modes and the
% blue to the LP31 modes. Good agreement is seen between the modes found in
% the mode finder of BPM-Matlab and the modes on the RP site. Also, the
% LPlm labels found by BPM-Matlab (shown in the figure titles of the
% plotted modes) are correct. Interestingly, we see the LP01 mode exhibit a
% local minimum at the core center, a feature that is never present in the
% LP0m modes of step-index core fibers.
% 
% Plotting of the refractive index is added to the RI function for
% visualization.
% 
% The beam propagator is not run in this example.

%% General and solver-related settings
P.name = mfilename;

%% Resolution-related parameters (check for convergence)
P.Lx_main = 50e-6;        % [m] x side length of main area
P.Ly_main = 50e-6;        % [m] y side length of main area
P.Nx_main = 150;          % x resolution of main area
P.Ny_main = 150;          % y resolution of main area
P.padfactor = 1.25;  % How much absorbing padding to add on the sides of the main area (1 means no padding, 2 means the absorbing padding on both sides is of thickness Lx_main/2)
P.dz_target = 1e-6; % [m] z step size to aim for
P.alpha = 3e14;             % [1/m^3] "Absorption coefficient" per squared unit length distance out from edge of main area

%% Problem definition
P.lambda = 800e-9; % [m] Wavelength
P.n_background = 1.44; % [] (may be complex) Background refractive index, (in this case, the cladding)
P.n_0 = 1.443; % [] reference refractive index

P.n.func = @calcRI;

%% Mode finder
nModes = 15; % For mode finding
plotModes = true; % If true, will plot the found modes
sortByLoss = false; % If true, sorts the list of found modes in order of ascending loss. If false, sorts in order of ascending imaginary part of eigenvalue (descending propagation constant)
singleCoreModes = false; % If true, finds modes for each core/shape individually. Note that the resulting "modes" will only be true modes of the entire structure if the core-to-core coupling is negligible.
P = findModes(P,nModes,singleCoreModes,sortByLoss,plotModes);

%% USER DEFINED RI FUNCTIONS
function n = calcRI(X,Y,n_background,nParameters)
% n may be complex
n = n_background*ones(size(X)); % Start by setting all pixels to n_background
R = sqrt(X.^2 + Y.^2);
n(R < 12e-6) = n_background + 0.0042*exp(-((R(R < 12e-6) - 3e-6)/3.5e-6).^2) + 0.001*exp(-((R(R < 12e-6) - 6e-6)/2e-6).^2);

figure(201);clf;
imagesc(X(1:end,1),Y(1,1:end),n.');
axis equal tight xy;colorbar;
xlabel('x [m]');ylabel('y [m]');
title('Radially symmetric refractive index');
drawnow;
end