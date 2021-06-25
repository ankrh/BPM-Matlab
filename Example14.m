clear P % Parameters struct

% This example is a reproduction of the RP Fiber Power demo file “Long
% period Bragg grating.fpw”, in which a fiber refractive index profile is
% defined with a super-Gaussian shape and a sinusoidal variation along z.
% This fiber functions as a mode converter from LP01 to LP03.




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

% To define the refractive index profile, you can use either P.shapes or
% P.n. Here we use P.n, which can be either a function that takes X, Y,
% n_background and nParameters as inputs and provides the complex (or real)
% refractive index profile as output, or it can be a struct with 3 fields:
% an 'n' field which is the complex refractive index matrix, and 'Lx' and
% 'Ly' fields that describe the side lengths of the provided n matrix. In
% the case of a struct, the provided n will be adapted to the new grid
% using the interpn function.
P.n.func = @calcRIP;
P.n.Nx = 200;
P.n.Ny = 200;
P.n.Nz = 1500;
P.n.Lx = 100e-6;
P.n.Ly = 100e-6;

P.nParameters = {0};

% Lx = 20e-6; Ly = 20e-6;
% Nx = 100; Ny = 80; Nz = 30;
% dx = Lx/Nx; dy = Ly/Ny; dz = P.Lz/(Nz-1);
% x = dx*(-(Nx-1)/2:(Nx-1)/2);
% y = dy*(-(Ny-1)/2:(Ny-1)/2);
% z = dz*(0:Nz-1);
% [X,Y,Z] = ndgrid(x,y,z);
% n = P.n_background + (1 + 0.001i)*max(0,0.02*(1 - ((X-3e-6)/6e-6).^2 - (Y/4e-6).^2 - (Z/1.8e-3).^2));
% P.n = struct('Lx',Lx,'Ly',Ly,'n',n);
% plotVolumetric(201,x,y,z,real(P.n.n));

nModes = 20; % For mode finding
plotModes = false; % If true, will plot the found modes
sortByLoss = false; % If true, sorts the list of found modes in order of ascending loss. If false, sorts in order of ascending imaginary part of eigenvalue (descending propagation constant)
singleCoreModes = false; % If true, finds modes for each core/shape individually. Note that the resulting "modes" will only be true modes of the entire structure if the core-to-core coupling is negligible.
P = findModes(P,nModes,singleCoreModes,sortByLoss,plotModes);
P.modes(1).label = 'LP01';
P.modes(6).label = 'LP02';
P.modes(15).label = 'LP03';
P.calcModeOverlaps = true;
Lambda = P.lambda/(P.modes(1).neff - P.modes(15).neff)
pointsPerPeriod = P.n.Nz/(P.Lz/Lambda)
P.nParameters = {2*pi/Lambda};

% P.E can be either a function that takes X, Y and Eparameters as inputs
% and provides the complex E field as output, or it can be a struct with 3
% fields: a 'field' field which is the complex E-field matrix, and 'Lx' and
% 'Ly' fields that describe the side lengths of the provided E matrix. In
% the case of a struct, the provided E field will be adapted to the new
% grid using the interpn function.
% P.E = @calcInitialE; % Defined at the end of this file
P.E = P.modes(1);

% Run solver
tic
P = FD_BPM(P);
toc

%% USER DEFINED RI INITIALIZATION FUNCTION
function n = calcRIP(X,Y,Z,n_background,nParameters) % This user-customizable function is the one that will be used to define the refractive index profile. You can either use analytical expressions or, e.g., load data from a file and use interpn to adapt the data to the simulation grid.
% n may be complex
n = n_background*ones(size(X)); % Start by setting all pixels to n_background
r_co = 30e-6;
r_cl = 40e-6;
n_cl = n_background;
n_co = 1.442;
dbeta = nParameters{1};
R = sqrt(X.^2 + Y.^2);
n(R <= r_cl) = n_cl + exp(-log(2)*(R(R <= r_cl)/r_co).^8).*(n_co - n_cl).*(1 + 0.1*sin(dbeta*Z(R <= r_cl))); % Equation for super-Gaussian modulating sinusoidally in z
end
