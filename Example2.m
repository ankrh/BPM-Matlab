clear P % Parameters struct

%% General and solver-related settings
P.name = mfilename;
P.useAllCPUs = true;
P.useGPU = false;
P.intNorm = false; % Choose true for field to be normalized w.r.t. max intensity, false to normalize such that total power is 1

%% Visualization parameters
P.saveVideo = false; % To save the field intensity and phase profiles at different transverse planes
P.updates = 30;            % Number of times to update plot. Must be at least 1, showing the final state.
P.downsampleImages = false; % Due to a weird MATLAB bug, MATLAB may crash when having created imagesc (or image) plots with dimensions larger than roughly 2500x2500 and then calling mex functions repeatedly. This flag will enable downsampling to 500x500 of all data before plotting, hopefully avoiding the issue.
P.displayScaling = 1;  % Zooms in on figures 1 & 3a,b. Set to 1 for no zooming.  

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
P.n_cladding = 1.45; % [] Cladding refractive index
P.n_0 = 1.46;
P.Lz = 2e-3; % [m] z propagation distances for this segment

% In the shapes 2D array, each row is a shape such as a core in a fiber.
% Column 1 are the x coordinates, column 2 are the y coordinates, column 3
% are radii, column 4 are the types of the shapes, column 5 are the peak
% refractive indices and column 6 is the g parameter, only needed if any of
% the shapes are GRIN lenses.

% Shape types are 1: Circular step-index disk, 2: Antialiased circular
% step-index disk, 3: Parabolic graded index disk, 4: GRIN lens focusing in
% both x and y, 5: GRIN lens focusing only in y.
P.shapes = [ -7e-6   -7e-6    10e-6  1  1.46;
             15e-6    0     1.25e-6  2  1.46;
              2e-6   12e-6    10e-6  3  1.465];

% P.E can be either a function that takes X, Y and Eparameters as inputs
% and provides the complex E field as output, or it can be a (2x1) cell array in
% which the first cell is a complex E field matrix and the second cell is a
% (2x1) array [Lx, Ly] that describe the side lengths of the provided E
% matrix. In the case of a cell array, the provided E field will be adapted
% to the new grid using the interp2 function.
P.E = @calcInitialE; % Defined at the end of this file

%% Run solver
[E_out,shapes_out] = FD_BPM(P);

%% Next segment
P.Lz = 5e-3;
P.taperScaling = 0.15;
P.twistRate = 2*pi/P.Lz;
P.shapes = shapes_out;
P.E = E_out;

%% Run solver
[E_out,shapes_out] = FD_BPM(P);

%% Next segment
P.Lz = 2e-3;
P.taperScaling = 1;
P.twistRate = 0;
P.shapes = shapes_out;
P.E = E_out;

%% Run solver
[E_out,shapes_out] = FD_BPM(P);

%% Next segment
P.Lz = 3e-3;
P.taperScaling = 1;
P.twistRate = 0;
P.shapes = shapes_out(3,:);
P.E = E_out;

%% Run solver
[E_out,shapes_out] = FD_BPM(P);

%% USER DEFINED E-FIELD INITIALIZATION FUNCTION
function E = calcInitialE(X,Y,Eparameters) % Function to determine the initial E field. Eparameters is a cell array of additional parameters such as beam size
w_0 = 5e-6;
amplitude1 = exp(-((X-1.5e-5).^2+Y.^2)/w_0^2);
amplitude2 = 2*exp(-((X+12e-6).^2+(Y+7e-6).^2)/w_0^2);
phase1 = zeros(size(X));
phase2 = 8e5*Y;
E = amplitude1.*exp(1i*phase1) + amplitude2.*exp(1i*phase2); % Electric field
end