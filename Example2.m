clear P % Parameters struct

% This example consists of multiple segments of fiber, showing that the
% output of one simulation can be used as input to the next. Segment 2 is
% twisting and tapering. Segments 1-3 have 3 cores, while segment 4 only
% has one of the cores. Gaussian beams are launched into two of the cores,
% one of those being in an off-center and tilted orientation. One can see
% that as the fiber tapers, the smaller core no longer guides the light to
% any significant extent and the larger cores are also shedding the light
% in a lot of their modes. In segment 3, periodic mode coupling between the
% two larger cores can be seen.

%% General and solver-related settings
P.name = mfilename;
P.useAllCPUs = true; % If false, BPM-Matlab will leave one processor unused. Useful for doing other work on the PC while simulations are running.
P.useGPU = false; % (Default: false) Use CUDA acceleration for NVIDIA GPUs

%% Visualization parameters
P.updates = 30;            % Number of times to update plot. Must be at least 1, showing the final state.
P.downsampleImages = false; % Due to a weird MATLAB bug, MATLAB may crash when having created imagesc (or image) plots with dimensions larger than roughly 2500x2500 and then calling mex functions repeatedly. This flag will enable downsampling to 500x500 of all data before plotting, hopefully avoiding the issue.
P.displayScaling = 1;  % Zooms in on figures. Set to 1 for no zooming.

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
P.n_background = 1.45; % [] Cladding refractive index
P.n_0 = 1.46; % [] reference refractive index
P.Lz = 2e-3; % [m] z propagation distances for this segment
P.taperScaling = 1; % [] the ratio of the width of the structure at the end of the segment to the width at the beginning of the segment.
P.twistRate = 0; % [rad/m] the rate of rotation in units of radians per meter.
P.figTitle = 'Segment 1';

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
% and provides the complex E field as output, or it can be a struct with 3
% fields: a 'field' field which is the complex E-field matrix, and 'Lx' and
% 'Ly' fields that describe the side lengths of the provided E matrix. In
% the case of a struct, the provided E field will be adapted to the new
% grid using the interpn function.
P.E = @calcInitialE; % Defined at the end of this file

% Run solver
P = FD_BPM(P);

%% Next segment
P.figTitle = 'Segment 2';
P.Lz = 5e-3;
P.taperScaling = 0.15;
P.twistRate = 2*pi/P.Lz;

% Run solver
P = FD_BPM(P);

%% Next segment
P.figTitle = 'Segment 3';
P.Lz = 2e-3;
P.taperScaling = 1;
P.twistRate = 0;

% Run solver
P = FD_BPM(P);

%% Next segment
P.figTitle = 'Segment 4';
P.Lz = 3e-3;
P.shapes = P.shapes(3,:);

% Run solver
P = FD_BPM(P);

%% USER DEFINED E-FIELD INITIALIZATION FUNCTION
function E = calcInitialE(X,Y,Eparameters) % Function to determine the initial E field. Eparameters is a cell array of additional parameters such as beam size
w_0 = 5e-6;
amplitude1 = exp(-((X-1.5e-5).^2+Y.^2)/w_0^2);
amplitude2 = 2*exp(-((X+12e-6).^2+(Y+7e-6).^2)/w_0^2);
phase1 = zeros(size(X));
phase2 = 8e5*Y;
E = amplitude1.*exp(1i*phase1) + amplitude2.*exp(1i*phase2); % Electric field
end