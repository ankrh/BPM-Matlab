clear P % Parameters struct

% This example consists of a straight Fermat's Golden Spiral multicore
% fiber. An input E-field has been pre-computed and is loaded from the file
% ExampleInputField.mat. The field is first propagated through the
% multicore fiber with FDBPM and then through air using FFTBPM. With this
% E-field input, a focal spot is formed at a distance of 5 mm from the
% distal end of the fiber in air. Compare to example 7, in which the fibre
% is bent.

% A video is saved, showing both simulation segments. Because we want both
% segments in one video file, we have to defer saving the video file efter
% the first segment by setting P.finalizeVideo = false in the first segment

%% General and solver-related settings
P.useAllCPUs = true;
P.useGPU = false;

%% Visualization parameters
P.figTitle = 'In straight fibre';
P.videoName = 'Example6.avi';
P.saveVideo = true; % To save the field intensity and phase profiles at different transverse planes
P.finalizeVideo = false;
P.updates = 30;            % Number of times to update plot. Must be at least 1, showing the final state.
P.downsampleImages = false; % Due to a weird MATLAB bug, MATLAB may crash when having created imagesc (or image) plots with dimensions larger than roughly 2500x2500 and then calling mex functions repeatedly. This flag will enable downsampling to 500x500 of all data before plotting, hopefully avoiding the issue.
P.displayScaling = 1;  % Zooms in on figures. Set to 1 for no zooming.  

%% Resolution-related parameters (check for convergence)
P.Lx_main = 200e-6;        % [m] x side length of main area
P.Ly_main = 200e-6;        % [m] y side length of main area
P.Nx_main = 600;          % x resolution of main area
P.Ny_main = 600;          % y resolution of main area
P.padfactor = 1.1;  % How much absorbing padding to add on the sides of the main area (1 means no padding, 2 means the absorbing padding on both sides is of thickness Lx_main/2)
P.dz_target = 0.4e-6; % [m] z step size to aim for
P.alpha = 3e14;             % [1/m^3] "Absorption coefficient" per squared unit length distance out from edge of main area

%% Problem definition - straight multicore fibre
P.lambda = 980e-9; % [m] Wavelength
P.n_cladding = 1.45; % [] Cladding refractive index
P.n_0 = 1.46;
P.Lz = 0.3e-3; % [m] z propagation distances for this segment
P.taperScaling = 1;
P.twistRate = 0;
P.bendDirection = 0;
P.bendingRoC = Inf;

nCores = 30;  %[] Number of cores in the multicore fibre
pitch = 15e-6; % [m] Intercore spacing
R = 2e-6; % [m] Core radius
n_core = 1.46; % Cores' refractive index 

% In the shapes 2D array, each row is a shape such as a core in a fiber.
% Column 1 are the x coordinates, column 2 are the y coordinates, column 3
% are radii, column 4 are the types of the shapes, column 5 are the peak
% refractive indices and column 6 is the g parameter, only needed if any of
% the shapes are GRIN lenses.

% Shape types are 1: Circular step-index disk, 2: Antialiased circular
% step-index disk, 3: Parabolic graded index disk, 4: GRIN lens focusing in
% both x and y, 5: GRIN lens focusing only in y.
P.shapes = NaN(nCores,5);
P.shapes(:,1) = pitch*sqrt(1:nCores).*cos((1:nCores)*pi*(3-sqrt(5))); % Fermat's Golden Spiral core x positions
P.shapes(:,2) = pitch*sqrt(1:nCores).*sin((1:nCores)*pi*(3-sqrt(5))); % Fermat's Golden Spiral core y positions
P.shapes(:,3) = R;
P.shapes(:,4) = 2; % Anti-aliased step index disks
P.shapes(:,5) = n_core;

% P.E can be either a function that takes X, Y and Eparameters as inputs
% and provides the complex E field as output, or it can be a struct with 3
% fields: a 'field' field which is the complex E-field matrix, and 'Lx' and
% 'Ly' fields that describe the side lengths of the provided E matrix. In
% the case of a struct, the provided E field will be adapted to the new
% grid using the interp2 function.
load('exampleInputField.mat','E');
P.E = E;

% Run solver
P = FD_BPM(P);

%% Free space FFTBPM propagation from fibre distal end
P.Lx_main = 1200e-6;        % [m] x side length of main area
P.Ly_main = 1200e-6;        % [m] y side length of main area
P.Nx_main = 2000;          % x resolution of main area
P.Ny_main = 2000;          % y resolution of main area
P.alpha = 8e13;             % [1/m^3] "Absorption coefficient" per squared unit length distance out from edge of main area

P.figTitle = 'In air';
P.Lz = 5e-3;

P.finalizeVideo = true;

% Run solver
[E_out_fft] = FFT_BPM(P);

