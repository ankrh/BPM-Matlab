P = BPMmatlab.model;

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
P.plotZoom = 1;  % Zooms in on figures. Set to 1 for no zooming.

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
P.n_background = 1.45; % [] (may be complex) Background refractive index, (in this case, the cladding)
P.n_0 = 1.46; % [] reference refractive index
P.Lz = 2e-3; % [m] z propagation distances for this segment
P.taperScaling = 1; % [] the ratio of the width of the structure at the end of the segment to the width at the beginning of the segment.
P.twistRate = 0; % [rad/m] the rate of rotation in units of radians per meter.
P.figTitle = 'Segment 1';

P = initializeRIfromFunction(P,@calcRI);
P = initializeEfromFunction(P,@calcInitialE);

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
P = initializeRIfromFunction(P,@calcRIsegment4);

% Run solver
P = FD_BPM(P);

%% USER DEFINED E-FIELD INITIALIZATION FUNCTION
function E = calcInitialE(X,Y,Eparameters) % Function to determine the initial E field. Eparameters is a cell array of additional parameters such as beam size
w_0 = 5e-6;
amplitude1 = exp(-((X-1.5e-5).^2+Y.^2)/w_0^2);
amplitude2 = 2*exp(-((X+12e-6).^2+(Y+7e-6).^2)/w_0^2);
phase1 = zeros(size(X));
phase2 = 8e5*Y; % This beam component is defined with a phase that increases linearly with the y-coordinate, corresponding to a beam that is tilted slightly in the -y direction
E = amplitude1.*exp(1i*phase1) + amplitude2.*exp(1i*phase2); % Electric field
end

%% USER DEFINED RI FUNCTIONS
function n = calcRI(X,Y,n_background,nParameters)
% n may be complex
n = n_background*ones(size(X)); % Start by setting all pixels to n_background
% Cores 1 and 2 are step index:
n((X + 7e-6).^2 + (Y + 7e-6).^2 < 10e-6^2) = 1.46;
n((X - 15e-6).^2 + Y.^2 < 1.25e-6^2) = 1.46;

% Core 3 is graded index:
corepos = [2e-6 12e-6];
r = 10e-6;
R3 = sqrt((X - corepos(1)).^2 + (Y - corepos(2)).^2);
n(R3 < r) = n_background + (1.465 - n_background)*(1 - (R3(R3 < r)/r).^2); % Equation for parabolic graded index core
end

function n = calcRIsegment4(X,Y,n_background,nParameters)
% n may be complex
n = n_background*ones(size(X)); % Start by setting all pixels to n_background

% Here we only want the third core, which has been tapered down and rotated
% by 2*pi:
corepos_initial = [2e-6 12e-6];
r_initial = 10e-6;
corepos_final = 0.15*corepos_initial; % The cores have rotated one whole revolution (2*pi), so we don't need to apply a rotation here.
r_final = 0.15*r_initial;
R3 = sqrt((X - corepos_final(1)).^2 + (Y - corepos_final(2)).^2);
n(R3 < r_final) = n_background + (1.465 - n_background)*(1 - (R3(R3 < r_final)/r_final).^2); % Equation for parabolic graded index core
end
