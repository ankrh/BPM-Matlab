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
P.name = mfilename;
P.useAllCPUs = true; % If false, BPM-Matlab will leave one processor unused. Useful for doing other work on the PC while simulations are running.
P.useGPU = false; % (Default: false) Use CUDA acceleration for NVIDIA GPUs

%% Visualization parameters
P.figTitle = 'In straight fibre';
P.videoName = 'Example6.avi';
P.saveVideo = true; % To save the field intensity and phase profiles at different transverse planes
P.finalizeVideo = false; % finalizeVideo should only be set to true in the last segment to be simulated.
P.updates = 30;            % Number of times to update plot. Must be at least 1, showing the final state.
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
P.n_background = 1.45; % [] (may be complex) Background refractive index, (in this case, the cladding)
P.n_0 = 1.46; % [] reference refractive index
P.Lz = 0.3e-3; % [m] z propagation distances for this segment
P.bendDirection = 0; % [degrees] direction of the bending, in a polar coordinate system with 0° to the right (towards positive x) and increasing angles in counterclockwise direction
P.bendingRoC = Inf; % [m] radius of curvature of the bend

P.n.func = @calcRI;

load('./ExampleData/ExampleInputField.mat','E');
P.E = E; % See the readme file for details

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

P.finalizeVideo = true; % finalizeVideo should only be set to true in the last segment to be simulated.

% Run solver
[E_out_fft] = FFT_BPM(P);

%% USER DEFINED RI FUNCTIONS
function n = calcRI(X,Y,n_background,nParameters)
% n may be complex
nCores = 30;
pitch = 15e-6;
R = 2e-6;
n_core = 1.46;

n = n_background*ones(size(X)); % Start by setting all pixels to n_background
for iCore=1:nCores
  xCore = pitch*sqrt(iCore).*cos(iCore*pi*(3-sqrt(5)));
  yCore = pitch*sqrt(iCore).*sin(iCore*pi*(3-sqrt(5)));
  n((X - xCore).^2 + (Y - yCore).^2 < R^2) = n_core;
end
end