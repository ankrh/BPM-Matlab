P = BPMmatlab.model;

% This example shows a basic case of launching a Gaussian beam
% into the GRINTECH lens GT-CFRL-100-025-20-CC (810) realized using FDBPM 
% and the subsequent propagation of the E field from GRIN lens's output in
% air using FFTBPM, showing that the resulting beam is collimated.

%% Part 1 run with FDBPM
%% General and solver-related settings
P.name = mfilename;
P.useAllCPUs = true; % If false, BPM-Matlab will leave one processor unused. Useful for doing other work on the PC while simulations are running.
P.useGPU = false; % (Default: false) Use CUDA acceleration for NVIDIA GPUs

%% Visualization parameters
P.figNum = 1;
P.figTitle = 'In GRIN Lens';
P.updates = 100;            % Number of times to update plot. Must be at least 1, showing the final state.
P.plotZoom = 1;  % Zooms in on figures. Set to 1 for no zooming.  

% The colormap options for the different subplots are 
% GPBGYR, HSV, parula, gray, cividis
P.intensityColormap = 'gray';

%% Resolution-related parameters (check for convergence)
P.Lx_main = 0.7e-3;        % [m] x side length of main area
P.Ly_main = 0.7e-3;        % [m] y side length of main area
P.Nx_main = 1000;          % x resolution of main area
P.Ny_main = 1000;          % y resolution of main area
P.padfactor = 1;  % How much absorbing padding to add on the sides of the main area (1 means no padding, 2 means the absorbing padding on both sides is of thickness Lx_main/2)
P.dz_target = 20e-6; % [m] z step size to aim for
P.alpha = 3e14;             % [1/m^3] "Absorption coefficient" per squared unit length distance out from edge of main area

%% Problem definition
P.lambda = 810e-9; % [m] Wavelength
P.n_background = 1.0; % [] (may be complex) Background refractive index, (in this case, air)
P.n_0 = 1.521; % [] reference refractive index
P.Lz = 6.05e-3; % [m] z propagation distances for this segment

P = initializeRIfromFunction(P,@calcRI);
P = initializeEfromFunction(P,@calcInitialE);

% Run solver
P = FD_BPM(P);

%% Output E from FDBPM propagating in air with FFTBPM
P.figNum = 2;
P.figTitle = 'In air';
P.n_0 = 1;

P.Lx_main = 700e-6;        % [m] x side length of main area
P.Ly_main = 700e-6;        % [m] y side length of main area
P.Nx_main = 500;          % x resolution of main area
P.Ny_main = 500;          % y resolution of main area
P.alpha = 8e13;             % [1/m^3] "Absorption coefficient" per squared unit length distance out from edge of main area
P.Lz = 1e-2; % [m] z propagation distances for this segment

% Run solver
FFT_BPM(P);

%% USER DEFINED E-FIELD INITIALIZATION FUNCTION
function E = calcInitialE(X,Y,Eparameters) % Function to determine the initial E field. Eparameters is a cell array of additional parameters such as beam size
w_0 = 5e-6;
offset = 0;
amplitude = exp(-((X-offset).^2+Y.^2)/w_0^2);
phase = zeros(size(X));
E = amplitude.*exp(1i*phase); % Electric field
end

%% USER DEFINED RI FUNCTIONS
function n = calcRI(X,Y,n_background,nParameters)
% n may be complex
n = n_background*ones(size(X)); % Start by setting all pixels to n_background
R = sqrt(X.^2 + Y.^2);
g = 260;
n(R < 500e-6) = 1.521*sech(g*R(R < 500e-6));
end

