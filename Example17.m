P = BPMmatlab.model;

% This example shows the speedup that can be obtained by launching
% simulations on the GPU compared to the CPU. The speedup is greatest for
% large Nx and Ny, such as in this example. The number of updates is set to
% only the mandatory 1, to minimize the time that needs to be spent on plot
% updates. The simulation can take around 60 seconds on CPU, so be patient.

% On a test PC with a 3.5 GHz Intel Xeon E5-1650 v3 CPU (launched in 2014)
% and an Nvidia GeForce GTX 970 GPU (also launched in 2014), the speedup
% for this example is about x14.

% On another test PC with a 2.9 GHz Intel Core i5-9400F CPU (launched 2019)
% and an Nvidia GeForce RTX 2070 GPU (launched 2018), the speedup is about
% x18.

% Higher-end GPUs are not likely to offer a much greater speedup than this,
% since the paraxial Helmholtz equation that BPM-Matlab is based on does
% not lend itself particularly well to parallelized solving.

%% General and solver-related settings
P.name = mfilename;
P.useAllCPUs = true; % If false, BPM-Matlab will leave one processor unused. Useful for doing other work on the PC while simulations are running.

%% Visualization parameters
P.updates = 1;            % Number of times to update plot. Must be at least 1, showing the final state.
P.plotEmax = 0.5; % Max of color scale in the intensity plot, relative to the peak of initial intensity

%% Resolution-related parameters (check for convergence)
P.Lx_main = 30e-6;        % [m] x side length of main area
P.Ly_main = 30e-6;        % [m] y side length of main area
P.Nx_main = 1500;          % x resolution of main area
P.Ny_main = 1500;          % y resolution of main area
P.padfactor = 1.5;  % How much absorbing padding to add on the sides of the main area (1 means no padding, 2 means the absorbing padding on both sides is of thickness Lx_main/2)
P.dz_target = 5e-7; % [m] z step size to aim for
P.alpha = 3e14;             % [1/m^3] "Absorption coefficient" per squared unit length distance out from edge of main area

%% Problem definition
P.lambda = 1000e-9; % [m] Wavelength
P.n_background = 1.45; % [] (may be complex) Background refractive index (in this case, the cladding)
P.n_0 = 1.46; % [] reference refractive index
P.Lz = 1e-4; % [m] z propagation distances for this segment

P = initializeRIfromFunction(P,@calcRI);
P = initializeEfromFunction(P,@calcInitialE);

% Run solver
fprintf('CPU:\n');
P.useGPU = false; % (Default: false) Use CUDA acceleration for NVIDIA GPUs
FD_BPM(P);
fprintf('GPU:\n');
P.useGPU = true; % (Default: false) Use CUDA acceleration for NVIDIA GPUs
P = FD_BPM(P);

%% USER DEFINED E-FIELD INITIALIZATION FUNCTION
function E = calcInitialE(X,Y,Eparameters) % Function to determine the initial E field. Eparameters is a cell array of additional parameters such as beam size
w_0 = 2.5e-6;
offset = 2.5e-6;
amplitude = exp(-((X-offset).^2+Y.^2)/w_0^2);
phase = zeros(size(X));
E = amplitude.*exp(1i*phase); % Electric field
end

%% USER DEFINED RI FUNCTIONS
function n = calcRI(X,Y,n_background,nParameters)
% n may be complex
n = n_background*ones(size(X)); % Start by setting all pixels to n_background
n(X.^2 + Y.^2 < 5e-6^2) = 1.46;
end