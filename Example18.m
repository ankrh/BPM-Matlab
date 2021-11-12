clear P % Parameters struct

% This example shows a x-symmetric RI and field distribution simulation

%% General and solver-related settings
P.name = mfilename;
P.useAllCPUs = false; % If false, BPM-Matlab will leave one processor unused. Useful for doing other work on the PC while simulations are running.
P.useGPU = false; % (Default: false) Use CUDA acceleration for NVIDIA GPUs
P.xSymmetry = 0; % Symmetry under mirroring in the x axis: 0: No symmetry, 1: Symmetric RIP and field, 2: Symmetric RIP and anti-symmetric field.
P.ySymmetry = 2; % Symmetry under mirroring in the y axis: 0: No symmetry, 1: Symmetric RIP and field, 2: Symmetric RIP and anti-symmetric field.

%% Visualization parameters
P.updates = 100;            % Number of times to update plot. Must be at least 1, showing the final state.
P.plotEmax = 0.5; % Max of color scale in the intensity plot, relative to the peak of initial intensity
P.figNum = 1;

%% Resolution-related parameters (check for convergence)
P.Lx_main = 20e-6/(1 + (P.ySymmetry ~= 0));        % [m] x side length of main area
P.Ly_main = 20e-6/(1 + (P.xSymmetry ~= 0));        % [m] y side length of main area
P.Nx_main = 200/(1 + (P.ySymmetry ~= 0));          % x resolution of main area
P.Ny_main = 200/(1 + (P.xSymmetry ~= 0));          % y resolution of main area
P.padfactor = 1.5;  % How much absorbing padding to add on the sides of the main area (1 means no padding, 2 means the absorbing padding on both sides is of thickness Lx_main/2)
P.dz_target = 1e-6; % [m] z step size to aim for
P.alpha = 3e14;             % [1/m^3] "Absorption coefficient" per squared unit length distance out from edge of main area

%% Problem definition
P.lambda = 1000e-9; % [m] Wavelength
P.n_background = 1.45; % [] (may be complex) Background refractive index (in this case, the cladding)
P.n_0 = 1.46; % [] reference refractive index
P.Lz = 2e-3; % [m] z propagation distances for this segment

P.n.func = @calcRI;

P.E = @calcInitialE; % Defined at the end of this file. See the readme file for details

P = findModes(P,5,'plotModes',true);

% % Run solver
% FD_BPM(P);
% 
% P.xSymmetry = 0;
% P.Lx_main = 20e-6/(1 + (P.ySymmetry ~= 0));        % [m] x side length of main area
% P.Ly_main = 20e-6/(1 + (P.xSymmetry ~= 0));        % [m] y side length of main area
% P.Nx_main = 200/(1 + (P.ySymmetry ~= 0));          % x resolution of main area
% P.Ny_main = 200/(1 + (P.xSymmetry ~= 0));          % y resolution of main area
% P.figNum = 4;
% FD_BPM(P);

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