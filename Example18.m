P = BPMmatlab.model;

% This example shows how to use the xSymmetry and ySymmetry settings to
% speed up calculations.
% 
% xSymmetry and ySymmetry describe which assumptions we make regarding 
% symmetries under mirroring in the x and y axes and can each be set to either
% 'NoSymmetry': (default), no symmetry assumption
% 'Symmetry': ordinary symmetry for both the refractive index and the electric field
% 'AntiSymmetry': ordinary symmetry for the refractive index and anti-symmetry for the electric field
% 
% In all simulations, FD_BPM is called without output variables, so that
% the P struct is not modified. This makes it easier for us to run the same
% simulation multiple times with different symmetry settings. Also in all
% simulations, the outer edge of the simulation area is the same, as well
% as the pixel sizes dx and dy. This is accomplished by halving Lx_main and
% Nx_main if there is y symmetry, and halving Ly_main and Ny_main if there
% is x symmetry.
% 
% In the first two simulations, a Gaussian beam is launched into a step
% index fiber, slightly offset in the x direction. First, the simulation is
% run with xSymmetry = 1, because the field has ordinary symmetry under
% mirroring in the x axis. Then the same simulation is run without symmetry
% assumptions. Observe how the result is the same, but the speed of the
% symmetry-assuming simulation is higher.
% 
% In the next two simulations, a mode superposition is created out of the
% four modes that are antisymmetric in x but symmetric in y (LP11o, LP31o,
% LP12o, LP51o). This superposition is launched first with xSymmetry = 2
% and ySymmetry = 1. Then the same simulation is run without symmetry
% assumptions. Again observe how the result is identical, but assuming
% symmetries will make the simulation faster.
% 
% When running a simulation for the first time, MATLAB may sometimes be
% slow, for some reason. Try running the simulation twice and/or setting
% updates = 1 to get the correct benchmarking speeds.

%% General and solver-related settings
P.name = mfilename;
P.useAllCPUs = false; % If false, BPM-Matlab will leave one processor unused. Useful for doing other work on the PC while simulations are running.
P.useGPU = false; % (Default: false) Use CUDA acceleration for NVIDIA GPUs

%% Visualization parameters
P.updates = 30;            % Number of times to update plot. Must be at least 1, showing the final state.
P.plotEmax = 1; % Max of color scale in the intensity plot, relative to the peak of initial intensity
P.figNum = 1;

%% Resolution-related parameters (check for convergence)
P.Lx_main = 20e-6;        % [m] x side length of main area
P.Ly_main = 10e-6;        % [m] y side length of main area
P.Nx_main = 200;          % x resolution of main area
P.Ny_main = 100;          % y resolution of main area
P.padfactor = 1.5;  % How much absorbing padding to add on the sides of the main area (1 means no padding, 2 means the absorbing padding on both sides is of thickness Lx_main/2)
P.dz_target = 1e-6; % [m] z step size to aim for
P.alpha = 3e14;             % [1/m^3] "Absorption coefficient" per squared unit length distance out from edge of main area

%% Problem definition
P.lambda = 1000e-9; % [m] Wavelength
P.n_background = 1.45; % [] (may be complex) Background refractive index (in this case, the cladding)
P.n_0 = 1.47; % [] reference refractive index
P.Lz = 2e-4; % [m] z propagation distances for this segment

P.xSymmetry = 'Symmetry'; % Symmetry under mirroring in the x axis
P.ySymmetry = 'NoSymmetry'; % Symmetry under mirroring in the y axis

P = initializeRIfromFunction(P,@calcRI);
P = initializeEfromFunction(P,@calcInitialE);

% Run solver
FD_BPM(P);

%% Same simulation, without use of symmetry
P.figNum = 4;
P.xSymmetry = 'NoSymmetry';
P.ySymmetry = 'NoSymmetry';
P.Ly_main = 20e-6;        % [m] y side length of main area
P.Ny_main = 200;          % y resolution of main area
FD_BPM(P);

%% Simulation with antisymmetry in x but ordinary symmetry in y
P.figNum = 8;
P.xSymmetry = 'AntiSymmetry';
P.ySymmetry = 'Symmetry';
P.Lx_main = 10e-6;        % [m] x side length of main area
P.Nx_main = 100;          % x resolution of main area
P.Ly_main = 10e-6;        % [m] y side length of main area
P.Ny_main = 100;          % y resolution of main area

P = findModes(P,10,'plotModes',true);
P.E = modeSuperposition(P,1:4);

FD_BPM(P);

%% Same as above, without use of symmetry
P.figNum = 12;
P.xSymmetry = 'NoSymmetry';
P.ySymmetry = 'NoSymmetry';
P.Lx_main = 20e-6;        % [m] x side length of main area
P.Nx_main = 200;          % x resolution of main area
P.Ly_main = 20e-6;        % [m] y side length of main area
P.Ny_main = 200;          % y resolution of main area

FD_BPM(P);

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
n(X.^2 + Y.^2 < 5e-6^2) = 1.47;
end