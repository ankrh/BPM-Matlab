clear P % Parameters struct

% This example consists of multiple segments of hexagonal multicore fiber.
% Segments 1-3 are straight fibres of equal lengths (Lz/3). This example is
% meant to propagate the field without any input phase applied to the proximal
% end so as to retrieve the acquired phase at the fiber distal end upon
% propagation in Example 8

%% General and solver-related settings
P.name = [mfilename, '_noInputPhase'];
P.useAllCPUs = true;
P.useGPU = false;

%% Visualization parameters
P.saveVideo = false; % To save the field intensity and phase profiles at different transverse planes
P.saveData = true; % To save the struct P and output E field
P.updates = 30;            % Number of times to update plot. Must be at least 1, showing the final state.
P.downsampleImages = false; % Due to a weird MATLAB bug, MATLAB may crash when having created imagesc (or image) plots with dimensions larger than roughly 2500x2500 and then calling mex functions repeatedly. This flag will enable downsampling to 500x500 of all data before plotting, hopefully avoiding the issue.
P.displayScaling = 1;  % Zooms in on figures 1 & 3a,b. Set to 1 for no zooming.  
dataName = [P.name '.mat'];
E_final = {};  % Initialising Eoutput array which is finally saved after all segment simulations 
powers_final = {}; 

%% Resolution-related parameters (check for convergence)
P.Lx_main = 150e-6;        % [m] x side length of main area
P.Ly_main = 150e-6;        % [m] y side length of main area
P.Nx_main = 500;          % x resolution of main area
P.Ny_main = 500;          % y resolution of main area
P.padfactor = 1.5;  % How much absorbing padding to add on the sides of the main area (1 means no padding, 2 means the absorbing padding on both sides is of thickness Lx_main/2)
P.dz_target = 0.5e-6; % [m] z step size to aim for
P.alpha = 3e14;             % [1/m^3] "Absorption coefficient" per squared unit length distance out from edge of main area

%% Problem definition - straight multicore fibre
P.figTitle = 'Segment 1';
P.lambda = 980e-9; % [m] Wavelength
P.n_cladding = 1.45; % [] Cladding refractive index
P.n_0 = 1.46;
P.Lz = 0.1e-3; % [m] z propagation distances for this segment
P.taperScaling = 1;
P.twistRate = 0; %2*pi/P.Lz;
P.bendingRoC = Inf;
P.bendDirection = 0;
numberOfCores = 37;  %[] Number of cores in the multicore fibre
pitch = 20e-6; % [m] Intercore spacing
R = 2e-6; % [m] Core radius
n_core = 1.46; % Cores' refractive index 
w_0 = 2.35e-6; % [m] Input Gaussian waist

% In the shapes 2D array, each row is a shape such as a core in a fiber.
% Column 1 are the x coordinates, column 2 are the y coordinates, column 3
% are radii, column 4 are the types of the shapes, column 5 are the peak
% refractive indices and column 6 is the g parameter, only needed if any of
% the shapes are GRIN lenses.

% Shape types are 1: Circular step-index disk, 2: Antialiased circular
% step-index disk, 3: Parabolic graded index disk, 4: GRIN lens focusing in
% both x and y, 5: GRIN lens focusing only in y.
shapeType = 2; % [] Shape type of all the cores

FibreParameters = {numberOfCores, pitch, R, n_core, shapeType}; 
P.shapes = getShapeParameters(FibreParameters);  % Defined at the end of this file

% P.E can be either a function that takes X, Y and Eparameters as inputs
% and provides the complex E field as output, or it can be a struct with 3
% fields: a 'field' field which is the complex E-field matrix, and 'Lx' and
% 'Ly' fields that describe the side lengths of the provided E matrix. In
% the case of a struct, the provided E field will be adapted to the new
% grid using the interp2 function.
P.Eparameters = {w_0, numberOfCores, P.shapes};
P.E = @calcInitialE; % Defined at the end of this file

% Run solver
[E_out,shapes_out,powers_out,P] = FD_BPM(P);
[E_final, powers_final] = addToSaveData(1, E_out, powers_out, E_final, powers_final);

%% Second segment - bent multicore fibre
P.figTitle = 'Segment 2';
P.Lz = 0.1e-3;
P.taperScaling = 1;
P.twistRate = 0; %2*pi/P.Lz;
P.bendingRoC = Inf;
P.bendDirection = 0;
P.shapes = shapes_out;
P.E = E_out;

% Run solver
[E_out,shapes_out,powers_out,P] = FD_BPM(P);
[E_final, powers_final] = addToSaveData(2, E_out, powers_out, E_final, powers_final);

%% Third segment - straight multicore fibre
P.figTitle = 'Segment 3';
P.Lz = 0.1e-3;
P.taperScaling = 1;
P.twistRate = 0; %2*pi/P.Lz;
P.bendingRoC = Inf;
P.bendDirection = 0;
P.shapes = shapes_out;
P.E = E_out;

% Run solver
[E_out,shapes_out,powers_out,P] = FD_BPM(P);
[E_final, powers_final] = addToSaveData(3, E_out, powers_out, E_final, powers_final);

if P.saveData 
    save(dataName, 'P','E_final','powers_final','shapes_out');
end

S = load('train');
sound(S.y.*0.3,S.Fs);

%% USER DEFINED E-FIELD INITIALIZATION FUNCTION
function E = calcInitialE(X,Y,Eparameters) % Function to determine the initial E field. Eparameters is a cell array of additional parameters such as beam size
w_0 = Eparameters{1};
numberOfCores = Eparameters{2};
shapeParameters = Eparameters{3};
amplitude = zeros(size(X));
for i = 1:numberOfCores
    amplitude = amplitude+exp(-((X-shapeParameters(i)).^2+(Y-shapeParameters(i+numberOfCores)).^2)/w_0^2);
end
phase = zeros(size(X));
E = amplitude.*exp(1i*phase);
end

%% USER DEFINED SHAPE-PARAMETERS INITIALIZATION FUNCTION FOR MULTICORE FIBRE
function shapeParameters = getShapeParameters(FibreParameters)
numberOfCores = FibreParameters{1};
pitch = FibreParameters{2};
R = FibreParameters{3};
n_core = FibreParameters{4};
shapeType = FibreParameters{5};
shapeParameters = NaN(5,numberOfCores); % Initialize output array
shapeParameters(3,:) = R; % All cores have radius R
shapeParameters(4,:) = shapeType; % All cores have the same shape type
shapeParameters(5,:) = n_core; % All cores have the same refractive index
shapeParameters(1:2,1) = [0; 0]; % x and y of the center core
shellSideIdx = 1; % Which side of the hex are we filling?
shellSideCoreIdx = 0; % How many cores on this side have been filled so far?
shellNum = 1; % Which shell are we in? The center core is not counted as a shell.
for coreIdx = 2:numberOfCores
  if shellSideCoreIdx == 0 % If this is the first core in this shell
    shapeParameters(1:2,coreIdx) = [shellNum*pitch; 0];
  else % Find new core position by adding onto the previous core's position
    shapeParameters(1:2,coreIdx) = shapeParameters(1:2,coreIdx-1) + [pitch*cos(shellSideIdx*pi/3 + pi/3); pitch*sin(shellSideIdx*pi/3 + pi/3)];
  end

  if shellSideCoreIdx == shellNum % If this side has been filled
    shellSideIdx = shellSideIdx + 1;
    shellSideCoreIdx = 1;
  else % Continue filling this side
    shellSideCoreIdx = shellSideCoreIdx + 1;
  end

  if shellSideCoreIdx == shellNum && shellSideIdx == 6 % Last core on last side would be a replicate of the first one drawn in this shell, so skip
    shellNum = shellNum + 1;
    shellSideIdx = 1;
    shellSideCoreIdx = 0;
  end
end
shapeParameters = shapeParameters.'; % Format: x y R shapeType n_core: in one row
end

function [E_output, powers_output] = addToSaveData(segment, E_out, powers_out, E_output, powers_output)
    E_output{segment} = E_out; 
    powers_output{segment} = powers_out; 
end