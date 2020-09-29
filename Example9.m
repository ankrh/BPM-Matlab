clear P % Parameters struct

% This example consists of multiple segments of hexagonal multicore fiber.
% Segments 1-3 are straight fibres of equal lengths (Lz/3). Bending could
% be introduced to segment 2 for verifying bending-related field variations
% in the multicore fibre. Running Example7.m is a prerequisite, from which
% the acquired phase upon propagation is obtained

%% General and solver-related settings
savedCaseName = 'Example7_noInputPhase';
currentTestCaseName = '_withInputPhase';
P.name = [mfilename, currentTestCaseName];
P.useAllCPUs = true;
P.useGPU = false;

%% Visualization parameters
P.saveVideo = true; % To save the field intensity and phase profiles at different transverse planes
P.saveData = true; % To save the struct P  
P.updates = 30;            % Number of times to update plot. Must be at least 1, showing the final state.
P.downsampleImages = false; % Due to a weird MATLAB bug, MATLAB may crash when having created imagesc (or image) plots with dimensions larger than roughly 2500x2500 and then calling mex functions repeatedly. This flag will enable downsampling to 500x500 of all data before plotting, hopefully avoiding the issue.
P.displayScaling = 1;  % Zooms in on figures 1 & 3a,b. Set to 1 for no zooming.  
dataName = [P.name '.mat'];
P.videoName = [P.name '.avi'];
E_final = {};  % Initialising Eoutput array which is finally saved after all segment simulations 
powers_final = {}; 
if P.saveVideo 
  P.videoHandle = VideoWriter(P.videoName);  % Create the video handle if you want to save video from all the frames
  open(P.videoHandle);
end

%% Resolution-related parameters (check for convergence)
P.Lx_main = 200e-6;        % [m] x side length of main area
P.Ly_main = 200e-6;        % [m] y side length of main area
P.Nx_main = 600;          % x resolution of main area
P.Ny_main = 600;          % y resolution of main area
P.padfactor = 1.5;  % How much absorbing padding to add on the sides of the main area (1 means no padding, 2 means the absorbing padding on both sides is of thickness Lx_main/2)
P.dz_target = 0.4e-6; % [m] z step size to aim for
P.alpha = 3e14;             % [1/m^3] "Absorption coefficient" per squared unit length distance out from edge of main area
tic
%% Problem definition - straight multicore fibre
P.figTitle = 'Segment 1';
P.lambda = 980e-9; % [m] Wavelength
P.n_cladding = 1.45; % [] Cladding refractive index
P.n_0 = 1.46;
P.Lz = 0.1e-3; % [m] z propagation distances for this segment
P.taperScaling = 1;
P.twistRate = 0;
P.bendingRoC = Inf;
P.bendDirection = 0;
numberOfCores = 30;  %[] Number of cores in the multicore fibre
pitch = 15e-6; % [m] Intercore spacing
R = 2e-6; % [m] Core radius
n_core = 1.46; % Cores' refractive index 
w_0 = 2.35e-6; % [m] Input Gaussian waist
focusType = 1; % [] 1: Point focus, 2: Cylindrical focus at fibre distal end
k_0 = 2*pi/P.lambda; 
focalLength = 5e-3; % [m] Focal length of applied phase

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
P.Eparameters = {w_0, numberOfCores, P.shapes, k_0, focusType, focalLength, pitch, savedCaseName};
P.E = @calcInitialE; % Defined at the end of this file

% Run solver
P = FD_BPM(P);
% [E_final, powers_final] = addToSaveData(1, E_out, powers_out, E_final, powers_final);

%% Second segment - bent multicore fibre
P.figTitle = 'Segment 2';
P.Lz = 0.1e-3;
P.taperScaling = 1;
P.twistRate = 0;
P.bendingRoC = Inf;
P.bendDirection = 0;
P.shapes = P.shapesFinal;
P.E = P.Efinal;

% Run solver
P = FD_BPM(P);
% [E_final, powers_final] = addToSaveData(2, E_out, powers_out, E_final, powers_final);

%% Third segment - straight multicore fibre
P.figTitle = 'Segment 3';
P.Lz = 0.1e-3;
P.taperScaling = 1;
P.twistRate = 0;
P.bendingRoC = Inf;
P.bendDirection = 0;
P.shapes = P.shapesFinal;
P.E = P.Efinal;

% Run solver
P = FD_BPM(P);
% [E_final, powers_final] = addToSaveData(3, E_out, powers_out, E_final, powers_final);

%% Fourth segment - Free space FFTBPM propagation from fibre distal end
P.Lx_main = 1200e-6;        % [m] x side length of main area
P.Ly_main = 1200e-6;        % [m] y side length of main area
P.Nx_main = 2000;          % x resolution of main area
P.Ny_main = 2000;          % y resolution of main area
P.alpha = 8e13;             % [1/m^3] "Absorption coefficient" per squared unit length distance out from edge of main area

P.figNum = 2;
P.figTitle = 'In air';
P.Lz = focalLength;
P.E = P.Efinal;

% Run solver
[E_out_fft] = FFT_BPM(P);
[E_final, powers_final] = addToSaveData(4, E_out_fft, [], E_final, powers_final);

if P.saveData 
%     save(dataName, 'P','E_final','powers_final','shapes_out');
end
if P.saveVideo
	close(P.videoHandle);
end
toc
S = load('train');
sound(S.y.*0.3,S.Fs);

%% USER DEFINED E-FIELD INITIALIZATION FUNCTION
function E = calcInitialE(X,Y,Eparameters) % Function to determine the initial E field. Eparameters is a cell array of additional parameters such as beam size
w_0 = Eparameters{1};
numberOfCores = Eparameters{2};
shapeParameters = Eparameters{3};
k_0 = Eparameters{4};
focusType = Eparameters{5};
focalLength = Eparameters{6};
pitch = Eparameters{7};
savedCaseName = Eparameters{8};
Nx= size(X,1); Ny= size(Y,1);
dx = X(2,1) - X(1,1);
dy = Y(1,2) - Y(1,1);
try
    load([savedCaseName]);  % Load the output E field data after propagating through a straight multicore fiber
catch
    error('You need to run Example7.m and save the data before Example9.m');
end
amplitude = zeros(size(X));
for i = 1:numberOfCores
    amplitude = amplitude+exp(-((X-shapeParameters(i)).^2+(Y-shapeParameters(i+numberOfCores)).^2)/w_0^2);
end
phase = zeros(size(X));
acquiredPhase = NaN(1,numberOfCores);    % Row vector (1D)
focusPhase = NaN(1,numberOfCores);

for idx = 1:numberOfCores
  acquiredPhase(idx) = angle(E_final{end}.field(Nx/2+ceil(shapeParameters(idx)/dx),Ny/2+ceil(shapeParameters(idx+numberOfCores)/dy)));  %Acquired phase of E field (Estruct.field) at distal end for previous travel through straight fibre
  switch focusType
      case 1
          focusPhase(idx) = k_0*((shapeParameters(idx))^2+(shapeParameters(idx+numberOfCores))^2)/(2*focalLength);  %Focusing phase for point focus
      case 2
          focusPhase(idx) = k_0*(shapeParameters(idx))^2/(2*focalLength);  %Focusing phase for horizontal line focus
  end           
  phase(sqrt((X-shapeParameters(idx)).^2+(Y-shapeParameters(idx+numberOfCores)).^2) < pitch/2) = focusPhase(idx)-acquiredPhase(idx);
end
E = amplitude.*exp(1i*phase);
end

%% USER DEFINED SHAPE-PARAMETERS INITIALIZATION FUNCTION FOR MULTICORE FIBRE
function shapeParameters = getShapeParameters(FibreParameters)
numberOfCores = FibreParameters{1};
pitch = FibreParameters{2};
R = FibreParameters{3};
n_core = FibreParameters{4};
shapeType = FibreParameters{5};
rho_n = pitch*sqrt(1:numberOfCores);
theta_n = (1:numberOfCores)*pi*(3-sqrt(5));
x_n = rho_n.*cos(theta_n);
y_n = rho_n.*sin(theta_n);
shapeParameters = NaN(5,numberOfCores); % Initialize output array
shapeParameters(3,:) = R; % All cores have radius R
shapeParameters(4,:) = shapeType; % All cores have the same shape type
shapeParameters(5,:) = n_core; % All cores have the same refractive index
for coreIdx = 1:numberOfCores
    shapeParameters(1:2,coreIdx) = [x_n(coreIdx); y_n(coreIdx)];
end
shapeParameters = shapeParameters.'; % Format: x y R shapeType n_core: in one row
end

function [E_output, powers_output] = addToSaveData(segment, E_out, powers_out, E_output, powers_output)
    E_output{segment} = E_out; 
    powers_output{segment} = powers_out; 
end