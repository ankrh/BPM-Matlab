clear P % Parameters struct

% This example starts shows the propagation of an LP mode in a multimode
% fiber. There are three sections in this MMF out of which the centre one
% could be bent to see the how the propagation of any mode is affected when the
% fiber is bent. 

updatestepsize = 1e-5;
%% General and solver-related settings
P.name = mfilename;
P.useAllCPUs = true;

%% Visualization parameters
P.saveVideo = false; % To save the field intensity and phase profiles at different transverse planes
P.saveData = false; % To save the struct P  
dataName = ['Data\Data_MMF\Data_',P.name '.mat'];
P.videoName = ['Video\Video_',P.name '.avi'];
P.calcModeOverlaps = true;  % Set it to true to calculate mode overlap integrals of propagating field with respect to different modes in the P.modes struct array
if P.saveVideo 
  P.videoHandle = VideoWriter(P.videoName);  % Create the video handle if you want to save video from all the frames
  open(P.videoHandle);
end

%% Resolution-related parameters (check for convergence)
P.Lx_main = 30e-6;        % [m] x side length of main area
P.Ly_main = 30e-6;        % [m] y side length of main area
P.Nx_main = 150;          % x resolution of main area
P.Ny_main = 150;          % y resolution of main area
P.padfactor = 1.5;  % How much absorbing padding to add on the sides of the main area (1 means no padding, 2 means the absorbing padding on both sides is of thickness Lx_main/2)
P.dz_target = 0.5e-6; % [m] z step size to aim for
P.alpha = 3e14;             % [1/m^3] "Absorption coefficient" per squared unit length distance out from edge of main area

%% Problem definition
P.lambda = 800e-9; % [m] Wavelength
P.n_cladding = 1.45; % [] Cladding refractive index
P.n_0 = 1.4666;

P.shapes = [ 0 0 6e-6  2  1.4666];

tic;
%% Segment 1
P.Lz = 5e-4; % [m] z propagation distances for this segment
P.updates = P.Lz/updatestepsize;
P.bendingRoC = Inf;
P.bendDirection = 0;

nModes = 30; % For mode finding
plotModes = true; % If true, will plot the found modes
sortByLoss = false; % If true, sorts the list of found modes in order of ascending loss. If false, sorts in order of ascending imaginary part of eigenvalue (descending propagation constant)
P = findModes(P,nModes,sortByLoss,plotModes);

% P.E = P.modes(1); % The LP01 mode, relatively sensitive to bending-induced mode coupling
P.E = P.modes(3); % The LP11e mode, about equally sensitive as LP01
% P.E = P.modes(15); % The LP04 mode, less sensitive

P = FD_BPM(P);

%% Segment 2
P.Lz = 1e-3; % [m] z propagation distances for this segment
P.updates = P.Lz/updatestepsize;
P.bendingRoC = 10e-3;
P.bendDirection = 0;
P.shapes = P.shapesFinal;
P.E = P.Efinal;

P = FD_BPM(P);

%% Segment 3
P.Lz = 1e-3; % [m] z propagation distances for this segment
P.updates = P.Lz/updatestepsize;
P.bendingRoC = 10e-3;
P.bendDirection = 90;
P.shapes = P.shapesFinal;
P.E = P.Efinal;

P = FD_BPM(P);

%% Segment 4
P.Lz = 5e-4; % [m] z propagation distances for this segment
P.updates = P.Lz/updatestepsize;
P.bendingRoC = Inf;
P.bendDirection = 0;
P.shapes = P.shapesFinal;
P.E = P.Efinal;

P = FD_BPM(P);

totalTime = toc;

%% Save output data and video 
if P.saveData 
  save(dataName,'P','totalTime');
end
if P.saveVideo
	close(P.videoHandle);
end

% S = load('train');
% sound(S.y.*0.3,S.Fs);

%% USER DEFINED E-FIELD INITIALIZATION FUNCTION
function E = calcInitialE(X,Y,Eparameters) % Function to determine the initial E field. Eparameters is a cell array of additional parameters such as beam size

end

function [E_output, powers_output, modeOverlap_output, P_output] = addToSaveData(segment, E_out, powers_out, modeOverlap_out, P_out, E_output, powers_output, modeOverlap_output, P_output)
    E_output{segment} = E_out; 
    powers_output{segment} = powers_out; 
    modeOverlap_output{segment} = modeOverlap_out;
    P_output{segment} = P_out;
end