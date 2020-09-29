clear P % Parameters struct

% This example simulates a photonic lantern with one few mode fiber (FMF) 
% as input and three output fibers: one standard single mode fiber (SSMF) 
% and two HI1016 fibers. The section with three output fibers are tapered to 
% increase the diameters of their cores and cladding uniformly. In this example,
% the LP11A mode of FMF is excited and the maximum output power will be
% coupled to the SSMF at the distal end. getLPmodes function is utilized to
% estimate the LP11A mode of FMF. 

%% Part 1 run with FDBPM
%% General and solver-related settings
P.name = mfilename;
P.useAllCPUs = true;
P.useGPU = false;

%% Visualization parameters
P.figNum = 1;
P.saveVideo = false; % To save the field intensity and phase profiles at different transverse planes
P.saveData = false; % To save the struct P  
P.updates = 50;            % Number of times to update plot. Must be at least 1, showing the final state.
P.downsampleImages = false; % Due to a weird MATLAB bug, MATLAB may crash when having created imagesc (or image) plots with dimensions larger than roughly 2500x2500 and then calling mex functions repeatedly. This flag will enable downsampling to 500x500 of all data before plotting, hopefully avoiding the issue.
P.displayScaling = 1;  % Zooms in on figures 1 & 3a,b. Set to 1 for no zooming.  
dataName = ['Data\Data_',P.name '_TMSI_LP11A.mat'];
P.videoName = ['Video\Video_',P.name '_TMSI_LP11A.avi'];
E_final = {};  % Initialising Eoutput array which is finally saved after all segment simulations 
powers_final = {}; 
modeOverlap_final = {};
P_final = {};
if P.saveVideo 
  P.videoHandle = VideoWriter(P.videoName);  % Create the video handle if you want to save video from all the frames
  open(P.videoHandle);
end

%The colormap options for the different subplots are 
%1: GBP, 2: HSV, 3:parula, 4: gray, 5: cividis
P.Intensity_colormap = 1; 
P.Phase_colormap = 2; 
P.n_colormap = 3; 

%% Resolution-related parameters (check for convergence)
P.Lx_main = 120e-6;        % [m] x side length of main area
P.Ly_main = 120e-6;        % [m] y side length of main area
P.Nx_main = 150;          % x resolution of main area
P.Ny_main = 150;          % y resolution of main area
P.padfactor = 1.5;  % How much absorbing padding to add on the sides of the main area (1 means no padding, 2 means the absorbing padding on both sides is of thickness Lx_main/2)
P.dz_target = 1e-6; % [m] z step size to aim for
P.alpha = 3e14;             % [1/m^3] "Absorption coefficient" per squared unit length distance out from edge of main area

%% Problem definition
P.lambda = 1550e-9; % [m] Wavelength
P.n_cladding = 1.4444; % [] Cladding refractive index - Air surroundings in the photonic lantern
P.n_0 = 1.4492;
P.Lz = 2e-2; % [m] z propagation distances for this segment
P.taperScaling = 1;  % [] Rate at which the segment profile is tapered
P.figTitle = 'TMSI Fiber';
P.n_colorlimits = [1.449 1.521];

% In the shapes 2D array, each row is a shape such as a core in a fiber.
% Column 1 are the x coordinates, column 2 are the y coordinates, column 3
% are radii, column 4 are the types of the shapes, column 5 are the peak
% refractive indices and column 6 is the g parameter, only needed if any of
% the shapes are GRIN lenses.

% Shape types are 1: Circular step-index disk, 2: Antialiased circular
% step-index disk, 3: Parabolic graded index disk, 4: GRIN lens focusing in
% both x and y, 5: GRIN lens focusing only in y.
P.shapes = [0 0 9e-6  2  1.4491];

% P.E can be either a function that takes X, Y and Eparameters as inputs
% and provides the complex E field as output, or it can be a struct with 3
% fields: a 'field' field which is the complex E-field matrix, and 'Lx' and
% 'Ly' fields that describe the side lengths of the provided E matrix. In
% the case of a struct, the provided E field will be adapted to the new
% grid using the interp2 function.
nModes = 5; % For mode finding
plotModes = false; % If true, will plot the found modes
sortByLoss = false; % If true, sorts the list of found modes in order of ascending loss. If false, sorts in order of ascending imaginary part of eigenvalue (descending propagation constant)
P = findModes(P,nModes,sortByLoss,plotModes);
P.E = P.modes(3);

% Set the Eparameters with LP mode type, core radius, n_core, n_cladding, 
% and wavelength to pass as parameters to the function getLPmodes inside 
% calcInitialE at the end of this script for calculating the LP modes this 
% MMF will support and set one mode as input E
P.Eparameters = {[2 1 'A'], P.shapes, P.n_cladding, P.lambda};  

% Run solver
tic
P = FD_BPM(P);
% [E_final, powers_final, modeOverlap_final, P_final] = addToSaveData(1, E_out, powers_out, modeOverlap, P, E_final, powers_final, modeOverlap_final, P_final);

%% Segment 2 run with FDBPM
P.Lx_main = 80e-6;        % [m] x side length of main area
P.Ly_main = 80e-6;        % [m] y side length of main area
P.Nx_main = 100;          % x resolution of main area
P.Ny_main = 100;          % y resolution of main area
P.figTitle = 'Gentle taper - 10 mm';
P.Lz = 10e-3; % [m] z propagation distances for this segment
P.n_cladding = 1.0; 
P.taperScaling = 50/23;
scaleFactor = 50/230*23/50;  % [] The factor by which the proximal end distances should be scaled down at TMSI interface
P.pitch = 3/4*125e-6*scaleFactor; % [m] Intercore separation in photonic lantern h=3/2 R_clad, R_clad is 125/2 um. This pitch is brought down by two tapers. 
P.shapes = [-P.pitch/2/sqrt(3)       P.pitch/2         62.5e-6*scaleFactor     1       1.4444;
                   P.pitch/sqrt(3)           0                     62.5e-6*scaleFactor     1       1.4444;
                   -P.pitch/2/sqrt(3)       -P.pitch/2        62.5e-6*scaleFactor     1       1.4444;
                   -P.pitch/2/sqrt(3)       P.pitch/2         4.5e-6*scaleFactor       1       1.4492;
                   P.pitch/sqrt(3)           0                     2.65e-6*scaleFactor     1       1.4511;
                   -P.pitch/2/sqrt(3)       -P.pitch/2        2.65e-6*scaleFactor     1       1.4511
                    ];
P.E = P.Efinal;
% P.n_colorlimits = [1.449 1.521];

% Run solver
P = FD_BPM(P);
% [E_final, powers_final, modeOverlap_final, P_final] = addToSaveData(2, E_out, powers_out, modeOverlap, P, E_final, powers_final, modeOverlap_final, P_final);

%% Segment 3 run with FDBPM
P.Lx_main = 250e-6;        % [m] x side length of main area
P.Ly_main = 250e-6;        % [m] y side length of main area
P.Nx_main = 200;          % x resolution of main area
P.Ny_main = 200;          % y resolution of main area
P.figTitle = 'Steep taper - 15 mm';
P.Lz = 15e-3; % [m] z propagation distances for this segment
P.shapes = P.shapesFinal;
P.taperScaling = 230/50;
P.E = P.Efinal;

% Run solver
P = FD_BPM(P);
% [E_final, powers_final, modeOverlap_final, P_final] = addToSaveData(3, E_out, powers_out, modeOverlap, P, E_final, powers_final, modeOverlap_final, P_final);

totalTime = toc; 
%% Save output data and video 
if P.saveData 
%     save(dataName, 'P_final','E_final','powers_final','shapes_out','modeOverlap_final','totalTime');
end
if P.saveVideo
	close(P.videoHandle);
end

S = load('train');
sound(S.y.*0.3,S.Fs);

%% USER DEFINED E-FIELD INITIALIZATION FUNCTION
function E = calcInitialE(X,Y,Eparameters) % Function to determine the initial E field. Eparameters is a cell array of additional parameters such as beam size
end

function [E_output, powers_output, modeOverlap_output, P_output] = addToSaveData(segment, E_out, powers_out, modeOverlap_out, P_out, E_output, powers_output, modeOverlap_output, P_output)
    E_output{segment} = E_out; 
    powers_output{segment} = powers_out; 
    modeOverlap_output{segment} = modeOverlap_out;
    P_output{segment} = P_out;
end