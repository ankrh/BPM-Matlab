clear P % Parameters struct

% This example illustrates the adiabatic mode conversion and shedding of
% unguided modes illustrated in figure 4 of the article "BPM-Matlab - An
% open-source optical propagation simulation tool in MATLAB". The model
% starts with a fiber with core radius 5 µm. Based on the results from
% findModes, we can see that it supports 8 modes:
% {'LP01';'LP11o';'LP11e';'LP21o';'LP21e';'LP02';'LP31o';'LP31e'}. The
% fiber tapers down to a final width of 1.5 µm, where only the LP01 mode is
% guided. Initially, we excite the fiber with a superposition of all 8
% guided modes.
% 
% From theory, we expect that
% At r = 4.80 µm, the LP31 modes are no longer guided (z = 0.57e-3 m)
% At r = 3.64 µm, the LP02 mode is no longer guided (z = 3.9e-3 m)
% At r = 3.58 µm, the LP21 modes are no longer guided (z = 4.1e-3 m)
% At r = 2.25 µm, the LP11 modes are no longer guided (z = 7.9e-3 m)
% 
% However, the light is not immediately lost when the fiber radius passes
% the above thresholds, it takes some distance before the light exits the
% simulation window and is absorbed in the surrounding absorber layer. At
% the end of the simulations, we see that the remaining light is that of
% the fundamental mode, and it contains 1/8th of the initial power (12.5
% %), which fits with adiabatic mode conversion of the initial fundamental
% mode to the final fundamental mode.
% 
% This example also illustrates the plotting of xz and yz slices. Each
% segment stores its slice data in separate elements of the xzSlice and
% yzSlice cell arrays.


%% General and solver-related settings
P.name = mfilename;
P.useAllCPUs = true;
P.useGPU = false;

%% Visualization parameters
updateFrequency = 5e4; % [1/m]
P.saveVideo = false; % To save the field intensity and phase profiles at different transverse planes
P.downsampleImages = false; % Due to a weird MATLAB bug, MATLAB may crash when having created imagesc (or image) plots with dimensions larger than roughly 2500x2500 and then calling mex functions repeatedly. This flag will enable downsampling to 500x500 of all data before plotting, hopefully avoiding the issue.
P.displayScaling = 1;  % Zooms in on figures 1 & 3a,b. Set to 1 for no zooming.  
P.plotEmax = 1; % Max of color scale in the intensity plot, relative to the peak of initial intensity

%% Resolution-related parameters (check for convergence)
P.Lx_main = 25e-6;        % [m] x side length of main area
P.Ly_main = 25e-6;        % [m] y side length of main area
P.Nx_main = 201;          % x resolution of main area
P.Ny_main = 201;          % y resolution of main area
P.padfactor = 1.5;  % How much absorbing padding to add on the sides of the main area (1 means no padding, 2 means the absorbing padding on both sides is of thickness Lx_main/2)
P.dz_target = 1e-6; % [m] z step size to aim for
P.alpha = 6e14;             % [1/m^3] "Absorption coefficient" per squared unit length distance out from edge of main area

%% Problem definition
P.lambda = 1000e-9; % [m] Wavelength
P.n_cladding = 1.45; % [] Cladding refractive index
P.n_0 = 1.46;
P.Lz = 10e-3; % [m] z propagation distances for this segment
P.updates = P.Lz*updateFrequency;            % Number of times to update plot. Must be at least 1, showing the final state.
P.taperScaling = 0.3;
P.figTitle = 'Segment 1';

% In the shapes 2D array, each row is a shape such as a core in a fiber.
% Column 1 are the x coordinates, column 2 are the y coordinates, column 3
% are radii, column 4 are the types of the shapes, column 5 are the peak
% refractive indices and column 6 is the g parameter, only needed if any of
% the shapes are GRIN lenses.

% Shape types are 1: Circular step-index disk, 2: Antialiased circular
% step-index disk, 3: Parabolic graded index disk, 4: GRIN lens focusing in
% both x and y, 5: GRIN lens focusing only in y.
P.shapes = [ 0 0  5e-6  2  1.46];

% P.E can be either a function that takes X, Y and Eparameters as inputs
% and provides the complex E field as output, or it can be a struct with 3
% fields: a 'field' field which is the complex E-field matrix, and 'Lx' and
% 'Ly' fields that describe the side lengths of the provided E matrix. In
% the case of a struct, the provided E field will be adapted to the new
% grid using the interp2 function.
nModes = 10; % For mode finding
plotModes = false; % If true, will plot the found modes
sortByLoss = false; % If true, sorts the list of found modes in order of ascending loss. If false, sorts in order of ascending imaginary part of eigenvalue (descending propagation constant)
P = findModes(P,nModes,sortByLoss,plotModes);
randvec = [0.974931557198556 + 0.222504963491603i ,  0.583350562881419 - 0.812220487790066i ,  0.915224078253766 - 0.402945264998293i , -0.432956662664553 - 0.901414737096290i ,  0.048613506502203 - 0.998817664534203i , -0.043136586317590 - 0.999069184251454i , -0.779356851920114 + 0.626580319963187i , -0.559597395571616 - 0.828764595569493i]; % Randomly generated array of mode coefficients for use in the superposition

P.E = P.modes(1); % We do this to get the correct values for the "Lx" and "Ly" fields.
P.E.field = randvec(1)*P.E.field; % The field needs to be multiplied by its coefficient
for iMode = 2:8
  P.E.field = P.E.field + randvec(iMode)*P.modes(iMode).field; % And here we add in the field contributions from the other modes
end
P.E.field = P.E.field*sqrt(sum(abs(P.E.field(:)).^2)/sum(abs(P.E.field(:)).^2)); % Normalize

% Run solver
P = FD_BPM(P);

%% Next segment
P.figTitle = 'Segment 2';
P.Lz = 5e-3;
P.updates = P.Lz*updateFrequency;            % Number of times to update plot. Must be at least 1, showing the final state.
P.taperScaling = 1;

% Run solver
P = FD_BPM(P);

%% Plot xz and yz slices
figure(2);clf;
subplot(2,1,1);
imagesc(P.z,P.x,abs([P.xzSlice{1} P.xzSlice{2}]).^2);
axis xy;
xlabel('z [m]');
ylabel('x [m]');
subplot(2,1,2);
imagesc(P.z,P.y,abs([P.yzSlice{1} P.yzSlice{2}]).^2);
axis xy;
xlabel('z [m]');
ylabel('y [m]');
