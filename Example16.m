clear P % Parameters struct

% This example shows offset and tilt of the E field in MMF splices.
% 
% First, the fundamental mode is injected into a step index MMF. Then,
% after the first segment, we simulate a bad splice by offsetting the beam
% by 4 µm in a 45 degrees direction (towards upper right) using the
% offsetField() function. The mode overlaps can be seen to change
% instantaneously at this splice. Next, after segment 2, another bad splice
% is simulated, in which segments 2 and 3 are at an angle of 2 degrees.
% Here we use the tiltField() function. Again, the mode overlaps change.
% 
% The phase gradient that tiltField() applies to the field depends on what
% the refractive index is in the wedge between the segments. tiltField
% assumes that the refractive index is uniformly equal to P.n_0, which
% should be a good approximation for weakly guiding fibers.

%% General and solver-related settings
P.name = mfilename;

%% Resolution-related parameters (check for convergence)
P.Lx_main = 30e-6;        % [m] x side length of main area
P.Ly_main = 30e-6;        % [m] y side length of main area
P.Nx_main = 200;          % x resolution of main area
P.Ny_main = 200;          % y resolution of main area
P.padfactor = 1.5;  % How much absorbing padding to add on the sides of the main area (1 means no padding, 2 means the absorbing padding on both sides is of thickness Lx_main/2)
P.dz_target = 1e-6; % [m] z step size to aim for
P.alpha = 5e14;             % [1/m^3] "Absorption coefficient" per squared unit length distance out from edge of main area

%% Problem definition, segment 1
P.lambda = 1000e-9; % [m] Wavelength
P.n_background = 1.45; % [] (may be complex) Background refractive index (in this case, the cladding)
P.n_0 = 1.46; % [] reference refractive index
P.n.func = @calcRI;

P.Lz = 2e-4; % [m] z propagation distances for this segment
P.updates = 20;            % Number of times to update plot. Must be at least 1, showing the final state.
P.plotEmax = 2;

P = findModes(P,10,'plotModes',false);

P.E = P.modes(1);

P.calcModeOverlaps = true;

% Run solver
P = FD_BPM(P);

%% Problem definition, segment 2
direction = 45; % [degrees] The angle of the direction in which the offset is performed
distance = 4e-6; % [m] The distancein the xy plane the field should be offset
P = offsetField(P,direction,distance);

P.Lz = 1e-3; % [m] z propagation distances for this segment
P.updates = 100;            % Number of times to update plot. Must be at least 1, showing the final state.

% Run solver
P = FD_BPM(P);

%% Problem definition, segment 3
direction = 90; % [degrees] The angle of the direction in the xy plane in which the tilt is performed
angle = 2; % [degrees] The angle of tilt imposed on the field
P = tiltField(P,direction,angle);

P.Lz = 1e-3; % [m] z propagation distances for this segment
P.updates = 100;            % Number of times to update plot. Must be at least 1, showing the final state.

% Run solver
P = FD_BPM(P);

%% USER DEFINED RI FUNCTIONS
function n = calcRI(X,Y,n_background,nParameters)
% n may be complex
n = n_background*ones(size(X)); % Start by setting all pixels to n_background
n(X.^2 + Y.^2 < 5e-6^2) = 1.46;
end