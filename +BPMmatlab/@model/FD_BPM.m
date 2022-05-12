function P = FD_BPM(P)
% Authors: Madhu Veetikazhy and Anders K. Hansen
% DTU Health and DTU Fotonik
% 
% Electric field propagation through fibre using Finite Difference Beam
% Propagation Method in 2D Reference: An assessment of FD BPM method -
% Youngchul Chung
%
% Douglas-Gunn Alternating Direction Implicit approach is based on
% me690-lctr-nts.pdf Eqs. 3.21, 3.22, 3.23a and 3.23b. Implicit steps are
% performed using the Thomson algorithm
% (https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm)
% ***********************************************************************

format long
format compact

if P.calcModeOverlaps && ~numel(P.modes)
  error('You have requested that mode overlaps be calculated, but you have not yet calculated any modes using findModes().');
end
if isempty(P.n.n)
  error('You have to initialize the refractive index object, for example by running initializeRIfromFunction().')
end
if isempty(P.E.field)
  error('You have to initialize the electric field object, either by passing a mode object into it or by running initializeEfromFunction().')
end

if (size(P.n.n,3) > 1) % A 3D array or 3D function has been provided
  if P.twistRate
    error('You cannot specify both twisting and a 3D refractive index profile');
  end
  if P.taperScaling ~= 1
    error('You cannot specify both tapering and a 3D refractive index profile');
  end
end

if P.xSymmetry ~= 0 && ~isinf(P.bendingRoC) && sind(P.bendDirection) || ...
   P.ySymmetry ~= 0 && ~isinf(P.bendingRoC) && cosd(P.bendDirection)
  error('The specified bending direction is inconsistent with the symmetry assumption');
end

if P.saveVideo
  if isempty(P.videoHandle)
    P.videoHandle = VideoWriter([P.name '.avi']);  % If videoHandle is empty, we create it
    P.videoHandle.FrameRate = 5;
    open(P.videoHandle);
  end
end

%% Check for GPU compatibility if needed
if P.useGPU
  if ~ispc
    error('GPU computing is only supported on Windows');
  end
  v = ver;
  if ~any(strcmp({v.Name},'Parallel Computing Toolbox'))
    error('You must have the Parallel Computing Toolbox installed to use GPU acceleration');
  end
  try
    GPUDev = gpuDevice;
    if(str2double(GPUDev.ComputeCapability) < 3)
      error('Your GPU is too old (CUDA compute capability < 3.0)')
    end
  catch
    error('No supported NVIDIA GPU found, or its driver is too old');
  end
end

%% Initialization of space and frequency grids
k_0 = 2*pi/P.lambda;  % [m^-1] Wavenumber

dx = P.dx;
dy = P.dy;
Nx = P.Nx;
Ny = P.Ny;
Lx = P.Lx;
Ly = P.Ly;

if P.calcModeOverlaps
  if P.modes(1).xSymmetry ~= P.xSymmetry || P.modes(1).ySymmetry ~= P.ySymmetry
    error('The pre-calculated mode fields were calculated with different symmetry assumptions. Please recalculate using findModes().');
  end
  if P.modes(1).Lx ~= Lx || P.modes(1).Ly ~= Ly || size(P.modes(1).field,1) ~= Nx || size(P.modes(1).field,2) ~= Ny
    error('The pre-calculated mode fields do not match the simulation Lx, Ly, Nx or Ny. Please recalculate using findModes().');
  end
end

x = P.x;
y = P.y;
[X,Y] = ndgrid(x,y);

% Check if the MATLAB version is one that has the mex/imagesc bug and, if
% so, get ready to do downsampling
hasMexImagescBug = verLessThan('matlab','9.8') && ~verLessThan('matlab','9.5');
if hasMexImagescBug && (Nx > 500 || Ny > 500) && ~P.disableDownsamplingWarning
  warning('MATLAB releases 2018b, 2019a and 2019b have a low-level bug that causes hard crashes when using mex functions and imagesc() with high resolutions (greater than roughly 500x500). For this reason, the plots will here be downsampled to <= 500x500. Avoid this by updating MATLAB to a version in which Mathworks have fixed the bug (R2020a or newer). You can disable this warning by setting P.disableDownsamplingWarning = true.');
end

if hasMexImagescBug && Nx > 500
  ix_plot = round(linspace(1,Nx,500));
else
  ix_plot = 1:Nx;
end
x_plot = x(ix_plot);

if hasMexImagescBug && Ny > 500
  iy_plot = round(linspace(1,Ny,500));
else
  iy_plot = 1:Ny;
end
y_plot = y(iy_plot);

%% Store inputs, if first segment
if ~P.priorData
  if ~isa(P.E,'function_handle')
    powerFraction_Einput = 1/(1 + (P.E.xSymmetry ~= 0))/(1 + (P.E.ySymmetry ~= 0)); % How large a fraction of the total power the input field contains
    P.E.field = P.E.field/sqrt(sum(abs(P.E.field(:)).^2)/powerFraction_Einput);
  end
end

%% Beam initialization
powerFraction = 1/(1 + (P.xSymmetry ~= 0))/(1 + (P.ySymmetry ~= 0)); % How large a fraction of the total power we are simulating

% Interpolate source E field to new grid
[Nx_source,Ny_source] = size(P.E.field);
dx_source = P.E.Lx/Nx_source;
dy_source = P.E.Ly/Ny_source;
x_source = getGridArray(Nx_source,dx_source,P.E.ySymmetry);
y_source = getGridArray(Ny_source,dy_source,P.E.xSymmetry);
[x_source,y_source,E_source] = calcFullField(x_source,y_source,P.E.field);
E = interpn(x_source,y_source,E_source,x,y.','linear',0);
E = E*sqrt(sum(abs(E_source(:)).^2)/sum(abs(E(:)).^2/powerFraction));

if P.ySymmetry == 2; E(X == 0) = 0; end
if P.xSymmetry == 2; E(Y == 0) = 0; end

E = complex(single(E)); % Force to be complex single

%% Refractive index initialization
dz_n = 0; % dz value for the refractive index profile, only used for 3D arrays

% Interpolate input n to the simulation xy grid
[Nx_source,Ny_source,Nz_source] = size(P.n.n);
dx_source = P.n.Lx/Nx_source;
dy_source = P.n.Ly/Ny_source;
x_source = getGridArray(Nx_source,dx_source,P.n.ySymmetry);
y_source = getGridArray(Ny_source,dy_source,P.n.xSymmetry);
[x_source,y_source,n_source] = calcFullRI(x_source,y_source,P.n.n);
if Nz_source == 1 % If P.n.n is 2D, 
  n = interpn(x_source,y_source,n_source,x,y.','linear',P.n_background);
  n_slice = n;
else % Otherwise it's 3D. We trim away any unnecessary repeated outer yz or xz slices that only contain n_background.
  if Nz_source > 1
    dz_n = P.Lz/(Nz_source-1);
  else
    dz_n = 0;
  end
  z_source = dz_n*(0:Nz_source-1);
  n = interpn(x_source,y_source,z_source,n_source,x,y.',z_source,'linear',P.n_background);
  n_slice = n(:,:,1);
  n = trimRI(n,P.n_background);
end

% If RI is 3D, plot it volumetrically
[Nx_n,Ny_n,Nz_n] = size(n);
if Nz_n > 1
  x_n = getGridArray(Nx_n,dx,P.ySymmetry);
  y_n = getGridArray(Ny_n,dy,P.xSymmetry);
  z_n = dz_n*(0:Nz_n-1);
  plotVolumetric(201,x_n,y_n,z_n,real(n),'BPM-Matlab_RI');
  title('Real part of refractive index');xlabel('x [m]');ylabel('y [m]');zlabel('z [m]');
end

%% Calculate z step size and positions
Nz = P.Nz; % Number of z steps in this segment
dz = P.dz;

if ~P.disableStepsizeWarning
  max_a = 5;
  max_d = 15;
  dz_max1 = max_a*4*dx^2*k_0*P.n_0;
  dz_max2 = max_a*4*dy^2*k_0*P.n_0;
  dz_max3 = max_d*2*P.n_0/k_0;
  if dz > min([dz_max1 dz_max2 dz_max3])
    warning('z step size is high (> %.1e m), which may introduce numerical artifacts. Alleviate this by decreasing dz_target, decreasing Nx_main and/or Ny_main, or increasing Lx_main and/or Ly_main. You can disable this warning by setting P.disableStepsizeWarning = true.',min([dz_max1 dz_max2 dz_max3]));
  end
end

zUpdateIdxs = round((1:P.updates)/P.updates*Nz); % Update indices relative to start of this segment
if P.priorData
  P.z = [P.z dz*zUpdateIdxs + P.z(end)];
else
  P.z = [0 dz*zUpdateIdxs];
end

%% Calculate proportionality factors for use in the mex function
ax = dz/(4i*dx^2*k_0*P.n_0);
ay = dz/(4i*dy^2*k_0*P.n_0);
d = -dz*k_0;

%% Calculate the edge absorber multiplier
xEdge = P.Lx_main*(1 + (P.ySymmetry ~= 0))/2;
yEdge = P.Ly_main*(1 + (P.xSymmetry ~= 0))/2;
multiplier = single(exp(-dz*max(0,max(abs(Y) - yEdge,abs(X) - xEdge)).^2*P.alpha)); % Is real

%% Figure initialization
h_f = figure(P.figNum);clf reset;
if strcmp(h_f.WindowStyle,'normal') 
  h_f.WindowState = 'maximized';
end

xlims = ([-1 1] + (P.ySymmetry ~= 0))*Lx/(2*P.plotZoom);
ylims = ([-1 1] + (P.xSymmetry ~= 0))*Ly/(2*P.plotZoom);

if P.xSymmetry ~= 0 && P.ySymmetry ~= 0
  redline_x = [0 P.Lx_main P.Lx_main];
  redline_y = [P.Ly_main P.Ly_main 0];
elseif P.xSymmetry ~= 0
  redline_x = [-P.Lx_main -P.Lx_main P.Lx_main P.Lx_main]/2;
  redline_y = [0 P.Ly_main P.Ly_main 0];
elseif P.ySymmetry ~= 0
  redline_x = [0 P.Lx_main P.Lx_main 0];
  redline_y = [-P.Ly_main -P.Ly_main P.Ly_main P.Ly_main]/2;
else
  redline_x = [-P.Lx_main P.Lx_main P.Lx_main -P.Lx_main -P.Lx_main]/2;
  redline_y = [P.Ly_main P.Ly_main -P.Ly_main -P.Ly_main P.Ly_main]/2;
end

h_axis1 = subplot(2,2,1);
h_im1 = imagesc(x_plot,y_plot,real(n_slice(ix_plot,iy_plot)).');
axis xy
axis equal
xlim(xlims);
ylim(ylims);
colorbar;
setColormap(gca,P.nColormap);
if P.n_colorlimits(2) > P.n_colorlimits(1) % Default is [0 0], so false
  h_axis1.CLim = P.n_colorlimits;
end
xlabel('x [m]');
ylabel('y [m]');
title('Real part of refractive index');

if P.xSymmetry ~= 0
  y0idx = 1;
else
  y0idx = round((Ny-1)/2+1);
end
if P.ySymmetry ~= 0
  x0idx = 1;
else
  x0idx = round((Nx-1)/2+1);
end
if P.priorData
  P.powers = [P.powers NaN(1,P.updates)];
  P.xzSlice{end+1} = NaN(Nx,P.updates);
  P.yzSlice{end+1} = NaN(Ny,P.updates);
else
  P.powers = NaN(1,P.updates+1);
  P.powers(1) = 1;
  P.xzSlice = {NaN(Nx,P.updates+1)};
  P.yzSlice = {NaN(Ny,P.updates+1)};
  P.xzSlice{1}(:,1) = E(:,y0idx);
  P.yzSlice{1}(:,1) = E(x0idx,:);
end
if P.storeE3D
  if P.priorData
    P.E3D{end+1} = complex(NaN(Nx,Ny,P.updates,'single'));
  else
    P.E3D{1} = complex(NaN(Nx,Ny,P.updates+1,'single'));
    P.E3D{1}(:,:,1) = E;
  end
end

h_ax2 = subplot(2,2,2);
h_plot2 = plot(P.z,P.powers,'linewidth',2);
xlim([0 P.z(end)]);
ylim([0 1.1]);
xlabel('Propagation distance [m]');
ylabel('Relative power remaining');
grid on; grid minor;

h_axis3a = subplot(2,2,3);
hold on;
box on;
h_im3a = imagesc(x_plot,y_plot,abs(E(ix_plot,iy_plot).').^2);
axis xy;
axis equal;
xlim(xlims);
ylim(ylims);
colorbar;
xlabel('x [m]');
ylabel('y [m]');
title('Intensity [W/m^2]');
line(redline_x,redline_y,'color','r','linestyle','--');
setColormap(gca,P.intensityColormap);
if P.plotEmax ~= 0
  caxis('manual');
  caxis([0 P.plotEmax*max(abs(E(:)).^2)]);
else
  caxis('auto');
end

h_axis3b = subplot(2,2,4);
hold on;
box on;
maxE0 = abs(max(E(:)));
h_im3b = imagesc(x_plot,y_plot,angle(E(ix_plot,iy_plot).'));
h_im3b.AlphaData = max(0,(1+log10(abs(E(ix_plot,iy_plot).'/maxE0).^2)/3));  %Logarithmic transparency in displaying phase outside cores
h_axis3b.Color = 0.7*[1 1 1];  % To set the color corresponding to phase outside the cores where there is no field at all
axis xy;
axis equal;
xlim(xlims);
ylim(ylims);
colorbar;
caxis([-pi pi]);
line(redline_x,redline_y,'color','r','linestyle','--');
xlabel('x [m]');
ylabel('y [m]');
title('Phase [rad]');
setColormap(gca,P.phaseColormap);

if ~verLessThan('matlab','9.5')
  sgtitle(P.figTitle,'FontSize',15,'FontWeight','bold');
end
drawnow;

if P.saveVideo
  frame = getframe(h_f);  %Get the frames
  writeVideo(P.videoHandle,frame);  %Stitch the frames to form a video and save
end

if P.calcModeOverlaps % Mode overlap figure
  nModes = length(P.modes);
  if P.priorData
    P.modeOverlaps = [P.modeOverlaps NaN(nModes,P.updates)];
  else
    P.modeOverlaps = NaN(nModes,P.updates+1);
    for iMode = 1:nModes
      P.modeOverlaps(iMode,1) = abs(sum(E(:).*conj(P.modes(iMode).field(:)))).^2/powerFraction; % The overlap integral calculation
    end
  end
  
  figure(P.figNum+1);clf reset;
  h_ax = axes;
  h_overlapplot = semilogy(P.z,P.modeOverlaps,'linewidth',2);
  xlim([0 P.z(end)]);
  ylim([1e-4 2]);
  xlabel('Propagation distance [m]');
  ylabel('Mode overlaps');
  legend(P.modes.label,'location','eastoutside','FontSize',6);
  grid on; grid minor;
  h_ax.LineStyleOrder = {'-','--',':','-.'};
end

%% Load variables into a parameters struct and start looping, one iteration per update
mexParameters = struct('dx',single(dx),'dy',single(dy),'dz',single(dz),'taperPerStep',single((1-P.taperScaling)/Nz),'twistPerStep',single(P.twistRate*P.Lz/Nz),...
  'multiplier',single(multiplier),'n_mat',complex(single(n)),'dz_n',single(dz_n),'d',single(d),'n_0',single(P.n_0),...
  'ax',single(ax),'ay',single(ay),'useAllCPUs',P.useAllCPUs,'RoC',single(P.bendingRoC),'rho_e',single(P.rho_e),'bendDirection',single(P.bendDirection),...
  'inputPrecisePower',P.powers(end-length(zUpdateIdxs))*powerFraction,'xSymmetry',uint8(P.xSymmetry),'ySymmetry',uint8(P.ySymmetry));

% fprintf("dz = %.2e, ax = %.2f i, ay = %.2f i, d = %.2f\n",dz,ax/1i,ay/1i,d);
mexParameters.iz_start = int32(0); % z index of the first z position to step from for the first call to FDBPMpropagator, in C indexing (starting from 0)
mexParameters.iz_end = int32(zUpdateIdxs(1)); % last z index to step into for the first call to FDBPMpropagator, in C indexing (starting from 0)
timeInMex = 0;
tic;
for updidx = 1:length(zUpdateIdxs)
  if updidx > 1
    mexParameters.iz_start = int32(zUpdateIdxs(updidx-1));
    mexParameters.iz_end   = int32(zUpdateIdxs(updidx));
  end
  checkMexInputs(E,mexParameters);
  beforeMex = toc;
  if P.useGPU
    [E,n_slice,precisePower] = FDBPMpropagator_CUDA(E,mexParameters);
  else
    [E,n_slice,precisePower] = FDBPMpropagator(E,mexParameters);
  end
  timeInMex = timeInMex + toc - beforeMex;
  if P.storeE3D
    P.E3D{end}(:,:,end-length(zUpdateIdxs)+updidx) = E;
  end
  
  %% Update figure contents
  h_im1.CData = real(n_slice(ix_plot,iy_plot)).'; % Refractive index at this update
  h_im3a.CData = abs(E(ix_plot,iy_plot).').^2; % Intensity at this update
  h_im3b.CData = angle(E(ix_plot,iy_plot).'); % Phase at this update
  h_im3b.AlphaData = max(0,(1+log10(abs(E(ix_plot,iy_plot).'/max(abs(E(:)))).^2)/3));  %Logarithmic transparency in displaying phase outside cores
  if updidx == 1
    caxis(h_axis1,'auto'); % To refresh the numbers on the color bar
  end
  if P.plotEmax == 0
    caxis(h_axis3a,'auto');
  end
  if P.n_colorlimits(2) > P.n_colorlimits(1) % Default is [0 0], so false
    h_axis1.CLim = P.n_colorlimits;
  end

  mexParameters.inputPrecisePower = precisePower;
  P.powers(end-length(zUpdateIdxs)+updidx) = precisePower/powerFraction;
  h_plot2.YData = P.powers;
  h_ax2.YLim = [0 1.1*max(P.powers)];

  P.xzSlice{end}(:,end-length(zUpdateIdxs)+updidx) = E(:,y0idx);
  P.yzSlice{end}(:,end-length(zUpdateIdxs)+updidx) = E(x0idx,:);

  if P.calcModeOverlaps
    for iMode = 1:nModes
      P.modeOverlaps(iMode,end-length(zUpdateIdxs)+updidx) = abs(sum(E(:).*conj(P.modes(iMode).field(:)))).^2/powerFraction; % The overlap integral calculation
      h_overlapplot(iMode).YData = P.modeOverlaps(iMode,:);
    end
  end
  drawnow;
  
  if P.saveVideo
    frame = getframe(h_f); 
    writeVideo(P.videoHandle,frame); 
  end
end
totalTime = toc;
fprintf('Segment done. Elapsed time is %.1f s. %.0f%% of the time was spent updating plots.\n',totalTime,100*(totalTime - timeInMex)/totalTime);
if ~P.disablePlotTimeWarning && timeInMex < 0.5*totalTime && totalTime > 20
  if P.saveVideo
    warning('More than half of the execution time is spent doing plot updates. You can speed up this simulation by decreasing the number of updates or setting P.saveVideo = false. You can disable this warning by setting P.disablePlotTimeWarning = true.');
  else
    warning('More than half of the execution time is spent doing plot updates. You can speed up this simulation by decreasing the number of updates. You can disable this warning by setting P.disablePlotTimeWarning = true.');
  end
end

%% If storing the 3D E data, plot it volumetrically
if P.storeE3D
  if numel(P.E3D) == 1
    z = P.z;
  else
    z = P.z(end-P.updates+1:end);
  end
  if numel(z) > 1
    plotVolumetric(201 + numel(P.E3D),x,y,z,abs(P.E3D{end}).^2,'BPM-Matlab_I');
  end
end

%% Calculate and plot the far field of the final E
figure(P.figNum+2);clf reset;
if hasMexImagescBug
  N_FFhalf = 250; % Points to have at negative theta_x and theta_y in the far field
else
  N_FFhalf = 1000; % Points to have at negative theta_x and theta_y in the far field
end
N_FF = 2*N_FFhalf + 1; % Total number of points in theta_x and theta_y
theta_max = 30; % [deg] max angle

theta_x = linspace(-theta_max,theta_max,N_FF);
kx = 2*pi/P.lambda*theta_x/180*pi;
dx_FF = 2*pi/(kx(2)-kx(1))/N_FF;
x_FF = (-(N_FF-1)/2:(N_FF-1)/2)*dx_FF;

theta_y = linspace(-theta_max,theta_max,N_FF);
ky = 2*pi/P.lambda*theta_y/180*pi;
dy_FF = 2*pi/(ky(2)-ky(1))/N_FF;
y_FF = (-(N_FF-1)/2:(N_FF-1)/2)*dy_FF;

[x_full,y_full,E_full] = calcFullField(x,y,E);

E_interp = interpn(x_full,y_full.',E_full,x_FF,y_FF.','linear',0);

E_FF = fftshift(fft2(ifftshift(conj(E_interp))));

subplot(1,2,1);
imagesc(theta_x,theta_y,abs(E_FF.').^2);
axis xy equal tight;
Theta_x = 4*std(theta_x,sum(abs(E_FF).^2,2));
Theta_y = 4*std(theta_y,sum(abs(E_FF).^2,1));
if ~verLessThan('matlab','9.5')
  sgtitle('Far field in air, in paraxial approximation');
end
title({'Intensity','Divergence 4\sigma full-angles:',['\Theta_x = ' num2str(Theta_x,3) ' deg, \Theta_y = ' num2str(Theta_y,3) ' deg']});
xlabel('\theta_x [deg]');
ylabel('\theta_y [deg]');
colorbar;
setColormap(gca,P.intensityColormap);

subplot(1,2,2);
imagesc(theta_x,theta_y,angle(E_FF.'),'AlphaData',max(0,(1+log10(abs(E_FF.'/max(abs(E_FF(:)))).^2)/3)));
set(gca,'Color',0.7*[1 1 1]);
axis xy equal tight;
title('Phase');
xlabel('\theta_x [deg]');
ylabel('\theta_y [deg]');
colorbar;
setColormap(gca,P.phaseColormap);

%% Store the final E field and n as the new inputs
P.E.Lx = P.Lx;
P.E.Ly = P.Ly;
P.E.field = E;
P.E.xSymmetry = P.xSymmetry;
P.E.ySymmetry = P.ySymmetry;

P.n.Lx = P.Lx;
P.n.Ly = P.Ly;
P.n.n = n_slice;
P.n.xSymmetry = P.xSymmetry;
P.n.ySymmetry = P.ySymmetry;

P.priorData = true;
end

function checkMexInputs(E,P)
typename = 'single';
assert(all(isfinite(E(:))));
assert(~isreal(E));
assert(isa(P.dx,typename));
assert(isreal(P.dx));
assert(isa(P.dy,typename));
assert(isreal(P.dy));
assert(isa(P.dz_n,typename));
assert(isreal(P.dz_n));
assert(isa(P.taperPerStep,typename));
assert(isreal(P.taperPerStep));
assert(isa(P.twistPerStep,typename));
assert(isreal(P.twistPerStep));
assert(isa(P.n_mat,typename));
assert(~isreal(P.n_mat));
assert(isa(P.d,typename));
assert(isreal(P.d));
assert(isa(P.n_0,typename));
assert(isreal(P.n_0));
assert(isa(P.ax,typename));
assert(~isreal(P.ax));
assert(isa(P.ay,typename));
assert(~isreal(P.ay));
assert(isa(P.RoC,typename));
assert(isreal(P.RoC));
assert(isa(P.rho_e,typename));
assert(isreal(P.rho_e));
assert(isa(P.bendDirection,typename));
assert(isreal(P.bendDirection));
assert(isa(P.xSymmetry,'uint8'));
assert(isreal(P.xSymmetry));
assert(isa(P.ySymmetry,'uint8'));
assert(isreal(P.ySymmetry));
end