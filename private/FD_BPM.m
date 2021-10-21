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

k_0 = 2*pi/P.lambda;  % [m^-1] Wavenumber

if isfield(P,'n_cladding')
  error('Error: n_cladding has been renamed n_background');
end
if isfield(P,'shapes')
  error('The P.shapes field has been deprecated. Use the P.n field to define the refractive index instead, as shown in the example files.');
end
if ~isfield(P,'n')
  error('You must specify the P.n field');
end
if isa(P.n,'function_handle')
  error('The method of defining a refractive index function has changed. The function handle must now be stored in P.n.func instead of P.n.');
end
if ~isfield(P.n,'n') && ~isfield(P.n,'func')
  error('P.n must contain either the field "func" or "n"');
end
if ~isfield(P.n,'func') && ~(isfield(P.n,'Lx') && isfield(P.n,'Ly'))
  error('You must specify the side lengths Lx and Ly if you provide an array (2D or 3D) for the refractive index');
end
if isfield(P.n,'Nx') || isfield(P.n,'Ny')
  warning('The P.n.Nx and P.n.Ny fields are no longer necessary. 3D RI functions are now simply evaluated with the same xy resolution as the simulation grid.');
end
if isfield(P.n,'func') && nargin(P.n.func) == 5 && ~isfield(P.n,'Nz')
  error('You must specify the refractive index array z resolution P.n.Nz if you provide a 3D refractive index function');
end
if ~isfield(P,'storeE3D')
  P.storeE3D = false;
end
if ~isfield(P,'figNum')
  P.figNum = 1;
end
if ~isfield(P,'figTitle')
  P.figTitle = '';
end
if ~isfield(P,'saveVideo')
  P.saveVideo = false;
end
if ~isfield(P,'finalizeVideo')
  P.finalizeVideo = true;
end
if ~isfield(P,'saveData')
  P.saveData = false;
end
if ~isfield(P,'useGPU')
  P.useGPU = false;
end
if ~isfield(P,'useAllCPUs')
  P.useAllCPUs = false;
end
if ~isfield(P,'downsampleImages')
  P.downsampleImages = false;
end
if ~isfield(P,'displayScaling')
  P.displayScaling = 1;
end
if ~isfield(P,'disableStepsizeWarning')
  P.disableStepsizeWarning = false;
end
if ~isfield(P,'Eparameters')
  P.Eparameters = {};
end
if ~isfield(P,'nParameters')
  P.nParameters = {};
end
if ~isfield(P,'taperScaling')
  P.taperScaling = 1;
end
if ~isfield(P,'twistRate')
  P.twistRate = 0;
end
if (isfield(P.n,'n') && size(P.n.n,3) > 1) || (isfield(P.n,'func') && nargin(P.n.func) == 5) % A 3D array or 3D function has been provided
  if P.twistRate
    error('You cannot specify both twisting and a 3D refractive index profile');
  end
  if P.taperScaling ~= 1
    error('You cannot specify both tapering and a 3D refractive index profile');
  end
end
if ~isfield(P,'rho_e')
  P.rho_e = 0.22;
end
if ~isfield(P,'bendingRoC')
  P.bendingRoC = Inf;
end
if ~isfield(P,'bendDirection')
  P.bendDirection = 0;
end
if P.saveVideo && ~isfield(P,'videoName')
  P.videoName = [P.name '.avi'];
end
if ~isfield(P,'Intensity_colormap')
  P.Intensity_colormap = 1;
end
if ~isfield(P,'Phase_colormap')
  P.Phase_colormap = 2;
end
if ~isfield(P,'n_colormap')
  P.n_colormap = 3;
end

if P.saveVideo
  if isfield(P,'videoHandle')
    video = P.videoHandle;
  else
    video = VideoWriter(P.videoName);  % If videoHandle is not passed from the model file, we create one
    video.FrameRate = 5;
    open(video);
  end
end

typename = 'single';

if ~isfield(P,'calcModeOverlaps')
  P.calcModeOverlaps = false; 
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
dx = P.Lx_main/P.Nx_main;
dy = P.Ly_main/P.Ny_main;

targetLx = P.padfactor*P.Lx_main;
targetLy = P.padfactor*P.Ly_main;

Nx = round(targetLx/dx);
if rem(Nx,2) ~= rem(P.Nx_main,2)
  Nx = Nx + 1; % Ensure that if Nx_main was set odd (to have a x slice at the center), Nx will also be odd
end
Ny = round(targetLy/dy);
if rem(Ny,2) ~= rem(P.Ny_main,2)
  Ny = Ny + 1; % Ensure that if Ny_main was set odd (to have a y slice at the center), Ny will also be odd
end
Lx = Nx*dx;
Ly = Ny*dy;

if P.calcModeOverlaps && (P.modes(1).Lx ~= Lx || P.modes(1).Ly ~= Ly || size(P.modes(1).field,1) ~= Nx || size(P.modes(1).field,2) ~= Ny)
  error('The pre-calculated mode fields do not match the simulation Lx, Ly, Nx or Ny');
end

x = dx*(-(Nx-1)/2:(Nx-1)/2);
y = dy*(-(Ny-1)/2:(Ny-1)/2);
[X,Y] = ndgrid(x,y);

if P.downsampleImages
  if Nx>500
    ix_plot = round(linspace(1,Nx,500));
  else
    ix_plot = 1:Nx;
  end
  x_plot = x(ix_plot);

  if Ny>500
    iy_plot = round(linspace(1,Ny,500));
  else
    iy_plot = 1:Ny;
  end
  y_plot = y(iy_plot);
end

%% Store inputs, if first segment
priorData = isfield(P,'originalEinput');
if ~priorData
  if ~isa(P.E,'function_handle')
    P.E.field = P.E.field/sqrt(sum(abs(P.E.field(:)).^2));
  end
  P.originalEinput = P.E;
end

%% Beam initialization
if isa(P.E,'function_handle')
  E = P.E(X,Y,P.Eparameters); % Call function to initialize E field
  E = E/sqrt(sum(abs(E(:)).^2));
else % Interpolate source E field to new grid
  [Nx_source,Ny_source] = size(P.E.field);
  dx_source = P.E.Lx/Nx_source;
  dy_source = P.E.Ly/Ny_source;
  x_source = dx_source*(-(Nx_source-1)/2:(Nx_source-1)/2);
  y_source = dy_source*(-(Ny_source-1)/2:(Ny_source-1)/2);
  E = interpn(x_source,y_source,P.E.field,x,y.','linear',0);
  E = E*sqrt(sum(abs(P.E.field(:)).^2)/sum(abs(E(:)).^2));
end

E = complex(single(E)); % Force to be complex single

if ~priorData
  P.Einitial = E;
end

%% Refractive index initialization
dz_n = 0; % dz value for the refractive index profile, only used for 3D arrays
if isfield(P.n,'func') % If P.n has a function field
  if nargin(P.n.func) == 4 % It's 2D
    n = single(P.n.func(X,Y,P.n_background,P.nParameters));
    n_slice = n;
  else % It's 3D. We evaluate it on the grid and trim the result
    dz_n = P.Lz/(P.n.Nz-1);
    z_n = dz_n*(0:P.n.Nz-1);
    [X_n,Y_n,Z_n] = ndgrid(single(x),single(y),single(z_n));
    n = single(P.n.func(X_n,Y_n,Z_n,P.n_background,P.nParameters));
    clear X_n Y_n Z_n
    n_slice = n(:,:,1);
    n = trimRI(n,P.n_background);
  end
else % Otherwise a P.n.n array must have been specified, so we will interpolate it to the simulation xy grid
  [Nx_source,Ny_source,Nz_source] = size(P.n.n);
  dx_source = P.n.Lx/Nx_source;
  dy_source = P.n.Ly/Ny_source;
  x_source = dx_source*(-(Nx_source-1)/2:(Nx_source-1)/2);
  y_source = dy_source*(-(Ny_source-1)/2:(Ny_source-1)/2);
  if Nz_source == 1 % If P.n.n is 2D, 
    n = interpn(x_source,y_source,P.n.n,x,y.','linear',P.n_background);
    n_slice = n;
  else % Otherwise it's 3D. We trim away any unnecessary repeated outer yz or xz slices that only contain n_background.
    dz_n = P.Lz/(Nz_source-1);
    z_source = dz_n*(0:Nz_source-1);
    n = interpn(x_source,y_source,z_source,P.n.n,x,y.',z_source,'linear',P.n_background);
    n_slice = n(:,:,1);
    n = trimRI(n,P.n_background);
  end
end

% If RI is 3D, plot it volumetrically
[Nx_n,Ny_n,Nz_n] = size(n);
if Nz_n > 1
  x_n = dx*(-(Nx_n-1)/2:(Nx_n-1)/2);
  y_n = dy*(-(Ny_n-1)/2:(Ny_n-1)/2);
  z_n = dz_n*(0:Nz_n-1);
  plotVolumetric(201,x_n,y_n,z_n,real(n),'BPM-Matlab_RI');
  title('Real part of refractive index');xlabel('x [m]');ylabel('y [m]');zlabel('z [m]');
end

%% Calculate z step size and positions
Nz = max(P.updates,round(P.Lz/P.dz_target)); % Number of z steps in this segment
dz = P.Lz/Nz;

if ~P.disableStepsizeWarning
  max_a = 5;
  max_d = 15;
  dz_max1 = max_a*4*dx^2*k_0*P.n_0;
  dz_max2 = max_a*4*dy^2*k_0*P.n_0;
  dz_max3 = max_d*2*P.n_0/k_0;
  if dz > min([dz_max1 dz_max2 dz_max3])
    warning('z step size is high (> %.1e m), which may introduce numerical artifacts. You can disable this warning by setting P.disableStepsizeWarning = true.',min([dz_max1 dz_max2 dz_max3]));
  end
end

zUpdateIdxs = round((1:P.updates)/P.updates*Nz); % Update indices relative to start of this segment
if priorData
  P.z = [P.z dz*zUpdateIdxs + P.z(end)];
else
  P.z = [0 dz*zUpdateIdxs];
end

%% Calculate proportionality factors for use in the mex function
ax = dz/(4i*dx^2*k_0*P.n_0);
ay = dz/(4i*dy^2*k_0*P.n_0);
d = -dz*k_0;

%% Calculate the edge absorber multiplier
multiplier = single(exp(-dz*max(0,max(abs(Y) - P.Ly_main/2,abs(X) - P.Lx_main/2)).^2*P.alpha)); % Is real

%% Figure initialization
h_f = figure(P.figNum);clf reset;
if strcmp(h_f.WindowStyle,'normal') 
  h_f.WindowState = 'maximized';
end

h_axis1 = subplot(2,2,1);
if P.downsampleImages
  h_im1 = imagesc(x_plot,y_plot,real(n_slice(ix_plot,iy_plot)).');
else
  h_im1 = imagesc(x,y,real(n_slice).');
end
axis xy
axis equal
xlim([-1 1]*Lx/(2*P.displayScaling));
ylim([-1 1]*Ly/(2*P.displayScaling));
colorbar;
setColormap(gca,P.n_colormap);
if isfield(P,'n_colorlimits')
  caxis(P.n_colorlimits);
end
xlabel('x [m]');
ylabel('y [m]');
title('Real part of refractive index');

if priorData
  P.powers = [P.powers NaN(1,P.updates)];
  P.xzSlice{end+1} = NaN(Nx,P.updates);
  P.yzSlice{end+1} = NaN(Ny,P.updates);
else
  P.powers = NaN(1,P.updates+1);
  P.powers(1) = 1;
  P.xzSlice = {NaN(Nx,P.updates+1)};
  P.yzSlice = {NaN(Ny,P.updates+1)};
  P.xzSlice{1}(:,1) = E(:,round((Nx-1)/2+1));
  P.yzSlice{1}(:,1) = E(round((Ny-1)/2+1),:);
end
if P.storeE3D
  P.E3D = complex(NaN(Nx,Ny,P.updates+1,'single'));
  P.E3D(:,:,1) = E;
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
if P.downsampleImages
  h_im3a = imagesc(x_plot,y_plot,abs(E(ix_plot,iy_plot).').^2);
else
  h_im3a = imagesc(x,y,abs(E.').^2);
end
axis xy;
axis equal;
xlim([-1 1]*Lx/(2*P.displayScaling));
ylim([-1 1]*Ly/(2*P.displayScaling));
colorbar;
xlabel('x [m]');
ylabel('y [m]');
title('Intensity [W/m^2]');
line([-P.Lx_main P.Lx_main P.Lx_main -P.Lx_main -P.Lx_main]/2,[P.Ly_main P.Ly_main -P.Ly_main -P.Ly_main P.Ly_main]/2,'color','r','linestyle','--');
setColormap(gca,P.Intensity_colormap);
if isfield(P,'plotEmax')
  caxis('manual');
  caxis([0 P.plotEmax*max(abs(E(:)).^2)]);
else
  caxis('auto');
end

h_axis3b = subplot(2,2,4);
hold on;
box on;
maxE0 = abs(max(E(:)));
if P.downsampleImages
  h_im3b = imagesc(x_plot,y_plot,angle(E(ix_plot,iy_plot).'));
  h_im3b.AlphaData = max(0,(1+log10(abs(E(ix_plot,iy_plot).'/maxE0).^2)/3));  %Logarithmic transparency in displaying phase outside cores
else
  h_im3b = imagesc(x,y,angle(E.'));
  h_im3b.AlphaData = max(0,(1+log10(abs(E.'/maxE0).^2)/3));  %Logarithmic transparency in displaying phase outside cores
end
h_axis3b.Color = 0.7*[1 1 1];  % To set the color corresponding to phase outside the cores where there is no field at all
axis xy;
axis equal;
xlim([-1 1]*Lx/(2*P.displayScaling));
ylim([-1 1]*Ly/(2*P.displayScaling));
colorbar;
caxis([-pi pi]);
line([-P.Lx_main P.Lx_main P.Lx_main -P.Lx_main -P.Lx_main]/2,[P.Ly_main P.Ly_main -P.Ly_main -P.Ly_main P.Ly_main]/2,'color','r','linestyle','--');
xlabel('x [m]');
ylabel('y [m]');
title('Phase [rad]');
setColormap(gca,P.Phase_colormap);

sgtitle(P.figTitle,'FontSize',15,'FontWeight','bold');
drawnow;

if P.saveVideo
  frame = getframe(h_f);  %Get the frames
  writeVideo(video,frame);  %Stitch the frames to form a video and save
end

if P.calcModeOverlaps % Mode overlap figure
  nModes = length(P.modes);
  if priorData
    P.modeOverlaps = [P.modeOverlaps NaN(nModes,P.updates)];
  else
    P.modeOverlaps = NaN(nModes,P.updates+1);
    for iMode = 1:nModes
      P.modeOverlaps(iMode,1) = abs(sum(E(:).*conj(P.modes(iMode).field(:)))).^2; % The overlap integral calculation
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

% tic;
%% Load variables into a parameters struct and start looping, one iteration per update
mexParameters = struct('dx',single(dx),'dy',single(dy),'dz',single(dz),'taperPerStep',single((1-P.taperScaling)/Nz),'twistPerStep',single(P.twistRate*P.Lz/Nz),...
  'multiplier',single(multiplier),'n_mat',complex(single(n)),'dz_n',single(dz_n),'d',single(d),'n_0',single(P.n_0),...
  'ax',single(ax),'ay',single(ay),'useAllCPUs',P.useAllCPUs,'RoC',single(P.bendingRoC),'rho_e',single(P.rho_e),'bendDirection',single(P.bendDirection),...
  'inputPrecisePower',P.powers(end-length(zUpdateIdxs)));

% fprintf("dz = %.2e, ax = %.2f i, ay = %.2f i, d = %.2f\n",dz,ax/1i,ay/1i,d);
mexParameters.iz_start = int32(0); % z index of the first z position to step from for the first call to FDBPMpropagator, in C indexing (starting from 0)
mexParameters.iz_end = int32(zUpdateIdxs(1)); % last z index to step into for the first call to FDBPMpropagator, in C indexing (starting from 0)
for updidx = 1:length(zUpdateIdxs)
  if updidx > 1
    mexParameters.iz_start = int32(zUpdateIdxs(updidx-1));
    mexParameters.iz_end   = int32(zUpdateIdxs(updidx));
  end
  checkMexInputs(E,mexParameters,typename);
  if P.useGPU
    [E,n_slice,precisePower] = FDBPMpropagator_CUDA(E,mexParameters);
  else
    [E,n_slice,precisePower] = FDBPMpropagator(E,mexParameters);
  end
  if P.storeE3D
    P.E3D(:,:,updidx + 1) = E;
  end
  
  %% Update figure contents
  if P.downsampleImages
    h_im1.CData = real(n_slice(ix_plot,iy_plot)).'; % Refractive index at this update
    h_im3a.CData = abs(E(ix_plot,iy_plot).').^2; % Intensity at this update
    h_im3b.CData = angle(E(ix_plot,iy_plot).'); % Phase at this update
    h_im3b.AlphaData = max(0,(1+log10(abs(E(ix_plot,iy_plot).'/max(abs(E(:)))).^2)/3));  %Logarithmic transparency in displaying phase outside cores
  else
    h_im1.CData = real(n_slice).'; % Refractive index at this update
    h_im3a.CData = abs(E.').^2; % Intensity at this update
    h_im3b.CData = angle(E.'); % Phase at this update
    h_im3b.AlphaData = max(0,(1+log10(abs(E.'/max(abs(E(:)))).^2)/3));  %Logarithmic transparency in displaying phase outside cores
  end
  if updidx == 1
    caxis(h_axis1,'auto'); % To refresh the numbers on the color bar
  end
  if ~isfield(P,'plotEmax')
    caxis(h_axis3a,'auto');
  end
  
  mexParameters.inputPrecisePower = precisePower;
  P.powers(end-length(zUpdateIdxs)+updidx) = precisePower;
  h_plot2.YData = P.powers;
  h_ax2.YLim = [0 1.1*max(P.powers)];

  P.xzSlice{end}(:,end-length(zUpdateIdxs)+updidx) = E(:,round((Nx-1)/2+1));
  P.yzSlice{end}(:,end-length(zUpdateIdxs)+updidx) = E(round((Ny-1)/2+1),:);

  if P.calcModeOverlaps
    for iMode = 1:nModes
      P.modeOverlaps(iMode,end-length(zUpdateIdxs)+updidx) = abs(sum(E(:).*conj(P.modes(iMode).field(:)))).^2; % The overlap integral calculation
      h_overlapplot(iMode).YData = P.modeOverlaps(iMode,:);
    end
  end
  drawnow;
  
  if P.saveVideo
    frame = getframe(h_f); 
    writeVideo(video,frame); 
  end
end
% toc

if P.saveVideo
  if P.finalizeVideo
  	close(video);
  else
    P.videoHandle = video;
  end
end

%% If storing the 3D E data, plot it volumetrically
if P.storeE3D
  plotVolumetric(202,x,y,P.z,abs(P.E3D).^2,'BPM-Matlab_I');
end

%% Calculate and plot the far field of the final E
figure(P.figNum+2);clf reset;
N_FFhalf = 1000; % Points to have at negative theta_x and theta_y in the far field
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

E_interp = interpn(x,y.',E,x_FF,y_FF.','linear',0);

E_FF = fftshift(fft2(ifftshift(conj(E_interp))));

subplot(1,2,1);
imagesc(theta_x,theta_y,abs(E_FF.').^2);
axis xy equal tight;
Theta_x = 4*std(theta_x,sum(abs(E_FF).^2,2));
Theta_y = 4*std(theta_y,sum(abs(E_FF).^2,1));
sgtitle('Far field in air, in paraxial approximation');
title({'Intensity','Divergence 4\sigma full-angles:',['\Theta_x = ' num2str(Theta_x,3) '°, \Theta_y = ' num2str(Theta_y,3) '°']});
xlabel('\theta_x [°]');
ylabel('\theta_y [°]');
colorbar;
setColormap(gca,P.Intensity_colormap);

subplot(1,2,2);
imagesc(theta_x,theta_y,angle(E_FF.'),'AlphaData',max(0,(1+log10(abs(E_FF.'/max(abs(E_FF(:)))).^2)/3)));
set(gca,'Color',0.7*[1 1 1]);
axis xy equal tight;
title('Phase');
xlabel('\theta_x [°]');
ylabel('\theta_y [°]');
colorbar;
setColormap(gca,P.Phase_colormap);

%% Store the final E field and n as the new inputs
P.E = struct('field',E,'Lx',Lx,'Ly',Ly);
P.n = struct('n',n_slice,'Lx',Lx,'Ly',Ly);

P.x = x;
P.y = y;

% S = load('train');
% sound(S.y.*0.1,S.Fs);
end

function checkMexInputs(E,P,typename)
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
end

function setColormap(gca,colormapType)
    switch colormapType
     case 1
        colormap(gca,GPBGYRcolormap);
     case 2
        colormap(gca,hsv/1.5);
     case 3
        colormap(gca,parula);
     case 4
        colormap(gca,gray);
     case 5
        colormap(gca,cividisColormap);
    end
end

function n = trimRI(n,n_background)
n = single(n);
n_background = single(n_background);
[Nx,Ny,~] = size(n);
xmin = find(any(n ~= n_background,[2 3]),1,'first');
xmax = find(any(n ~= n_background,[2 3]),1,'last');
xtrim = min(xmin-1,Nx - xmax);
ymin = find(any(n ~= n_background,[1 3]),1,'first');
ymax = find(any(n ~= n_background,[1 3]),1,'last');
ytrim = min(ymin-1,Ny - ymax);

n_temp = n(xtrim+1:Nx-xtrim,ytrim+1:Ny-ytrim,:);
n = n_background*ones(size(n_temp) + [2 2 0], class(n_temp));
n(2:end-1, 2:end-1, :) = n_temp;
end