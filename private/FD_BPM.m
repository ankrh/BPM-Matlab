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

if ~(ispc || ismac)
  error('Error: BPM-Matlab doesn''t support Linux systems yet');
end

format long
format compact

k_0 = 2*pi/P.lambda;  % [m^-1] Wavenumber

if isfield(P,'n_cladding')
  error('Error: n_cladding has been renamed n_background');
end
if isfield(P,'nFunc') && isfield(P,'shapes')
  error('Error: You must specify exactly one of the fields "shapes" and "nFunc"');
end
if isfield(P,'shapes') && isempty(P.shapes)
  P.shapes = [0 0 0 1 0];
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
if ~isfield(P,'rho_e')
  P.rho_e = 0.22;
end
if ~isfield(P,'bendingRoC')
  P.bendingRoC = Inf;
end
if ~isfield(P,'bendDirection')
  P.bendDirection = 0;
end
if isfield(P,'shapes') && size(P.shapes,2) == 5
  if any(P.shapes(:,4) == 4) || any(P.shapes(:,4) == 5)
    error('Since you have a GRIN lens, you must define the gradient constant g in the shapes array');
  else
    P.shapes(:,6) = NaN;
  end
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
  if isfield(P,'shapes')
    P.originalShapesInput = P.shapes;
  end
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
  [X_source,Y_source] = ndgrid(x_source,y_source);
  E = interpn(X_source,Y_source,P.E.field,X,Y,'linear',0);
  E = E*sqrt(sum(abs(P.E.field(:)).^2)/sum(abs(E(:)).^2));
end

E = complex(single(E)); % Force to be complex single

if ~priorData
  P.Einitial = E;
end

%% Refractiv index initialization
% If shapes is defined, calculate refractive index profile
if isfield(P,'shapes')
  n = P.n_background*ones(Nx,Ny,'single'); % May be complex
  for iShape = 1:size(P.shapes,1)
    switch P.shapes(iShape,4)
      case 1
        n((X - P.shapes(iShape,1)).^2 + (Y - P.shapes(iShape,2)).^2 < P.shapes(iShape,3)^2) = P.shapes(iShape,5);
      case 2
        delta = max(dx,dy);
        r_diff = sqrt((X - P.shapes(iShape,1)).^2 + (Y - P.shapes(iShape,2)).^2) - P.shapes(iShape,3) + delta/2;
        n(r_diff < 0) = P.shapes(iShape,5);
        n(r_diff >= 0 & r_diff < delta) = r_diff(r_diff >= 0 & r_diff < delta)/delta*(P.n_background - P.shapes(iShape,5)) + P.shapes(iShape,5);
        clear r_diff
      case 3
        r_ratio_sqr = ((X - P.shapes(iShape,1)).^2 + (Y - P.shapes(iShape,2)).^2)/P.shapes(iShape,3)^2;
        n(r_ratio_sqr < 1) = r_ratio_sqr(r_ratio_sqr < 1)*(P.n_background - P.shapes(iShape,5)) + P.shapes(iShape,5);
        clear r_ratio_sqr
      case 4
        r = sqrt((X - P.shapes(iShape,1)).^2 + (Y - P.shapes(iShape,2)).^2);
        n(r < P.shapes(iShape,3)) = P.shapes(iShape,5)*sech(P.shapes(iShape,6)*r(r < P.shapes(iShape,3)));
        clear r
      case 5
        r = abs(Y - P.shapes(iShape,2));
        n(r < P.shapes(iShape,3)) = P.shapes(iShape,5)*sech(P.shapes(iShape,6)*r(r < P.shapes(iShape,3)));
        clear r
    end
  end
elseif isa(P.n,'function_handle') % Otherwise if P.n is a function
  n = single(P.n(X,Y,P.n_background,P.nParameters));
else % Otherwise P.n must be a struct
  [Nx_source,Ny_source] = size(P.n.n);
  dx_source = P.n.Lx/Nx_source;
  dy_source = P.n.Ly/Ny_source;
  x_source = dx_source*(-(Nx_source-1)/2:(Nx_source-1)/2);
  y_source = dy_source*(-(Ny_source-1)/2:(Ny_source-1)/2);
  [X_source,Y_source] = ndgrid(x_source,y_source);
  n = single(interpn(X_source,Y_source,double(P.n.n),X,Y,'linear',P.n_background));
end

%% Calculate z step size and positions
Nz = max(P.updates,round(P.Lz/P.dz_target)); % Number of z steps in this segment
dz = P.Lz/Nz;

if ~P.disableStepsizeWarning
  max_a = 5;
  max_d = 2.5;
  dz_max1 = max_a*4*dx^2*k_0*P.n_0;
  dz_max2 = max_a*4*dy^2*k_0*P.n_0;
  dz_max3 = max_d*2*P.n_0/k_0;
  if isfield(P,'shapes') && (any(P.shapes(:,4) == 1) || any(P.shapes(:,4) == 2))
    if dz > min([dz_max1 dz_max2 dz_max3])
      warning('z step size is high (> %.1e m), which may introduce numerical artifacts. You can disable this warning by setting P.disableStepsizeWarning = true.',min([dz_max1 dz_max2 dz_max3]));
    end
  else
    if dz > min([dz_max1 dz_max2])
      warning('z step size is high (> %.1e m), which may introduce numerical artifacts. You can disable this warning by setting P.disableStepsizeWarning = true.',min(dz_max1,dz_max2));
    end
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

%% Calculate the static and dynamic multipliers
multiplier = single(exp(-dz*max(0,max(abs(Y) - P.Ly_main/2,abs(X) - P.Lx_main/2)).^2*P.alpha)); % Is real

%% Figure initialization
h_f = figure(P.figNum);clf;
h_f.WindowState = 'maximized';

h_axis1 = subplot(2,2,1);
if P.downsampleImages
  h_im1 = imagesc(x_plot,y_plot,real(n(ix_plot,iy_plot)).');
else
  h_im1 = imagesc(x,y,real(n).');
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

subplot(2,2,2);
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
if isfield(P,'shapes') && (P.taperScaling == 1 && P.twistRate == 0)
  for iShape = 1:size(P.shapes,1)
    if P.shapes(iShape,4) <=3
      line(P.shapes(iShape,1) + P.shapes(iShape,3)*cos(linspace(0,2*pi,100)),P.shapes(iShape,2) + P.shapes(iShape,3)*sin(linspace(0,2*pi,100)),'color','w','linestyle','--');
    end
  end
end
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
if isfield(P,'shapes') && (P.taperScaling == 1 && P.twistRate == 0)
  for iShape = 1:size(P.shapes,1)
    if P.shapes(iShape,4) <=3
      line(P.shapes(iShape,1) + P.shapes(iShape,3)*cos(linspace(0,2*pi,100)),P.shapes(iShape,2) + P.shapes(iShape,3)*sin(linspace(0,2*pi,100)),'color','w','linestyle','--');
    end
  end
end
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
  
  figure(P.figNum+1);clf;
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
% mexParameters = struct('dx',single(dx),'dy',single(dy),'taperPerStep',single((1-P.taperScaling)/Nz),'twistPerStep',single(P.twistRate*P.Lz/Nz),...
%   'shapes',single(real(P.shapes)),'shapeAbsorptions',single(shapeAbsorptions),'n_background',single(real(P.n_background)),'backgroundAbsorption',single(backgroundAbsorption),'multiplier',complex(single(multiplier)),'d',single(d),'n_0',single(P.n_0),...
%   'ax',single(ax),'ay',single(ay),'useAllCPUs',P.useAllCPUs,'RoC',single(P.bendingRoC),'rho_e',single(P.rho_e),'bendDirection',single(P.bendDirection),...
%   'inputPrecisePower',P.powers(end-length(zUpdateIdxs)));
mexParameters = struct('dx',single(dx),'dy',single(dy),'taperPerStep',single((1-P.taperScaling)/Nz),'twistPerStep',single(P.twistRate*P.Lz/Nz),...
  'multiplier',single(multiplier),'n_mat',complex(single(n)),'d',single(d),'n_0',single(P.n_0),...
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
    [E,n,precisePower] = FDBPMpropagator_CUDA(E,mexParameters);
  else
    [E,n,precisePower] = FDBPMpropagator(E,mexParameters);
  end

  %% Update figure contents
  if P.downsampleImages
    h_im1.CData = real(n(ix_plot,iy_plot)).'; % Refractive index at this update
    h_im3a.CData = abs(E(ix_plot,iy_plot).').^2; % Intensity at this update
    h_im3b.CData = angle(E(ix_plot,iy_plot).'); % Phase at this update
    h_im3b.AlphaData = max(0,(1+log10(abs(E(ix_plot,iy_plot).'/max(abs(E(:)))).^2)/3));  %Logarithmic transparency in displaying phase outside cores
  else
    h_im1.CData = real(n).'; % Refractive index at this update
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

%% Calculate the output shapes and store the final E field as the new input field
if isfield(P,'shapes')
  shapesFinal = P.shapes;
  shapesFinal(:,1) = P.taperScaling*(cos(P.twistRate*P.Lz)*P.shapes(:,1) - sin(P.twistRate*P.Lz)*P.shapes(:,2));
  shapesFinal(:,2) = P.taperScaling*(sin(P.twistRate*P.Lz)*P.shapes(:,1) + cos(P.twistRate*P.Lz)*P.shapes(:,2));
  shapesFinal(:,3) = P.taperScaling*P.shapes(:,3);
  P.shapes = shapesFinal;
end

P.E = struct('field',E,'Lx',Lx,'Ly',Ly);
P.n = struct('n',n,'Lx',Lx,'Ly',Ly);

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
