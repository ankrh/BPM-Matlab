function [Estruct,shapes_out,powers,P] = FD_BPM(P)
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

if isempty(P.shapes)
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
if ~isfield(P,'Eparameters')
  P.Eparameters = {};
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
if size(P.shapes,2) == 5
  if any(P.shapes(:,4) == 4) || any(P.shapes(:,4) == 5)
    error('Since you have a GRIN lens, you must define the gradient constant g in the shapes array');
  else
    P.shapes(:,6) = NaN;
  end
end
if ~isfield(P,'videoName')
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

if P.saveVideo && ~isfield(P,'videoHandle')
  video = VideoWriter(P.videoName);  % If videoHandle is not passed from Example.m file, video of only the last segment will be saved
  open(video);
elseif P.saveVideo
  video = P.videoHandle;  %videoHandle is passed from Example.m file
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

%% Calculate the output shapes
shapes_out = P.shapes;
shapes_out(:,1) = P.taperScaling*(cos(P.twistRate*P.Lz)*P.shapes(:,1) - sin(P.twistRate*P.Lz)*P.shapes(:,2));
shapes_out(:,2) = P.taperScaling*(sin(P.twistRate*P.Lz)*P.shapes(:,1) + cos(P.twistRate*P.Lz)*P.shapes(:,2));
shapes_out(:,3) = P.taperScaling*P.shapes(:,3);

%% Beam initialization
if isa(P.E,'function_handle')
  E = P.E(X,Y,P.Eparameters); % Call function to initialize E field
  P.E_0 = E; % E_0 variable is only for powers measurement, when FDBPM is called from multiple segments
else % Interpolate source E field to new grid
  [Nx_Esource,Ny_Esource] = size(P.E.field);
  dx_Esource = P.E.Lx/Nx_Esource;
  dy_Esource = P.E.Ly/Ny_Esource;
  x_Esource = dx_Esource*(-(Nx_Esource-1)/2:(Nx_Esource-1)/2);
  y_Esource = dy_Esource*(-(Ny_Esource-1)/2:(Ny_Esource-1)/2);
  [X_source,Y_source] = ndgrid(x_Esource,y_Esource);
  E = interp2(X_source.',Y_source.',P.E.field.',X.',Y.','linear',0).';
end

E = complex(single(E)); % Force to be complex single precision
if isfield(P,'E_0')
    P_0 = sum(abs(P.E_0(:)).^2);  % For powers w.r.t. the input E field from segment 1
else
    P_0 = sum(abs(E(:)).^2);  % If the input E from segment 1 is missing, use the current E as initial field
end

%% Calculate z step size and positions
Nz = max(P.updates,round(P.Lz/P.dz_target)); % Number of z steps in this segment
dz = P.Lz/Nz;

zUpdateIdxs = round((1:P.updates)/P.updates*Nz); % Update indices relative to start of this segment
zUpdates = [0 dz*zUpdateIdxs];

%% Calculate proportionality factors for use in the mex function
ax = dz/(4i*dx^2*k_0*P.n_0);
ay = dz/(4i*dy^2*k_0*P.n_0);
d = -dz*k_0/(2*P.n_0); % defined such that in each step in the mex function, E = E*multiplier*exp(1i*d*(n^2-n_0^2))

%% Define the multiplier
absorber = exp(-dz*max(0,max(abs(Y) - P.Ly_main/2,abs(X) - P.Lx_main/2)).^2*P.alpha);
multiplier = absorber; % This could also include a phase gradient due to bending

%% Figure initialization
h_f = figure(P.figNum);clf;
h_f.WindowState = 'maximized';

h_axis1 = subplot(2,2,1);
if P.downsampleImages
  h_im1 = imagesc(x_plot,y_plot,zeros(min(500,Ny),min(500,Nx),'single'));
else
  h_im1 = imagesc(x,y,zeros(Ny,Nx,'single'));
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
title('Refractive index');

powers = NaN(1,P.updates+1);
powers(1) = sum(abs(E(:)).^2)/P_0;
subplot(2,2,2);
h_plot2 = plot(zUpdates,powers,'linewidth',2);
xlim([0 P.Lz]);
xlabel('Propagation distance [m]');
ylabel('Relative power remaining');

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
  caxis([0 P.plotEmax]);
else
  caxis('auto');
end

h_axis3b = subplot(2,2,4);
hold on;
box on;
maxE0 = max(E(:));
if P.downsampleImages
  h_im3b = imagesc(x_plot,y_plot,angle(E(ix_plot,iy_plot).'/maxE0)); 
  h_im3b.AlphaData = max(0,(1+log10(abs(E(ix_plot,iy_plot).'/maxE0).^2)/3));  %Logarithmic transparency in displaying phase outside cores
else
  h_im3b = imagesc(x,y,angle(E.'/maxE0)); 
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

% tic;
%% Load variables into a parameters struct and start looping, one iteration per update
mexParameters = struct('dx',single(dx),'dy',single(dy),'taperPerStep',single((1-P.taperScaling)/Nz),'twistPerStep',single(P.twistRate*P.Lz/Nz),...
  'shapes',single(P.shapes),'n_cladding',single(P.n_cladding),'multiplier',complex(single(multiplier)),'d',single(d),'n_0',single(P.n_0),...
  'ax',single(ax),'ay',single(ay),'useAllCPUs',P.useAllCPUs,'RoC',single(P.bendingRoC),'rho_e',single(P.rho_e),'bendDirection',single(P.bendDirection));

mexParameters.iz_start = int32(0); % z index of the first z position to step from for the first call to FDBPMpropagator, in C indexing (starting from 0)
mexParameters.iz_end = int32(zUpdateIdxs(1)); % last z index to step into for the first call to FDBPMpropagator, in C indexing (starting from 0)
for updidx = 1:length(zUpdateIdxs)
  if updidx > 1
    mexParameters.iz_start = int32(zUpdateIdxs(updidx-1));
    mexParameters.iz_end   = int32(zUpdateIdxs(updidx));
  end
  checkMexInputs(E,mexParameters);
  if P.useGPU
    [E,n] = FDBPMpropagator_CUDA(E,mexParameters);
  else
    [E,n] = FDBPMpropagator(E,mexParameters);
  end

  %% Update figure contents
  if P.downsampleImages
    h_im1.CData = n(ix_plot,iy_plot).'; % Refractive index at this update
    h_im3a.CData = abs(E(ix_plot,iy_plot).').^2; % Intensity at this update
    h_im3b.CData = angle(E(ix_plot,iy_plot).'/maxE0); % Phase at this update
    h_im3b.AlphaData = max(0,(1+log10(abs(E(ix_plot,iy_plot).'/max(abs(E(:)))).^2)/3));  %Logarithmic transparency in displaying phase outside cores
  else
    h_im1.CData = n.'; % Refractive index at this update
    h_im3a.CData = abs(E.').^2; % Intensity at this update
    h_im3b.CData = angle(E.'/maxE0); % Phase at this update
    h_im3b.AlphaData = max(0,(1+log10(abs(E.'/max(abs(E(:)))).^2)/3));  %Logarithmic transparency in displaying phase outside cores
  end
  if updidx == 1
    caxis(h_axis1,'auto'); % To refresh the numbers on the color bar
  end
  if ~isfield(P,'plotEmax')
    caxis(h_axis3a,'auto');
  end
  
  powers(updidx+1) = sum(abs(E(:)).^2)/P_0; 
  h_plot2.YData = powers;
  drawnow;
  
  if P.saveVideo
    frame = getframe(h_f); 
    writeVideo(video,frame); 
  end
end
% toc

if P.saveVideo && ~isfield(P,'videoHandle')
	close(video);
end

Estruct = struct('field',E,'Lx',Lx,'Ly',Ly,'x',x,'y',y);

% S = load('train');
% sound(S.y.*0.1,S.Fs);
end

function checkMexInputs(E,parameters)
assert(all(isfinite(E(:))));
assert(~isreal(E));
assert(isa(E,'single'));
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