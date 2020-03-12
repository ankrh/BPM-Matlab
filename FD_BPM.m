% Author: Madhu Veettikazhy
% Date:12 November 2018
% Gaussian beam (CW) electric field propagation from the
% through fibre using Finite Difference Beam Propagation Method
% in 2D
% Reference: An assessment of FD BPM method - Youngchul Chung
%
% Douglas-Gunn Alternating Direction Implicit approach is based on me690-lctr-nts.pdf Eqs. 3.21, 3.22, 3.23a and 3.23b.
% Implicit steps are performed using the Thomson algorithm (https://en.wikipedia.org/wiki/Tridiagonal_matrix_algorithm)
% ***************************************************************************************************
format long
format compact

% Test of execution speed, 10 steps through a 2000x2000 grid, only the one
% mandatory update:
% 
% calctype 1: 657 seconds
% calctype 2: 48.7 seconds
% calctype 3: 1.56 seconds
% 
% That means method 3 is 421 times faster than method 1
% 
% GPU speedup comparison on 2020-02-02, with 1000 steps through a 2000x2000 grid, one update:
% calctype 3: 114 seconds
% calctype 4: 83.0 seconds when compiling with GCC, slightly slower when compiling with MSVC
% calctype 5: 12.2 seconds
% The speedup from calctype 3 to 5 is 9x
% Thus, the speedup from calctype 1 to 5 is ~3800x

%% General parameters
lambda = 1e-6;          % [m] Wavelength
w_0 = 3e-6;             % [m] Initial waist plane 1/e^2 radius of the gaussian beam
Lz = 1e-2;              % [m] z propagation distance

taperScaling = 0.5; % Specifies how much the refractive index profile of the last z slice should be scaled relative to the first z slice, linearly scaling in between
twistRate = 0; % Specifies how rapidly the fiber twists, measured in radians per metre

n_cladding = 1.45;
n_core = 1.46;

%% Resolution-related parameters
targetzstepsize = 1e-6; % [m] z step size to aim for
Lx_main = 50e-6;        % [m] x side length of main area
Ly_main = 50e-6;        % [m] y side length of main area
Nx_main = 200;          % x resolution of main area
Ny_main = 200;          % y resolution of main area

%% Solver-related parameters
useAllCPUs = true;
useGPU = false;

targetLx = 1.5*Lx_main;   % [m] Full area x side length, including absorber layer
targetLy = 1.5*Ly_main;   % [m] Full area y side length, including absorber layer
alpha = 3e14;             % [1/m^3] "Absorption coefficient" per unit length distance out from edge of main area, squared

%% Visualization parameters
updates = 100;            % Number of times to update plot. Must be at least 1, showing the final state. If set higher than Nz (such as Inf), script will simply update every step.
colormax = 1e10;          % Maximum to use for the color scale in figure 3a

%% Check for GPU compatibility if needed
if useGPU
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
Nz = round(Lz/targetzstepsize) % Number of z steps

dx = Lx_main/Nx_main;
dy = Ly_main/Ny_main;
dz = Lz/Nz;

Nx = round(targetLx/dx);
if rem(Nx,2) ~= rem(Nx_main,2)
  Nx = Nx + 1; % Ensure that if Nx_main was set odd (to have a x slice at the center), Nx will also be odd
end
Ny = round(targetLy/dy);
if rem(Ny,2) ~= rem(Ny_main,2)
  Ny = Ny + 1; % Ensure that if Ny_main was set odd (to have a y slice at the center), Ny will also be odd
end
Lx = Nx*dx;
Ly = Ny*dy;

x = dx*(-(Nx-1)/2:(Nx-1)/2);
y = dy*(-(Ny-1)/2:(Ny-1)/2);
z = dz*(0:Nz); % Has Nz + 1 elements
[X,Y] = ndgrid(x,y);

k_0 = 2*pi/lambda; % [m^-1] Wavenumber

%% Initialisation of optical fibre parameters
n_0 = n_core;

shapeTypes = [1 ; 2];
shapeParameters = [1e-5  2e-5 1.5e-5 ; -1e-5 2e-5 1.5e-5];
shapeRIs = [n_core ; n_core];

%% Beam initialization
% amplitude = exp(-((X-Lx_main/4).^2+Y.^2)/w_0^2) - exp(-((X+Lx_main/4).^2+Y.^2)/w_0^2); % Gaussian field amplitude
amplitude = exp(-((X-Lx_main/10).^2+Y.^2)/w_0^2); % Gaussian field amplitude
% amplitude = exp(-(X.^2+Y.^2)/w_0^2); % Gaussian field amplitude
phase = 8e5*Y;
E = amplitude.*exp(1i*phase); % Electric field
E = complex(single(E/sqrt(dx*dy*sum(abs(E(:)).^2))));

%% Fibre propagation and plotting
updatezindices = round((1:min(Nz,updates))/min(Nz,updates)*Nz);
z_updates = dz*updatezindices;
nextupdatenumber = 1;

ax = dz/(4i*dx^2*k_0*n_0);
ay = dz/(4i*dy^2*k_0*n_0);

d = -dz*k_0/(2*n_0); % E = E*multiplier*exp(1i*d*(n^2-n_0^2))

absorber = exp(-dz*max(0,max(abs(Y) - Ly_main/2,abs(X) - Lx_main/2)).^2*alpha);
multiplier = absorber; % This could also include a phase gradient due to bending

parameters = struct('dx',single(dx),'dy',single(dy),'taperPerStep',single(1-(1-taperScaling)/Nz),'twistPerStep',single(twistRate*Lz/Nz),...
  'shapeTypes',uint8(shapeTypes),'shapeParameters',single(shapeParameters),'shapeRIs',single(shapeRIs),'n_cladding',single(n_cladding),'multiplier',complex(single(multiplier)),...
  'd',single(d),'n_0',single(n_0),'ax',single(ax),'ay',single(ay),'useAllCPUs',useAllCPUs);

%% Figure initialization
figure(1);clf;
h_im1 = imagesc(x,y,zeros(Ny,Nx));
axis xy
axis equal
xlim([-Lx/2 Lx/2]);
ylim([-Ly/2 Ly/2]);
colorbar;
xlabel('x [m]');
ylabel('y [m]');
title('Refractive index');

powers = NaN(1,min(Nz,updates)+1);
powers(1) = sum(dx*dy*abs(E(:)).^2);
figure(2);clf;
plot([0 z_updates],powers,'YDataSource','powers','linewidth',2);
xlim([0 Lz]);
xlabel('Propagation distance [m]');
ylabel('Relative power remaining');

figure(3);clf;
subplot(2,1,1);
hold on;
h_im3a = imagesc(x,y,abs(E.').^2);
axis xy;
axis equal;
xlim([-Lx/2 Lx/2]);
ylim([-Ly/2 Ly/2]);
colorbar;
caxis('manual');
xlabel('x [m]');
ylabel('y [m]');
title('Intensity [W/m^2]');
% if max(n_mat(:) > min(n_mat(:))); contour(X,Y,n_mat,(n_cladding+eps(n_cladding))*[1 1],'color','w','linestyle','--'); end
line([-Lx_main Lx_main Lx_main -Lx_main -Lx_main]/2,[Ly_main Ly_main -Ly_main -Ly_main Ly_main]/2,'color','r','linestyle','--');
caxis([0 colormax]);
subplot(2,1,2);
hold on;
h_im3b = imagesc(x,y,angle(E.'));
axis xy;
axis equal;
xlim([-Lx/2 Lx/2]);
ylim([-Ly/2 Ly/2]);
colorbar;
caxis([-pi pi]);
% if max(n_mat(:) > min(n_mat(:))); contour(X,Y,n_mat,(n_cladding+eps(n_cladding))*[1 1],'color','w','linestyle','--'); end
line([-Lx_main Lx_main Lx_main -Lx_main -Lx_main]/2,[Ly_main Ly_main -Ly_main -Ly_main Ly_main]/2,'color','r','linestyle','--');
xlabel('x [m]');
ylabel('y [m]');
title('Phase [rad]');
colormap(gca,hsv/1.5);

parameters.iz_start = int32(0); % z index of the first z position to step from for the first call to FDBPMpropagator, in C indexing (starting from 0)
parameters.iz_end = int32(updatezindices(1)); % last z index to step into for the first call to FDBPMpropagator, in C indexing (starting from 0)
tic;
for updidx = 1:length(updatezindices)
  if updidx > 1
    parameters.iz_start = int32(updatezindices(updidx-1));
    parameters.iz_end   = int32(updatezindices(updidx));
  end
  if useGPU
    [E,d] = FDBPMpropagator_CUDA(E,parameters); % d is a parameter from which we can retrieve the refractive index profile that the mex function calculated
  else
    [E,d] = FDBPMpropagator(E,parameters); % d is a parameter from which we can retrieve the refractive index profile that the mex function calculated
  end

  h_im1.CData = sqrt(log(d)*n_0*2i/dz/k_0 + n_0^2); % Refractive index at this update
  h_im3a.CData = abs(E.').^2; % Intensity at this update
  h_im3b.CData = angle(E.'); % Phase at this update
  powers(updidx+1) = dx*dy*sum(abs(E(:)).^2);
  refreshdata(2);
  drawnow;
end
toc

function RIprofiles = RIprofilesFunc(z,n_core)
% In this function, the refractive index profile at each z position is
% defined in terms of fundamental mathematical shapes. Currently only
% circular disks are supported.

% Shapes is a cell array with one cell for each z position in the z array
% passed in. Each cell contains a 2D array in which each column corresponds
% to a shape, and the rows correspond to different parameters of the shape.
% Each column consists of these rows: [shape type; x center coord; y center
% coord; radius; refractive index]. shape type can be 1: Disk, 2:
% Rectangle, 3: Triangle. Currently, only shape type 1 is implemented.

for iz = 1:length(z)
  RIprofiles{iz} =  [1 1; 1e-5 -1e-5; 2e-5 2e-5; 1.5e-5 1.5e-5; n_core n_core];
end
end