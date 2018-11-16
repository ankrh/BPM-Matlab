%%
% Author: Madhu Veettikazhy
% Date: 16 November 2018
% Gaussian beam (CW) electric field 2D propagation in air
% ***************************************************************************************************

%% User-specified general parameters
Nx = 2200;                                                        % x resolution
Ny = 1;                                                             % y resolution
Nz = 1000;                                                        % Number of z steps

Lx = 400e-6;                                                     % [m] x side length
Lz = 1e-3;                                                          % [m] z propagation distance

lambda = 1e-6;                                                % [m] Wavelength 
w = 10e-6;                                                       % [m] Initial waist plane 1/e^2 radius of the gaussian beam

%% Initialization of space and frequency grids
dx = Lx/Nx;
dz = Lz/Nz;

x = linspace(-Lx/2,Lx/2,Nx);
y = linspace(-Ly/2,Ly/2,Ny);
[X,Y] = ndgrid(x,y);

kx = 2*pi/Lx*(-Nx/2:Nx/2-1);
ky = 2*pi/Ly*(-Ny/2:Ny/2-1);
[kX,kY] = ndgrid(kx,ky);

%% Beam initialization
amplitude = exp(-(X.^2)/w^2);                             % Gaussian field amplitude
phase = zeros(Nx,Ny);                                         % Phase
E_ref = amplitude.*exp(1i*phase);                      % Electric field
       
%% Fresnel Propagation and plotting
prop_kernel = ifftshift(exp(-1i*dz*(kX.^2)*lambda/(4*pi))); % Fresnel propagation kernel

for zidx = 1:Nz
    E_ref = ifft2(fft2(E_ref).*prop_kernel);
    figure(1); plot(x,abs(E_ref.').^2);
    title(['Intensity at z = ' num2str(zidx*dz,'%.1e') ' m']);
    axis([-Lx/2 Lx/2 0 1]);
    drawnow;
end


