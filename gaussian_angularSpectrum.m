%%
% Author: Madhu Veettikazhy
% Date: 22 March 2018
% Gaussian beam (CW) electric field propagation from the 
% focal point (initial waist plane width is defined here) 
% back to the distal end of the optical fiber. 
% ***************************************************************************************************

%% User-specified general parameters
Nx = 600;                                                        % x resolution
Ny = 300;                                                         % y resolution
Nz = 100;                                                        % Number of z steps

Lx = 250e-6;                                                     % [m] x side length
Ly = 300e-6;                                                         % [m] y side length
Lz = 2e-3;                                                       % [m] z propagation distance

lambda = 632.8e-9;                                               % [m] Wavelength 
w = 10e-6;                                                       % [m] Initial waist plane 1/e^2 radius of the gaussian beam

%% Initialization of space and frequency grids
dx = Lx/Nx;
dy = Ly/Ny;
dz = Lz/Nz;

x = dx*(-Nx/2:Nx/2-1);
y = dy*(-Ny/2:Ny/2-1);
[X,Y] = ndgrid(x,y);

kx = 2*pi/Lx*(-Nx/2:Nx/2-1);
ky = 2*pi/Ly*(-Ny/2:Ny/2-1);
[kX,kY] = ndgrid(kx,ky);

%% Beam initialization
amplitude = exp(-2*(X.^2+Y.^2)/w^2);                             % Gaussian field amplitude
phase = zeros(Nx,Ny);                                            % Phase
E = amplitude.*exp(1i*phase);                                    % Electric field

%% Figure initialization
h_fig1 = figure(1);
set(h_fig1,'Position',[50 50 Nx*2+150 Ny+150]); clf;
h_ax1 = axes('Units','pixels','Position', [50 75 Nx Ny]);
h_im1 = imagesc(x,y,abs(E.').^2);
axis xy
set(h_ax1, 'xlimmode','manual', 'ylimmode','manual');

h_ax2 = axes('Units','pixels','Position', [Nx+100 75 Nx Ny]);
h_im2 = imagesc(x,y,angle(E.'));
axis xy
set(h_ax2, 'xlimmode','manual', 'ylimmode','manual');
       
%% Fresnel Propagation and plotting
prop_kernel = ifftshift(exp(-1i*dz*(kX.^2 + kY.^2)*lambda/(4*pi))); % Fresnel propagation kernel

for zidx = 1:Nz
    E = ifft2(fft2(E).*prop_kernel);
    
    h_ax1.Title.String = ['Intensity at z = ' num2str(zidx*dz,'%.1e') ' m'];
    h_im1.CData = abs(E.').^2;
    h_ax2.Title.String = ['Phase at z = ' num2str(zidx*dz,'%.1e') ' m'];
    h_im2.CData = angle(E.');
    drawnow;
end






%%  Angular Spectrum Propagation
% 
% dist_z = 1e-6;                       % altitude (meters) 
% phy_x = w;               % physical width (meters) 
% phy_y = w;               % physical length (meters) 
% Fs_x = N/phy_x; 
% Fs_y = N/phy_y; 
% dx2 = Fs_x^(-1); 
% dy2 = Fs_y^(-1); 
% dFx = Fs_x/N;
% dFy = Fs_y/N; 
% Fx = (-Fs_x/2:dFx:(Fs_x/2 - dFx)); 
% Fy = (-Fs_y/2:dFy:(Fs_y/2 - dFy)); 
% 
% % alpha and beta (wavenumber components) 
% alpha = lambda.*Fx; 
% beta = lambda.*Fy;
% 
% % gamma_cust 
% gamma_cust = zeros(length(beta), length(alpha)); 
% for j = 1:length(beta)     
%     for i = 1:length(alpha)         
%         if (alpha(i)^2 + beta(j)^2) > 1             
%             gamma_cust(j, i) = 0;         
%         else
%             gamma_cust(j, i) = sqrt(1 - alpha(i)^2 - beta(j)^2);         
%         end
%     end
% end
% 
% % angular spectrm based formula 
% U1=E_gauss;
% for iter=1:100
% U1 = ifft2(ifftshift(fftshift(fft2(U1)).*exp(1i*k.*gamma_cust.*dist_z))); 
% I1 = (1/(16*pi)).*(U1.*conj(U1)); 
% imagesc(I1),title(['Propagated step: ', num2str(iter)]); 
% drawnow;
% % pause(1);
% end
