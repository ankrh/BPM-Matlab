%%
% Author: Madhu Veettikazhy
% Date:12 November 2018
% Gaussian beam (CW) electric field propagation from the 
% distal end to the proximal end of the optical fibre. 
% ***************************************************************************************************

%% User-specified general parameters
Nx = 200;                                                            % x resolution
Ny = 200;                                                            % y resolution
Nz_air = 100;                                                      % Number of z steps in air
Nz_fibre = 100;                                                   % Number of z steps inside fibre

Lx = 20e-6;                                                         % [m] x side length
Ly = 20e-6;                                                         % [m] y side length
Lz_air = 100e-6;                                                  % [m] z propagation distance in air
Lz_fibre = 1e-3;                                                  % [m] z propagation distance in fibre

lambda = 1e-6;                                                    % [m] Wavelength 
k_0 = 2*pi/lambda;                                              % [m^-1] Wavenumber
w = 2.33e-6;                                                        % [m] Initial waist plane 1/e^2 radius of the gaussian beam

x_offset = 0;                                                        % [m] Offset in the gaussian beam in x direction
y_offset = 0;                                                        % [m] Offset in the gaussian beam in y direction
%% Initialisation of space and frequency grids
dx = Lx/Nx;
dy = Ly/Ny;
dz_air = Lz_air/Nz_air;
dz_fibre = Lz_fibre/Nz_fibre;

x = dx*(-Nx/2:Nx/2-1);
y = dy*(-Ny/2:Ny/2-1);
[X,Y] = ndgrid(x,y);
RHO = sqrt(X.^2+Y.^2);

kx = 2*pi/Lx*(-Nx/2:Nx/2-1);
ky = 2*pi/Ly*(-Ny/2:Ny/2-1);
[kX,kY] = ndgrid(kx,ky);

%% Initialisation of optical fibre parameters
n_core =1.48;                                                      % [] Core refractive index
n_cladding = 1.47;                                               % [] Cladding refractive index
n_0 = n_cladding;                                                % [] Reference refractive index

core_radius=2e-6;                                                % [m] Radius of fibre core
n_fibre=ones(Nx,Ny)*n_core;                              % Refractive index matrix of the complete fibre
n_fibre(RHO> core_radius) = n_cladding;


%% Beam initialization
amplitude = exp(-((X-x_offset).^2+(Y-y_offset).^2)/w^2);      % Gaussian field amplitude
phase = zeros(Nx,Ny);                                                           % Phase
E = amplitude.*exp(1i*phase);                                               % Electric field
I_in = max(max(abs(E.').^2)); 
% mesh(X,Y,abs(E.').^2);

%% Figure initialization
h_fig1 = figure(1);
set(h_fig1,'Position',[50 50 4*Nx+125 2*Ny+75]); clf;
h_ax1 = axes('Units','pixels','Position', [50 50 2*Nx 2*Ny]);
h_im1 = imagesc(x,y,angle(E.'));
axis tight
set(h_ax1, 'xlimmode','manual', 'ylimmode','manual');

h_ax2 = axes('Units','pixels','Position', [2*Nx+100 50 2*Nx 2*Ny]);
h_im2 = imagesc(x,y,abs(E.').^2);
axis xy
set(h_ax2, 'xlimmode','manual', 'ylimmode','manual');


%% Fibre propagation and plotting
delta_n_2 = n_fibre.^2-n_0^2;                              %delta_n^2 in the expression for FD BPM
delta_n_2_column = delta_n_2(:);                         %Converting from NxN square matrix to N^2x1 column matrix

% Matrices a, b, and c in the expression for FD BPM
a = dz_fibre/(2*dx^2);
b = 2*dz_fibre/(dx^2) - (k_0^2*dz_fibre/2).*delta_n_2_column + (2*1i*k_0*n_0);
c = (-2*dz_fibre/(dx^2)) + (k_0^2*dz_fibre/2).*delta_n_2_column+ (2*1i*k_0*n_0);

%Sparse matrices
N=Nx*Ny;
M_i_lhs = sparse(1:N,1:N,b(1:N),N,N);
M_iMinus1 = sparse(2:N,1:N-1,a,N,N);
M_iPlus1 = sparse(1:N-1,2:N,a,N,N);

for i=Nx:Nx:N-Nx                                                %Dirichlet's boundary condition
    M_iMinus1(i+1,i)=0;                                        %Setting the boundary values to zero                     
    M_iPlus1(i,i+1)=0;
end

M_jMinus1 = sparse(Nx+1:N,1:N-Nx,a,N,N);
M_jPlus1 = sparse(1:N-Nx,Nx+1:N,a,N,N);
M_lhs = -M_iMinus1-M_iPlus1+M_i_lhs-M_jMinus1-M_jPlus1;

M_i_rhs = sparse(1:N,1:N,c(1:N),N,N);
M_rhs = M_iMinus1+M_iPlus1+M_i_rhs+M_jMinus1+M_jPlus1;

%FD BPM
E_column=E(:);                                                   %Converting 2D matrix into a column matrix

theta = 0 : 0.01 : 2*pi;                                          %For plotting core circle
x_circle = core_radius * cos(theta) ;
y_circle = core_radius * sin(theta) ;

video = VideoWriter('fibre_field.avi');                   %For saving the propagation frames as a video
open(video);
tic
for zidx = 1:Nz_fibre
    RHS = M_rhs*E_column;
    E_column = M_lhs\RHS;
    E_2D = reshape(E_column, [Ny Nx]);             %Reshaping column matrix to 2D for plotting
    
    h_ax1.Title.String = ['Phase at z = ' num2str(zidx*dz_fibre,'%.1e') ' m'];
    h_im1.CData = angle(E_2D.');
    h_ax2.Title.String = ['Fibre intensity at z = ' num2str(zidx*dz_fibre,'%.1e') ' m'];
    h_im2.CData = abs(E_2D.').^2;
    hold on;
    plot(x_circle, y_circle,'white', 'LineWidth', 1);    %Plotting the core circle  
    hold off;
    drawnow;
    frame = getframe(gcf);                                     %Get the frames 
    writeVideo(video,frame);                                  %Stitch the frames to form a video and save
end
toc
close(video);

