%%
% Author: Madhu Veettikazhy
% Date:12 November 2018
% Gaussian beam (CW) electric field propagation from the 
% distal end to the proximal end of the optical fibre. 
% ***************************************************************************************************
saveVideo = false;

%% User-specified general parameters
Nx = 400;                                                            % x resolution
Ny = 400;                                                            % y resolution
Nz_fibre = 100;                                                   % Number of z steps inside fibre

Lx = 40e-6;                                                         % [m] x side length
Ly = 40e-6;                                                         % [m] y side length
Lz_fibre = 1e-3;                                                  % [m] z propagation distance in fibre

lambda = 1e-6;                                                    % [m] Wavelength 
w = 2.38e-6;                                                        % [m] Initial waist plane 1/e^2 radius of the gaussian beam

x_offset = 0;                                                        % [m] Offset in the gaussian beam in x direction
y_offset = 0;                                                        % [m] Offset in the gaussian beam in y direction

n_core =1.48;                                                      % [] Core refractive index
n_cladding = 1.47;                                               % [] Cladding refractive index
core_radius=2e-6;                                                % [m] Radius of fibre core

%% Initialisation of space and frequency grids
dx = Lx/Nx;
dy = Ly/Ny;
dz_fibre = Lz_fibre/Nz_fibre;

x = linspace(-Lx/2,Lx/2,Nx);
y = linspace(-Ly/2,Ly/2,Ny);
[X,Y] = ndgrid(x,y);
RHO = sqrt(X.^2+Y.^2);

k_0 = 2*pi/lambda;                                              % [m^-1] Wavenumber

%% Initialisation of optical fibre parameters
n_0 = n_cladding;                                                % [] Reference refractive index
n_fibre=ones(Nx,Ny)*n_core;                              % Refractive index matrix of the complete fibre
n_fibre(RHO> core_radius) = n_cladding;

%% Beam initialization
amplitude = exp(-((X-x_offset).^2+(Y-y_offset).^2)/w^2);      % Gaussian field amplitude
phase = zeros(Nx,Ny);                                                           % Phase
E = amplitude.*exp(1i*phase);                                               % Electric field
% I_in = max(max(abs(E.').^2)); 
% mesh(X,Y,abs(E.').^2);

%% Figure initialization
theta = 0 : 0.01 : 2*pi;                                          %For plotting core circle
x_circle = core_radius * cos(theta) ;
y_circle = core_radius * sin(theta) ;

h_fig1 = figure(1);
% set(h_fig1,'Position',[50 50 1000 500]); clf;
h_ax1 = subplot(1,2,1);
% h_ax1 = axes('Units','pixels','Position', [50 50 400 400]);
h_im1 = imagesc(x,y,angle(E.'));
axis xy
axis equal
xlim(Lx*[-1/2 1/2]);
ylim(Ly*[-1/2 1/2]);
set(h_ax1, 'xlimmode','manual', 'ylimmode','manual');
line(x_circle, y_circle,'Color','white', 'LineWidth', 1);    %Plotting the core circle  

h_ax2 = subplot(1,2,2);
% h_ax2 = axes('Units','pixels','Position', [550 50 400 400]);
h_im2 = imagesc(x,y,abs(E.').^2);
axis xy
axis equal
xlim(Lx*[-1/2 1/2]);
ylim(Ly*[-1/2 1/2]);
set(h_ax2, 'xlimmode','manual', 'ylimmode','manual');
line(x_circle, y_circle,'Color','white', 'LineWidth', 1);    %Plotting the core circle  

% h_ax3 = subplot(1,3,3);
% h_pl3 = plot(x,abs(E(:,round(Ny/2))).^2);
% line(x,abs(E(:,round(Ny/2))).^2,'Linewidth',2,'Color','red');
% xlim([-0.5 0.5]*1e-5);

%% Fibre propagation and plotting
delta_n_2 = n_fibre.^2-n_0^2;                              %delta_n^2 in the expression for FD BPM
delta_n_2_column = delta_n_2(:);                         %Converting from NxN square matrix to N^2x1 column matrix

% Matrices a, b, and c in the expression for FD BPM
ax = dz_fibre/(2*dx^2);
ay = dz_fibre/(2*dy^2);
b = dz_fibre*(1/dx^2+1/dy^2) - (k_0^2*dz_fibre/2).*delta_n_2_column + (2*1i*k_0*n_0);
c = -dz_fibre*(1/dx^2+1/dy^2) + (k_0^2*dz_fibre/2).*delta_n_2_column+ (2*1i*k_0*n_0);

%Sparse matrices
N=Nx*Ny;
M_i_lhs = sparse(1:N,1:N,b(1:N),N,N);
M_iMinus1 = sparse(2:N,1:N-1,ax,N,N);
M_iPlus1 = sparse(1:N-1,2:N,ax,N,N);

for i=Nx:Nx:N-Nx                                                %Dirichlet's boundary condition
    M_iMinus1(i+1,i)=0;                                        %Setting the boundary values to zero                     
    M_iPlus1(i,i+1)=0;
end

M_jMinus1 = sparse(Nx+1:N,1:N-Nx,ay,N,N);
M_jPlus1 = sparse(1:N-Nx,Nx+1:N,ay,N,N);
M_lhs = -M_iMinus1-M_iPlus1+M_i_lhs-M_jMinus1-M_jPlus1;

M_i_rhs = sparse(1:N,1:N,c(1:N),N,N);
M_rhs = M_iMinus1+M_iPlus1+M_i_rhs+M_jMinus1+M_jPlus1;

%FD BPM
E_column=E(:);                                                   %Converting 2D matrix into a column matrix

if saveVideo
	video = VideoWriter('fibre_field.avi');                   %For saving the propagation frames as a video
	open(video);
end
tic
for zidx = 1:Nz_fibre
	RHS = M_rhs*E_column;
	E_column = M_lhs\RHS;
	E_2D = reshape(E_column, [Nx Ny]);             %Reshaping column matrix to 2D for plotting

	h_ax1.Title.String = ['Phase at z = ' num2str(zidx*dz_fibre,'%.1e') ' m'];
	h_im1.CData = angle(E_2D.');
	h_ax2.Title.String = ['Fibre intensity at z = ' num2str(zidx*dz_fibre,'%.1e') ' m'];
	h_im2.CData = abs(E_2D.').^2;
% 	h_pl3.YData = abs(E_2D(:,round(Ny/2))).^2;

	drawnow;
	if saveVideo
		frame = getframe(gcf);                                     %Get the frames 
		writeVideo(video,frame);                                  %Stitch the frames to form a video and save
	end
end
toc
if saveVideo
	close(video);
end
