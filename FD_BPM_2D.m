% Author: Madhu Veettikazhy
% Date:12 November 2018
% Gaussian beam (CW) electric field propagation from the 
% through fibre using Finite Difference Beam Propagation Method
% in 2D
% Reference: An assessment of FD BPM method - Youngchul Chung
% ***************************************************************************************************

%% User-specified general parameters
Nx = 2000;                                                          % x resolution
Ny = 1;                                                                % y resolution - set to 1 since 2D
Nz = 1000;                                                           % Number of z steps inside fibre
Lx = 200e-6;                                                        % [m] x side length
Lz = 1e-3;                                                            % [m] z propagation distance in fibre

lambda = 1e-6;                                                    % [m] Wavelength
w = 10e-6;                                                           % [m] Initial waist plane 1/e^2 radius of the gaussian beam

x_offset = 0;                                                        % [m] Offset in the gaussian beam in x direction

n_core =1.47;                                                      % [] Core refractive index
n_cladding = 1.47;                                              % [] Cladding refractive index
core_radius=10e-6;                                             % [m] Radius of fibre core

%% Initialisation of space and frequency grids
dx = Lx/Nx;
dz = Lz/Nz;

x = linspace(-Lx/2,Lx/2,Nx);
k_0 = 2*pi/lambda;                                              % [m^-1] Wavenumber

%% Initialisation of optical fibre parameters
n_0 = n_cladding;
n_fibre=ones(Nx,Ny)*n_core;
n_fibre(abs(x)> core_radius) = n_cladding;

%% Beam initialization
amplitude = exp(-(x-x_offset).^2/w^2);                       % Gaussian field amplitude
phase =  zeros(Ny,Nx);      %2.5e5.*(x);                      % Phase
E = amplitude.*exp(1i*phase);                                     % Electric field at z=0
% plot(x,abs(E).^2);

%% Absorber layer at boundaries
% Reference: Beam-propagation analysis of loss in bent optical waveguides and fibers
% J. Saijonmaa, D. Yevick
pix_x_a=round(Nx/2.5);
pix_x_b=round(Nx/30);
gamma=1/10;
absorb(pix_x_a:Nx-pix_x_a) = 1;
absorb(pix_x_b:pix_x_a) = (0.5*(1-cos((pi*((pix_x_b:pix_x_a)-pix_x_b)./(pix_x_a-pix_x_b))))).^gamma;
absorb(Nx-pix_x_a:Nx-pix_x_b) = (0.5*(1-cos(pi*((Nx-pix_x_a:Nx-pix_x_b)-(Nx-pix_x_b))./(pix_x_a-pix_x_b)))).^gamma;
absorb(1:pix_x_b) = 0;
absorb(Nx-pix_x_b:Nx) = 0;
% line(x, absorb.','Color','red', 'LineWidth', 1);    
absorb_column=absorb(:);

absorb_pseudo=absorb;                                      %absorb_pseudo is only for plotting purposes 
absorb_pseudo(pix_x_b:pix_x_a)=0;                  %to see the computational window  of interest in a rectangle
absorb_pseudo(Nx-pix_x_a:Nx-pix_x_b)=0;

%% Fibre propagation and plotting
delta_n_2 = n_fibre.^2-n_0^2;                              %delta_n^2 in the expression for FD BPM
delta_n_2_column = delta_n_2(:);                         %Converting from NxN square matrix to N^2x1 column matrix
E_column=E(:);                                                   %Converting 2D matrix into a column matrix
N=Nx*Ny;                                                            %N*N - size of sparse matrices

% Sparse matrix
a = dz/(2*dx^2);
b = dz/(dx^2) - (k_0^2*dz/2).*delta_n_2_column + (2*1i*k_0*n_0);
c = (-dz/(dx^2)) + (k_0^2*dz/2).*delta_n_2_column+ (2*1i*k_0*n_0);

M_i_lhs = sparse(1:N,1:N,b(1:N),N,N);
M_iMinus1 = sparse(2:N,1:N-1,a,N,N);
M_iPlus1 = sparse(1:N-1,2:N,a,N,N);
for i=Nx:Nx:N-Nx                                               %Dirichlet's boundary condition
    M_iMinus1(i+1,i)=0;
    M_iPlus1(i,i+1)=0;
end
M_lhs = -M_iMinus1+M_i_lhs-M_iPlus1;

M_i_rhs = sparse(1:N,1:N,c(1:N),N,N);
M_rhs = M_iMinus1+M_i_rhs+M_iPlus1;

figure(1);plot(x,abs(E').^2); 

%For plotting the diverged gaussian field from gaussian_Prop_2D.m
% NxNy_ref = size(E_ref);
% Nx_ref=NxNy_ref(1);
% Nx_ref_start=round(Nx_ref/2)-round(Nx/2);
% Nx_ref_end=round(Nx_ref/2)+round(Nx/2)-1;

for zidx = 1:Nz
    RHS = M_rhs*E_column;
    E_column = M_lhs\RHS;
    E_column=E_column.*absorb_column;
    E_2D = reshape(E_column, [Ny Nx]);
    
    plot(x,abs(E_2D').^2);                                                                                                      %Current field in blue
    title(['Intensity at z = ' num2str(zidx*dz,'%.1e') ' m']);
    axis([-Lx/2 Lx/2 0 1]);
%     hold on;
%     plot(x,abs(E_ref(Nx_ref_start:Nx_ref_end,:).').^2,'Color','green', 'LineWidth', 0.25);        %Reference field in green
%     line(x, absorb.','Color','red', 'LineWidth', 1);                                                                      %Absorbed profile in red
%     line(x, absorb_pseudo.','Color','yellow', 'LineWidth', 0.25);                                               %Computational window of interest in yellow
%     hold off;
    drawnow;
end