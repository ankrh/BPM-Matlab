%% Function for estimating all the LP modes supported by an optical fiber

% 1) Definition of the values considered for 
% azimuthal mode number, the normalized propagation constant, for which
% roots is found and V-number, including approximate number of modes supported.  
%
% 2) Function call of Transcendentaleq_fiber.m which calculates the
% transcendental equation for each azimuthal mode number supported by the
% fibre parameters and its roots. See Transcendentaleq_fiber.m for more information. 
%
% 3) Function call of LPmodeCalc.m which calculates all the LPlm (A) modes 
% (cosine) and LPlm (B) modes (sine) for each azimuthal mode number and
% determined roots. See LPmodeCalc.m for more information. 

%% 
function LP_modes_AB = getLPmodes(X,Y,core_radius,n_core,n_cladding,lambda)

%% Parameter section %%
N=size(X,2); %Resolution.

% Beam parameters
k0=2*pi/lambda; % Wavenumber for excitaiton wavelength.
E0=1; % Initial amplitude. 

% Initialization of the b-V curves and subsequently field calculation. 
N_b = 10^6; % Resolution of normalized propagation constant axis.  
bv = (1:N_b)./(N_b+1); % Normalized propagation constant
V = k0.*core_radius.*sqrt(n_core^2-n_cladding^2); % Normalized frequency
N_mode=V^2/2; % approximate number of modes supported by the fiber. 


%% The transcendentaleq_fiber function %%
% In this section we calculate the trancendental equations for each
% specified azimuthal mode number (L) and find the proper roots, 
% corresponding to the radial solutions for each azimuthal mode number.
[~,~,~,broots_mat] = Transcendentaleq_fiber(V,bv);

%% The LPmodeCalc function %%
% In this section the total number of modes is calculated on the basis of
% fibre parameters, beam parameters and determined roots from
% Transcendentaleq_fiber.m. 
intNorm=false; 
[LP_modes_AB] = LPmodeCalc(core_radius,n_core,n_cladding,N,X,Y,k0,E0,broots_mat,intNorm);
% % Allocating calculated LPlm (A) modes (cosine) and LPlm (B) modes (sine)
% % to different cell structures for mode overlap calculations. 
% Ecell_modeA=LP_modesAB{1,1};
% Ecell_modeB=LP_modesAB{2,1};

end

