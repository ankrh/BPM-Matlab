%% Calculation of transverse LP mode profiles%%

% References:
% http://m2matlabdb.ma.tum.de/download.jsp?MC_ID=12&SC_ID=18&MP_ID=108
% https://www.rp-photonics.com/lp_modes.html
% Foundation of Guided-Wave Optics,Chen-Lin Chen, 2007 by Wiley & Sons. 
% Gloge D, Weakly Guiding Fibers, Appl. Opt. 1971, Vol 10, No. 10. 
% Fundamentals of Optical Fibers, John Buck, 2004 by Wiley and Sons. 
% Numerical Methods for Engineers and Scientists Using MATLAB, 2013 CRC
% Press.

% updated 02/04-2020.

% Who to run the code
% 1) Input appropriate parameters for given fibre of interest, such as grid
%    parameters, wavelength, n_core and n_clad etc. Here set the Azimuthal 
%    mode number L=0,1,2,3 etc.
%   
% 2) Set Rootfinding to true. The rootfinding algorithm will find the roots
%    for f=LHS-RHS on the basis of an "change in sign approach". The roots
%    is stored in broots. Since the trancendental equation is inherently 
%    singular, the rootfinding algorithm will find roots at the points of 
%    singularity. These "roots" are disregarded, because they do not 
%    correspond to actual guided modes.  The trancendental equation and 
%    subsequently the stored roots can be depicted in figure (1), by 
%    setting LHSRHSplot to true. 
%    
% 3) Set LPmodeCalc to true. This will calculate and store the modal
%    fields, including the abs.^2 of the modal field, in a cell array for 
%    each root found in step 2. A particular array in the cell structure
%    for plotting of the modal field can be accessed by changing idx_root.
%    The radial mode number M should be changed accordingly for correct
%    titling of the figures.
%
% 5) Set Modaloverlap to true. The modaloverlap calculates the 
%    overlap of the calculated LPmodes with an incident Gaussian beam. 
%
% 4) Set SaveData to true. When satisfied with the calculations, the
%    matrices containing relevant modal field data is saved.  

% Preamble 
clear all; close all; clc; tic;
Rootfinding=true; LHSRHSplot=false; LPmodeCalc=true; 
Modaloverlap=false; SaveData=false;

FileName = 'LPmodes_25um_980nm_L3.mat'; %added 02-05-2020
dataName =[FileName '.mat']; % 02-05-2020

% Initialization of the fiber
N=600; %Resolution.
Lx = 375e-6; % Grid length in x-direction.
Ly = 375e-6; % Grid length in y-direction.
r_core = 2.65e-6; % Core radius.
d_fiber = 2.*r_core; % Fiber diameter.
n_core = 1.4511; % Core index.
n_clad = 1.4444; % Cladding index.
Delta=(n_core-n_clad)/n_clad; % Index constrast, check for weakly guiding approximation.
xv=linspace(-Lx/2,Lx/2,N); % Defining the x-direction;
yv=linspace(-Ly/2,Ly/2,N); % Defining the y-direction;
[X,Y]=ndgrid(xv,yv); % generating the computational mesh (x,y)-grid. 

% Initialization of the b-V curves and subsequently field calculation. 
L=0; % The azimuthal mode number.
lambda=1.550e-6; % Excitation wavelength. 
k_wave = 2*pi/lambda; % Wavenumber.
E_0=1; % initial amplitude. 
N_b=10^6; % Resolution of normalized propagation constant axis.  
bv = (1:N_b)/(N_b+1); % Normalized propagation constant
V=k_wave.*r_core.*sqrt(n_core^2-n_clad^2); % Normalized frequency
W=V.*sqrt(bv); % Coefficient, Modified bessel function of second kind.
U=V.*sqrt(1-bv); % Coefficient, Bessel function of first kind.  

if Rootfinding
  % Firstly calculation of the left and right handed side of
  % transcendental equation.
  if L==0
    % Calculation of the left hand side at rho=r_core.
    LHS=U.*besselj(1,U)./besselj(0,U);
    % Calculation of the right hand side at rho=r_core.
    RHS=W.*besselk(1,W)./besselk(0,W);
  elseif L>=1
    % Calculation of the left hand side at rho=r_core.
    LHS=U.*besselj(L-1,U)./besselj(L,U);
    % Calculation of the right hand side at rho=r_core.
    RHS=-W.*besselk(L-1,W)./besselk(L,W);
  end
  
  f = LHS - RHS; % Definition of the trancendental equation. 
  signchangepositions = find(abs(diff(f>0))); % finding changes of sign. 
  % finding root values of f = LHS - RHS at the change of sign positions. 
  for i=length(signchangepositions):-1:1
    if abs(f(signchangepositions(i))-f(signchangepositions(i) + 1)) > abs(f(signchangepositions(i) - 1)-f(signchangepositions(i) + 2))
      signchangepositions(i) = []; % Remove this element because it's a singularity
    end
  end
  broots = fliplr(bv(signchangepositions)) % Storing of root values from root finding algorithm.
  
  if LHSRHSplot
    % plotting the transcendental equation to depict intersections of LHS and
    % RHS.
    plotylims = [-10 10];
    figure(1)
    plot(bv,LHS,bv,RHS,bv,LHS-RHS,[0 1],[0 0],'k','LineWidth',2)
    grid on; grid minor; box on;
    ylim(plotylims)
    xlabel('Normalized propagation constant b')
    legend('LHS','RHS','LHS - RHS','Location','northeast','AutoUpdate','off')
    title('Transcendental equation')
    for i=1:length(broots)
      line([broots(i) broots(i)],[plotylims(1) 0],'linestyle','--','color','k')
    end
  end
end

if LPmodeCalc
  Emat_field = cell(length(broots),1); % initializing E matrix.
  for M = 1:length(broots)
    Emat_field{M} = NaN(N,N);
    b = broots(M); % Use the root as normalized propagation constant.
    % Calculation of beta from normalized propagation constant.
    beta = k_wave.*sqrt((n_clad^2+(b.*(n_core^2-n_clad^2)))); % Propagation constant.
    % Calculation of normalized phase and attenuation constants.
    u = r_core*sqrt(n_core^2*k_wave^2-beta^2);
    w = abs(r_core*sqrt(beta^2-n_clad^2*k_wave^2));
    
    R_ratio = sqrt(X.^2 + Y.^2)/r_core;
    
    Emat_field{M}(R_ratio <  1) = E_0.*besselj(L,u*R_ratio(R_ratio < 1)).*cos(L*atan2(Y(R_ratio < 1),X(R_ratio < 1)));
    Emat_field{M}(R_ratio >= 1) = E_0.*besselj(L,u)./besselk(L,w).*besselk(L,w*R_ratio(R_ratio >= 1)).*cos(L*atan2(Y(R_ratio >= 1),X(R_ratio >= 1)));
  end
  
  % plotting one of the calculated LP mode on the basis of the electric field.
  M = 1; % Choose a particular root from vector, which contains the roots, broots. The m'th root in broots
  % corresponds to the m'th radial mode, m=1,2,3,4...
  figure(2)
  meshz(X,Y,abs(Emat_field{M}.').^2)
  colormap(jet((2^8)))
  axis tight
  title(['Intensity plot for LP_',num2str(L),'_',num2str(M),'mode',' (b=',num2str(broots(M)),')'],'fontsize',15)
  
  figure(3);
  subplot(2,1,1);
  box on;
  imagesc(xv,yv,abs(Emat_field{M}.').^2);
  colormap(jet((2^8)));
  axis equal tight
  c = colorbar;
  c.Label.String = 'Intensity [Arb.]';
  set(c,'fontsize',15)
  xlabel('x [m]','fontsize',12)
  ylabel('y [m]','fontsize',12)
  title(['Intensity plot for LP_',num2str(L),'_',num2str(M),'mode',' (b=',num2str(broots(M)),')'],'fontsize',15)
  
  subplot(2,1,2);
  box on;
  imagesc(xv,yv,angle(Emat_field{M}.'));
  colormap(gca,hsv);
  axis equal tight
  caxis([-pi pi])
  c = colorbar;
  c.Label.String = 'Phase [rad]';
  set(c,'fontsize',15)
  xlabel('x [m]','fontsize',12)
  ylabel('y [m]','fontsize',12)
  title(['Phase plot for LP_',num2str(L),'_',num2str(M),'mode',' (b=',num2str(broots(M)),')'],'fontsize',15)
end

if Modaloverlap
  % Calculation of the incident field.
  % In the simplest case, the incident optical field is approximated as an
  % Gaussian function at z=0 (proximal end of optical fiber).
  w0 = r_core.*(0.65+(1.619./V^(3/2))+(2.879/V^6)); % beam waist of the incident Gaussian optical field,
  % optimized by Marcuse's formula.
  
  phase = 0;
  Emat_gauss = exp(-(X.^2 + Y.^2)/w0^2).*exp(1i*phase);
  Enorm_gauss = Emat_gauss./sqrt(sum(abs(Emat_gauss(:)).^2));
  
  % plotting the incident intensity across the grid.
  figure(4)
  imagesc(xv,yv,abs(Enorm_gauss.').^2);
  axis equal tight
  colormap(jet((2^8)));
  c = colorbar;
  c.Label.String = 'Intensity [Arb.]';
  set(c,'fontsize',15)
  xlabel('x [m]','fontsize',12)
  ylabel('y [m]','fontsize',12)
  title('Incident optical intensity','fontsize',15)
  
  % Calculating the overlap integral between the incident field and the
  % selected modal field.
  Cmat = cell(length(broots),1);
  Cmat_abs = cell(length(broots),1);
  for M=1:length(broots)
    Enorm_mode = Emat_field{M}./sqrt(sum(abs(Emat_field{M}(:)).^2));
    Cmat{M} = sum(Enorm_gauss(:).*conj(Enorm_mode(:)));
    Cmat_abs{M} = abs(Cmat{M}).^2;
  end
  Cmat
  Cmat_abs  

% Bar plot of overlap coefficients and fractional power
figure(5) % values of overlap coefficients 
bar([1:length(broots)],[Cmat{:}])
xlabel('Modal number m','fontsize',12)
ylabel('Overlap coefficient','fontsize',12)

figure(6) %  values of absolute squared coefficients (power)
bar([1:length(broots)],[Cmat_abs{:}])
xlabel('Modal number m','fontsize',12)
ylabel('Fractional power','fontsize',12)
title('Launching efficiency \eta','fontsize',15)

end



if SaveData
  save(dataName,'Emat_field','Emat_gauss','Cmat','Cmat_abs')
end

toc;
%% END %%