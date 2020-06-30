%% Calculation of mode decomposition and tranverse profiles %%
% This file is working example of how to use the function package developed
% during a special course project " Ultrafast pulse fibre delivery for
% Lightsheet microscopy".

%% Overview of the script and function calls %%
% This Matlab script is written in 6 separate sections. 
% 1) First a short preamble section, handling clearing of workspace and command window etc.
% Here visualization of loaded data and calculations is also enabled.
%
% 2) The parameter section, definition of all the fiber parameters for the
% calculations, definition of grid parameters and beam parameters such as
% wavelength etc. Finally the definition of the values considered for 
% azimuthal mode number, the normalized propagation constant, for which
% roots is found and V-number. Including approximate number of modes supported.  
%
% 3) Function call of Transcendentaleq_fiber.m which calculates the
% transcendental equation for each azimuthal mode number supported by the
% fibre parameters and its roots. See Transcendentaleq_fiber.m for more information. 
%
% 4) Function call of LPmodeCalc.m which calculates all the LPlm (A) modes 
% (cosine) and LPlm (B) modes (sine) for each azimuthal mode number and
% determined roots. See LPmodeCalc.m for more information. 
% 
% 5) Function call of ModeOverlap function which calculates the mode
% decomposition of loaded electric field data into all input modes. See 
% ModeOverlap.m for more information. A validation test is displayed in the
% command window of the mode overlap calculation against the total power.
% 
% 6) The final section, is a plotting section which is enabled in line 38.
% See line 108 for more information about the particular plot outputs. 
% 
%% Updated: 21-06-2020

%% Preamble %%
clear all; close all; clc; tic;

LHSRHSplot=true; Modeplot=true;  Dataplot=true; Overlapplot=true;

%% Parameter section %%
% Initialization of the fiber
r = 25e-6; % Core radius.
n_core = 1.4666; % Core index.
n_clad = 1.45; % Cladding index.
Delta = (n_core-n_clad)/n_clad; % Index constrast, check for weakly guiding approximation.

% Initialization of the computational grid. 
N=600; %Resolution.
Lx = 135e-6; % Grid length in x-direction.
Ly = 135e-6; % Grid length in y-direction.
xv=linspace(-Lx/2,Lx/2,N); % Defining the x-direction;
yv=linspace(-Ly/2,Ly/2,N); % Defining the y-direction;
[X,Y]=ndgrid(xv,yv); % generating the computational mesh (x,y)-grid. 

% Beam parameters
lambda=1050e-9; % Excitation wavelength. 
k0=2*pi/lambda; % Wavenumber for excitaiton wavelength.
E0=1; % Initial amplitude. 

% Initialization of the b-V curves and subsequently field calculation. 
N_b = 10^6; % Resolution of normalized propagation constant axis.  
bv = (1:N_b)./(N_b+1); % Normalized propagation constant
V = k0.*r.*sqrt(n_core^2-n_clad^2); % Normalized frequency
N_mode=V^2/2; % approximate number of modes supported by the fiber. 


%% The transcendentaleq_fiber function %%
% In this section we calculate the trancendental equations for each
% specified azimuthal mode number (L) and find the proper roots, 
% corresponding to the radial solutions for each azimuthal mode number.
% The transcendental equation can be plotted, set LHSRHSplot to true, in
% line 16. See line 60-79 for the lines of code. 
[f,LHS_mat,RHS_mat,broots_mat] = Transcendentaleq_fiber(V,bv);
NL = length(broots_mat)
Nm = length(broots_mat{1})

%% The LPmodeCalc function %%
% In this section the total number of modes is calculated on the basis of
% fibre parameters, beam parameters and determined roots from
% Transcendentaleq_fiber.m. 
intNorm=false; 
[Ecell_modesAB] = LPmodeCalc(r,n_core,n_clad,N,X,Y,k0,E0,broots_mat,intNorm);
% Allocating calculated LPlm (A) modes (cosine) and LPlm (B) modes (sine)
% to different cell structures for mode overlap calculations. 
Ecell_modeA=Ecell_modesAB{1,1};
Ecell_modeB=Ecell_modesAB{2,1};

%% The ModeOverlap function %%
% In this section the mode overlap is calculated for each mode determined
% by LPmodeCalc.m. 
% Loading calculated distal end electric field data from a FD_BPM script.
load('Video_FG050LGA_1050nm_1cm_3cmy_LP09','E') % Loading data set for mode overlap calculation. 
Ecell_data=E; % 
[Ccell_outA] = ModeOverlap(Ecell_data,Ecell_modeA);
[Ccell_outB] = ModeOverlap(Ecell_data,Ecell_modeB);
% Allocating calculated mode overlap integrals, for LPlm (A) modes (cosine) and LPlm
% (B) modes (sine) respectively, to different cell structures for visualization. 
Cmat_absA=Ccell_outA{2,1};
Cmat_absB=Ccell_outB{2,1};
% Testing valditity of the mode overlap calculation against the total power. 
fprintf('The sum of all the modes'' fractional powers (total power) is %f\n',sum([Cmat_absA(~isnan(Cmat_absA)) ; Cmat_absB(~isnan(Cmat_absB))]));

%% Plotting section
% In this section plots of determined roots and transcendental equations
% from Transcendentaleq_fiber.m is selected by LHSRHSplot=true in line 38.
% See line 88-106.
% Furthermore plots of calculated LPlm (A) modes (cosine) and LPlm (B) modes (sine)
% from LPmodeCalc.m is enabled in line 16 by Modeplot=true. See line
% 108-133. In addition the loaded data can be displayed by Dataplot=true. 
% See line 138-157. 
% Finally the mode decomposition of loaded data into LPlm (A) modes 
% (cosine) and LPlm (B) modes (sine) respectively, can be visualized by 
%  Overlapplot=true. See line 161-179. 

if LHSRHSplot
  idx_fbplot=2; % If LHSRHSplot is true, choose a specific trancendental equation to plot.
  % plotting the transcendental equation to depict intersections of LHS and
  % RHS.
  plotylims = [-50 50];
  LHSplot=LHS_mat{idx_fbplot}; RHSplot=RHS_mat{idx_fbplot};
  fplot=LHSplot-RHSplot;
  figure(1)
  plot(bv,LHSplot,bv,RHSplot,bv,fplot,[0 1],[0 0],'k','LineWidth',2)
  grid on; grid minor; box on;
  ylim(plotylims)
  xlabel('Normalized propagation constant b')
  legend('LHS','RHS','LHS - RHS','Location','northeast','AutoUpdate','off')
  title(['Transcendental equation for L = ' num2str(idx_fbplot-1)]);
  brootplot=broots_mat{idx_fbplot};
  for j6=1:length(brootplot)
    line([brootplot(j6) brootplot(j6)],[plotylims(1) 0],'linestyle','--','color','k')
  end
end

if Modeplot
  idx_L=1; % If Modeplot is true, choose a specific azimuthal mode number to plot
  idx_m=9; % If Modeplot is true, choose a specific m mode number to plot.
  % In line 96, call either Ecell_modeA or Ecell_modeB for cosine (A) or
  % sin-dependent modes (B).
  mode_plot=Ecell_modeA{idx_L,idx_m}; bvplot=broots_mat{idx_L,1}; b=bvplot(idx_m);
  
  figure(2)
  subplot(1,2,1)
  imagesc(xv,yv,abs(mode_plot.').^2);
  axis equal tight
  colormap(jet((2^8)));
  c = colorbar;
  c.Label.String = 'Intensity [Arb.]';
  set(c,'fontsize',15)
  xlabel('x [m]','fontsize',12)
  ylabel('y [m]','fontsize',12)
  subplot(1,2,2)
  box on;
  imagesc(xv,yv,angle(mode_plot.'));
  colormap(gca,hsv);
  axis equal tight
  caxis([-pi pi])
  c = colorbar;
  c.Label.String = 'Phase [rad]';
  set(c,'fontsize',15)
  xlabel('x [m]','fontsize',12)
  ylabel('y [m]','fontsize',12)
  sgt = sgtitle(['Intensity plot for LP_',num2str(idx_L-1),'_',num2str(idx_m),'mode','(b=',num2str(b),')'],'fontsize',20);
end

if Dataplot
  figure(3)
  subplot(1,2,1)
  imagesc(xv,yv,abs(E.').^2);
  axis equal tight
  colormap(jet((2^8)));
  c = colorbar;
  c.Label.String = 'Intensity [Arb.]';
  set(c,'fontsize',15)
  xlabel('x [m]','fontsize',12)
  ylabel('y [m]','fontsize',12)
  subplot(1,2,2)
  meshz(X,Y,abs(E.').^2)
  axis tight
  colormap(jet((2^8)))
  c = colorbar;
  c.Label.String = 'Intensity [Arb.]';
  set(c,'fontsize',15)
  sgt = sgtitle('Loaded distal end electric field','fontsize',20);
end

if Overlapplot
  figure(4);clf;
  [Lmat,Mmat] = meshgrid(-length(broots_mat)+1:length(broots_mat)-1,1:length(broots_mat{1}));
  data = [flipud(Cmat_absB(2:end,:)) ; Cmat_absA].';
  scatterbar3(Mmat,Lmat,data,1);
  axis tight;
  h = get(gca,'DataAspectRatio');
  if h(3)==1
    set(gca,'DataAspectRatio',[1 1 1/max(h(1:2))]);
  else
    set(gca,'DataAspectRatio',[1 1 h(3)]);
  end
  view([-90 90]); grid on;
  colormap(jet((2^8)));
  set(gca,'CLim',[0 1]);
  h_colbar = colorbar('southoutside');
  h_colbar.Label.String = 'Fractional power';
  title('LP mode decomposition - 1050nm - 1cm - 3cm(y) - LP09')
  xlabel('m')
  ylabel('L (negative values are the odd modes)')
  zlabel('Fractional power')
end

%%
toc;

