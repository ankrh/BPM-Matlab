%%
% Author: Madhu Veettikazhy
% Date: 12 July 2019
% To propagate the beam from fibre distal end to the focus point 
% ***************************************************************************************************

% clear all; 
close all;
saveFocusVideo = false; 
saveFocusData = false; 
% FileName = 'MCF6_twist_8,2mm_1000nm_APC';
% focusVideoName = ['Data_' today '/' FileName '.avi'];
% load([FileName '.mat']');
focusVideoName = 'Focus5mm_MM_fewModes15_05mm_noBend_APC.avi';
% load('Data_MC_Fermats40_testPointMove_noBend_phased.mat');
% E=E_MCfermats_40cores_noBend;

%% User-specified general parameters
Lx_main = 800e-6;                                                     % [m] x side length of main area
Ly_main = 800e-6;                                                    % [m] y side length of main area
focus = 5e-3;
Lz = focus;                                                          % [m] z propagation distance
targetzstepsize = 10e-6;                                          % [m] z step size to aim for

Nx_main = 400;                                                         % x resolution of main area
Ny_main = 400;                                                             % y resolution of main area

lambda = 0.98e-6;                                                % [m] Wavelength 
w = 30e-6;                                                       % [m] Initial waist plane 1/e^2 radius of the gaussian beam

updates = 200;          % Number of times to update plot. If set higher than Nz (such as Inf), script will simply update every step.
displayScaling = 2;

absorbertype = 4; % 1: No absorber, 2: constant absorber, 3: linear absorber, 4: quadratic absorber
targetLx = 3*Lx_main; % [m] Full area x side length
targetLy = 3*Ly_main;    % [m] Full area y side length
alpha = 10e1; % [1/m] "Absorption coefficient", used for absorbertype 2
beta = 10e4; % [1/m^2] "Absorption coefficient" per unit length distance out from edge of main area, used for absorbertype 3
gamma = 3e8; % [1/m^3] "Absorption coefficient" per unit length distance out from edge of main area, squared, used for absorbertype 4

%% Initialization of space and frequency grids
Nz = round(Lz/targetzstepsize);                                       % Number of z steps

dx_air = Lx_main/Nx_main;
dy_air = Ly_main/Ny_main;
dz_air = Lz/Nz;

Nx = round(targetLx/dx_air);
if rem(Nx,2) ~= rem(Nx_main,2)
    Nx = Nx + 1; % Ensure that if Nx_main was set odd (to have a x slice at the center), Nx will also be odd
end
Ny = round(targetLy/dy_air);
if rem(Ny,2) ~= rem(Ny_main,2)
    Ny = Ny + 1; % Ensure that if Ny_main was set odd (to have a y slice at the center), Ny will also be odd
end
Lx = Nx*dx_air;
Ly = Ny*dy_air;

x = dx_air*(-(Nx-1)/2:(Nx-1)/2);
y = dy_air*(-(Ny-1)/2:(Ny-1)/2);
[X_fft,Y_fft] = ndgrid(x,y);

kx = 2*pi/Lx*([0:floor((Nx-1)/2) floor(-(Nx-1)/2):-1]);
ky = 2*pi/Ly*([0:floor((Ny-1)/2) floor(-(Ny-1)/2):-1]);
[kX,kY] = ndgrid(kx,ky);

switch absorbertype
    case 1
        absorber = ones(Nx,Ny);
    case 2
        absorber = ones(Nx,Ny);
        absorber(abs(X_fft) > Lx_main/2 | abs(Y_fft) > Ly_main/2) = exp(-dz_air*alpha);
    case 3
        absorber = exp(-dz_air*max(0,max(abs(Y_fft) - Ly_main/2,abs(X_fft) - Lx_main/2))*beta);
    case 4
        absorber = exp(-dz_air*max(0,max(abs(Y_fft) - Ly_main/2,abs(X_fft) - Lx_main/2)).^2*gamma);
end

%% Beam initialization
% amplitude = exp(-(X.^2+Y.^2)/w^2);                             % Gaussian field amplitude
% phase = 0*X;                                         % Phase
% E_ref = interp2(X.',Y.',amplitude.*exp(1i*phase),X_fft.',Y_fft.','linear',0);                      % Electric field

E_ref = interp2(X.',Y.',E,X_fft.',Y_fft.','linear',0);                      % Electric field

if saveFocusVideo
    video = VideoWriter(focusVideoName);                   %For saving the propagation frames as a video
    open(video);
end
%% Fresnel Propagation and plotting
prop_kernel = exp(-1i*dz_air*(kX.^2+kY.^2)*lambda/(4*pi)); % Fresnel propagation kernel

figure(1);
figure1_Settings = [];
screenSize = get(0,'MonitorPositions');
figure1_Settings.monitorNumber = 1; % Main monitor number
figure1_Settings.size.figureHeight = screenSize(figure1_Settings.monitorNumber,4); % Figure height (pixels)
figure1_Settings.size.figureWidth = screenSize(figure1_Settings.monitorNumber,3); % Figure width (pixels)
set(gcf, 'Position',  [0, 0, figure1_Settings.size.figureWidth, figure1_Settings.size.figureHeight]);
h_im = imagesc(x,y,abs(E_ref.').^2);
colormap(jet);
% h_im1_xticklabel = get(gca,'XTickLabel');
% set(gca,'XTickLabel',h_im1_xticklabel,'FontSize',10,'fontweight','bold')
h_imtitle = title({'Intensity profile',' at z = 0 m'});
h_imtitle.FontSize = 20;
axis xy
axis equal
xlim([-Lx/displayScaling Lx/displayScaling]);
ylim([-Ly/displayScaling Ly/displayScaling]);
colorbar;

if saveFocusVideo
    frame = getframe(gcf);                                     %Get the frames 
    writeVideo(video,frame);                                  %Stitch the frames to form a video and save
end

updatesliceindices = unique(round(linspace(1,Nz,min(Nz,updates))));
nextupdatesliceindicesindex = 1;

tic
for zidx = 1:Nz
%     pause(0.05);
    E_ref = absorber.*ifft2(fft2(E_ref).*prop_kernel);
    
    if zidx == updatesliceindices(nextupdatesliceindicesindex)
        h_im.CData = abs(E_ref.').^2;
        h_imtitle.String = {['Intensity profile'];['at z = ' num2str(zidx*dz_air,'%.1e') ' m']};

        nextupdatesliceindicesindex = nextupdatesliceindicesindex + 1;
        drawnow;
        if saveFocusVideo
            frame = getframe(gcf);                                     %Get the frames 
            writeVideo(video,frame);                                  %Stitch the frames to form a video and save
        end
    end
end
toc

if saveFocusVideo
	close(video);
end
% figure('Renderer', 'painters', 'Position', [0 0 400 400]); imagesc(x,y,abs((E_ref(100:400,100:400)).').^2); axis tight equal; set(gca,'xtick',[]); set(gca,'ytick',[]);
if saveFocusData
    E_air_focus = E_ref; 
    save('Data_MM_fewModes15_05mm_noBend_APC_focus.mat','powers','Lz','z_updates',...
            'E','n_mat','nx','ny','dx','dy','Lambda','x','y','X','Y','RHO','w_0','k_0','lambda','E_0','num_cores','E_air_focus');
end
    
WarnWave = [sin(1:.6:400), sin(1:.7:400), sin(1:.4:400)];
Audio = audioplayer(WarnWave, 22050);
play(Audio);