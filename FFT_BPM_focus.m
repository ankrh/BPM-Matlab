%%
% Author: Madhu Veettikazhy
% Date: 12 July 2019
% To propagate the beam from fibre distal end to the focus point 
% ***************************************************************************************************

close all;
saveFocusVideo = false; 
saveFocusData = false; 
focusVideoName = 'GRIN_collimatedBeam_air_2D'; 

%% User-specified general parameters
Lx_main = 300e-6;                                                     % [m] x side length of main area
Ly_main = 300e-6;                                                    % [m] y side length of main area
focus = 8e-3;
Lz = focus;                                                          % [m] z propagation distance
targetzstepsize = 400e-6;                                          % [m] z step size to aim for

Nx_main = 1000;                                                         % x resolution of main area
Ny_main = 1000;                                                             % y resolution of main area

lambda = 0.98e-6;                                                % [m] Wavelength 
w = 30e-6;                                                       % [m] Initial waist plane 1/e^2 radius of the gaussian beam

updates = 200;          % Number of times to update plot. If set higher than Nz (such as Inf), script will simply update every step.
displayScaling = 3;

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

E_ref = interp2(X.',Y.',E,X_fft.',Y_fft.','linear',0);                      % Electric field, E from the FD_BPM.m output

if saveFocusVideo
    video = VideoWriter(focusVideoName);                   %For saving the propagation frames as a video
    open(video);
end
%% Fresnel Propagation and plotting
prop_kernel = exp(1i*dz_air*(kX.^2+kY.^2)*lambda/(4*pi)); % Fresnel propagation kernel

figure(1);
figure1_Settings = [];
screenSize = get(0,'MonitorPositions');
figure1_Settings.monitorNumber = 1; % Main monitor number
figure1_Settings.size.figureHeight = screenSize(figure1_Settings.monitorNumber,4); % Figure height (pixels)
figure1_Settings.size.figureWidth = screenSize(figure1_Settings.monitorNumber,3); % Figure width (pixels)
set(gcf, 'Position',  [0, 0, figure1_Settings.size.figureWidth, figure1_Settings.size.figureHeight]);

subplot(1,2,1);
h_im_I = imagesc(x,y,abs(E_ref.').^2);
colormap(jet);
h_imItitle = title({'Intensity profile',' at z = 0 m'});
h_imItitle.FontSize = 20;
axis xy
axis equal
xlim([-Lx/displayScaling Lx/displayScaling]);
ylim([-Ly/displayScaling Ly/displayScaling]);
colorbar;

subplot(1,2,2);
h_im_phi = imagesc(x,y,angle(E_ref.'));
colormap(gca,hsv/1.5);
h_imPhiTitle = title({'Phase profile',' at z = 0 m'});
h_imPhiTitle.FontSize = 20;
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
    E_ref = absorber.*ifft2(fft2(E_ref).*prop_kernel);
    
    if zidx == updatesliceindices(nextupdatesliceindicesindex)
        h_im_I.CData = abs(E_ref.').^2;
        h_imItitle.String = {['Intensity profile'];['at z = ' num2str(zidx*dz_air,'%.1e') ' m']};
        h_im_phi.CData = angle(E_ref.');
        h_imPhiTitle.String = {['Phase profile'];['at z = ' num2str(zidx*dz_air,'%.1e') ' m']};

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

if saveFocusData
    E_air_focus = E_ref; 
%     save('FileName.mat','powers','Lz','z_updates',...
%             'E','n_mat','nx','ny','dx','dy','Lambda','x','y','X','Y','RHO','w_0','k_0','lambda','E_0','num_cores','E_air_focus');
save('Data_GRIN_lens_collimated_air.mat',...
            'E','Nx','Ny','dx','dy','x','y','E_air_focus');
end
    
S = load('train');
sound(S.y.*0.1,S.Fs);