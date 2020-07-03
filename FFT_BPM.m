function [Estruct] = FFT_BPM(P)
%%
% Author: Madhu Veettikazhy
% Date: 12 July 2019
% ***************************************************************************************************

videoName = [P.name '.avi'];

%% Initialization of space and frequency grids
Nz = round(Lz/dz_target);  % [] Number of z steps

targetLx = P.padfactor*P.Lx_main;
targetLy = P.padfactor*P.Ly_main;

dx = Lx_main/Nx_main;
dy = Ly_main/Ny_main;
dz = Lz/Nz;

Nx = round(targetLx/dx);
if rem(Nx,2) ~= rem(Nx_main,2)
  Nx = Nx + 1; % Ensure that if Nx_main was set odd (to have a x slice at the center), Nx will also be odd
end
Ny = round(targetLy/dy);
if rem(Ny,2) ~= rem(Ny_main,2)
  Ny = Ny + 1; % Ensure that if Ny_main was set odd (to have a y slice at the center), Ny will also be odd
end
Lx = Nx*dx;
Ly = Ny*dy;

x = dx*(-(Nx-1)/2:(Nx-1)/2);
y = dy*(-(Ny-1)/2:(Ny-1)/2);
[X_fft,Y_fft] = ndgrid(x,y);

kx = 2*pi/Lx*([0:floor((Nx-1)/2) floor(-(Nx-1)/2):-1]);
ky = 2*pi/Ly*([0:floor((Ny-1)/2) floor(-(Ny-1)/2):-1]);
[kX,kY] = ndgrid(kx,ky);

switch absorbertype
  case 1
    absorber = ones(Nx,Ny);
  case 2
    absorber = ones(Nx,Ny);
    absorber(abs(X_fft) > Lx_main/2 | abs(Y_fft) > Ly_main/2) = exp(-dz*alpha);
  case 3
    absorber = exp(-dz*max(0,max(abs(Y_fft) - Ly_main/2,abs(X_fft) - Lx_main/2))*beta);
  case 4
    absorber = exp(-dz*max(0,max(abs(Y_fft) - Ly_main/2,abs(X_fft) - Lx_main/2)).^2*gamma);
end

%% Beam initialization
% amplitude = exp(-(X.^2+Y.^2)/w^2);                             % Gaussian field amplitude
% phase = 0*X;                                         % Phase
% E_ref = interp2(X.',Y.',amplitude.*exp(1i*phase),X_fft.',Y_fft.','linear',0);                      % Electric field

E_ref = interp2(X.',Y.',E,X_fft.',Y_fft.','linear',0);                      % Electric field, E from the FD_BPM.m output

if P.saveVideo
  video = VideoWriter(videoName);                   %For saving the propagation frames as a video
  open(video);
end
%% Fresnel Propagation and plotting
prop_kernel = exp(1i*dz*(kX.^2+kY.^2)*lambda/(4*pi)); % Fresnel propagation kernel

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
xlim([-Lx/(2*displayScaling) Lx/(2*displayScaling)]);
ylim([-Ly/(2*displayScaling) Ly/(2*displayScaling)]);
colorbar;

subplot(1,2,2);
h_im_phi = imagesc(x,y,angle(E_ref.'));
colormap(gca,hsv/1.5);
h_imPhiTitle = title({'Phase profile',' at z = 0 m'});
h_imPhiTitle.FontSize = 20;
axis xy
axis equal
xlim([-Lx/(2*displayScaling) Lx/(2*displayScaling)]);
ylim([-Ly/(2*displayScaling) Ly/(2*displayScaling)]);
colorbar;

if P.saveVideo
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
    h_imItitle.String = {['Intensity profile'];['at z = ' num2str(zidx*dz,'%.1e') ' m']};
    h_im_phi.CData = angle(E_ref.');
    h_imPhiTitle.String = {['Phase profile'];['at z = ' num2str(zidx*dz,'%.1e') ' m']};
    
    nextupdatesliceindicesindex = nextupdatesliceindicesindex + 1;
    drawnow;
    if P.saveVideo
      frame = getframe(gcf);                                     %Get the frames
      writeVideo(video,frame);                                  %Stitch the frames to form a video and save
    end
  end
end
toc

if P.saveVideo
  close(video);
end

S = load('train');
sound(S.y.*0.1,S.Fs);