function [Estruct, P] = FFT_BPM(P)
%%
% Author: Madhu Veettikazhy
% Date: 12 July 2019
% ***************************************************************************************************

if ~isfield(P,'figNum')
  P.figNum = 1;
end
if ~isfield(P,'figTitle')
  P.figTitle = '';
end
if ~isfield(P,'Eparameters')
  P.Eparameters = {};
end
if ~isfield(P,'videoName')
  P.videoName = [P.name '.avi'];
end
if ~isfield(P,'Intensity_colormap')
  P.Intensity_colormap = 1;
end
if ~isfield(P,'Phase_colormap')
  P.Phase_colormap = 2;
end

%% Initialization of space and frequency grids
Nz = P.updates; % Number of z steps in this segment

targetLx = P.padfactor*P.Lx_main;
targetLy = P.padfactor*P.Ly_main;

dx = P.Lx_main/P.Nx_main;
dy = P.Ly_main/P.Ny_main;
dz = P.Lz/Nz;

Nx = round(targetLx/dx);
if rem(Nx,2) ~= rem(P.Nx_main,2)
  Nx = Nx + 1; % Ensure that if Nx_main was set odd (to have a x slice at the center), Nx will also be odd
end
Ny = round(targetLy/dy);
if rem(Ny,2) ~= rem(P.Ny_main,2)
  Ny = Ny + 1; % Ensure that if Ny_main was set odd (to have a y slice at the center), Ny will also be odd
end
Lx = Nx*dx;
Ly = Ny*dy;

x = dx*(-(Nx-1)/2:(Nx-1)/2);
y = dy*(-(Ny-1)/2:(Ny-1)/2);
[X,Y] = ndgrid(x,y);

kx = 2*pi/Lx*([0:floor((Nx-1)/2) floor(-(Nx-1)/2):-1]);
ky = 2*pi/Ly*([0:floor((Ny-1)/2) floor(-(Ny-1)/2):-1]);
[kX,kY] = ndgrid(kx,ky);

absorber = exp(-dz*max(0,max(abs(Y) - P.Ly_main/2,abs(X) - P.Lx_main/2)).^2*P.alpha);
prop_kernel = exp(1i*dz*(kX.^2+kY.^2)*P.lambda/(4*pi*P.n_0)); % Fresnel propagation kernel

%% Beam initialization
if isa(P.E,'function_handle')
  E = P.E(X,Y,P.Eparameters); % Call function to initialize E field
  P.E_0 = E; % To save the initial electric field
else % Interpolate source E field to new grid
  [Nx_Esource,Ny_Esource] = size(P.E.field);
  dx_Esource = P.E.Lx/Nx_Esource;
  dy_Esource = P.E.Ly/Ny_Esource;
  x_Esource = dx_Esource*(-(Nx_Esource-1)/2:(Nx_Esource-1)/2);
  y_Esource = dy_Esource*(-(Ny_Esource-1)/2:(Ny_Esource-1)/2);
  [X_source,Y_source] = ndgrid(x_Esource,y_Esource);
  E = interp2(X_source.',Y_source.',P.E.field.',X.',Y.','linear',0).';
end

%% Fresnel Propagation and plotting
if P.saveVideo && ~isfield(P,'videoHandle')
  video = VideoWriter(P.videoName);  %For saving the propagation frames as a video
  open(video);
elseif P.saveVideo
  video = P.videoHandle;  %videoHandle is passed from Example.m file
end

h_f = figure(P.figNum);clf;
h_f.WindowState = 'maximized';

h_ax1 = subplot(1,2,1);
h_im_I = imagesc(x,y,abs(E.').^2);
if isfield(P,'plotEmax')
  caxis('manual');
  caxis([0 P.plotEmax]);
else
  caxis('auto');
end
setColormap(gca,P.Intensity_colormap);
h_imItitle = title({'Intensity profile',' at z = 0 m'});
h_imItitle.FontSize = 20;
axis xy
axis equal
xlim([-1 1]*Lx/(2*P.displayScaling));
ylim([-1 1]*Ly/(2*P.displayScaling));
line([-P.Lx_main P.Lx_main P.Lx_main -P.Lx_main -P.Lx_main]/2,[P.Ly_main P.Ly_main -P.Ly_main -P.Ly_main P.Ly_main]/2,'color','r','linestyle','--');
colorbar;

h_ax2 = subplot(1,2,2);
h_im_phi = imagesc(x,y,angle(E.'));
h_im_phi.AlphaData = max(0,(1+log10(abs(E.'/max(abs(E(:)))).^2)/3));  %Logarithmic transparency in displaying phase outside cores
h_ax2.Color = 0.7*[1 1 1];
setColormap(gca,P.Phase_colormap);
h_imPhiTitle = title({'Phase profile',' at z = 0 m'});
h_imPhiTitle.FontSize = 20;
axis xy
axis equal
xlim([-1 1]*Lx/(2*P.displayScaling));
ylim([-1 1]*Ly/(2*P.displayScaling));
line([-P.Lx_main P.Lx_main P.Lx_main -P.Lx_main -P.Lx_main]/2,[P.Ly_main P.Ly_main -P.Ly_main -P.Ly_main P.Ly_main]/2,'color','r','linestyle','--');
colorbar;

if P.saveVideo
  frame = getframe(h_f);                                     %Get the frames
  writeVideo(video,frame);                                  %Stitch the frames to form a video and save
end

sgtitle(P.figTitle,'FontSize',20,'FontWeight','bold');

updatesliceindices = unique(round(linspace(1,Nz,min(Nz,P.updates))));
nextupdatesliceindicesindex = 1;

% tic
for zidx = 1:Nz
  E = absorber.*ifft2(fft2(E).*prop_kernel);
  
  if zidx == updatesliceindices(nextupdatesliceindicesindex)
    h_im_I.CData = abs(E.').^2;
    if ~isfield(P,'plotEmax')
      caxis(h_ax1,'auto'); %  To refresh the numbers on the color bar
    end
    h_imItitle.String = {'Intensity profile';['at z = ' num2str(zidx*dz,'%.1e') ' m']};
    h_im_phi.CData = angle(E.');
    h_im_phi.AlphaData = max(0,(1+log10(abs(E.'/max(abs(E(:)))).^2)/3));  %Logarithmic transparency in displaying phase outside cores
    h_imPhiTitle.String = {'Phase profile';['at z = ' num2str(zidx*dz,'%.1e') ' m']};
    
    nextupdatesliceindicesindex = nextupdatesliceindicesindex + 1;
    drawnow;
    if P.saveVideo
      frame = getframe(h_f);                                     %Get the frames
      writeVideo(video,frame);                                  %Stitch the frames to form a video and save
    end
  end
end
% toc

if P.saveVideo && ~isfield(P,'videoHandle')
  close(video);
end

Estruct = struct('field',E,'Lx',Lx,'Ly',Ly,'x',x,'y',y,'dx',dx,'dy',dy,'dz',dz);

% S = load('train');
% sound(S.y.*0.1,S.Fs);

function setColormap(gca,colormapType)
    switch colormapType
     case 1
        colormap(gca,GPBGYRcolormap);
     case 2
        colormap(gca,hsv/1.5);
     case 3
        colormap(gca,parula);
     case 4
        colormap(gca,gray);
     case 5
        colormap(gca,cividisColormap);
    end
end
end