function P = FFT_BPM(P)
%% Initialization of space and frequency grids
P.dz_target = P.Lz/P.updates; % Set the step size so that we get the requested number of updates, one per step
Nz = P.Nz; % Number of z steps in this segment

dx = P.dx;
dy = P.dy;
dz = P.dz;
Nx = P.Nx;
Ny = P.Ny;
Lx = P.Lx;
Ly = P.Ly;

kx = 2*pi/Lx*([0:floor((Nx-1)/2) floor(-(Nx-1)/2):-1]);
ky = 2*pi/Ly*([0:floor((Ny-1)/2) floor(-(Ny-1)/2):-1]);
[kX,kY] = ndgrid(single(kx),single(ky));

prop_kernel = exp(1i*dz*(kX.^2+kY.^2)*P.lambda/(4*pi*P.n_0)); % Fresnel propagation kernel
clear kX kY

x = dx*(-(Nx-1)/2:(Nx-1)/2);
y = dy*(-(Ny-1)/2:(Ny-1)/2);
[X,Y] = ndgrid(single(x),single(y));

absorber = exp(-dz*max(0,max(abs(Y) - P.Ly_main/2,abs(X) - P.Lx_main/2)).^2*P.alpha);

if P.priorData
  P.z = [P.z dz*(1:Nz) + P.z(end)];
else
  P.z = [0 dz*(1:Nz)];
end

%% Beam initialization
[Nx_source,Ny_source] = size(P.E.field);
dx_source = P.E.Lx/Nx_source;
dy_source = P.E.Ly/Ny_source;
x_source = getGridArray(Nx_source,dx_source,P.E.ySymmetry);
y_source = getGridArray(Ny_source,dy_source,P.E.xSymmetry);
[x_source,y_source,E_source] = calcFullField(x_source,y_source,P.E.field);
E = interpn(x_source,y_source,E_source,x,y.','linear',0);
E = E*sqrt(sum(abs(E_source(:)).^2)/sum(abs(E(:)).^2));

clear X Y X_source Y_source

%% Fresnel Propagation and plotting
if P.saveVideo
  if isempty(P.videoHandle)
    P.videoHandle = VideoWriter([P.name '.avi']);  % If videoHandle is empty, we create it
    P.videoHandle.FrameRate = 5;
    open(P.videoHandle);
  end
end

h_f = figure(P.figNum);clf;
h_f.WindowState = 'maximized';

if P.xSymmetry ~= 0
  y0idx = 1;
else
  y0idx = round((Ny-1)/2+1);
end
if P.ySymmetry ~= 0
  x0idx = 1;
else
  x0idx = round((Nx-1)/2+1);
end
if P.priorData
  P.powers = [P.powers NaN(1,Nz)];
  P.xzSlice{end+1} = NaN(Nx,Nz);
  P.yzSlice{end+1} = NaN(Ny,Nz);
else
  P.powers = NaN(1,Nz+1);
  P.powers(1) = 1;
  P.xzSlice = {NaN(Nx,Nz+1)};
  P.yzSlice = {NaN(Ny,Nz+1)};
  P.xzSlice{1}(:,1) = E(:,y0idx);
  P.yzSlice{1}(:,1) = E(x0idx,:);
end
if P.storeE3D
  if P.priorData
    P.E3D{end+1} = complex(NaN(Nx,Ny,Nz,'single'));
  else
    P.E3D{1} = complex(NaN(Nx,Ny,Nz+1,'single'));
    P.E3D{1}(:,:,1) = E;
  end
end

xlims = [-1 1]*Lx/(2*P.plotZoom);
ylims = [-1 1]*Ly/(2*P.plotZoom);

if P.xSymmetry ~= 0 && P.ySymmetry ~= 0
  redline_x = [0 P.Lx_main P.Lx_main];
  redline_y = [P.Ly_main P.Ly_main 0];
elseif P.xSymmetry ~= 0
  redline_x = [-P.Lx_main -P.Lx_main P.Lx_main P.Lx_main]/2;
  redline_y = [0 P.Ly_main P.Ly_main 0];
elseif P.ySymmetry ~= 0
  redline_x = [0 P.Lx_main P.Lx_main 0];
  redline_y = [-P.Ly_main -P.Ly_main P.Ly_main P.Ly_main]/2;
else
  redline_x = [-P.Lx_main P.Lx_main P.Lx_main -P.Lx_main -P.Lx_main]/2;
  redline_y = [P.Ly_main P.Ly_main -P.Ly_main -P.Ly_main P.Ly_main]/2;
end

h_ax2 = subplot(2,2,2);
h_plot2 = plot(P.z,P.powers,'linewidth',2);
xlim([0 P.z(end)]);
ylim([0 1.1]);
xlabel('Propagation distance [m]');
ylabel('Relative power remaining');
grid on; grid minor;

h_axis3a = subplot(2,2,3);
hold on;
box on;
h_im3a = imagesc(x,y,abs(E.').^2);
axis xy;
axis equal;
xlim(xlims);
ylim(ylims);
colorbar;
xlabel('x [m]');
ylabel('y [m]');
title('Intensity [W/m^2]');
line(redline_x,redline_y,'color','r','linestyle','--');
setColormap(gca,P.intensityColormap);
if isfield(P,'plotEmax')
  caxis('manual');
  caxis([0 P.plotEmax*max(abs(E(:)).^2)]);
else
  caxis('auto');
end

h_axis3b = subplot(2,2,4);
hold on;
box on;
maxE0 = abs(max(E(:)));
h_im3b = imagesc(x,y,angle(E.'));
h_im3b.AlphaData = max(0,(1+log10(abs(E.'/maxE0).^2)/3));  %Logarithmic transparency in displaying phase outside cores
h_axis3b.Color = 0.7*[1 1 1];  % To set the color corresponding to phase outside the cores where there is no field at all
axis xy;
axis equal;
xlim(xlims);
ylim(ylims);
colorbar;
caxis([-pi pi]);
line(redline_x,redline_y,'color','r','linestyle','--');
xlabel('x [m]');
ylabel('y [m]');
title('Phase [rad]');
setColormap(gca,P.phaseColormap);

if P.saveVideo
  frame = getframe(h_f);                                     %Get the frames
  writeVideo(P.videoHandle,frame);                                  %Stitch the frames to form a video and save
end

if ~verLessThan('matlab','9.5')
  sgtitle(P.figTitle,'FontSize',20,'FontWeight','bold');
end

for zidx = 1:Nz
  E = absorber.*ifft2(fft2(E).*prop_kernel);

  if P.storeE3D
    P.E3D{end}(:,:,end-Nz+zidx) = E;
  end
  P.powers(end-Nz+zidx) = sum(abs(E(:).^2));
  h_plot2.YData = P.powers;
  h_ax2.YLim = [0 1.1*max(P.powers)];
  P.xzSlice{end}(:,end-Nz+zidx) = E(:,y0idx);
  P.yzSlice{end}(:,end-Nz+zidx) = E(x0idx,:);

  h_im3a.CData = abs(E.').^2;
  if ~isfield(P,'plotEmax')
    caxis(h_axis3a,'auto'); %  To refresh the numbers on the color bar
  end
  h_im3b.CData = angle(E.');
  h_im3b.AlphaData = max(0,(1+log10(abs(E.'/max(abs(E(:)))).^2)/3));  %Logarithmic transparency in displaying phase outside cores

  drawnow;
  if P.saveVideo
    frame = getframe(h_f);                                     %Get the frames
    writeVideo(P.videoHandle,frame);                                  %Stitch the frames to form a video and save
  end
end

Estruct = struct('field',E,'Lx',Lx,'Ly',Ly,'x',x,'y',y,'dx',dx,'dy',dy,'dz',dz);
end