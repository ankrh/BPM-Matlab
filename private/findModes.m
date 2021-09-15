function P = findModes(P,nModes,varargin)
singleCoreModes = false;
sortByLoss = false;
plotModes = true;
for i=1:2:numel(varargin)
  switch varargin{i}
    case 'singleCoreModes'
      singleCoreModes = varargin{i+1};
    case 'sortByLoss'
      sortByLoss = varargin{i+1};
    case 'plotModes'
      plotModes = varargin{i+1};
    otherwise
      error('Argument #%d for findModes not recognized. The syntax for findModes has recently changed. See the examples and the readme.',i+2);
  end
end

if isfield(P,'n_cladding')
  error('Error: n_cladding has been renamed n_background');
end
if isfield(P,'shapes')
  error('The P.shapes field has been deprecated. Use the P.n field to define the refractive index instead, as shown in the example files.');
end
if ~isfield(P,'n')
  error('You must specify the P.n field');
end
if isa(P.n,'function_handle')
  error('The method of defining a refractive index function has changed. The function handle must now be stored in P.n.func instead of P.n.');
end
if ~isfield(P.n,'n') && ~isfield(P.n,'func')
  error('P.n must contain either the field "func" or "n"');
end
if ~isfield(P.n,'func') && ~(isfield(P.n,'Lx') && isfield(P.n,'Ly'))
  error('You must specify the side lengths Lx and Ly if you provide an array (2D or 3D) for the refractive index');
end
if isfield(P.n,'func') && nargin(P.n.func) == 5 && ~isfield(P.n,'Nz')
  error('You must specify the refractive index array z resolution P.n.Nz if you provide a 3D refractive index function');
end
if ~isfield(P,'nParameters')
  P.nParameters = {};
end
if ~isfield(P,'rho_e')
  P.rho_e = 0.22;
end
if ~isfield(P,'bendingRoC')
  P.bendingRoC = Inf;
end
if ~isfield(P,'bendDirection')
  P.bendDirection = 0;
end
if ~isfield(P,'Intensity_colormap')
  P.Intensity_colormap = 1;
end
if ~isfield(P,'Phase_colormap')
  P.Phase_colormap = 2;
end

if isfield(P,'modes')
  P = rmfield(P,'modes');
end

k_0 = 2*pi/P.lambda;  % [m^-1] Wavenumber
dx = P.Lx_main/P.Nx_main;
dy = P.Ly_main/P.Ny_main;

targetLx = P.padfactor*P.Lx_main;
targetLy = P.padfactor*P.Ly_main;

Nx = round(targetLx/dx);
if rem(Nx,2) ~= rem(P.Nx_main,2)
  Nx = Nx + 1; % Ensure that if Nx_main was set odd (to have a x slice at the center), Nx will also be odd
end
Ny = round(targetLy/dy);
if rem(Ny,2) ~= rem(P.Ny_main,2)
  Ny = Ny + 1; % Ensure that if Ny_main was set odd (to have a y slice at the center), Ny will also be odd
end
N = Nx*Ny;                                                            %N*N - size of sparse matrices
if nModes >= N - 1
  error('Error: The number of modes requested must be less than the pixels in the full simulation window minus one (roughly Nx_main*padfactor*Ny_main*padfactor - 1)');
end
Lx = Nx*dx;
Ly = Ny*dy;

x = dx*(-(Nx-1)/2:(Nx-1)/2);
y = dy*(-(Ny-1)/2:(Ny-1)/2);
[X,Y] = ndgrid(x,y);

V = [];
D = [];

fprintf('Finding modes...\n');
tic

if isfield(P.n,'func') % If P.n has a function field
  if nargin(P.n.func) == 4 % It's 2D
    n = single(P.n.func(X,Y,P.n_background,P.nParameters));
  else % It's 3D. We evaluate it on the grid for z = 0
    n = single(P.n.func(X,Y,zeros(size(X)),P.n_background,P.nParameters));
  end
else % Otherwise a P.n.n array must have been specified
  [Nx_source,Ny_source,Nz_source] = size(P.n.n);
  dx_source = P.n.Lx/Nx_source;
  dy_source = P.n.Ly/Ny_source;
  x_source = dx_source*(-(Nx_source-1)/2:(Nx_source-1)/2);
  y_source = dy_source*(-(Ny_source-1)/2:(Ny_source-1)/2);
  if Nz_source == 1 % If P.n.n is 2D, interpolate it to the simulation grid
    n = interpn(x_source,y_source,P.n.n,x,y.','linear',P.n_background);
  else % Otherwise it's 3D so we take the first slice and interpolate to the simulation grid
    n = interpn(x_source,y_source,P.n.n(:,:,1),x,y.','linear',P.n_background);
  end
end

anycomplex = ~isreal(n);
  
if singleCoreModes
  coreIdxs = findCores(n,P.n_background);
else
  coreIdxs = ones(size(n));
end
nCores = max(coreIdxs(:));

iMode = 0;
for iCore = 1:nCores
  n_core = n;
  n_core(coreIdxs ~= iCore) = P.n_background;
  
  if ~isfinite(P.bendingRoC)
    [radiallySymmetric,xC,yC] = testRadialSymmetry(X,Y,n_core,P.n_background); % Numerically estimates whether this core is radially symmetric, and finds the centroid coordinates (xC, yC)
  else
    radiallySymmetric = false;
  end
  
  n_core = double(n_core); % Needed for sparse operations
  
  n_core_bent = real(n_core).*(1-(real(n_core).^2.*(X*cosd(P.bendDirection) + Y*sind(P.bendDirection))*P.rho_e/(2*P.bendingRoC))).*exp((X*cosd(P.bendDirection) + Y*sind(P.bendDirection))/P.bendingRoC);
  
  delta_n_2 = real(n_core_bent).^2 - P.n_0^2;                              %delta_n^2 in the expression for FD BPM

  dz = 1e-10;
  absorber = exp(-dz*(max(0,max(abs(Y) - P.Ly_main/2,abs(X) - P.Lx_main/2)).^2*P.alpha + 2*pi*imag(n_core)/P.lambda)); % First part is edge absorber, second is absorption from the imaginary part of the refractive indicres
  ax = 1.00001*dz/(dx^2*2i*k_0*P.n_0);
%   ax = dz/(dx^2*2i*k_0*P.n_0);
  ay = dz/(dy^2*2i*k_0*P.n_0);


  M_rhs = sparse(1:N,1:N,absorber(1:N) + delta_n_2(1:N)*dz*k_0/(2i*P.n_0),N,N) + ...
    sparse(1:N-1,2:N,[repmat([repmat(ax,1,Nx-1) 0],1,Ny-1) repmat(ax,1,Nx-1)],N,N) + ...
    sparse(2:N,1:N-1,[repmat([repmat(ax,1,Nx-1) 0],1,Ny-1) repmat(ax,1,Nx-1)],N,N) + ...
    sparse(1:N-Nx,1+Nx:N,ay,N,N) + ...
    sparse(1+Nx:N,1:N-Nx,ay,N,N);
  M_rhs(1:N+1:N*N) = M_rhs(1:N+1:N*N) - repmat([ax repmat(2*ax,1,Nx-2) ax],1,Ny);
  M_rhs(1:N+1:N*N) = M_rhs(1:N+1:N*N) - [repmat(ay,1,Nx) repmat(2*ay,1,Nx*(Ny-2)) repmat(ay,1,Nx)];
  absorberdiag = sparse(1:N,1:N,absorber(1:N),N,N);
  M_rhs = M_rhs*absorberdiag;

  [V,D] = eigs(M_rhs,ceil(nModes/nCores),1,'Display',false,'SubspaceDimension',min(N,ceil(nModes/nCores)*10));
  D = diag(D);
  
  kappa = (1-real(D))/dz*P.lambda/(4*pi);
  realn = sqrt(P.n_0^2 - 2*P.n_0*imag(D/(dz*k_0)));

  neff = realn(realn > P.n_background) + 1i*kappa((realn > P.n_background));
  V = V(:,realn > P.n_background);
  D = D(realn > P.n_background);
  
  for iCoreMode = 1:numel(D)
    iMode = iMode + 1;
    P.modes(iMode).Lx = Lx;
    P.modes(iMode).Ly = Ly;
    E = reshape(V(:,iCoreMode),[Nx Ny]);
    E = E.*exp(-1i*angle(max(E(:)))); % Shift phase arbitrarily so that the resulting modes are (nearly) real
    P.modes(iMode).field = E;
    if anycomplex
      P.modes(iMode).neff = neff(iCoreMode);
    else
      P.modes(iMode).neff = real(neff(iCoreMode)); % If the user did not specify any complex refractive indices, then the only contribution to kappa would be from the edge absorber, and showing a complex neff would just confuse people and not be very physically meaningful
    end
    if radiallySymmetric
      [~,iMax] = max(E(:));
      xMax = X(iMax);
      yMax = Y(iMax);
      theta = atan2(yMax - yC,xMax - xC);
      radialE = interpn(X,Y,E,xC + linspace(0,max(Lx,Ly),1000)*cos(theta),yC + linspace(0,max(Lx,Ly),1000)*sin(theta));
      radialEpruned = radialE(abs(radialE) > 0.1*max(abs(radialE)));
      m = sum(abs(diff(angle(radialEpruned*exp(1i*pi/2)) > 0))) + 1;

      R = sqrt((xMax - xC)^2 + (yMax - yC)^2);
      azimuthalE = interpn(X,Y,E,xC + R*cos(theta + linspace(0,2*pi,1000)),yC + R*sin(theta + linspace(0,2*pi,1000)));
      azimuthalEpruned = azimuthalE(abs(azimuthalE) > 0.1*max(abs(azimuthalE)));
      l = sum(abs(diff(angle(azimuthalEpruned*exp(1i*pi/2)) > 0)))/2;

      if l > 0
        Emaxmirrored = interpn(X,Y,E,xMax,2*yC - yMax);
        if real(E(iMax)/Emaxmirrored) < 0
          parity = 'o';
        else
          parity = 'e';
        end
      else
        parity = '';
      end
      P.modes(iMode).label = [', LP' num2str(l) num2str(m) parity];
    else
      P.modes(iMode).label = '';
    end
  end
end

if ~isfield(P,'modes')
  fprintf('\b Done, %.1f seconds elapsed.\nNo guided modes found.\n',toc);
  return
end

if sortByLoss
  [~,sortedidxs] = sort(imag([P.modes.neff]),'ascend');
else
  [~,sortedidxs] = sort(real([P.modes.neff]),'descend');
end
P.modes = P.modes(sortedidxs(1:min(numel(P.modes),nModes)));

fprintf('\b Done, %.1f seconds elapsed.\n%d guided modes found.\n',toc,numel(P.modes));

for iMode = 1:numel(P.modes)
  P.modes(iMode).label = ['Mode ' num2str(iMode) P.modes(iMode).label];
  if plotModes
    E = P.modes(iMode).field;
    h_f = figure(100+iMode);
    h_f.WindowStyle = 'docked';
    subplot(1,2,1);
    imagesc(x,y,abs(E.').^2);
    axis equal; axis tight; axis xy;
    setColormap(gca,P.Intensity_colormap);
    subplot(1,2,2);
    imagesc(x,y,angle(E.'),'AlphaData',max(0,(1+log10(abs(E.'/max(E(:))).^2)/3)));
    set(gca,'Color',0.7*[1 1 1]);  % To set the color corresponding to phase outside the cores where there is no field at all
    caxis([-pi pi]);
    axis equal; axis tight; axis xy;
    setColormap(gca,P.Phase_colormap);
    neff = P.modes(iMode).neff;
    if anycomplex
      sgtitle({[P.modes(iMode).label ', n_{eff} = ' num2str(real(neff),'%.6g') ' + ' num2str(imag(neff),'%.3g') 'i'],['rough loss estimate: ' num2str(imag(neff)*4*pi/P.lambda,'%.3g') ' m^{-1} (' num2str(-10*log10(exp(-1))*imag(neff)*4*pi/P.lambda,'%.3g') ' dB/m)']});
    else
      sgtitle({[P.modes(iMode).label ', n_{eff} = ' num2str(real(neff),'%.6g')],['rough loss estimate: ' num2str(imag(neff)*4*pi/P.lambda,'%.3g') ' m^{-1} (' num2str(-10*log10(exp(-1))*imag(neff)*4*pi/P.lambda,'%.3g') ' dB/m)']});
    end
  end
end
drawnow;
end

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