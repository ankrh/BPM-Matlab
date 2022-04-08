function P = initializeRIfromFunction(P,hFunc,varargin)
if ~isa(hFunc,"function_handle")
  error('Argument 2 is wrong type. Must be function handle.');
end
if numel(varargin)
  nParameters = varargin{1};
  if ~isa(nParameters,"cell")
    error('Argument 3 is wrong type. Must be cell array.');
  end
else
  nParameters = {};
end

if nargin(hFunc) == 4 % It's 2D
  [X,Y] = ndgrid(single(P.x),single(P.y));
  P.n.n = single(hFunc(X,Y,P.n_background,nParameters));
else % It's 3D. We evaluate it on the grid and trim the result
  dz_n = P.Lz/(P.n.Nz-1);
  z_n = dz_n*(0:P.n.Nz-1);
  [X,Y,Z_n] = ndgrid(single(x),single(y),single(z_n));
  n = single(hFunc(X,Y,Z_n,P.n_background,nParameters));
  clear X Y Z_n
  P.n.n = trimRI(n,P.n_background);
end

P.n.Lx = P.Lx;
P.n.Ly = P.Ly;
P.n.xSymmetry = P.xSymmetry;
P.n.ySymmetry = P.ySymmetry;
end

function n = trimRI(n,n_background)
n_background = single(n_background);
[Nx,Ny,~] = size(n);
xmin = find(any(n ~= n_background,[2 3]),1,'first');
xmax = find(any(n ~= n_background,[2 3]),1,'last');
xtrim = min(xmin-1,Nx - xmax);
ymin = find(any(n ~= n_background,[1 3]),1,'first');
ymax = find(any(n ~= n_background,[1 3]),1,'last');
ytrim = min(ymin-1,Ny - ymax);

n_temp = n(xtrim+1:Nx-xtrim,ytrim+1:Ny-ytrim,:);
n = n_background*ones(size(n_temp) + [2 2 0], class(n_temp));
n(2:end-1, 2:end-1, :) = n_temp;
end