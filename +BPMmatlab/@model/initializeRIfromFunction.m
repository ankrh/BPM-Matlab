function P = initializeRIfromFunction(P,hFunc,varargin)
if ~isa(hFunc,"function_handle")
  error('Argument 2 is the wrong type. It must be a function handle.');
end
if nargin(hFunc) ~= 4 && nargin(hFunc) ~= 5
  error('The function handle you provided takes neither 4 arguments (X,Y,n_background,nParameters) nor 5 arguments (X,Y,Z,n_background,nParameters). It must conform to one of these syntaxes.');
end

if numel(varargin) >= 1
  nParameters = varargin{1};
  if ~isa(nParameters,"cell")
    error('Argument 3 is the wrong type. It must be a cell array.');
  end
else
  nParameters = {};
end
if numel(varargin) >= 2
  Nz = varargin{2};
  if mod(Nz,1) || Nz <= 0
    error('Argument 4 is the wrong type. It must be a positive integer.');
  end
else
  Nz = 1;
end
if Nz > 1 && nargin(hFunc) == 4
  error('You have specified an Nz > 1, but the function handle provided only takes 4 arguments (X,Y,n_background,nParameters). It must take 5 arguments (X,Y,Z,n_background,nParameters).');
end

if nargin(hFunc) == 4 % It's 2D
  [X,Y] = ndgrid(single(P.x),single(P.y));
  P.n.n = single(hFunc(X,Y,P.n_background,nParameters));
else % It's 3D. We evaluate it on the grid and trim the result
  if Nz > 1
    dz_n = P.Lz/(Nz-1);
  else
    dz_n = 0;
  end
  z_n = dz_n*(0:Nz-1);
  [X,Y,Z_n] = ndgrid(single(P.x),single(P.y),single(z_n));
  n = single(hFunc(X,Y,Z_n,P.n_background,nParameters));
  clear X Y Z_n
  P.n.n = trimRI(n,P.n_background);
end

P.n.Lx = P.dx*size(P.n.n,1);
P.n.Ly = P.dy*size(P.n.n,2);
P.n.xSymmetry = P.xSymmetry;
P.n.ySymmetry = P.ySymmetry;
end