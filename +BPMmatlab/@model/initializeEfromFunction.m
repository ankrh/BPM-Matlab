function P = initializeEfromFunction(P,hFunc,varargin)
if ~isa(hFunc,"function_handle")
  error('Argument 2 is of the wrong type. It must be a function handle.');
end
if numel(varargin)
  Eparameters = varargin{1};
  if ~isa(Eparameters,"cell")
    error('Argument 3 is wrong type. Must be cell array.');
  end
else
  Eparameters = {};
end

[X,Y] = ndgrid(single(P.x),single(P.y));

E = hFunc(X,Y,Eparameters); % Call function to initialize E field

powerFraction = 1/(1 + (P.xSymmetry ~= 0))/(1 + (P.ySymmetry ~= 0)); % How large a fraction of the total power we are simulating

P.E.field = E/sqrt(sum(abs(E(:)).^2)/powerFraction);

P.E.Lx = P.Lx;
P.E.Ly = P.Ly;
P.E.xSymmetry = P.xSymmetry;
P.E.ySymmetry = P.ySymmetry;
end