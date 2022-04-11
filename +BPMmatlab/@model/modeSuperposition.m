function result = modeSuperposition(P,modeIdxs,varargin)
if numel(varargin) == 0
  coeffs = ones(numel(modeIdxs),1);
else
  coeffs = varargin{1};
end
result = BPMmatlab.electricFieldProfile;
result.Lx = P.modes(modeIdxs(1)).Lx;
result.Ly = P.modes(modeIdxs(1)).Ly;
result.field = zeros(size(P.modes(1).field));
result.xSymmetry = P.modes(modeIdxs(1)).xSymmetry;
result.ySymmetry = P.modes(modeIdxs(1)).ySymmetry;
for modeIdx = 1:numel(modeIdxs)
  result.field = result.field + coeffs(modeIdx)*P.modes(modeIdxs(modeIdx)).field;
end
end