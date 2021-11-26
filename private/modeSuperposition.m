function result = modeSuperposition(P,modeIdxs,varargin)
if numel(varargin) == 0
  coeffs = ones(numel(modeIdxs),1);
else
  coeffs = varargin{1};
end
result = struct('Lx',P.modes(modeIdxs(1)).Lx,'Ly',P.modes(modeIdxs(1)).Ly,'field',zeros(size(P.modes(1).field)),'xSymmetry',P.modes(modeIdxs(1)).xSymmetry,'ySymmetry',P.modes(modeIdxs(1)).ySymmetry);
for modeIdx = 1:numel(modeIdxs)
  result.field = result.field + coeffs(modeIdx)*P.modes(modeIdxs(modeIdx)).field;
end
end