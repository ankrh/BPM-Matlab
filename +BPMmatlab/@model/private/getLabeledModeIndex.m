function idx = getLabeledModeIndex(P,label)
idx = find(contains({P.modes.label},label),1);
if isempty(idx)
  error('Mode %s not found',label);
end
end