function P = clearData(P)
if isfield(P,'z')
  P = rmfield(P,'z');
end
if isfield(P,'powers')
  P = rmfield(P,'powers');
end
if isfield(P,'modeOverlaps')
  P = rmfield(P,'modeOverlaps');
end
if isfield(P,'originalEinput')
  P.E = P.originalEinput;
  P = rmfield(P,'originalEinput');
end
if isfield(P,'originalShapesInput')
  P.shapes = P.originalShapesInput;
  P = rmfield(P,'originalShapesInput');
end
end