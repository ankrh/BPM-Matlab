function P = offsetField(P,direction,distance)
if isa(P.E,'function_handle')
  error('Error: BPM-Matlab currently doesn''t support field offsets in segments where the field is specified as a function handle. You can avoid this by adding a short initial propagation segment and applying the field offset between the first and second segments.');
end
[X,Y] = ndgrid(P.x,P.y);
P.E.field = interpn(X + distance*cosd(direction),Y + distance*sind(direction),P.E.field,X,Y,'linear',0);
end