function P = tiltField(P,direction,angle)
if isa(P.E,'function_handle')
  error('Error: BPM-Matlab currently doesn''t support field tilts in segments where the field is specified as a function handle. You can avoid this by adding a short initial propagation segment and applying the field tilt between the first and second segments.');
end
[X,Y] = ndgrid(P.x,P.y);
k = P.n_0*2*pi/P.lambda*angle/180*pi; % Refractive index in wedge in between fiber segments assumed to be uniformly P.n_0
P.E.field = P.E.field.*exp(-1i*k*(cosd(direction).*X + sind(direction).*Y));
end