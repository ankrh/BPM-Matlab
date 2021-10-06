function P = tiltField(P,direction,angle)
[X,Y] = ndgrid(P.x,P.y);
k = P.n_0*2*pi/P.lambda*angle/180*pi; % Refractive index in wedge in between fiber segments assumed to be uniformly P.n_0
P.E.field = P.E.field.*exp(-1i*k*(cosd(direction).*X + sind(direction).*Y));
end