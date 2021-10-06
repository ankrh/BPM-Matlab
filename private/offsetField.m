function P = offsetField(P,direction,distance)
[X,Y] = ndgrid(P.x,P.y);
P.E.field = interpn(X + distance*cosd(direction),Y + distance*sind(direction),P.E.field,X,Y,'linear',0);
end