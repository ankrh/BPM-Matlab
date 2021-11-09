function [x,y,E] = calcFullField(x,y,E)
if min(x) > 0 % ySymmetric
  x = cat(2,-flip(x,2),x);
  E = cat(1, flip(E,1),E);
end
if min(y) > 0 % xSymmetric
  y = cat(2,-flip(y,2),y);
  E = cat(2, flip(E,2),E);
end
end