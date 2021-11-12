function [x,y,E] = calcFullField(x,y,E)
switch sign(min(x))
  case -1 % ySymmetry == 0, no symmetry
    % Do nothing
  case 0 % ySymmetry == 2, anti-symmetric
    x = cat(2,-flip(x(2:end  ),2),x);
    E = cat(1,-flip(E(2:end,:),1),E);
  case 1 % ySymmetry == 1, symmetric
    x = cat(2,-flip(x,2),x);
    E = cat(1, flip(E,1),E);
end
switch sign(min(y))
  case -1 % xSymmetry == 0, no symmetry
    % Do nothing
  case 0 % xSymmetry == 2, anti-symmetric
    y = cat(2,-flip(y(  2:end),2),y);
    E = cat(2,-flip(E(:,2:end),1),E);
  case 1 % xSymmetry == 1, symmetric
    y = cat(2,-flip(y,2),y);
    E = cat(2, flip(E,2),E);
end
end