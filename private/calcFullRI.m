function [x,y,n] = calcFullRI(x,y,n)
switch sign(min(x))
  case -1 % ySymmetry == 0, no symmetry
    % Do nothing
  case 0 % ySymmetry == 2, anti-symmetric
    x = cat(2,-flip(x(2:end  ),2),x);
    n = cat(1, flip(n(2:end,:),1),n);
  case 1 % ySymmetry == 1, symmetric
    x = cat(2,-flip(x         ,2),x);
    n = cat(1, flip(n         ,1),n);
end
switch sign(min(y))
  case -1 % xSymmetry == 0, no symmetry
    % Do nothing
  case 0 % xSymmetry == 2, anti-symmetric
    y = cat(2,-flip(y(  2:end),2),y);
    n = cat(2, flip(n(:,2:end),2),n);
  case 1 % xSymmetry == 1, symmetric
    y = cat(2,-flip(y         ,2),y);
    n = cat(2, flip(n         ,2),n);
end
end