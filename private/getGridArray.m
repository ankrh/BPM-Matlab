function x = getGridArray(Nx,dx,symmetry)
switch symmetry
  case 0
    x = dx*(-Nx/2+1/2:Nx/2-1/2);
  case 1
    x = dx*(1/2:Nx-1/2);
  case 2
    x = dx*(0:Nx-1);
end
end