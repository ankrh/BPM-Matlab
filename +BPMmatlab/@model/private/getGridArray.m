function x = getGridArray(Nx,dx,symmetry)
switch symmetry
  case 'NoSymmetry'
    x = dx*(-Nx/2+1/2:Nx/2-1/2);
  case 'Symmetry'
    x = dx*(1/2:Nx-1/2);
  case 'AntiSymmetry'
    x = dx*(0:Nx-1);
end
end