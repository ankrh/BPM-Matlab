function coreIdxs = findCores(n,n_background)
A = single(n) ~= single(n_background);

[Nx,Ny] = size(A);

B = zeros(Nx,Ny);
n = 0;
for ix=1:Nx
  on = false;
  for iy=1:Ny
    if ~on && A(ix,iy)
      n = n + 1;
      B(ix,iy) = n;
      on = true;
    elseif on && A(ix,iy)
      B(ix,iy) = n;
    elseif on && ~A(ix,iy)
      on = false;
    end
  end
end

for iy=1:Ny
  for ix=2:Nx
    if B(ix,iy) && B(ix-1,iy)
      B(B == B(ix,iy)) = B(ix-1,iy);
    end
  end
end

[~,~,ic] = unique(B);

coreIdxs = reshape(ic-1,size(B));
end
