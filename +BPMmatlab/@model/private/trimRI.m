function n = trimRI(n,n_background)
n_background = single(n_background);
[Nx,Ny,~] = size(n);
xmin = find(any(n ~= n_background,[2 3]),1,'first');
xmax = find(any(n ~= n_background,[2 3]),1,'last');
xtrim = min(xmin-1,Nx - xmax);
ymin = find(any(n ~= n_background,[1 3]),1,'first');
ymax = find(any(n ~= n_background,[1 3]),1,'last');
ytrim = min(ymin-1,Ny - ymax);

n_temp = n(xtrim+1:Nx-xtrim,ytrim+1:Ny-ytrim,:);
if ndims(n_temp) == 3
  n = n_background*ones(size(n_temp) + [2 2 0], class(n_temp));
else
  n = n_background*ones(size(n_temp) + [2 2  ], class(n_temp));
end
n(2:end-1, 2:end-1, :) = n_temp;
end