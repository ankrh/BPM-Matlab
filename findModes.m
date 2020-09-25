function P = findModes(P,nModes,plotModes)
k_0 = 2*pi/P.lambda;  % [m^-1] Wavenumber
Nx = P.Nx_main;
Ny = P.Ny_main;

dx = P.Lx_main/Nx;
dy = P.Ly_main/Ny;

x = dx*(-(Nx-1)/2:(Nx-1)/2);
y = dy*(-(Ny-1)/2:(Ny-1)/2);
[X,Y] = ndgrid(x,y);

n_mat = P.n_cladding*ones(Nx,Ny);
for iShape = 1:size(P.shapes,1)
  switch P.shapes(iShape,4)
    case 1
      n_mat((X-P.shapes(iShape,1)).^2 + (Y-P.shapes(iShape,2)).^2 < P.shapes(iShape,3)^2) = P.shapes(iShape,5);
    case 2
      delta = max(dx,dy);
      r_diff = sqrt((X-P.shapes(iShape,1)).^2 + (Y-P.shapes(iShape,2)).^2) - P.shapes(iShape,3) + delta/2;
      n_mat(r_diff < delta) = min(r_diff(r_diff < delta)/delta*(P.n_cladding - P.shapes(iShape,5)) + P.shapes(iShape,5),P.shapes(iShape,5));
    case 3
      r_ratio_sqr = ((X-P.shapes(iShape,1)).^2 + (Y-P.shapes(iShape,2)).^2)/P.shapes(iShape,3)^2;
      n_mat(r_ratio_sqr < 1) = r_ratio_sqr(r_ratio_sqr < 1)*(P.n_cladding - P.shapes(iShape,5)) + P.shapes(iShape,5);
    case 4
      r_ratio_sqr = ((X-P.shapes(iShape,1)).^2 + (Y-P.shapes(iShape,2)).^2)/P.shapes(iShape,3)^2;
      r_abs = sqrt((X-P.shapes(iShape,1)).^2 + (Y-P.shapes(iShape,2)).^2);
      n_mat(r_ratio_sqr < 1) = 2*P.shapes(iShape,5)*exp(P.shapes(iShape,6)*r_abs(r_ratio_sqr < 1))./(exp(2*P.shapes(iShape,6)*r_abs(r_ratio_sqr < 1)) + 1);
    case 5
      r_ratio_sqr = (Y-P.shapes(iShape,2)).^2/P.shapes(iShape,3)^2;
      r_abs = Y - P.shapes(iShape,2);
      n_mat(r_ratio_sqr < 1) = 2*P.shapes(iShape,5)*exp(P.shapes(iShape,6)*r_abs(r_ratio_sqr < 1))./(exp(2*P.shapes(iShape,6)*r_abs(r_ratio_sqr < 1)) + 1);
  end
end
delta_n_2 = n_mat.^2 - P.n_0^2;                              %delta_n^2 in the expression for FD BPM

dz = 1e-10;
ax = 1.00001*dz/(dx^2*2i*k_0*P.n_0);
% ax = dz/(dx^2*2i*k_0*P.n_0);
ay = dz/(dy^2*2i*k_0*P.n_0);

N = Nx*Ny;                                                            %N*N - size of sparse matrices
M_rhs = sparse(1:N,1:N,1 + delta_n_2(:).'*dz*k_0/(2i*P.n_0),N,N) + ...
        sparse(1:N-1,2:N,[repmat([repmat(ax,1,Nx-1) 0],1,Ny-1) repmat(ax,1,Nx-1)],N,N) + ...
        sparse(2:N,1:N-1,[repmat([repmat(ax,1,Nx-1) 0],1,Ny-1) repmat(ax,1,Nx-1)],N,N) + ...
        sparse(1:N-Nx,1+Nx:N,ay,N,N) + ...
        sparse(1+Nx:N,1:N-Nx,ay,N,N);
M_rhs(1:N+1:N*N) = M_rhs(1:N+1:N*N) - repmat([ax repmat(2*ax,1,Nx-2) ax],1,Ny);
M_rhs(1:N+1:N*N) = M_rhs(1:N+1:N*N) - [repmat(ay,1,Nx) repmat(2*ay,1,Nx*(Ny-2)) repmat(ay,1,Nx)];

fprintf('Finding modes...\n');
tic
% [V,D] = eigs(M_rhs,nModes,'smallestimag','Display',true,'SubspaceDimension',nModes*10);
[V,D] = eigs(M_rhs,nModes,'smallestimag','Display',false,'SubspaceDimension',nModes*10);
fprintf('\b Done, %.1f seconds elapsed\n',toc);

for iMode = nModes:-1:1
  P.modes(iMode).Lx = P.Lx_main;
  P.modes(iMode).Ly = P.Ly_main;
  P.modes(iMode).field = reshape(V(:,iMode),[Nx Ny]);
  if plotModes
    h_f = figure(100+iMode);
    h_f.WindowStyle = 'docked';
    subplot(1,2,1);
    imagesc(abs(P.modes(iMode).field.').^2);
    axis equal; axis tight;
    subplot(1,2,2);
    imagesc(angle(P.modes(iMode).field).');
    axis equal; axis tight;
    sgtitle(['Mode ' num2str(iMode) ', eigenvalue ' num2str(D(iMode,iMode))]);
  end
end
end