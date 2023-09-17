function [radiallySymmetric,xC,yC] = testRadialSymmetry(X,Y,n,n_background,xSymmetry,ySymmetry)
n = double(n) - double(n_background);

if ySymmetry
  xC = 0;
else
  xC = sum(X(:).*abs(n(:)).^2)/sum(abs(n(:)).^2); % x centroid
end
if xSymmetry
  yC = 0;
else
  yC = sum(Y(:).*abs(n(:)).^2)/sum(abs(n(:)).^2); % y centroid
end

R = sqrt((X-xC).^2 + (Y-yC).^2); % Distances of all pixels from the centroid
[~,sortIdxs] = sort(R(:));
nsortedRounded = round(n(sortIdxs),5);

monotonicity_real = sign(diff(real(nsortedRounded)));
reversals_real = sum(abs(diff(monotonicity_real(monotonicity_real ~= 0))/2)); % Number of times the monotonicity of the real part changes (increasing/decreasing) as a function of the radial distance
monotonicity_imag = sign(diff(imag(nsortedRounded)));
reversals_imag = sum(abs(diff(monotonicity_imag(monotonicity_imag ~= 0))/2)); % Number of times the monotonicity of the imaginary part changes (increasing/decreasing) as a function of the radial distance

Rsorted = R(sortIdxs);
% figure(201);clf reset;plot(Rsorted,nsortedRounded);grid on;grid minor;
% figure(202);clf reset;imagesc(X(:,1),Y(1,:),n.'); axis equal tight;
% figure(203);clf reset;plot(Rsorted(1:end-1),diff(nsortedRounded));grid on;grid minor;

radiallySymmetric = reversals_real < 5 && reversals_imag < 5; % We arbitrarily choose that five reversals is the max we will allow. Non-radially-symmetric RI distributions will have many reversals.
end