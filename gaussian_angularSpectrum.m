%%
% MADHU V
% 22 March 2018
% Gaussian continous wave propagation

%%
clc;
clear all;
close all;
%%

N=600;       %Resolution 
L=250e-6;         %Active area 
delta=L/N;     %Pixel pitch 
lambda=632.8e-9; 
lambdabar=lambda/(2*pi);
w=10e-6;         %Initial waist plane width
k=1/lambdabar;
format LONGENG

%% Calculate values for the LG mode
gauss=zeros(N);          %Gaussian 
phase=zeros(N);           %Phase

for x=-N/2:N/2-1
	for y=-N/2:N/2-1
        xp=(x+0.0000000001);
        yp=(y+0.0000000001);
        r=sqrt(xp^2+yp^2)*delta;
        gauss(x+N/2+1,y+N/2+1)=exp(-(r^2)/w^2);          % exp(-(x^2+y^2)/w^2)
        phase(x+N/2+1,y+N/2+1)=1; %exp(1i*(r^2));               
	end
end

%% Change the name of PSI function accordingly
E_gauss=gauss.*phase;
E_gauss=E_gauss/sqrt(max(max(abs(E_gauss).^2)));
mesh(angle(E_gauss));

clearvars gaus Phs x y r 

%% Fresnel Propagation

dist=10e-6;           %CHANGE THIS VALUE FOR VARYING RESOLUTION OF DISTANCE

prop_kernel=zeros(N);        %Kernel for propagation through a distance 'dist'
for im=-N/2:N/2-1
    for jn=-N/2:N/2-1
        prop_kernel(jn+N/2+1,im+N/2+1)=exp(-1i*dist*lambda*pi*((im/(N*delta))^2+(jn/(N*delta))^2));
    end
end

Qf_field = E_gauss;
for iter=1:100 
    %% Fresnel propagation
    tmp=fftshift(fft2(ifftshift(Qf_field)));
    Qf_field=fftshift(ifft2(ifftshift(tmp.*prop_kernel)));
    
  subplot(2,1,1) , mesh(abs(Qf_field)),title(['Propagated step: ', num2str(iter)]);
  subplot(2,1,2) , imagesc(angle(Qf_field)),title(['Propagated step: ', num2str(iter)]);
    drawnow;
end

% % % Angular Spectrum Propagation
% dist_z = 10e-6;                       % altitude (meters) 
% phy_x = w;               % physical width (meters) 
% phy_y = w;               % physical length (meters) 
% Fs_x = N/phy_x; 
% Fs_y = N/phy_y; 
% dx2 = Fs_x^(-1); 
% dy2 = Fs_y^(-1); 
% dFx = Fs_x/N;
% dFy = Fs_y/N; 
% Fx = (-Fs_x/2:dFx:(Fs_x/2 - dFx)); 
% Fy = (-Fs_y/2:dFy:(Fs_y/2 - dFy)); 
% 
% % alpha and beta (wavenumber components) 
% alpha = lambda.*Fx; 
% beta = lambda.*Fy;
% 
% % gamma_cust 
% gamma_cust = zeros(length(beta), length(alpha)); 
% for j = 1:length(beta)     
%     for i = 1:length(alpha)         
%         if (alpha(i)^2 + beta(j)^2) > 1             
%             gamma_cust(j, i) = 0;         
%         else
%             gamma_cust(j, i) = sqrt(1 - alpha(i)^2 - beta(j)^2);         
%         end
%     end
% end
% 
% % angular spectrm based formula 
% U1=E_gauss;
% for iter=1:100
% U1 = ifft2(ifftshift(fftshift(fft2(U1)).*exp(1i*k.*gamma_cust.*dist_z))); 
% I1 = (1/(16*pi)).*(U1.*conj(U1)); 
% imagesc(I1),title(['Propagated step: ', num2str(iter)]); 
% drawnow;
% % pause(1);
% end
