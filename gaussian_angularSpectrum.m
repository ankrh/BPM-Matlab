%%
% Author: Madhu Veettikazhy
% Date: 22 March 2018
% Gaussian beam (CW) electric field propagation from the 
% focal point (initial waist plane width is defined here) 
% back to the distal end of the optical fiber. 
% ***************************************************************************************************

clc;
clear all;
close all;
%% General parameters

N=600;                                                                          % Resolution 
L=250e-6;                                                                      % Active area 
delta=L/N;                                                                      % Pixel pitch 
lambda=632.8e-9;                                                           % Wavelength 
lambdabar=lambda/(2*pi);
w=10e-6;                                                                        % Initial waist plane width of the gaussian beam
k=1/lambdabar;                                                               % Wavenumber
format LONGENG

%% Gaussian beam generation

x=(-N/2:N/2-1);
y=x;
[X,Y]=meshgrid(x,y);  

R = sqrt(X.^2+Y.^2);                                                        % r and phi co-ordinates
PHI = atan2(Y,X);
gauss = exp(-(R.^2).*(delta/w)^2);                                      % Gaussian field amplitude
phase = 1;                                                                          % Phase
E_gauss=gauss.*phase;                                                       % Electric field
E_gauss=E_gauss/sqrt(max(max(abs(E_gauss).^2)));            % Normalization of E field
% mesh(abs(E_gauss));

clearvars gauss phase x y                                                      % Clear unwanted variables

%% Fresnel Propagation

dist_z=1e-6;                                                                              % Steps in distance for every iteration (in meters)
prop_kernel= exp((-1i*dist_z*lambda*pi/(L^2)).*(R.^2));           % Fresnel propagation kernel
total_steps = 100;

E_field_prop = E_gauss;

%Propagating the E_gauss field in 'total_steps' steps of 'dist_z' meters
for iter = 1:total_steps
        E_field_fourier = fftshift(fft2(ifftshift(E_field_prop)));
        E_field_prop = fftshift(ifft2(ifftshift(E_field_fourier.*prop_kernel)));
        
        prop_distance = dist_z*iter/1e-6;                                       % In microns
        subplot(1,2,1) , mesh(abs(E_field_prop)),title('Field amplitude');
        subplot(1,2,2) , imagesc(angle(E_field_prop)),title('Phase');
        suptitle(['Propagated distance: ', num2str(prop_distance),' microns'])
        drawnow;
end

%%  Angular Spectrum Propagation
% 
% dist_z = 1e-6;                       % altitude (meters) 
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
