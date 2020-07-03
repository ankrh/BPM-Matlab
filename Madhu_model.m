



% P.n_colorlimits = [1.507 1.521]; % Minimum and maximum caxis for refractive index plot (figure 1)

% if fibreType == 4 % Photonic Lantern
%   Eparameters = {w_0,fibreType,shapeParameters,numberOfCores,pitch,k_0,lambda,Ecell_modeA{1,9},photonicLanternInput,SavedFileName};    % Cell array of parameters that the E field initialization function (defined at the end of this file) will need
% elseif exist('Ecell_modeA','var') % For LP mode propagation through MMF or SMF
%   Eparameters = {w_0,fibreType,shapeParameters,numberOfCores,pitch,k_0,lambda,Ecell_modeA{1,1},SavedFileName};  %9 input fields
% else
%   Eparameters = {w_0,fibreType,shapeParameters,numberOfCores,pitch,k_0,lambda,SavedFileName};  %8 input fields
% end

% pitch = 20e-6;  % 3/4*125e-6; % [m] Intercore separation in multicore fibre, in photonic lantern h=3/2 R_clad, R_clad is 125/2 um
% numberOfCores =37; % [] Numer of cores in the multicore fibre, not used for photonic lantern
% coreRadius = 0.5e-3;  %[m] Core radius of the SMF/MMF/MCF/GRIN lens
% fibreType = 5;  % Type of fibre for E field initialization - 1: Single/Multimode, 2: Hex multicore, 3: Fermat's multicore 4: Photonic Lantern
%                          %  21/31: Hex/Fermat's multicore with phase programming, 5: GRIN Lens
% FibreParameters = {fibreType,numberOfCores,pitch,coreRadius}; 
% photoelasticCoeff = 0.22;  %[] coefficient depending on Poisson’s ratio and componentsof the photoelastic tensor - in bending expression
% photonicLanternInput = 2; %[] 1: LP01 input to SSMF, 2: LP01 input to HI1060 on the right, 3: LP01 input to HI1060 below

% bendingRoC{1} = Inf;  %[m] Bending radius of curvature for the fibre section
% bendDirection{1} = 0;  % [deg] The angle of bending direction, 0: bending in +x, 90: bending in +y

%% USER DEFINED E-FIELD INITIALIZATION FUNCTION
% function E = calcInitialE(X,Y,Eparameters) % Function to determine the initial E field. Eparameters is a cell array of additional parameters such as beam size
% [Nx, Ny] = size(X);
% dx = X(2,1) - X(1,1);
% dy = Y(1,2) - Y(1,1);
% % amplitude = exp(-((X-Lx_main/4).^2+Y.^2)/w_0^2) - exp(-((X+Lx_main/4).^2+Y.^2)/w_0^2); % Gaussian field amplitude
% % amplitude = exp(-((X-Lx_main/10).^2+Y.^2)/w_0^2); % Gaussian field amplitude
% % amplitude = exp(-(X.^2+Y.^2)/w_0^2); % Gaussian field amplitude
% w_0 = Eparameters{1};
% fibreType = Eparameters{2};
% shapeParameters = Eparameters{3};
% numberOfCores = Eparameters{4};
% pitch = Eparameters{5};
% k_0 = Eparameters{6};
% lambda = Eparameters{7}; 
% focus = 5e-3; 
% SavedFileName = Eparameters{end}; 
% 
% switch fibreType
%   case 1
%     E = Eparameters{8}; % LP mode for SMF and MMF
%     if size(E,1) == 1 % when it is savedFileName
%       disp('New field is calculated to propagate through the fibre with the given w_0, not the LP mode. If you need to input LP mode, please run Example_LPmodes.m');
%       amplitude = exp(-((X-shapeParameters{1}(1)).^2+(Y-shapeParameters{1}(2)).^2)/w_0^2);
%       phase = zeros(size(X));
%       E = amplitude.*exp(1i*phase);
%     end
%   case {2, 3}
%     amplitude = zeros(size(X));
%     for i = 1:3:numel(shapeParameters{1})
%       amplitude = amplitude+exp(-((X-shapeParameters{1}(i)).^2+(Y-shapeParameters{1}(i+1)).^2)/w_0^2);
%     end
%     phase = zeros(size(X));
%     E = amplitude.*exp(1i*phase);
%   case 4
%     LPmode = Eparameters{8};
%     photonicLanternInput = Eparameters{9};
%     E = zeros(size(X));
%     pixelsX = Nx/10; pixelsY = Ny/10;
%     x_coord_pixel_1  = -ceil(pitch/2/sqrt(3)/dx);           y_coord_pixel_1 = (pitch/2/dy);
%     x_coord_pixel_2  = ceil(pitch/sqrt(3)/dx);               y_coord_pixel_2 = 0;
%     x_coord_pixel_3  = -ceil(pitch/2/sqrt(3)/dx);           y_coord_pixel_3 = -(pitch/2/dy);
%     switch photonicLanternInput
%       case 1
%         E(Nx/2+1+x_coord_pixel_1-pixelsX:Nx/2+1+x_coord_pixel_1+pixelsX,Ny/2+y_coord_pixel_1-pixelsY:Ny/2+y_coord_pixel_1+pixelsY) ...
%           =LPmode(Nx/2-pixelsX:Nx/2+pixelsX,Nx/2-pixelsY:Nx/2+pixelsY);
%       case 2
%         E(Nx/2+1+x_coord_pixel_2-pixelsX:Nx/2+1+x_coord_pixel_2+pixelsX,Ny/2+y_coord_pixel_2-pixelsY:Ny/2+y_coord_pixel_2+pixelsY) ...
%           =LPmode(Nx/2-pixelsX:Nx/2+pixelsX,Nx/2-pixelsY:Nx/2+pixelsY);
%       case 3
%         E(Nx/2+1+x_coord_pixel_3-pixelsX:Nx/2+1+x_coord_pixel_3+pixelsX,Ny/2+y_coord_pixel_3-pixelsY:Ny/2+y_coord_pixel_3+pixelsY) ...
%           =LPmode(Nx/2-pixelsX:Nx/2+pixelsX,Nx/2-pixelsY:Nx/2+pixelsY);
%       case 23     %Simulataneously exciting 2 and 3 input fibres
%         E(Nx/2+1+x_coord_pixel_2-pixelsX:Nx/2+1+x_coord_pixel_2+pixelsX,Ny/2+y_coord_pixel_2-pixelsY:Ny/2+y_coord_pixel_2+pixelsY) ...
%           =LPmode(Nx/2-pixelsX:Nx/2+pixelsX,Nx/2-pixelsY:Nx/2+pixelsY);
%         E(Nx/2+1+x_coord_pixel_3-pixelsX:Nx/2+1+x_coord_pixel_3+pixelsX,Ny/2+y_coord_pixel_3-pixelsY:Ny/2+y_coord_pixel_3+pixelsY) ...
%           =LPmode(Nx/2-pixelsX:Nx/2+pixelsX,Nx/2-pixelsY:Nx/2+pixelsY);
%     end
%     
%   case {21, 31}   % Should have run case 2 or 3 before and saved the E data
%     load([SavedFileName,'.mat']);
%     amplitude = zeros(size(X));
%     for i = 1:3:numel(shapeParameters{1})
%       amplitude = amplitude+exp(-((X-shapeParameters{1}(i)).^2+(Y-shapeParameters{1}(i+1)).^2)/w_0^2);
%     end
%     phase = zeros(size(X));
%     acquiredPhase = NaN(1,numberOfCores);    % Row vector (1D)
%     focusPhase = NaN(1,numberOfCores);
%     
%     for idx = 1:numberOfCores
%       acquiredPhase(idx) = angle(E(Nx/2+ceil(shapeParameters{1}(idx*2+idx-2)/dx),Ny/2+ceil(shapeParameters{1}(idx*2+idx-1)/dy)));  %Acquired phase of E field at distal end for previous travel through fibre
%       focusPhase(idx) = -k_0*((shapeParameters{1}(idx*2+idx-2))^2+(shapeParameters{1}(idx*2+idx-1))^2)/(2*focus);  %Focusing phase for point focus
%       %             focusPhase(idx) = -k_0*(shapeParameters{1}(idx*2+idx-1))^2/(2*focus);  %Focusing phase for horizontal line focus
%       phase(sqrt((X-shapeParameters{1}(idx*2+idx-2)).^2+(Y-shapeParameters{1}(idx*2+idx-1)).^2) < pitch/2) = focusPhase(idx)-acquiredPhase(idx);
%     end
%     E = amplitude.*exp(1i*phase);
%     case 5
%       amplitude = exp(-((X-shapeParameters{1}(1)).^2+(Y-shapeParameters{1}(2)).^2)/w_0^2);
%       phase = zeros(size(X));
%       E = amplitude.*exp(1i*phase);
% end
% 
% % amplitude2 = 2*exp(-((X+12e-6).^2+(Y+7e-6).^2)/w_0^2);
% % phase2 = 8e5*Y;
% % if ~E
% % E = amplitude.*exp(1i*phase);% + amplitude2.*exp(1i*phase2); % Electric field
% % end
% end


%% USER DEFINED SHAPE-PARAMETERS INITIALIZATION FUNCTION FOR MULTICORE FIBRE
% function shapeParameters = getShapeParameters(segment,FibreParameters)
% fibreType = FibreParameters{1};
% 
% switch fibreType
%   case 1
%       R = FibreParameters{4};
%     shapeParameters{segment} = [0; % x values
%       0; % y values
%       R]; % r values
%   case {2, 21}
% 	numberOfCores = FibreParameters{2};
%     pitch = FibreParameters{3};
%     R = FibreParameters{4};
%     shapeParameters{segment} = NaN(3,numberOfCores); % Initialize output array
%     shapeParameters{segment}(3,:) = R; % All cores have radius R
%     shapeParameters{segment}(1:2,1) = [0; 0]; % x and y of the center core
%     shellSideIdx = 1; % Which side of the hex are we filling?
%     shellSideCoreIdx = 0; % How many cores on this side have been filled so far?
%     shellNum = 1; % Which shell are we in? The center core is not counted as a shell.
%     for coreIdx = 2:numberOfCores
%       if shellSideCoreIdx == 0 % If this is the first core in this shell
%         shapeParameters{segment}(1:2,coreIdx) = [shellNum*pitch; 0];
%       else % Find new core position by adding onto the previous core's position
%         shapeParameters{segment}(1:2,coreIdx) = shapeParameters{segment}(1:2,coreIdx-1) + [pitch*cos(shellSideIdx*pi/3 + pi/3); pitch*sin(shellSideIdx*pi/3 + pi/3)];
%       end
%       
%       if shellSideCoreIdx == shellNum % If this side has been filled
%         shellSideIdx = shellSideIdx + 1;
%         shellSideCoreIdx = 1;
%       else % Continue filling this side
%         shellSideCoreIdx = shellSideCoreIdx + 1;
%       end
%       
%       if shellSideCoreIdx == shellNum && shellSideIdx == 6 % Last core on last side would be a replicate of the first one drawn in this shell, so skip
%         shellNum = shellNum + 1;
%         shellSideIdx = 1;
%         shellSideCoreIdx = 0;
%       end
%     end
%   case {3, 31}
% 	numberOfCores = FibreParameters{2};
%     pitch = FibreParameters{3};
%     R = FibreParameters{4};     
%     rho_n = pitch*sqrt(1:numberOfCores);
%     theta_n = (1:numberOfCores)*pi*(3-sqrt(5));
%     x_n = rho_n.*cos(theta_n);
%     y_n = rho_n.*sin(theta_n);
%     shapeParameters{segment} = NaN(3,numberOfCores); % Initialize output array
%     shapeParameters{segment}(3,:) = R; % All cores have radius R
%     for coreIdx = 1:numberOfCores
%         shapeParameters{segment}(1:2,coreIdx) = [x_n(coreIdx); y_n(coreIdx)];
%     end
%   case 4
%     pitch = FibreParameters{3};
%     shapeParameters{segment} = [-pitch/2/sqrt(3)    pitch/sqrt(3)     -pitch/2/sqrt(3)    -pitch/2/sqrt(3)    pitch/sqrt(3)     -pitch/2/sqrt(3);  % x values
%       pitch/2     0    -pitch/2     pitch/2    0     -pitch/2; % y values
%       62.5e-6   62.5e-6   62.5e-6   4.5e-6      2.65e-6     2.65e-6]; % r values
%   case 5
%     R = FibreParameters{4};
%     shapeParameters{segment} = [0; % x values
%       0; % y values
%       R]; % r values
%   otherwise
%     disp('This fibre type is not supported');
%     return;
% end
% end
% 
% function shapeRIs = getShapeRIs(segment,fibreType,shapeParameters,n_core)
% switch fibreType
%   case 4
%     shapeRIs{segment} = [1.4444 1.4444 1.4444 1.4492 1.4511 1.4511]; %
%   otherwise
%     shapeRIs{segment} = n_core*ones(1,size(shapeParameters{1},2));
% end
% end