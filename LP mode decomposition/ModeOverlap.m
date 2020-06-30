%% Function ModeOverlap -  Overlap integral calculation %%
% This function calculates the overlap integral for two input fields and
% deliver an output cell containing the coefficients and the absolute
% squared values (fractional power). 

% Input parameters:
% Emat_data: A NxN complex single valued matrix, containing the electric field
% for the mode overlap is being considered. 
% Ecell_mode: A cell structure or matrix containing the electric fields for
% which the mode overlap integral is being considered, with respect to
% Emat_data.

% Output parameters:
% Ccell_out: A 2x1 cell structure containing the overlap coefficients in 
% entry 1 and the absolute squared coefficient value (fractional power) in
% entry 2. 


% Note that Emat_data must be an NxN complex single valued matrix. Ecell_mode
% must be either a cell structure containing containing a number of complex
% single valued matrices or a NxN complex singled valued matrix. 

% See: Foundation of Guided-Wave Optics,Chen-Lin Chen, 2007 by Wiley & Sons. 
% See: Fundamentals of Optical Fibers, John Buck, 2004 by Wiley and Sons.
% See: J. Demas, L. Rishøj and S. Ramachandran, Optics Express, Volume 23
% Issue 22, Free-space beam shaping for precise control and conversion of modes
% in optical fiber. 

%% Updated: 19-06-2020

%%
function [Ccell_out] = ModeOverlap(Emat_data,Ecell_mode)
%
[NL,Nm] = size(Ecell_mode); % Initialization of NL and Nm. 
Ccell_coeff = cell(NL,Nm); % Initializing overlap coefficient cell structure. 
Ccell_abs = cell(NL,Nm); % Initializing Overlap coefficient.^2 cell structure.  
% Calculation of the overlap integral and subsequently storrage of the
% numbers in the initialzied cell structures. 
for j1=1:NL
  for j2=1:Nm
    Emat_mode=Ecell_mode{j1,j2};  % Initialization the overlap integral one cell entry at a time.
    if isempty(Emat_mode)== 1 % If Emat_mode is empty, no overlap integral calculation is preformed
      % NaN is passed to C_coeff cell structure. 
      C_coef = NaN;
    else
      C_coef = sum(Emat_data(:).*conj(Emat_mode(:))); % The overlap integral calculation.
    end
    C_abs = abs(C_coef).^2;
    Ccell_coeff{j1,j2}=C_coef;
    Ccell_abs{j1,j2}=C_abs;
  end
end
% Post processing of the mode overlap cells for proper visualization.  
% Initialization of the output cell structures. 
Ccell_out = cell(2,1);
Ccellout_abs = cell(NL,Nm);
% An single valued number is entered in Ccellout_abs, for all empty entries
% in Ecell_mode, yielding no overlap integral calculation. 
for j3=1:NL
  for j4=1:Nm
    Ccellout_abs{j3,j4}=single(Ccell_abs{j3,j4});
  end
end
% Conversion of Ccellout_abs cell structure into a NxN complex single
% valued matrix for proper visualization. 
Cmatout_abs=cell2mat(Ccellout_abs);
% Allocating calculated coefficients and absolute squared coefficient
% values to specific output cell entries. 
Ccell_out{1,1}=Ccell_coeff; Ccell_out{2,1}=Cmatout_abs; 
end  