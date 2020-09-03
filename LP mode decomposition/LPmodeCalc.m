%% Function LPmodeCalc -  Hybrid LPlm modes - Step index fibre %%
% For a given set of fiber parameters, beam parameters, mode numbers and grid
% definition, this function calculates the hybrid LPlm modes, both cosine 
% and sine phi-dependence of a step index fibre. 

% Input parameters:
% r: Radius of core.
% n_core: Core index. 
% n_clad: Cladding index.
% N: Resolution of computational grid.
% X: Mesh grid of x coordinates. 
% Y: Mesh grid of y coordinates. 
% k0: wave number of incident beam.
% E0: Amplitude of incident beam.
% broots_mat: A cell structure containing the determined roots for each
% azimuthal mode number suppported by the fiber parameters.
% intNorm: Choose true for field to be normalized w.r.t. max intensity, false to normalize such that total power is 1

% Output parameters:
% Ecell_modesAB: Cell structure containing calculated LPlm modal fields for
% each broots_mat value passed to function. The modes is stored in a cell 
% NL x length(broots_mat) structure containing N x N matrices, one
% for each mode, both cosine (A) and sine (B) dependence. The cosine
% phi-dependent modes is in {1,1} an sine phi-dependent modes is in {1,2}. 
% Note that the modes are complex double precisions.


% See: https://www.rp-photonics.com/lp_modes.html
% See: Foundation of Guided-Wave Optics,Chen-Lin Chen, 2007 by Wiley & Sons. 
% See: Gloge D, Weakly Guiding Fibers, Appl. Opt. 1971, Vol 10, No. 10. 
% See: Fundamentals of Optical Fibers, John Buck, 2004 by Wiley and Sons.

%% Updated: 19-06-2020

%%
function [Ecell_modesAB] = LPmodeCalc(r,n_core,n_clad,N,X,Y,k0,E0,broots_mat,intNorm)
NL = length(broots_mat);
% The normalized radius of the fiber. 
R= sqrt(X.^2 + Y.^2)/r;
% Initializing the cell structure to contain LPlmA and LPlmB modes. 
Ecell_modesAB= cell(2,1); 
% Initializing the cell structure for that E matrix containing LPlm A (cosine) and LPlm B(sine).
Ecell_fieldA = cell(NL,length(broots_mat{1}));
Ecell_fieldB = cell(NL,length(broots_mat{1}));
for j1 = 1:NL
  L = j1-1; % The azimuthal mode number.
  brootv=broots_mat{j1}; % The determined roots for the particular azimuthal mode number.
  for j2 = 1:length(brootv)
    % Initialization of the electric field matrices.
    Emat_fieldA = NaN(N,N);
    Emat_fieldB = NaN(N,N);
    b = brootv(j2); % Use the root as normalized propagation constant.
    % Calculation of beta from normalized propagation constant.
    beta = k0.*sqrt((n_clad.^2+(b.*(n_core.^2-n_clad.^2))));
    % Calculation of normalized phase and attenuation constants.
    U = r*sqrt(n_core.^2*k0.^2-beta.^2);
    W = abs(r*sqrt(beta.^2-n_clad.^2*k0.^2));
    % Calculation of the modal fields for cosine(phi) dependence (A)
    Emat_fieldA(R <  1) = E0.*besselj(L,U*R(R < 1)).*cos(L*atan2(Y(R < 1),X(R < 1)));
    Emat_fieldA(R >= 1) = E0.*besselj(L,U)./besselk(L,W).*besselk(L,W*R(R >= 1)).*cos(L*atan2(Y(R >= 1),X(R >= 1)));
    if L > 0
      % Calculation of the modal fields for sine(phi) dependence (B)
      Emat_fieldB(R <  1) = E0.*besselj(L,U*R(R < 1)).*sin(L*atan2(Y(R < 1),X(R < 1)));
      Emat_fieldB(R >= 1) = E0.*besselj(L,U)./besselk(L,W).*besselk(L,W*R(R >= 1)).*sin(L*atan2(Y(R >= 1),X(R >= 1)));
    else
      Emat_fieldB = [];
    end
    Ecell_fieldA{j1,j2} = Emat_fieldA;  Ecell_fieldB{j1,j2} = Emat_fieldB;
    % Normalization of calculated modal fields, depending on intNorm
    % input.
    if L==0
      if intNorm
        Ecell_fieldA{j1,j2} = complex(Emat_fieldA./sqrt(sum(max(Emat_fieldA(:)).^2)));
      else
        Ecell_fieldA{j1,j2} = complex(Emat_fieldA./sqrt(sum(abs(Emat_fieldA(:)).^2)));
      end
    elseif L>=1
      if intNorm
        Ecell_fieldA{j1,j2} = complex(Emat_fieldA./sqrt(sum(max(Emat_fieldA(:)).^2)));
        Ecell_fieldB{j1,j2} = complex(Emat_fieldB./sqrt(sum(max(Emat_fieldB(:)).^2)));
      else
        Ecell_fieldA{j1,j2} = complex(Emat_fieldA./sqrt(sum(abs(Emat_fieldA(:)).^2)));
        Ecell_fieldB{j1,j2} = complex(Emat_fieldB./sqrt(sum(abs(Emat_fieldB(:)).^2)));
      end
    end
  end
end
% Allocation of calculated modal fields into separate cell entries, with
% respect to LPlm (A) modes (cosine) and LPlm (B) modes (sine).
Ecell_modesAB{1,1}=Ecell_fieldA; Ecell_modesAB{2,1}=Ecell_fieldB;
end
    
 

