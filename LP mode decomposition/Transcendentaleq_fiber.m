%% Function - Trancendental equation - Step index fibre %%
% For a step index fibre, this function calculates the trancendental
% equation for each supported azimuthal mode number (L). The transcendental
% equation is written in terms of a normalized propagation constant. Then
% this function determines the roots by a bracket approach, looking for 
% change in sign within a bracket. All singular roots is disregarded. 

% Input parameters:
% V: The normalized frequency V=k0.*r_core.*sqrt(n_core^2-n_clad^2).
% bv: A vector that contains the normalized propagation constant
% elements (b-axis).
 

% Ouput parameters: 
% fv: The calculated trancendenetal equations. The fv vector is (1 x resolution of bv)
% double input vector. det(f)=0 for LHS=RHS. fv=LHS-RHS. 
% LHS_mat: Left-hand side of fv expression. This part describes the
% core region. 
% RHS_mat: Right-hand side of fv expression. This part describes the
% cladding region.
% broots_mat: A matrix structure containing the determined roots for each
% azimuthal mode number. broots_mat is a L x 1 cell structure, with each row in 
% column containing multiple entries. L being the number of supported
% azimuthal mode numbers. 

% See: Foundation of Guided-Wave Optics,Chen-Lin Chen, 2007 by Wiley & Sons.
% See: Numerical Methods for Engineers and Scientists Using MATLAB, 2013 CRC Press.

%% Updated: 19-06-2020

%%
function [fv,LHS_mat,RHS_mat,broots_mat] = Transcendentaleq_fiber(V,bv)
% Initialization of LHS and RHS cells structures for calculation storaged. 
LHS_mat=cell(length(bv),1); RHS_mat=cell(length(bv),1); 
W=V.*sqrt(bv); % Normalized coefficient , Modified bessel function of second kind.
U=V.*sqrt(1-bv); % Normalized coefficient, Bessel function of first kind.  

% Initialization of the while loop, which calculates the trancendental
% equation and subsequently find all the roots for each transcendental equation 
% calculated this way. The done = true, when no further indices is returned empty
% before the allocation to the output cell structure. Rootfinding section, see line 62-68. 
L = 0;
done = false;
while ~done
  j1 = L+1;
  % Firstly calculation of the left and right handed side of
  % transcendental equation.
  if L==0
    % Calculation of the left hand side at rho=r_core.
    LHS=U.*besselj(1,U)./besselj(0,U);
    LHS_mat{j1}=LHS;
    % Calculation of the right hand side at rho=r_core.
    RHS=W.*besselk(1,W)./besselk(0,W);
    RHS_mat{j1}=RHS;
  elseif L>=1
    % Calculation of the left hand side at rho=r_core.
    LHS=U.*besselj(L-1,U)./besselj(L,U);
    LHS_mat{j1}=LHS;
    % Calculation of the right hand side at rho=r_core.
    RHS=-W.*besselk(L-1,W)./besselk(L,W);
    RHS_mat{j1}=RHS;
  end
  fv = LHS - RHS; % Definition of the trancendental equation.
  signchangepositions = find(abs(diff(fv>0))); % finding changes of sign indices.
  % Finding root values of f = LHS - RHS at the change of sign indices.
  for j2=length(signchangepositions):-1:1
    if abs(fv(signchangepositions(j2))-fv(signchangepositions(j2) + 1)) > abs(fv(signchangepositions(j2) - 1)-fv(signchangepositions(j2) + 2))
       signchangepositions(j2) = []; % Remove this element if it's a singular. 
    end
  end
  % Displaying the roots in the correct order, with respect to LPlm -
  % mode designation.
  broots = fliplr(bv(signchangepositions));
  if ~isempty(broots)
    % Allocation of root values from root finding algorithm to the cell structure.
    broots_mat{j1,:} = broots; %#ok<AGROW>
    L = L+1;
  else
    done = true;
  end
end

end