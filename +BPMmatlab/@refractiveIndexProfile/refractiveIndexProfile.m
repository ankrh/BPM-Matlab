classdef refractiveIndexProfile
  properties
    Lx(1,1) double {mustBePositive} = 1
    Ly(1,1) double {mustBePositive} = 1
    n(:,:,:) single {mustBeFinite} = [] % May be complex and may be either 2D or 3D
    xSymmetry(1,1) BPMmatlab.symmetry = BPMmatlab.symmetry.NoSymmetry
    ySymmetry(1,1) BPMmatlab.symmetry = BPMmatlab.symmetry.NoSymmetry
  end
end