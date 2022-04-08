classdef electricFieldProfile
  properties
    Lx(1,1) double {mustBePositive} = 1
    Ly(1,1) double {mustBePositive} = 1
    field(:,:) single {mustBeFinite} = [] % Is usually complex
    xSymmetry(1,1) BPMmatlab.symmetry = BPMmatlab.symmetry.NoSymmetry
    ySymmetry(1,1) BPMmatlab.symmetry = BPMmatlab.symmetry.NoSymmetry
    label (1,:) char = ''
    neff (1,1) double = NaN
  end
end