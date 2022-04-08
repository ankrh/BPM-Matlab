classdef model
  % MODEL Top-level BPM-Matlab model class
  %
  %   BPMmatlab.model collects all objects (as properties) and methods
  %   required to define a full BPM-Matlab model.

  properties
    % Visualization parameters
    name (1,:) char
    figTitle (1,:) char = ''
    figNum (1,1) double {mustBeInteger, mustBePositive} = 1
    updates (1,1) double {mustBeInteger, mustBePositive} = 50
    plotEmax (1,1) double {mustBeNonnegative} = 0
    plotZoom (1,1) double {mustBeGreaterThanOrEqual(plotZoom,1)} = 1
    storeE3D (1,1) logical = false
    saveVideo (1,1) logical = false
    intensityColormap (1,1) BPMmatlab.colormap = 'GPBGYR'
    phaseColormap (1,1) BPMmatlab.colormap = 'HSV'
    nColormap (1,1) BPMmatlab.colormap = 'Parula'
    calcModeOverlaps (1,1) logical = false
    disableStepsizeWarning (1,1) logical = false
    disablePlotTimeWarning (1,1) logical = false

    % Solver parameters
    useAllCPUs (1,1) logical = false
    useGPU (1,1) logical = false
    Nx_main (1,1) double {mustBeInteger, mustBePositive} = 2
    Ny_main (1,1) double {mustBeInteger, mustBePositive} = 2
    xSymmetry (1,1) BPMmatlab.symmetry = 'NoSymmetry'
    ySymmetry (1,1) BPMmatlab.symmetry = 'NoSymmetry'
    dz_target (1,1) double {mustBePositive} = 1e-6
    padfactor (1,1) double {mustBeGreaterThanOrEqual(padfactor,1)} = 1.5

    % Geometry parameters
    Lx_main (1,1) double {mustBePositive} = 1
    Ly_main (1,1) double {mustBePositive} = 1
    Lz (1,1) double {mustBePositive} = 1
    taperScaling (1,1) double {mustBePositive} = 1
    twistRate (1,1) double {mustBeFinite, mustBeReal} = 0
    bendingRoC (1,1) double {mustBePositive} = Inf
    bendDirection (1,1) double {mustBeFinite, mustBeReal} = 0

    % Optical and material parameters
    alpha (1,1) double {mustBePositive} = 3e14
    lambda (1,1) double {mustBePositive} = 1
    n_background (1,1) double {mustBePositive} = 1
    n_0 (1,1) double {mustBePositive} = 1
    rho_e (1,1) double {mustBePositive} = 0.22

    % Refractive index profile
    n (1,1) BPMmatlab.refractiveIndexProfile

    % Electric field to propagate
    E (1,1) BPMmatlab.electricFieldProfile
  end
  
  properties (Dependent)
    Nx
    Ny
    Nz
    Lx
    Ly
    dx
    dy
    dz
    x
    y
  end
  
  properties (SetAccess = private)
    powers
    xzSlice
    yzSlice
    videoHandle
    modes BPMmatlab.electricFieldProfile
    E3D single
    z
    priorData = false
  end
  
  methods
    function P = model() % Constructor
      P.name = ['BPM-Matlab model ' char(datetime('now','Format','d-MMM-y HH.mm.ss'))];
    end

    function dx = get.dx(P)
      dx = P.Lx_main/P.Nx_main;
    end

    function dy = get.dy(P)
      dy = P.Ly_main/P.Ny_main;
    end

    function dz = get.dz(P)
      dz = P.Lz/P.Nz;
    end

    function Nx = get.Nx(P)
      targetLx = P.padfactor*P.Lx_main;
      Nx = round(targetLx/P.dx);
      if P.ySymmetry
        Nx = Nx + (rem(Nx,2) ~= rem(P.Nx_main,2)); % Ensure that if Nx_main was set odd (to have a x slice at the center), Nx will also be odd
      end
    end

    function Ny = get.Ny(P)
      targetLy = P.padfactor*P.Ly_main;
      Ny = round(targetLy/P.dy);
      if P.xSymmetry
        Ny = Ny + (rem(Ny,2) ~= rem(P.Ny_main,2)); % Ensure that if Ny_main was set odd (to have a y slice at the center), Ny will also be odd
      end
    end

    function Nz = get.Nz(P)
      Nz = max(P.updates,round(P.Lz/P.dz_target));
    end

    function Lx = get.Lx(P)
      Lx = P.dx*P.Nx;
    end

    function Ly = get.Ly(P)
      Ly = P.dy*P.Ny;
    end

    function x = get.x(P)
      x = getGridArray(P.Nx,P.dx,P.ySymmetry);
    end

    function y = get.y(P)
      y = getGridArray(P.Ny,P.dy,P.xSymmetry);
    end

    function finalizeVideo(P)
  	  close(P.videoHandle);
    end

    P = FD_BPM(P)
    P = FFT_BPM(P)
    P = initializeRIfromFunction(P,hFunc,nParameters)
    P = initializeEfromFunction(P,hFunc,Eparameters)

    [E,n_slice,precisePower] = FDBPMpropagator(E,mexParameters);
    [E,n_slice,precisePower] = FDBPMpropagator_CUDA(E,mexParameters);
  end
end