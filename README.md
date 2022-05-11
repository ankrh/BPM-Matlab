# BPM-Matlab

## What is this repository for?
BPM-Matlab is a MATLAB-based numerical simulation tool for solving the paraxial Helmholtz equation using the Douglas-Gunn Alternating Direction Implicit (DG-ADI) method to efficiently model the electric field propagation using a Finite-Difference Beam Propagation Method (FD-BPM) in a wide variety of optical fiber geometries with arbitrary refractive index profiles.

Included is also a solver for the Fast Fourier Transform Beam Propagation Method, FFT_BPM.

You can find the latest version of BPM-Matlab at https://github.com/ankrh/BPM-Matlab.


## LICENSE
This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>.


## How do I get set up?
### Requirements
- Windows 7, macOS 10.12 (Sierra) or newer or Linux
- MATLAB R2018a or newer
- (For GPU accelerated computation) A Windows PC with a CUDA-enabled graphics card and the MATLAB Parallel Computing Toolbox

### BPM-Matlab function files
All the functions needed for running BPM-Matlab are located in subfolders of the folder that the examples are in. If you set the folder containing the examples to your MATLAB working direction or you add it to the MATLAB path, all the functions will be automatically located. If you choose to keep your model files (see below) in a different folder to the example files, you will need to manually add the folder with the example files to your MATLAB path.

## How do I use BPM-Matlab?
### Compilation
The folders include all the executables necessary, so you don't need to compile anything. If, however, you want to change the routine in the FDBPMpropagator.c source code located in the folder "src"), you will need to recompile the mex-file. Check out the header of the source-file on how to do so. The source code is written in such a way that it can be compiled as either C or C++ code using either GCC, MSVC, clang or as CUDA source code using NVCC.

### Model files
In BPM-Matlab, you set up your model in a single m-file. You can find a few examples to get you started. Once you're familiar with those, you can start to write your own model files based on those example files.

#### Segments
If the model contains multiple segments (for example, a stretch of straight fiber followed by a stretch of bent fiber), BPM-Matlab will add a new segment with each call to the FD_BPM() or FFT_BPM() solver functions, using as the input E-field the simulation output of the previous segment (see below). As such, you can stack segments one after the other by only re-defining the parameters that change. Check Example2.m for an example.

### List and explanation of parameters
The below parameters apply to the FD_BPM solver. The FFT_BPM solver uses a subset of the below parameters.

#### Visualization parameters
These parameters affect how the visualizations will look. They will not affect the actual calculation of the E-field result.

- `name`
(Default: 'BPM-Matlab model ' followed by a date/time stamp)
A descriptive name of the model. This will only be used if you choose to save a video of the simulation.

- `figTitle`
(Default: Empty char array)
A char array of text to use as the title of the various figures generated. This could, for example, be "Segment 1" or "Tapered segment" etc.

- `figNum`
(Default: 1)
The figure number offset to use. Change this if you want to avoid having one segment overwrite the figure contents of a previous segment.

- `updates`
(Default: 50)
The number of times during the simulation the plot in the frontend should update. This is useful for following the E-field evolution, but adds overhead to the calculation, as for each update, the currently calculated full simulation window has to be extracted and displayed. If you're only interested in the E-field at the end of the waveguide, set this to 1.

- `plotEmax`
(Default: 0)
Here you can set the maximum of the color scale in the intensity plot, relative to the peak of initial intensity. If set to 0, the intensity plot autoscales at each update.

- `plotZoom`
(Default: 1)
How far to zoom in in the plots. This does not affect the solver calculations. To zoom to twice the size, use plotZoom = 2, etc.

- `storeE3D`
(Default: false) If this is set to true, FD_BPM will store all the update slices of the propagated E-field in the 3D array P.E3D, and will also plot it volumetrically once the simulation is complete. Each segment will append its 3D array of E-values to a new cell in P.E3D. The reason for this is that the solver window parameters (Lx, Ly, Nx, Ny) may be different from one segment to the next.

- `saveVideo`  
(Default: false) Set this to true to record the field intensity and phase profiles at different transverse planes in preparation of saving to a video. To save the video, you must call finalizeVideo\(P\) after having run all the segments.

- `Intensity_colormap`, `Phase_colormap`, `n_colormap`  
(Defaults: 'GPBGYR', 'HSV' and 'parula', respectively)
An object that designates the colormap of the respective subplot. The options are  
  - 'GPBGYR'
  - 'HSV'
  - 'parula'
  - 'gray'
  - 'cividis'

- `calcModeOverlaps`
(Default: false)
If true, this will make the solver calculate the overlap of the E-field with any precalculated modes in each update. A plot will be generated to show these overlaps.

- `disableStepsizeWarning`
(Default: false)
If true, the FD_BPM solver will not warn you if the step sizes you have specified in the dz_target property are potentially so large that they will generate numerical artifacts.

- `disablePlotTimeWarning`
(Default: false)
If true, the FD_BPM solver will not warn you if the updates take more than 50% of the total execution time of the solver.

#### Solver parameters
These parameters will affect how the solver will numerically describe the problem. You have to ensure yourself that the pixel size and z step size are small enough for the simulation to converge. This is typically done by doing a manual parameter scan of the resolution parameters.

- `useAllCPUs`  
(Default: false)
On Windows and Linux, BPM-Matlab will by default use multithreading but leave one processor core unused, which is useful for doing other work on the PC while simulations are running. If you are in need for speed and don't plan to use your PC while doing the simulation, set this parameter to true. Multithreading on Mac is not supported.

- `useGPU`  
(Default: false)
This allows to use CUDA acceleration for NVIDIA GPUs. Set this to true if you have a CUDA enabled NVIDIA GPU.

- `Nx_main`, `Ny_main`  
The number of pixels in the main simulation window.

- `xSymmetry`, `ySymmetry`
(Default: 'NoSymmetry')
An object that designates whether the solver should use symmetry assumptions during the solving process. The possible values are:
  - 'NoSymmetry' - No symmetry assumptions should be made.
  - 'Symmetry' - Assume that both the refractive index and the electic field have ordinary symmetry under mirroring in the axis.
  - 'AntiSymmetry' - Assume that the refractive index exhibits ordinary symmetry and the electric field exhibits antisymmetry (sign inversion) under mirroring in the axis.

- `dz_target`
[m]
(Default: 1e-6)
The solver will divide the segment into steps that are as close as possible to this value. The finer the steps, the fewer numerical artifacts will be present in the result but the execution will be slower.

- `padfactor`
(Default: 1.5)
How much absorbing padding to add on the sides of main simulation window. The full simulation window then consists of the main simulation window in the middle and the absorbing padding around it. `P.padfactor` is the ratio of the widths (`Lx`, `Ly`) of the full simulation window to the widths (`Lx_main`, `Ly_main`) of the main simulation window. 1 means no padding, 2 means the full simulation window is twice as wide as the main simulation window, i.e., the absorbing padding is of thickness `Lx_main/2` and `Ly_main/2` on the sides.
The absorbing padding is sometimes required to deal with the non-guided energy: if this hits the boundary of the calculated frame, it would numerically get reflected and propagate back into the main area, thus leading to erroneous results. The absorbing padding deals with this problem by removing energy that has propagated outside the main area from the calculated frame.

- `alpha`
[1/m^3]
(Default: 3e14)
Absorption strength of the absorbing padding. The absorption of the padding is implemented as an absorption coefficient, the magnitude of which is proportional to the distance out from the edge of the main simulation window, squared, such that the energy propagating further out is attenuated more strongly. In other words, the field intensity at a distance d out from the main window is multiplied in each step of length dz by a `exp(-d^2×alpha×dz)`.

#### Geometry parameters
These parameters describe the spatial layout of the structure you want to model.

- `Lx_main`, `Ly_main`  
[m]
The physical size of the main simulation window. If you are simulating a fiber in which all the light you're interested in is in the core, then you could set this to be a bit larger than the core diameter. If you are simulating a fiber in which you are also interested in cladding-guided light, then you need to set these to be a bit larger than the cladding diameter.

- `Lz`  
[m]
Length of the segment propagated through.

- `taperScaling`
(Default: 1)
The ratio of the width of the structure at the end of the segment to the width at the beginning of the segment. The solver will assume a linear tapering from the start until the end.

- `twistRate`
[radians/m]
(Default: 0)
The rate of twisting rotation of the fiber.

- `bendingRoC`
[m]
(Default: Inf)
Radius of curvature of the fiber bend.

- `bendDirection`
[degrees]
(Default: 0)
Direction of the bending, in a polar coordinate system with 0° to the right (towards positive x) and increasing angles in counterclockwise direction.

#### Optical and material properties
- `lambda`
[m]
The wavelength in vacuum. If you want to simulate broadband light propagation, you have to do a simulation for each involved wavelength separately.

- `n_background`
(Default: 1)
Refractive index of the background, typically corresponding to the cladding.

- `n_0`
(Default: 1)
The reference refractive index, a parameter that is introduced during the derivation of the particular form of the paraxial Helmholtz equation used in BPM-Matlab. The refractive indices defined through `n_background` and in `n` are the actual indices used in the simulation. `n_0` (the reference refractive index) should be chosen such that the paraxial Helmholtz equation remains a good approximation, i.e., close or equal to the refractive indices where the main part of the energy is propagating in the simulation. For core-guided light, it would usually be a good idea to set `n_0 = n_core`.

- `rho_e`
(Default: 0.22)
A material contant used in the formula for the effective refractive index in a bent fiber. It comprises the Poisson ratio and the photo-elastic tensor components. For silica, it is 0.22.

#### Refractive index
The refractive index may be defined as either a 2D or 3D refractive index.

In the case of a 2D refractive index, it is assumed that it is unchanged throughout the segment, with the exception of any tapering and twisting you may have defined. For a 3D refractive index, the third dimension of the array corresponds to the different z coordinates along the segment.

Also be aware that after a segment has been run, the final refractive index profile is stored in the model (`P.n`), ready to be used automatically as the input to the next segment.

You are free to redefine the simulation grid from one segment to the next (for example, if one segment requires finer resolution than the others). In that case, the refractive index will automatically be interpolated from the old grid to the new grid using the interpn function without the need for user interaction.

There are two ways to define the refractive index:

##### Method 1, using initializeRIfromFunction()
This method is probably the most common and consists of defining a MATLAB function at the end of the model file that provides the 2D or 3D (possibly complex) refractive index profile as an output.

You run the function initializeRIfromFunction() as shown in example 1 (for 2D) or example 14 (for 3D) to make BPM-Matlab apply your function to the simulation grid in preparation for running the solver.

If you want to work with a 2D refractive index profile, the function you define must take 4 arguments as inputs: `X`, `Y`, `n_background` and `nParameters`. nParameters is a cell array that is often empty but is available for the user to use in whatever way they want. For example, a user might choose to pass in a core width, core refractive index etc. as different elements in the cell array nParameters.

If you want to work with a 3D refractive index, the function must take 5 arguments as inputs: `X`, `Y`, `Z`, `n_background` and `nParameters`.

The syntax of the initializeRIfromFunction() function is:
`P = initializeRIfromFunction(P,hFunc,nParameters,Nz_n);`
Here, `hFunc` is a MATLAB function handle, typically generated with the @ operator (see example 1).
The parameter `Nz_n` is only required when running initializeRIfromFunction() to initialize a 3D refractive index. It describes how many z slices the function should be evaluated at (see example 14). `Nz_n` should be high enough to resolve all the structure in the z direction, but the higher it is the more memory the solver will need and the slower it will run.

##### Method 2, defining `P.n` manually
This method can be used if your refractive index is non-trivial and you want to load it from a data file, for example. In this method you define the refractive index object `P.n` manually. See example 12.
`P.n` is a "BPMmatlab.refractiveIndexProfile" object with five properties:
- `n`
- `Lx`
- `Ly`
- `xSymmetry`
- `ySymmetry`

The `n` property is a 2D or 3D array of values corresponding to the (possibly complex) refractive index. For 3D arrays, the first slice along the third dimension corresponds to the refractive index at z = 0, and the last slice corresponds to that at z = P.Lz.

`Lx` and `Ly` are the side lengths that correpond to the first and second dimensions of the `n` array.

`xSymmetry` and `ySymmetry` have the same meaning as the ones described above for the solver itself, but in this object they describe which symmetry assumptions `n` should be interpreted with. For example, if `xSymmetry` = `ySymmetry` = 'Symmetry', the `P.n` array will be interpreted as being just the part of the refractive index that is in the first quadrant.

#### Initial electric field
Defining the initial electric field is very similar to defining the refractive index, see above. 

As with the refractive index, the solver will store the final electric field of a segment in `P.E`, ready to be automatically used as the input to the next segment.

You are free to redefine the simulation grid from one segment to the next (for example, if one segment requires finer resolution than the others). In that case, the E-field will automatically be interpolated from the old grid to the new grid using the interpn function without the need for user interaction.

There are three ways to define the initial E-field:

##### Method 1, using initializeEfromFunction()
You define a MATLAB function at the end of the model file that provides the 2D (possibly complex) electric field profile as an output. You run the function initializeEfromFunction() as shown in example 1 to make BPM-Matlab apply your function to the simulation grid in preparation for running the solver.

The function must always take 3 arguments: X, Y and Eparameters. It must return a 2D array of electric field values corresponding to the X and Y coordinate locations.

Similar to `nParameters`, `Eparameters` is a cell array that is at the user's disposal to optionally pass additional arguments to the function.

##### Method 2, defining `P.E` manually
As with the refractive index, if your initial electric field is non-trivial and doesn't have a simple analytical expression that you can put into a function, you can define `P.E` manually by setting the properties yourself.

`P.E` is a "BPMmatlab.electricFieldProfile" object with seven properties:
- `field`
- `Lx`
- `Ly`
- `xSymmetry`
- `ySymmetry`
- `label` (optional)
- `neff` (optional)

The `field` property is the complex E-field 2D array in the initial z slice of the segment.
`Lx` and `Ly` are the side lengths that correpond to the first and second dimensions of the `field` array.
`xSymmetry` and `ySymmetry` have the same meaning as the ones described above for the solver itself, but here describes which symmetry assumptions `field` should be interpreted with. For example, if `xSymmetry` = `ySymmetry` = 'Symmetry', the `field` array will be interpreted as being just the part of the electric field that is in the first quadrant.

`label` and `neff` are used only for electric fields found by the findModes() function, see below. They are not necessary for the user to specify.
`label` is a char array that contains a text description of the field, such as 'LP21' for a mode.
`neff` is the effective refractive index of the mode. It isn't used by the solver but provided for the user's convenience.

##### Method 3, using findModes()
This method is perhaps the most common.
BPM-Matlab includes a mode solver that can calculate the supported modes of the waveguide. The mode finder starts its search at effective refractive indices equal to `n_0`, so if the mode finder fails to identify the lower-order modes, try increasing `n_0`. The modes are *not* used for the propagation simulation itself. The mode finder will not return unguided modes. The modes will be labeled with LPlm designations if the refractive index is assessed to be radially symmetric, otherwise the modes will simply be labeled sequentially (Mode 1, 2, 3 etc.).
To use the mode solver, find the modes supported by the waveguide (after you set all the other parameters of `P`, especially `P.n`) by using
`P = findModes(P,nModes);` where `nModes` is the number of modes to try to find.

After `nModes`, three different Name,Value pairs can be added to the list of arguments. They can be any of the following:
- `plotModes`
(Default: `true`)
Set to `false` to tell the solver to not plot the found modes at the end of findModes.
- `sortByLoss`
(Default: `false`)
Set to `true` to sort the list of found modes in order of ascending loss. If `false`, sorts in order of ascending real part of the effective refractive index.
- `singleCoreModes`
(Default: `false`)
If `true`, finds modes for each core individually. Note that the resulting "modes" will only be true modes of the entire structure if the core-to-core coupling is negligible.

findModes will upon completion set `P.modes`, an array of "BPMmatlab.electricFieldProfile" objects with each element correspondng to one mode.

After calling findModes(), you can set `P.E` equal to one of the elements of `P.modes`. For example `P.E = P.modes(1)` will set the initial electric field to the first found field, which is often the fundamental mode, depending on the choice of `n_0`.

A function is provided to make it easy to inject a field that is a superposition of modes. After running findModes(), you can run `P.E = modeSuperposition(P,modeIdxs,coefficients)` with the arguments
- `modeIdxs`
An array containing the indices of the modes you want to superpose, for example [1 4 5] for superposing the first, fourth and fifth modes.
- `coefficients`
An optional array of the same length as `modeIdxs` containing the complex coefficients for each mode, thereby specifying both the amplitude and phase you want for each mode. If not specified, all coefficients will be assumed equal to 1.

### Invoking the solver
After you have defined your model, invoke `P = FD_BPM(P);` to call the actual simulation tool. During the simulation, the display will be updated according to `P.updates`.

#### Optional FFT-BPM solver
In addition to the FD_BPM solver, BPM-Matlab comes with an FFT-BPM solver. FD-BPM should be used for propagation through media with a non-uniform refractive index such as guiding structures, while FFT-BPM can be used for propagation through media with uniform refractive index. Check Example3.m for a comparison.

Unlike FD_BPM, FFT_BPM can step arbitrarily large steps in z without accumulating numerical artifacts, assuming that the beam does not diverge to the edge of the simulation window. If that is the case, however, you will need to put the step size low enough that the absorber layer around the window gets enough steps to get rid of the escpaing light without itself introducing numerical artifacts in the form of "reflections" going back into the simulation window.

Note that many, but not all the parameters listed above for FD_BPM are supported in FFT_BPM.

#### Multiple segments
Each call to `P = FD_BPM(P);` or `P = FFT_BPM(P);` will add a new segment to the simulation that uses as an input E-field and refractive index profiles the outputs of the previous segment, such that you can stack segments one after the other by only re-defining the parameters that changed. Check Example2.m for an example.

## Contribution guidelines
If you want to report a bug or are missing a feature, shoot us an email, see below.

### Who do I talk to?
This repository is maintained by the Technical University of Denmark "Biophotonics Imaging Group" as well as the "Diode Lasers and LED Systems" group.
The main responsibles are Anders Kragh Hansen: ankrh@fotonik.dtu.dk and Madhu Veettikazhy: madve@dtu.dk
