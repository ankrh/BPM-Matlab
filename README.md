# BPM-Matlab

## What is this repository for?

BPM-Matlab is a numerical simulation tool in which the Douglas-Gunn Alternating Direction Implicit (DG-ADI) method is used to efficiently model the electric field propagation using a Finite-Difference Beam Propagation Method (FD-BPM) in a wide variety of optical fiber geometries with arbitrary refractive index profiles.

You can find the latest version of BPM-Matlab and more information in the wiki on https://gitlab.gbar.dtu.dk/biophotonics/BPM-Matlab.


## LICENSE

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>.


## How do I get set up?

### Requirements:
- This software was tested on Windows 10 and macOS 11
- MATLAB R2018a or newer


## How do I use BPM-Matlab?

### BPM-Matlab function files:
All the functions needed for running BPM-Matlab are located in the folder "BPM-Matlab", and this folder therefore has to be on your MATLAB path when trying to run any BPM-Matlab model file. You will not need to modify any of the distributed files. We recommend keeping the model files in the "BPM-Matlab" folder itself, such that you don't need to manually add "BPM-Matlab" to your MATLAB path.

### Model files:
In BPM-Matlab, you set up your model in a single m-file. You can find a few examples to get you started in the BPM-Matlab folder. Once you're familiar with those, you can start to write your own model files based on those example files.

#### Segments
If the model contains multiple segments, BPM-Matlab will add a new segment with each call to the FD_BPM() solver function, using as the input E-field the simulation output of the previous segment, see towards the end of this file. As such, you can stack segments one after the other by only re-defining the parameters that change. Check Example2.m for an example.

#### General parameters
- `P.name`  
A descriptive name of the model. This will only be used if you choose to automatically save the output of the simulation. Good practice is to use the same name as the model filename.
- `P.useAllCPUs`  
BPM-Matlab will by default leave one processor core unused, which is useful for doing other work on the PC while simulations are running. If you are in need for speed and don't plan to use your PC while doing the simulation, set this parameter to true.
- `P.useGPU`  
This allows to use CUDA acceleration for NVIDIA GPUs. The default is false. Only set this to true if you have a CUDA enabled NVIDIA GPU.

##### Visualization parameters
- `P.updates`  
The number of times during the simulation the plot in the frontend should update. This is useful for following the E-field evolution, but adds overhead to the calculation, as for each update, the currently calculated full simulation window has to be extracted and displayed. If you're only interested in the E-field at the end of the waveguide, set this to 1.
- `P.plotEmax`  
Here you can set the maximum of the color scale in the intensity plot, relative to the peak of initial intensity. If left unset, the intensity plot autoscales.

##### Resolution related paramaters
Be aware that you have to ensure yourself that the pixel size and z step size are small enough for the simulation to converge. This is typically done by doing a manual parameter scan of the resolution parameters.

- `P.Lx_main`, `P.Ly_main`  
The physical size of the main simulation window, in meters.
- `P.Nx_main`, `P.Ny_main`  
The number of pixels in the main simulation window.
- `P.dz_target`  
The z step size between two calculated simulation windows, in meters.
- `P.padfactor`  
How much absorbing padding to add on the sides of main simulation window. The full simulation window then consists of the main simulation window in the middle and the absorbing padding around it. `P.padfactor` is the ratio of the widths (`Lx`, `Ly`) of the full simulation window to the widths (`Lx_main`, `Ly_main`) of the main simulation window. 1 means no padding, 2 means the full simulation window is twice as wide as the main simulation window, i.e., the absorbing padding is of thickness `Lx_main/2` on both sides.  
The absorbing padding is sometimes required to deal with the non-guided energy: if this hits the boundary of the calculated frame, it would numerically get reflected and propagate back into the main area, thus leading to erroneous results. The absorbing padding deals with this problem by removing energy that has propagated outside the main area from the calculated frame.
- `P.alpha`  
Absorption strength of the absorbing padding. The absorption of the padding is implemented as an absorption coefficient, the magnitude of which is proportional to the distance out from the edge of the main simulation window, squared, such that the energy propagating further out is attenuated more strongly. In other words, the field intensity at a distance d out from the main window is multiplied in each step of length dz by a `exp(-d^2×alpha×dz)`.

##### Physical properties
- `P.lambda`  
Wavelength, in metres.
- `P.n_background`  
Refractive index of the background, typically corresponding to the cladding.
- `P.n_0`  
The reference refractive index, introduced during the derivation of the particular form of the paraxial Helmholtz equation used in BPM-Matlab. The refractive indices defined through `P.n_clad` and in `P.shapes` are the actual indices used in the simulation. `P.n_0` (the reference refractive index) has to be chosen such that the paraxial Helmholtz equation remains a good approximation, i.e., close or equal to the refractive indices where the main part of the energy is propagating in the simulation.
- `P.Lz`  
Length of the segment propagated through, in metres.

There are five ways to define the refractive index profile:
1) Defining the P.shapes field in which different rows correspond to different geometric shapes, such as circular fiber cores.
2) Defining a function in P.n.func that takes X, Y, n_background and nParameters as inputs and provides the 2D refractive index as an output. You'd typically put this function definition at the end of the model file.
3) Defining a 2D array of refractive indices in P.n.n with x and y side lengths specified in P.n.Lx and P.n.Ly. This array could, e.g., be imported from a data file.
4) Defining a function in P.n.func that takes X, Y, Z, n_background and nParameters as inputs and provides the 3D refractive index as an output. P.n.Lx, P.n.Ly, P.n.Nx, P.n.Ny and P.n.Nz are the side lengths and resolutions of the window in which the refractive index function will be evaluated. The z side length will always be P.Lz. As with method 2, you'd usually put this function definition at the end of the model file.
5) Defining a 3D array of refractive indices in P.n.n with x and y side lengths specified in P.n.Lx and P.n.Ly. The first xy slice corresponds to the refractive index exactly at z = 0, and the last slice to z = P.Lz. As with method 3, this is useful for refractive index data imported from a file.

- `P.shapes`  
Shapes is a 2D array that describes the refractive index distribution in terms of geometric shapes. Each row is a shape such as a circular core in a fiber. Column 1 are the x coordinates, column 2 are the y coordinates, column 3 are radii, column 4 are the types of the shapes, column 5 are the peak refractive indices and column 6 is the g parameter, only needed if any of the shapes are GRIN lenses.  
Shape types in column 4 is an integer and can be one of the following:  
 1. Circular step-index disk
 2. Antialiased circular step-index disk
 3. Parabolic graded index disk
 4. GRIN lens focusing in both x and y
 5. GRIN lens focusing only in y.
- `P.n`  
This is an alternative way of defining the refractive index to `P.shapes` that allows more control of the refractive index. If `P.shapes` isn't defined, `P.n` must be. When using `P.n`, the user must define either a `P.n.n` field or a `P.n.func` field.
- `P.n.n`
A 2D or 3D array containing the (complex) refractive index values. For a 3D array, it is assumed to stretch in the z direction from `z = 0` to `z = P.Lz`. The spacing between points in the x and y directions is `P.n.Lx/P.n.Nx` and `P.n.Ly/P.n.Ny`, but the spacing in the z direction is `P.Lz/(P.n.Nz - 1)` since the first and last slices are taken to be exactly at `z = 0` and `z = P.Lz`.
- `P.n.Lx` and `P.n.Ly`
The x and y widths of the provided `P.n.n` data OR in the case of a 3D function in `P.n.func`, the side lengths of the window in which the refractive index function will be evaluated.
- `P.n.Nx`, `P.n.Ny` and `P.n.Nz`
The x, y and z resolutions that a 3D `P.n.func` function will be evaluated with. High values of P.n.Nx*P.n.Ny*P.n.Nz will require large amounts of memory (CPU or GPU).

- `P.E`  
There are three ways to define the initial E-field:
1) Defining P.E as a function that takes X, Y and Eparameters as inputs and provides the complex E-field as output. You'd typically define this function at the end of the model file. `Eparameters` is a cell array that can optionally pass additional arguments to the function.
2) Defining P.E as a struct with 3 fields: a 'field' field which is the complex E-field matrix, and 'Lx' and 'Ly' fields that describe the side lengths of the provided E matrix. In the case of a struct, the provided E-field will be adapted to the new grid using the interpn function.
3) After calling findModes(), setting P.E equal to one of the elements of P.modes, for example P.E = P.modes(1) for the fundamental mode. A superposition field can be created using the modeSuperposition() function.

#### Advanced parameters
##### Bending, Tapering, Twisting
The fibres can be bent, tapered and twisted as follows:

- `P.taperScaling`  
The ratio of the width of the structure at the end of the segment to the width at the beginning of the segment.
- `P.twistRate`  
The rate of rotation in units of radians per meter.
- `P.bendDirection`  
Direction of the bending in degrees, in a polar coordinate system with 0° to the right (towards positive x) and increasing angles in counterclockwise direction.
- `P.bendingRoC`  
Radius of curvature of the bend in meters.

##### Video saving
Optionally, the figure can be saved to a video during simulation, check Example6.m on how to do that. The following parameters can be used:

- `P.saveVideo`  
Set this to `true` to record the field intensity and phase profiles at different transverse planes in preparation of saving to a video. Default is `false`.
- `P.videoName`  
Set the name of the video to be saved.
- `P.finalizeVideo`  
This controls whether the video should be finalised and saved after the current segment. Set this to `true` if the current segment is the last segment in your simulation, and to `false` otherwise.

##### Colormaps of the graphs

- `P.Intensity_colormap`, `P.Phase_colormap`, `P.n_colormap`  
The integer assigned defines the colormap of the respective subplot. The options are  
 1. GPBGYR
 2. HSV
 3. parula
 4. gray
 5. cividis

### Invoking the solver
After you have defined your model, invoke `P = FD_BPM(P);` to call the actual simulation tool. During the simulation, the display will be updated according to `P.updates`.

#### Multiple segments
Each call to `P = FD_BPM(P);` will add a new segment to the simulation that uses as an input E-field the output of the previous segment, such that you can stack segments one after the other by only re-defining the parameters that changed. Check Example2.m for an example.

#### Optional FFT-BPM solver
BPM-Matlab comes with an FD-BPM and an FFT-BPM solver. FD-BPM should be used for propagation through media with a non-uniform refractive index, while FFT-BPM should be used for propagation through media with uniform refractive index. Check Example3.m for a comparison.

#### Optional Mode solver
BPM-Matlab includes a mode solver that can calculate the supported modes of the waveguide and calculate the overlap of the simulated E-field with these modes. Check Example9.m on how to do so. The modes are *not* used for the propagation simulation itself.
To use the mode solver, find the modes supported by the waveguide (after you set all the other parameters of P, especially P.shapes) by using
`P = findModes(P,nModes,sortByLoss,plotModes);` with

- `nModes`  
number of modes to find
- `plotModes`  
set to `true` if the solver should plot the found modes at the end of findModes
- `sortByLoss`
Set to `true` to sort the list of found modes in order of ascending loss. If `false`, sorts in order of ascending imaginary part of eigenvalue (descending propagation constant)
- `singleCoreModes`
If `true` and if `P.shapes` is defined, finds modes for each core/shape individually. Note that the resulting "modes" will only be true modes of the entire structure if the core-to-core coupling is negligible.

This will set `P.modes`, a struct array with each elemnt correspondng to one mode. Some of these may be non-guided (radiating). You may then set `P.calcModeOverlaps` to `true` to calculate mode overlap integrals of propagating field with respect to the different modes that were set in the `P.modes` struct array.

A function is provided to make it easy to inject a field that is a superposition of modes. After running `findModes`, you can run `P.E = modeSuperposition(P,modeIdxs,coefficients)` with the arguments
- `modeIdxs`
An array containing the indices of the modes you want to superpose, for example [1 4 5] for superposing the first, fourth and fifth modes.
- `coefficients`
An optional array of the same length as `modeIdxs` containing the complex coefficients for each mode, thereby specifying both the amplitude and phase you want for each mode. If not specified, all coefficients will be assumed equal to 1.


### Compilation
The folder includes all the executables necessary, so you don't need to compile anything. If, however, you want to change the routine in the FD-BPM.c source code located in the folder "src", you will need to recompile the respective mex-files. Check out the source-file on how to do so.


## Contribution guidelines

If you want to report a bug or are missing a feature, shoot us an email, see below.

### Who do I talk to?
This repository is part of the Technical University of Denmark "biophotonics" team.
The main responsibles are Anders Kragh Hansen: ankrh@fotonik.dtu.dk and Madhu Veettikazhy: madve@dtu.dk
