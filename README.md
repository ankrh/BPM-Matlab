# BPM-Matlab

## What is this repository for?

BPM-Matlab is a numerical simulation tool in which the Douglas-Gunn Alternating Direction Implicit (DG-ADI) method is used to efficiently model the electric field propagation in a wide variety of optical fiber geometries with arbitrary refractive index profiles.

You can find the latest version of BPM-Matlab and more information in the wiki on https://gitlab.gbar.dtu.dk/biophotonics/BPM-Matlab.


## LICENSE

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>.


## How do I get set up?

### Requirements:
- Windows 10
- MATLAB R2018a or newer


## How do I use BPM-Matlab?

### BPM-Matlab function files:
All the functions needed for running BPM-Matlab are located in the folder "BPM-Matlab", and this folder therefore has to be on your MATLAB path when trying to run any BPM-Matlab model file. You will not need to modify any of the distributed files. We recommend keeping the model files in the "BPM-Matlab" folder itself, such that you don't need to manually add "BPM-Matlab" to your MATLAB path.

### Model files:
In BPM-Matlab, you set up your model in a single m-file. You can find a few examples to get you started in the BPM-Matlab folder. Once you're familiar with those, you can start to write your own model files based on those example files.

The following description is based on Example1.m, which provides a minimal working example. There are more advanced parameters introduced in the later example model files. Check those for further explanations on the advanced parameters.
#### General parameters
- `P.name`  
A descriptive name of the model. This will only be used if you choose to automatically save the output of the simulation. Good practice is to use the same name as the model filename.
- `P.useAllCPUs`  
BPM-Matlab will by default leave one processor core unused, which is useful for doing other work on the PC while simulations are running. If you are in need for speed and don't plan to use your PC while doing the simulation, set this parameter to true.
- `P.useGPU`  
This allows to use CUDA acceleration for NVIDIA GPUs. The default is false. Only set this to true if you have a CUDA enabled NVIDIA GPU.

#### Visualization parameters
- `P.updates`  
The number of times during the simulation the plot in the frontend should update. This is useful for following the E-field evolution, but adds overhead to the calculation, as for each update, the currently calculated full simulation window has to be extracted and displayed. If you're only interested in the E-field at the end of the waveguide, set this to 1.
- `P.plotEmax`  
Here you can set the maximung of the color scale in the intensity plot, relative to the peak of initial intensity. If left unset, the intensity plot autoscales.


#### Resolution related paramaters
Be aware that you have to ensure yourself that the pixel size and step size are small enough for the simulation to converge. This is typically done by doing a manual parameter scan of the resolution parameters.

- `P.Lx_main`, `P.Ly_main`  
The physical size of the main simulation window, in meters.
- `P.Nx_main`, `P.Ny_main`  
The number of pixels in the main simulation window.
- `P.dz_target`  
The step size between two calculated simulation windows, in meters.
- `P.padfactor`  
How much absorbing padding to add on the sides of main simulation window. The full simulation window then consists of the main simulation window in the middle and the absorbing padding around it. `P.padfactor` is the ratio of the widths (`Lx`, `Ly`) of the full simulation window to the widths (`Lx_main`, `Ly_main`) of the main simulation window. 1 means no padding, 2 means the full simulation window is twice as wide as the main simulation window, i.e., the absorbing padding is of thickness `Lx_main/2` on both sides.  
The absorbing padding is sometimes required to deal with the non-guided energy: if this hits the boundary of the calculated frame, it would numerically get reflected and propagate back into the main area, thus leading to erroneous results. The absorbing padding deals with this problem by removing energy that has propagated outside the main area from the calculated frame.
- `P.alpha`  
Absorption strength of the absorbing padding. The absorption of the padding is implemented as an absorption coefficient, the magnitude of which is proportional to the distance out from the edge of the main simulation window, squared, such that the energy propagating further out is attenuated stronger in each step. In other words, the field intensity at a distance d out from the main window is multiplied in each step by a `exp(-d^2×alpha×dz)`.

#### Physical properties
- `P.lambda`  
Wavelength, in metres.
- `P.n_cladding`  
Refractive index of the cladding.
- `P.n_0`  
The reference refractive index, introduced during the derivation of the particular form of the paraxial Helmholtz equation used in BPM-Matlab. The refractive indices definede through `P.n_clad` and in `P.shapes` are the actual indices used in the simulation. `P.n_0` (the reference refractive index) has to be chosen such that the paraxial Helmholtz equation remains a good approximation, i.e., close or equal to the refractive indices where the main part of the energy is propagating in the simulation.
- `P.Lz`  
Length of the segment propagated through, in metres.
- `P.shapes`  
In the shapes 2D array, each row is a shape such as a core in a fiber. Column 1 are the x coordinates, column 2 are the y coordinates, column 3 are radii, column 4 are the types of the shapes, column 5 are the peak refractive indices and column 6 is the g parameter, only needed if any of the shapes are GRIN lenses.  
Shape types in column 4 can be one of the following:  
 1. Circular step-index disk
 2. Antialiased circular step-index disk
 3. Parabolic graded index disk
 4. GRIN lens focusing in both x and y
 5. GRIN lens focusing only in y.
- `P.E`  
Struct containing 3 fields definining the initial E-field, or a function handle to calculate this initial E-field.  
In the struct version, the 'field' field is the complex E-field matrix, and the `Lx` and `Ly` fields describe the side lengths of the provided E matrix. In the case of a struct, the provided E-field will be adapted to the new grid using the `interp2` function.  
In the function version, the function should be defined at the end of the model file and taks `X`, `Y` and `Eparameters` as inputs and provide the complex E field as output. `X` and `Y` are the physical locations of the pixels, while `Eparameters` is a cell array that can optionally pass additional arguments to the function.

#### Invoking the solver
`P = FD_BPM(P);`  
This calls the actual simulation tool.

### Compilation
The folder includes all the executables necessary, so you don't need to compile anything. If, however, you want to change the routine in either the FD-BPM.c source code located in the folder "src", you will need to recompile the respective mex-files. Check out the source-file on how to do so.


## Contribution guidelines

If you want to report a bug or are missing a feature, shoot us an email, see below.

### Who do I talk to?
This repository is part of the DTU "biophotonics" team.
The main responsible is Anders Kragh Hansen: ankrh@fotonik.dtu.dk
