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
- Windows 7 or newer
- MATLAB R2018a or newer


## How do I use BPM-Matlab?

### BPM-Matlab function files:
All the functions needed for running BPM-Matlab are located in the folder "BPM-Matlab", and this folder therefore has to be on your MATLAB path when trying to run any BPM-Matlab model file. You will not need to modify any of the distributed files. We recommend keeping the model files in the "BPM-Matlab" folder itself, such that you don't need to manually add "BPM-Matlab" to your MATLAB path.

### Model files:
In BPM-Matlab, you set up your model in a single m-file. You can find a few examples to get you started in the BPM-Matlab folder. Once you're familiar with those, you can start to write your own model files based on those example files.

Each model file must include at least the following definitions:
#### General parameters
- P.name  
A descriptive name of the model. Good practice is to use the same name as the model filename.
- P.updates  
The number of times during the simulation the plot in the frontend should update. This is useful for following the E-field evolution, but adds overhead to the calculation, as for each update, the currently calculated simulation frame has to be extracted and displayed. If you're only interested in the E-filed at the end of the waveguide, set this to 1.

#### Resolution related paramaters
You have to ensure that the pixel size and step size are small enough for the simulation to converge.
- P.Lx_main, P.Ly_main  
The physical size of the main area in the calculated frame, in meters.
- P.Nx_main, P.Ny_main  
The number of pixels in the main area in the calculated frame.
- P.dz_target  
The step size between two calculated frames, in meters.
- P.padfactor  
How much absorbing padding to add on the sides of main area. The calculated frame then consists of the main area in the middle and the absorbing padding around it. 1 means no padding, 2 means the absorbing padding on both sides is of thickness Lx_main/2.
- P.alpha  
Absorption strength of the absorbing padding, in squared unit length distance out from edge of the main area.

#### Physical properties
- P.lambda  
Wavelength, in metres.
- P.n_cladding  
Refractive index of the cladding.
- P.n_0  
"Base" refractive index.
- P.Lz  
Length of this segment, in metres.
- P.shapes  
In the shapes 2D array, each row is a shape such as a core in a fiber. Column 1 are the x coordinates, column 2 are the y coordinates, column 3 are radii, column 4 are the types of the shapes, column 5 are the peak refractive indices and column 6 is the g parameter, only needed if any of the shapes are GRIN lenses.  
Shape types are 1: Circular step-index disk, 2: Antialiased circular step-index disk, 3: Parabolic graded index disk, 4: GRIN lens focusing in both x and y, 5: GRIN lens focusing only in y.
- P.E  
Struct containing 3 fields definining the initial E-field, or a function handle to calculate this initial E-field.  
In the struct version, the 'field' field is the complex E-field matrix, and the 'Lx' and 'Ly' fields describe the side lengths of the provided E matrix. In the case of a struct, the provided E field will be adapted to the new grid using the interp2 function.  
In the function version, the function should be defined at the end of the model file and taks X, Y and Eparameters as inputs and provide the complex E field as output. X and Y are the physical locations of the pixel, while Eparameters is a cell array that can optionally pass additional arguments to the function.

#### Invoking the solver
P = FD_BPM(P);  
This calls the actual simulation tool.

### Compilation
The folder includes all the executables necessary, so you don't need to compile anything. If, however, you want to change the routine in either the FD-BPM.c source code located in the folder "src", you will need to recompile the respective mex-files. Check out the source-file on how to do so.


## Contribution guidelines

If you want to report a bug or are missing a feature, shoot us an email, see below.

### Who do I talk to?
This repository is part of the DTU "biophotonics" team.
The main responsible is Anders Kragh Hansen: ankrh@fotonik.dtu.dk
