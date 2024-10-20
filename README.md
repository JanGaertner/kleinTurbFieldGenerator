> [!NOTE]
> When using this software please cite this software package with:
> Gärtner, J. W. Generate turbulent fields with a digital filter (Version 1.0.0) (2024).
> doi: https://doi.org/10.5281/zenodo.13956967 Available: https://github.com/JanGaertner/kleinTurbFieldGenerator.git

[![DOI](https://zenodo.org/badge/875547596.svg)](https://doi.org/10.5281/zenodo.13956967)


# Turublent Fluctuations using Digital Filters

This tool generates turbulent fluctuations with the proposed model of 
Klein et al. [1]. However, not for an inlet boundary condition but for a 
given field. 

It is particularly useful for cases such as homogeneous isotropic turbulence (HIST) or shear layer simulations where turbulence needs to be introduced throughout the domain.

The tool also supports defining specific regions or blending functions, allowing turbulence to be applied selectively within the computational domain. These fluctuations are then superimposed onto the base velocity field, U, to simulate the desired turbulence.

## Download and Install

This tool is tested for OpenFOAM v2312, however it should also work with the 
earlier versions of OpenFOAM. You can download and install the tool in two
easy steps:

 1. Clone the repository from GitHub with: 
    ```bash
    git clone https://github.com/JanGaertner/kleinTurbFieldGenerator.git
    ```
 2. Compile the code by executing wmake within the repository:
    ```bash
    wmake src/
    ```

## Usage 

All settings of the tool are set in the kleinTurbFieldGeneratorDict dictionary
in the case system/ folder. An example file is given below:

```
/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  2312                                  |
|   \\  /    A nd           | Website:  www.openfoam.com                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/
FoamFile
{
    version     2.0;
    format      ascii;
    arch        "LSB;label=32;scalar=64";
    class       dictionary;
    location    "system";
    object      kleinTurbFieldGeneratorDict;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// Set the turbulent properties
Lt              6E-4;   // Turbulent length scale
turbIntensity   0.05;   // Turbulent intensity
nFilterCells       8;   // number of neighboring filter cells -- see [1]
UMean       (30 0 0);   // Mean velocity of the field

// Optional Parameters
// ===================

// Provide a cellSet, e.g., generated with topoSet, in which the 
// fluctuations are applied
cellSet shearLayer;

// Alternatively one can provide a filtering function which will be multiplied
// with the generated turb. fluctuations
// UPrime = blendingFunction * UPrime

filteringFunction
{
    type coded;         // All OpenFOAM Function1 types are supported
    filterDir 1;        // Direction of the filter, e.g., a blending in y-axis
                        // 0 : x-axis
                        // 1 : y-axis
                        // 2 : z-axis
    
    // Example for a blending in y direction provided by the coded
    // environment. Note: x variable is always the input to the function
    // Only 1-D inputs are allowed.
    code
    #{
        const scalar H = 2.88E-3;
        const scalar deltaU = 30;
        const scalar delta = 5;//(deltaU/H/10.25);
        scalar f =  0.5*(std::tanh(((x+0.5*H)/H)*delta) - std::tanh(((x-0.5*H)/H)*delta));
        if (f<1E-6)
            f = 0;
        return f;
    #};
}
```
## License
This OpenFOAM library is under the GNU General Public Licensce.

## References
[1] M. Klein, A. Sadiki, J. Janicka, "A digital filter based generation of 
        inflow data for spatially developing direct numerical or large eddy 
        simulations", Journal of Computational Physics, 2003, vol. 186