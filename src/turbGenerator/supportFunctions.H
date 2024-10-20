/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Description
    Support functions for the turbulent generator tool.

Author
    Jan Wilhelm Gärtner <jan.gaertner@outlook.de>
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include <random>
#include <set>
#include <stack>
#include "aveBox.H"
#include "cellSet.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    // Get all cell IDs in the given box, defined by its length, and the central
    // cell
    List<label> getCellsInBox
    (
        const fvMesh& mesh,
        const aveBox& box,
        const label cellID
    );

    // Collect all neighboring cells in given dimension/axis, e.g., all 
    // neighboring cells on the x-axis
    List<label> getNeighborCellsInAxis
    (
        const fvMesh& mesh,         // base mesh
        const label nNeighborCells, // How many neighbor cells to find
        const vector dir,           // Direction given as a vector
        const label cellID          // Current cell ID to start from
    );

    // Calculate the turbulent length scale based on the auto-correlation
    // S. B. Pope. Turbulent flows. Cambridge University Press, Cambridge, 2000.
    // ISBN 0521598869.
    vector calcTurbLengthScale
    (
        const volVectorField& U,
        const label n,
        const scalar deltaX,
        const fvMesh& mesh
    );

} // namespace Foam


// ************************************************************************* //
