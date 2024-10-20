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

\*---------------------------------------------------------------------------*/

#include "supportFunctions.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::List<Foam::label> Foam::getCellsInBox
(
    const fvMesh& mesh,
    const aveBox& box, 
    const label cellID
)
{
    // Found neighboring cell IDs
    List<label> cellList;
    std::set<label> visitedCells;

    // List of all cell centers
    const volVectorField& cellCenters = mesh.C(); 
    // Loop over all neighbors
    const labelListList& neighborIDList = mesh.cellCells();

    // Center point of the cellID
    const point& centerPoint = cellCenters[cellID];

    std::stack<label> cellsToProcess;

    cellsToProcess.push(cellID);
    visitedCells.insert(cellID);

    while (!cellsToProcess.empty())
    {
        const label i = cellsToProcess.top();
        cellsToProcess.pop();

        // Get neighboring cells for the current cell
        const List<label>& neighbors = neighborIDList[i];
        for (const label& nID : neighbors)
        {
            // Is the neighbor cell in the box
            if (box.containsPoint(centerPoint,cellCenters[nID]))
            {
                // if it is, has it already been visited?
                // Returns a pair, if second entry is true, it was not yet
                // in the visitedCells set.
                auto res = visitedCells.insert(nID);

                // if it has not been visited yet add it to the cells to be
                // processed
                if (res.second)
                {
                    cellsToProcess.push(nID);
                }
            }
        }
    }

    // Do not add the cellID to the cellList, only neighbors
    cellList.resize(visitedCells.size()-1);
    label i = 0;
    for (const label& ind : visitedCells)
    {
        if (ind != cellID)
            cellList[i++] = ind;
    }

    return cellList;
}


Foam::List<Foam::label> Foam::getNeighborCellsInAxis
(
    const fvMesh& mesh,         // base mesh
    const label nNeighborCells, // How many neighbor cells to find
    const vector dir,           // Direction given as a vector
    const label cellID          // Current cell ID to start from
)
{
    std::stack<label> foundCells;

    // Get cell centers
    const auto&  cellCenters = mesh.C();

    // Get list of faces
    const auto& faces = mesh.faces();

    // Add the current cell to foundCells
    foundCells.push(cellID);

    // Boolean flag to check if boundary patch has been detected
    bool boundaryDetected = false;

    while (foundCells.size() <= nNeighborCells && !boundaryDetected)
    {
        // Current cell to process
        const label currCellID = foundCells.top();
        
        // Get current cell center
        const vector cellCenter = cellCenters[currCellID];

        // Get the current cell
        const auto& cell = mesh.cells()[currCellID];
        
        const auto oldSize = foundCells.size();

        // Loop over the cell faces and find the face in the direction
        for (const label& fi : cell)
        {
            const face& f = faces[fi];

            // Get the unit normal vector
            const vector fCentre = f.centre(mesh.points());
            vector dist = fCentre - cellCenter;

            // normalize vector
            dist = dist/mag(dist);

            // If more than 90% is pointing in the dir, it is the searched face
            if ((dist & dir) > 0.9)
            {
                // get the neighboring cell
                label ownerCell = -1;
                if (fi < mesh.owner().size())
                    ownerCell = mesh.owner()[fi];

                label neighborCell = -1;
                if (fi < mesh.neighbour().size())
                    neighborCell = mesh.neighbour()[fi];
                
                if (neighborCell != -1 && neighborCell != currCellID)
                    foundCells.push(neighborCell);
                else if (ownerCell != -1 && ownerCell != currCellID)
                    foundCells.push(ownerCell);
                // If there is no neighbor cell the boundary patch is reached
                // and the search-loop breaks
                else
                    boundaryDetected = true;
                
                break;
            }
        }

        // If no more cells can be added
        if (oldSize == foundCells.size())
            break;
    }

    List<label> foundCellsList(foundCells.size()-1);
    // Add in reverse order so that the closest cell is in the first entry
    for (label i=foundCellsList.size()-1; i >= 0; i--)
    {
        foundCellsList[i] = foundCells.top();
        foundCells.pop();
    }

    return foundCellsList;
}


Foam::vector Foam::calcTurbLengthScale
(
    const volVectorField& U,
    const label n,
    const scalar deltaX,
    const fvMesh& mesh
)
{
    // Calculate the two point correlation and fluctuation for each point 
    // and n cells
    
    // Integral length scale L
    vector L(0,0,0);

    // Loop over each velocity component
    for (label cmptI=0; cmptI < 3; cmptI++)
    {
        // It is quite memory intensive, but in a first step we get all 
        // neighboring cells of each cell and store them in a list
        List<Pair<List<label>>> cellNeighbors(U.size());

        forAll(U,i)
        {
            // Get the cells in that direction
            // Get the neighboring cells in the direction
            vector dir(0,0,0);
            dir[cmptI] = 1;

            Pair<List<label>>& posAndNegNeighbors = cellNeighbors[i];

            posAndNegNeighbors.first() = getNeighborCellsInAxis
            (
                mesh,
                n,
                dir,
                i
            );

            dir[cmptI] = -1;
            posAndNegNeighbors.second() = getNeighborCellsInAxis
            (
                mesh,
                n,
                dir,
                i
            );
        }

        // Calculate the USqrMean for the comptI        
        scalar USqrMean = 0;
        forAll(U,i)
        {
            USqrMean += U[i][cmptI]*U[i][cmptI];
        }
        USqrMean /= U.size();

        // Loop over all cells to calculate the correlations
        for (label k=0; k < n; k++)
        {
            scalar UCorrMean = 0;
            label nEntries = 0;
            forAll(U,i)
            {
                const scalar u0 = U[i][cmptI];

                // Get the positive neighbor
                if (k < cellNeighbors[i].first().size())
                {
                    const scalar u1 = U[cellNeighbors[i].first()[k]][cmptI];
                    UCorrMean += u0*u1;
                    nEntries++;
                }


                // Get the positive neighbor
                if (k < cellNeighbors[i].second().size())
                {
                    const scalar u1 = U[cellNeighbors[i].second()[k]][cmptI];
                    UCorrMean += u0*u1;
                    nEntries++;
                }
            }
            UCorrMean /= nEntries;

            L[cmptI] += UCorrMean/USqrMean * deltaX;
        }
    }
    
    return L;
}



// ************************************************************************* //
