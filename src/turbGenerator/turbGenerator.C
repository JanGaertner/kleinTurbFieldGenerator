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

Application
    turbGenerator

Group
    grpIncompressibleSolvers

Description
    Add Turbulent Pertubations to a given velocity field
    Based on the Klein turbulence generator,
    [1] M. Klein, A. Sadiki, J. Janicka, "A digital filter based generation of 
        inflow data for spatially developing direct numerical or large eddy 
        simulations", Journal of Computational Physics, 2003, vol. 186
    
    Important: This tool assumes a uniform spacing of the grid
    
    Parameters
    ----------
    uMean : np.array
        The mean velocity field
    Lt : float
        Turbulent integral length scale
    turbIntensity : float
        Turbulent intensity, e.g., 0.05 for 5% turb. intensity
    mesh : fvMesh
        The fvMesh object to operate on
    deltaX : float
        The cell width, assuming a uniform mesh
    filteWidth : int
        Sets the filter width, by default it is set to,
        nFilterCells = (2*Lt)/delta
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



int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Generate turbulent fluctuations on a given mean velocity field"
    );

    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"

    // Read the Settings from dicionary 
    // ================================

    IOdictionary settings
    (
        IOobject
        (
            "turbGenerator",
            runTime.system(),
            runTime,
            IOobject::MUST_READ,
            IOobject::NO_WRITE
        )
    );

    // Get the turbulent length scale
    const scalar Lt = settings.get<scalar>("Lt");
    const scalar turbIntensity  = settings.get<scalar>("turbIntensity");

    // Given number of cells for the filter width, or multiples of deltaX
    const label nFilterCells = settings.get<label>("nFilterCells");

    // Read mean velocity
    const vector UMean = settings.get<vector>("UMean");

    // Check if a cellZone is defined
    const word cellSetName = settings.getOrDefault<word>("cellSet","none");

    List<label> activeCells;

    if (cellSetName != "none")
    {
        Info << "Cell Set active: "<<cellSetName << endl;
        // Try to load the cellSet
        cellSet set(mesh,cellSetName);

        activeCells.resize(set.size());
        label i=0;
        for (auto& cellI : set)
        {
            activeCells[i++] = cellI;
        }
        Info << "Adding turb fluctuations to " << activeCells.size() << " cells" <<endl;
    }

    // Add a filtering function
    // ========================
    // Provide a custom function which returns for each cell a scalar filter value
    // such that: U = U + filter * UPrime
    // The filter function is evaluated only for a single coordinate direction,
    // thus the direction has to be speciefied with the direction component
    // filterDir    0; // x-dir
    // filterDir    1; // y-dir
    // filterDir    2; // z-dir
    autoPtr<Function1<scalar>> filteringFunction;
    scalar filterDir = -1;
    if (settings.findDict("filteringFunction") != nullptr)
    {
        auto filterDict = settings.subDict("filteringFunction");
        filterDir = filterDict.get<scalar>("filterDir");
        filteringFunction = 
            Function1<scalar>::New
            (
                "filteringFunction",
                settings,
                &(mesh.thisDb())
            );
    }




    // Start Computation
    // =================

    // Estimate the cell width from the cubic root of the volume
    const auto& V = mesh.V();
    const scalar deltaX  = std::cbrt(V[0]);
    Info << "cell width calculated to be: "<<deltaX << endl;

    // Check that nFilterCells is set large enough
    if (nFilterCells * deltaX < 2 * Lt)
        Info << "Warning: nFilterCells*deltaX < 2*Lt (" 
             << nFilterCells * deltaX << " < "<<  2 * Lt <<")"<< endl;

    // Read the velocity field
    // =======================

    Info<< "Reading field U\n" << endl;
    const volVectorField U
    (
        IOobject
        (
            "U",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );


    // Calculate the filter coefficients based on Eq. (14) and Eq. (13) of [1]
    // Calculate for twice the filter width
    List<scalar> bkCoeff(2*nFilterCells+1);
    scalar sumCoeffSqr = 0;
    // Loop over -nFilterCells : nFilterCells
    forAll(bkCoeff,i)
    {
        bkCoeff[i] = 
            std::exp
            (
                -(M_PI*std::pow((i-nFilterCells)*deltaX,2))
                /(4.0*std::pow(Lt,2))
            );
        sumCoeffSqr += std::pow(bkCoeff[i],2);
    }

    for (scalar& coeff :  bkCoeff)
        coeff /= std::sqrt(sumCoeffSqr);

    // Initiate a field of random numbers to store the UPrime fluctuation for each 
    // cell. We substract the rnd numbers in the range -0.5:0.5
    std::random_device rnd_dev;
    std::mt19937 generator(rnd_dev());
    std::uniform_real_distribution<scalar> distr(-0.5,0.5);

    
    volVectorField UPrime
    (
        IOobject
        (
            "UPrime",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
        dimensionedVector("U",dimVelocity,vector(0,0,0))
    );

    // Initialize UPrime with random numbers -- See Eq. (9) and Eq. (16) of [1]
    forAll(UPrime,i)
    {
        UPrime[i].x() = distr(generator);
        UPrime[i].y() = distr(generator);
        UPrime[i].z() = distr(generator);
    }


    Info << "Calc UPrime with neighboring cells" << endl;

    // Loop over all three directions, x,y, and z
    // Represents the i,j, and k of Eq. (16)
    for (label dirI=0; dirI < 3; dirI++)
    {
        const volVectorField UPrimeBkp("UPrimeBkp",UPrime);

        // Loop over all cells
        forAll(UPrime,i)
        {
            // Get the neighboring cells in the direction
            vector dir(0,0,0);
            dir[dirI] = 1;

            List<label> posNeighborCells = getNeighborCellsInAxis
            (
                mesh,
                nFilterCells,
                dir,
                i
            );

            dir[dirI] = -1;
            List<label> negNeighborCells = getNeighborCellsInAxis
            (
                mesh,
                nFilterCells,
                dir,
                i
            );

            // Loop over all velocity components
            for (label cmptI=0; cmptI < 3; cmptI++)
            {
                UPrime[i][cmptI] = UPrimeBkp[i][cmptI] * bkCoeff[nFilterCells];

                forAll(negNeighborCells,k)
                {
                    const label& cellI = negNeighborCells[k];
                    // Index of the bkCoeff array
                    const label bkCoeffInd = nFilterCells-k-1;
                    UPrime[i][cmptI] += UPrimeBkp[cellI][cmptI]*bkCoeff[bkCoeffInd];
                }

                forAll(posNeighborCells,k)
                {
                    const label& cellI = posNeighborCells[k];
                    // Index of the bkCoeff array
                    const label bkCoeffInd = nFilterCells+k+1;
                    UPrime[i][cmptI] += UPrimeBkp[cellI][cmptI]*bkCoeff[bkCoeffInd];
                }
            }
        }
    }

    // Calculate the mean and standard deviation from UPrime
    vector meanUPrime(0,0,0);
    forAll(UPrime,i)
    {
        meanUPrime.x() += UPrime[i].x();
        meanUPrime.y() += UPrime[i].y();
        meanUPrime.z() += UPrime[i].z();
    }
    meanUPrime /= UPrime.size();

    vector sigma(0,0,0);
    forAll(UPrime,i)
    {
        sigma.x() += std::pow(UPrime[i].x()-meanUPrime.x(),2);
        sigma.y() += std::pow(UPrime[i].y()-meanUPrime.y(),2);
        sigma.z() += std::pow(UPrime[i].z()-meanUPrime.z(),2);
    }
    sigma.x() = std::sqrt(sigma.x()/UPrime.size());
    sigma.y() = std::sqrt(sigma.y()/UPrime.size());
    sigma.z() = std::sqrt(sigma.z()/UPrime.size());

    // Make the mean and standard deviation to 0 and 1
    forAll(UPrime,i)
    {
        UPrime[i].x() = (UPrime[i].x()-meanUPrime.x())/sigma.x();
        UPrime[i].y() = (UPrime[i].y()-meanUPrime.y())/sigma.y();
        UPrime[i].z() = (UPrime[i].z()-meanUPrime.z())/sigma.z();
    }

    Info << "Adding U" << endl;

    // Apply Lund transform to cofGorm the field to the specified energy
    // i.e. Reynolds stress
    symmTensor R;
    R.xx() = std::pow(turbIntensity*mag(UMean),2);
    R.yy() = std::pow(turbIntensity*mag(UMean),2);
    R.zz() = std::pow(turbIntensity*mag(UMean),2);

    Tensor<scalar> a;
    a.xx() = std::sqrt(R.xx());
    a.xy() = 0;
    a.xz() = 0;
    a.yx() = R.yx()/a.xx();
    a.yy() = std::sqrt(R.yy()-std::pow(a.yx(),2));
    a.yz() = 0;
    a.zx() = R.zx()/a.xx();
    a.zy() = (R.zy()-a.yx()*a.zx())/a.yy();
    a.zz() = std::sqrt(R.zz() - std::pow(a.zx(),2) - std::pow(a.zy(),2));

    forAll(U,i)
    {
        UPrime[i] = a & UPrime[i];
    }


    // Check the computed fields
    // =========================

    // Compute the turbulent intensity
    scalar calcTurbIntensity = 0;
    forAll(UPrime,i)
    {
        calcTurbIntensity += std::sqrt
        (
            (1.0/3.0)*
            (
                std::pow(UPrime[i].x(),2)
              + std::pow(UPrime[i].y(),2)
              + std::pow(UPrime[i].z(),2)
            )
        );
    }
    calcTurbIntensity /= UPrime.size();
    calcTurbIntensity /= mag(UMean);
    Info << "Set turb. intensity: " << turbIntensity*100.0 << "% calculated: " 
         << calcTurbIntensity*100.0 << endl;

    // Calculate the turblent length scale 
    const vector L = calcTurbLengthScale(UPrime,nFilterCells,deltaX,mesh);
    Info << "Set turb. length scale: " << Lt << " calculated: " << L << endl; 

    // Write out the solution
    // ======================

    if (filteringFunction.valid())
    {
        forAll(UPrime,i)
        {
            const point pos = mesh.C()[i];
            UPrime[i] = UPrime[i] * filteringFunction->value(pos[filterDir]);
        }
    }

    // Add the UPrime only to the defined cells in activeCells list -- if set
    if (activeCells.size() > 0)
    {
        volVectorField UNew("UNew",U);
        for (const label& cellI : activeCells)
            UNew[cellI] += UPrime[cellI];
        
        UNew.write();
    }
    else
    {
        volVectorField UNew("UNew",U + UPrime);
        UNew.write();
    }
    UPrime.write();

    return 0;
}


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
