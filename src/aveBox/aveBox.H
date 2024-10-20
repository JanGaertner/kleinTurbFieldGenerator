/*---------------------------------------------------------------------------*\
    888b     d888                                  d8888                   
    8888b   d8888                                 d88888                   
    88888b.d88888                                d88P888                   
    888Y88888P888  .d88b.  888  888  .d88b.     d88P 888 888  888  .d88b.  
    888 Y888P 888 d88""88b 888  888 d8P  Y8b   d88P  888 888  888 d8P  Y8b 
    888  Y8P  888 888  888 Y88  88P 88888888  d88P   888 Y88  88P 88888888 
    888   "   888 Y88..88P  Y8bd8P  Y8b.     d8888888888  Y8bd8P  Y8b.     
    888       888  "Y88P"    Y88P    "Y8888 d88P     888   Y88P    "Y8888   
-------------------------------------------------------------------------------                                                                                                                                                    
License
    This file is part of movingAverage.

    movingAverage is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    movingAverage is distributed in the hope that it will be useful, but 
    WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY 
    or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with movingAverage.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::aveBox

Description
    Provide information of the LES filter cell

SourceFiles
    aveBox.C

Author
    Jan Wilhelm Gärtner <jan.gaertner@outlook.de> Copyright (C) 2020>

\*---------------------------------------------------------------------------*/

#ifndef aveBox_H 
#define aveBox_H



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#include "volFields.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                          Class aveBox Declaration
\*---------------------------------------------------------------------------*/

class aveBox
{
    private:
    //- Private Data
        
        //- Dimensions of the box
        scalar width_;    // Always x coordinate
        scalar length_;   // Always y coordinate
        scalar height_;   // Always z coordinate

public:

    // Constructor
       
       //- Construct with three different entries
       aveBox(const scalar width, const scalar length, const scalar height) :
            width_(width),
            length_(length),
            height_(height)
            {};

        //- Construct with all the same entry for an equal length box
        aveBox(const scalar length) : width_(length),length_(length),height_(length) {};

        //- Construct from dictionary 
        //  Length must be read. If height and width are not specified, 
        //  construct equal length box
        aveBox(const dictionary& dict);

        virtual ~aveBox() {};

    // Accessor functions for member variables as const reference
        
        //- Return dimensions
        inline scalar width() const {return width_;}
        inline scalar length() const {return length_;}
        inline scalar height() const {return height_;}

    // Helper functions

        //- Return true if the point p is in the box with the center point
        //- pCenter
        virtual bool containsPoint(const point& pCenter, const point& p) const;

    // I-O

        //- Print box dimensions to screen
        virtual void printDim() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
