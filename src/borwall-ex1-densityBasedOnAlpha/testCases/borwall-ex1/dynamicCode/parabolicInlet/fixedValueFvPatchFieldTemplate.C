/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) YEAR AUTHOR,AFFILIATION
     \\/     M anipulation  |
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

#include "fixedValueFvPatchFieldTemplate.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "unitConversion.H"
//{{{ begin codeInclude

//}}} end codeInclude


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

//{{{ begin localCode

//}}} end localCode


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

extern "C"
{
    // dynamicCode:
    // SHA1 = 55f9e0c59db90711cbd8fd8a33cb506e13a7643e
    //
    // unique function name that can be checked if the correct library version
    // has been loaded
    void parabolicInlet_55f9e0c59db90711cbd8fd8a33cb506e13a7643e(bool load)
    {
        if (load)
        {
            // code that can be explicitly executed after loading
        }
        else
        {
            // code that can be explicitly executed before unloading
        }
    }
}

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

makeRemovablePatchTypeField
(
    fvPatchVectorField,
    parabolicInletFixedValueFvPatchVectorField
);


const char* const parabolicInletFixedValueFvPatchVectorField::SHA1sum =
    "55f9e0c59db90711cbd8fd8a33cb506e13a7643e";


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

parabolicInletFixedValueFvPatchVectorField::
parabolicInletFixedValueFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(p, iF)
{
    if (false)
    {
        Info<<"construct parabolicInlet sha1: 55f9e0c59db90711cbd8fd8a33cb506e13a7643e"
            " from patch/DimensionedField\n";
    }
}


parabolicInletFixedValueFvPatchVectorField::
parabolicInletFixedValueFvPatchVectorField
(
    const parabolicInletFixedValueFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<vector>(ptf, p, iF, mapper)
{
    if (false)
    {
        Info<<"construct parabolicInlet sha1: 55f9e0c59db90711cbd8fd8a33cb506e13a7643e"
            " from patch/DimensionedField/mapper\n";
    }
}


parabolicInletFixedValueFvPatchVectorField::
parabolicInletFixedValueFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<vector>(p, iF, dict)
{
    if (false)
    {
        Info<<"construct parabolicInlet sha1: 55f9e0c59db90711cbd8fd8a33cb506e13a7643e"
            " from patch/dictionary\n";
    }
}


parabolicInletFixedValueFvPatchVectorField::
parabolicInletFixedValueFvPatchVectorField
(
    const parabolicInletFixedValueFvPatchVectorField& ptf
)
:
    fixedValueFvPatchField<vector>(ptf)
{
    if (false)
    {
        Info<<"construct parabolicInlet sha1: 55f9e0c59db90711cbd8fd8a33cb506e13a7643e"
            " as copy\n";
    }
}


parabolicInletFixedValueFvPatchVectorField::
parabolicInletFixedValueFvPatchVectorField
(
    const parabolicInletFixedValueFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(ptf, iF)
{
    if (false)
    {
        Info<<"construct parabolicInlet sha1: 55f9e0c59db90711cbd8fd8a33cb506e13a7643e "
            "as copy/DimensionedField\n";
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

parabolicInletFixedValueFvPatchVectorField::
~parabolicInletFixedValueFvPatchVectorField()
{
    if (false)
    {
        Info<<"destroy parabolicInlet sha1: 55f9e0c59db90711cbd8fd8a33cb506e13a7643e\n";
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void parabolicInletFixedValueFvPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (false)
    {
        Info<<"updateCoeffs parabolicInlet sha1: 55f9e0c59db90711cbd8fd8a33cb506e13a7643e\n";
    }

//{{{ begin code
    #line 0 "/home/ajit/Desktop/borwall-example1/src/borwall-ex1-densityBasedOnAlpha/testCases/borwall-ex1/0/U.boundaryField.inlet"
vectorField faceCenters(this->patch().Cf());
            vectorField inletVal(faceCenters.size(), Zero);

            forAll(inletVal, faceI)
            {
                scalar yi = faceCenters[faceI].y();

                inletVal[faceI].x() = 0.1*4*yi*(1-yi);
            }
            operator==(inletVal);
//}}} end code

    this->fixedValueFvPatchField<vector>::updateCoeffs();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

