/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2021 OpenFOAM Foundation
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

Class
    Foam::twoPhaseChangeElsaModels::noPhaseElsaChange

Description
    Dummy phase change model for 'no phase change'

SourceFiles
    noPhaseElsaChange.C

\*---------------------------------------------------------------------------*/

#ifndef noPhaseElsaChange_H
#define noPhaseElsaChange_H

#include "twoPhaseChangeElsaModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace twoPhaseChangeElsaModels
{

/*---------------------------------------------------------------------------*\
                                 Class noPhaseElsaChange
\*---------------------------------------------------------------------------*/

class noPhaseElsaChange
:
    public twoPhaseChangeElsaModel
{

public:

    //- Runtime type information
    TypeName("none");


    // Constructors

        //- Construct for mixture
        noPhaseElsaChange
        (
            const immiscibleIncompressibleTwoPhaseElsaMixture& mixture
        );


    //- Destructor
    virtual ~noPhaseElsaChange()
    {}


    // Member Functions

        //- Return the mass condensation and vaporisation rates as a
        //  coefficient to multiply (1 - alphal) for the condensation rate
        //  and a coefficient to multiply  alphal for the vaporisation rate
        virtual Pair<tmp<volScalarField>> mDotAlphal() const;

        //- Return the mass condensation and vaporisation rates as coefficients
        //  to multiply (p - pSat)
        virtual Pair<tmp<volScalarField>> mDotP() const;

        //- Return unallocated explicit and implicit sources
        //  for the phase-fraction equation
        virtual Pair<tmp<volScalarField::Internal>> Salpha
        (
            volScalarField&
        ) const;

        //- Return the zero source matrix
        //  for the p_rgh pressure equation
        virtual tmp<fvScalarMatrix> Sp_rgh
        (
            const volScalarField& rho,
            const volScalarField& gh,
            volScalarField& p_rgh
        ) const;

        //- Return the zero source matrix for the momentum equation
        virtual tmp<fvVectorMatrix> SU
        (
            const volScalarField& rho,
            const surfaceScalarField& rhoPhi,
            volVectorField& U
        ) const;

        //- Correct the noPhaseElsaChange phaseChange model
        virtual void correct();

        //- Read the transportProperties dictionary and update
        virtual bool read();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace twoPhaseChangeElsaModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
