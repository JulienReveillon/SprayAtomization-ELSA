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

\*---------------------------------------------------------------------------*/

#include "noPhaseElsaChange.H"
#include "fvScalarMatrix.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace twoPhaseChangeElsaModels
{
    defineTypeNameAndDebug(noPhaseElsaChange, 0);
    addToRunTimeSelectionTable(twoPhaseChangeElsaModel, noPhaseElsaChange, dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::twoPhaseChangeElsaModels::noPhaseElsaChange::noPhaseElsaChange
(
    const immiscibleIncompressibleTwoPhaseElsaMixture& mixture
)
:
    twoPhaseChangeElsaModel(typeName, mixture)
{}


// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

Foam::Pair<Foam::tmp<Foam::volScalarField>>
Foam::twoPhaseChangeElsaModels::noPhaseElsaChange::mDotAlphal() const
{
    return Pair<tmp<volScalarField>>
    (
        volScalarField::null(),
        volScalarField::null()
    );
}


Foam::Pair<Foam::tmp<Foam::volScalarField>>
Foam::twoPhaseChangeElsaModels::noPhaseElsaChange::mDotP() const
{
    return Pair<tmp<volScalarField>>
    (
        volScalarField::null(),
        volScalarField::null()
    );
}


Foam::Pair<Foam::tmp<Foam::volScalarField::Internal>>
Foam::twoPhaseChangeElsaModels::noPhaseElsaChange::Salpha
(
    volScalarField& alpha
) const
{
    return Pair<tmp<volScalarField::Internal>>
    (
        tmp<volScalarField::Internal>(nullptr),
        tmp<volScalarField::Internal>(nullptr)
    );
}


Foam::tmp<Foam::fvScalarMatrix>
Foam::twoPhaseChangeElsaModels::noPhaseElsaChange::Sp_rgh
(
    const volScalarField& rho,
    const volScalarField& gh,
    volScalarField& p_rgh
) const
{
    return tmp<fvScalarMatrix>(new fvScalarMatrix(p_rgh, dimVolume/dimTime));
}


Foam::tmp<Foam::fvVectorMatrix>
Foam::twoPhaseChangeElsaModels::noPhaseElsaChange::SU
(
    const volScalarField& rho,
    const surfaceScalarField& rhoPhi,
    volVectorField& U
) const
{
    return tmp<fvVectorMatrix>
    (
        new fvVectorMatrix(U, dimMass*dimVelocity/dimTime)
    );
}


void Foam::twoPhaseChangeElsaModels::noPhaseElsaChange::correct()
{
    twoPhaseChangeElsaModel::correct();
}


bool Foam::twoPhaseChangeElsaModels::noPhaseElsaChange::read()
{
    return twoPhaseChangeElsaModel::read();
}


// ************************************************************************* //
