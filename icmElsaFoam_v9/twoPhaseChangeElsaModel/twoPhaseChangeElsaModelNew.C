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

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::twoPhaseChangeElsaModel> Foam::twoPhaseChangeElsaModel::New
(
    const immiscibleIncompressibleTwoPhaseElsaMixture& mixture
)
{
    IOobject twoPhaseChangeElsaModelIO
    (
        IOobject
        (
            phaseChangePropertiesName,
            mixture.U().time().constant(),
            mixture.U().db(),
            IOobject::MUST_READ,
            IOobject::NO_WRITE,
            false
        )
    );

    word modelType(twoPhaseChangeElsaModels::noPhaseElsaChange::typeName);

    if (twoPhaseChangeElsaModelIO.typeHeaderOk<IOdictionary>(false))
    {
        IOdictionary(twoPhaseChangeElsaModelIO).lookup
        (
            twoPhaseChangeElsaModel::typeName
        ) >> modelType;
    }
    else
    {
        Info<< "No phase change: "
            << twoPhaseChangeElsaModelIO.name()
            << " not found" << endl;
    }

    Info<< "Selecting phaseChange model " << modelType << endl;

    dictionaryConstructorTable::iterator cstrIter =
        dictionaryConstructorTablePtr_->find(modelType);

    if (cstrIter == dictionaryConstructorTablePtr_->end())
    {
        FatalErrorInFunction
            << "Unknown " << twoPhaseChangeElsaModel::typeName<< " type "
            << modelType << nl << nl
            << "Valid  twoPhaseChangeElsaModels are : " << endl
            << dictionaryConstructorTablePtr_->sortedToc()
            << exit(FatalError);
    }

    return autoPtr<twoPhaseChangeElsaModel>(cstrIter()(mixture));
}


// ************************************************************************* //
