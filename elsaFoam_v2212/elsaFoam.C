/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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
    elsaFoam

Group
    grpMultiphaseSolvers

Description
    Solver for mixing 2 incompressible fluids.                                  
                                                                                
    Turbulence modelling is generic, i.e. laminar, RAS or LES may be selected.  
                                                                                
Developing Team                                                                 
                                                                                
------------------------------------------------------------------------------- 
                                                                                
    |CORIA UMR 6614, CNRS, University and INSA of Rouen                         
                                                                                
    |Contributors: F.X. Demoulin*,J. Reveillon*, R. Lebas, B. Duret, J. Anez,   
     D. Ferrando                                                                
                                                                                
    |*Corresponding authors: demoulin@coria.fr  - Julien.Reveillon@coria.fr        
                                                                                
------------------------------------------------------------------------------- 
                                                                                
    Reference:                                                                  
                                                                                
    1) Lebas, Romain, et al. "Numerical simulation of primary break-up and         
                                                                                
       atomization: DNS and modelling study."                                   
                                                                                
       International Journal of Multiphase Flow 35.3 (2009): 247-260.           
                                                                                
    2) Demoulin, Francois-Xavier, et al. "Toward using direct numerical         
                                                                                
       simulation to improve primary break-up modeling."                        
                                                                                
       Atomization and Sprays 23.11 (2013).                                     
                                                                                
                                                                                
       and references Therein 
\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "MULES.H"
#include "subCycle.H"
#include "incompressibleTwoPhaseMixture.H"
#include "turbulentTransportModel.H"
#include "pimpleControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    argList::addNote
    (
        "Solver for mixing two incompressible fluids"
    );

    #include "postProcess.H"

    #include "addCheckCaseOptions.H"
    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "initContinuityErrs.H"
    #include "createFields.H"
    #include "createElsa.H"
    #include "createTimeControls.H"
    #include "CourantNo.H"
    #include "setInitialDeltaT.H"

    turbulence->validate();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readTimeControls.H"
        #include "CourantNo.H"
        #include "alphaCourantNo.H"
        #include "setDeltaT.H"

        ++runTime;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        mixture.correct();

        #include "alphaEqnSubCycle.H"
        #include "alphaDiffusionEqn.H"
        #include "sigmaPrimeElsaEqn.H"

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            #include "UEqn.H"

            // --- Pressure corrector loop
            while (pimple.correct())
            {
                #include "pEqn.H"
            }

            if (pimple.turbCorr())
            {
                turbulence->correct();
            }
        }

        runTime.write();

        runTime.printExecutionTime(Info);
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
