/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

Application
    icmElsaFoam

Description
    Solver for 2 incompressible, isothermal immiscible fluids using a VOF
    (volume of fluid) phase-fraction based interface capturing approach,
    with optional mesh motion and mesh topology changes including adaptive
    re-meshing. 
    The solver is changed considering the ELSA formalism. 
    The adaptation to the version 6 of OpenFOAM has been made by A.Remigi and L.Palanti

    Developing Team
   -------------------------------------------------------------------------
   |									   |
   |CORIA UMR 6614, CNRS, University and INSA of Rouen			   |
   |Contributors: F.X. Demoulin*,J. Reveillon*, R. Lebas, B. Duret, J. Anez|
   |*Corresponding authors: demoulin@coria.fr - Julien.Reveillon@coria.fr  |
   |									   |
   -------------------------------------------------------------------------
   
   Reference:

   1) Lebas, Romain, et al. "Numerical simulation of primary break-up and
   atomization: DNS and modelling study."
   International Journal of Multiphase Flow 35.3 (2009): 247-260.

   2) Demoulin, Francois-Xavier, et al. "Toward using direct numerical
   simulation to improve primary break-up modeling."
   Atomization and Sprays 23.11 (2013).

   3)Andreini, Antonio, et al. "Development of a turbulent liquid 
   flux model for Eulerian–Eulerian multiphase flow simulations."
   International Journal of Multiphase Flow 81 (2016): 88-103.
  
   and references Therein


\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "CMULES.H"
#include "EulerDdtScheme.H"
#include "localEulerDdtScheme.H"
#include "CrankNicolsonDdtScheme.H"
#include "subCycle.H"
#include "immiscibleIncompressibleTwoPhaseElsaMixture.H"
#include "noPhaseElsaChange.H"
#include "kinematicMomentumTransportModel.H" //v9, on v7 is turbulentTransportModel
#include "pimpleControl.H"
#include "pressureReference.H"
#include "fvModels.H"  //v9, on v7 is fvOptions
#include "fvConstraints.H" //v9, on v7 does not exist
#include "CorrectPhi.H"
#include "fvcSmooth.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
   #include "postProcess.H"

    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "initContinuityErrs.H"
    #include "createDyMControls.H"
    #include "createFields.H"
    #include "createFieldRefs.H"
    #include "createElsa.H"
    #include "createAlphaFluxes.H"
    #include "initCorrectPhi.H"
    #include "createUfIfPresent.H"

    turbulence->validate();
    
    if (!LTS)
    {
    #include "CourantNo.H"
    #include "setInitialDeltaT.H"
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
    Info<< "\nStarting time loop\n" << endl;

    while (runTime.run())
    {
        #include "readDyMControls.H"

        if (LTS)
        {
            #include "setRDeltaT.H"
        }
        else
        {
        #include "CourantNo.H"
        #include "alphaCourantNo.H"
        #include "setDeltaT.H"    
        }

        runTime++;

        Info<< "Time = " << runTime.timeName() << nl << endl;

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
            if (pimple.firstPimpleIter() || moveMeshOuterCorrectors)
            {
                tmp<volScalarField> divU;                                       
                                                                                
                if                                                              
                (                                                               
                    correctPhi                                                  
                 && !isType<twoPhaseChangeElsaModels::noPhaseElsaChange>(phaseChange)   
                )                                                               
                {                                                               
                    // Construct and register divU for mapping                  
                    divU = new volScalarField                                   
                    (                                                           
                        "divU0",                                                
                        fvc::div(fvc::absolute(phi, U))                         
                    );                                                          
                }                                                               
                                                                                
                fvModels.preUpdateMesh();      
                mesh.update();

                if (mesh.changing())
                {
                    // Do not apply previous time-step mesh compression flux
                    // if the mesh topology changed
                    if (mesh.topoChanging())
                    {
                        talphaPhi1Corr0.clear();
                    }

                    gh = (g & mesh.C()) - ghRef;
                    ghf = (g & mesh.Cf()) - ghRef;

                    MRF.update();

                    if (correctPhi)
                    {
                        // Calculate absolute flux
                        // from the mapped surface velocity
                        phi = mesh.Sf() & Uf();

                        #include "correctPhi.H"

                        // Make the flux relative to the mesh motion
                        fvc::makeRelative(phi, U);
                    }

                    mixture.correct();

                    if (checkMeshCourantNo)
                    {
                        #include "meshCourantNo.H"
                    }
                }
                divU.clear();
            }
            fvModels.correct();                                                 
                                                                                
            surfaceScalarField rhoPhi                                           
            (                                                                   
                IOobject                                                        
                (                                                               
                    "rhoPhi",                                                   
                    runTime.timeName(),                                         
                    mesh                                                        
                ),                                                              
                mesh,                                                           
                dimensionedScalar(dimMass/dimTime, 0)                           
            );   

            #include "alphaControls.H"
            #include "alphaElsaEqnSubCycle.H"
            #include "alphaElsaDiffusionEqn.H"

            mixture.correct();

            #include "sigmaPrimeElsaEqn.H"
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

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
