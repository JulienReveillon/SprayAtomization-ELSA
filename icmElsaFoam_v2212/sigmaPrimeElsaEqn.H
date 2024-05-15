{
    alpha1.max(0.0);

    alpha1.min(1.0);

    Info << "Computing sigma prime equation" << endl;

    sigmaMin.primitiveFieldRef() = 
    (
        cstElsaSigmaMin * sqrt(alpha1*(1.-alpha1)) *pow(mesh.V(),-1.0/3.0)
    );

    if ( webEquilMod.value() < 0.5)
    {
    // WeberEquil Lebas 2009 IJMF
        sigmaEquil.primitiveFieldRef() = 
        (
            rho * turbulence->k() / elsaWeberEquil / elsaSurfaceTensionCoef * 
            alpha1 * (1.-alpha1)
        );
        
        Info << "mod Lebas" << endl;
    }
     
    else
    {
     // Modified WeberEquil Benjamin Duret et al IJMF 2013
        sigmaEquil.primitiveFieldRef() = 
        (
            2. * (rho1 + rho2) * turbulence->k() / elsaWeberEquil /
            elsaSurfaceTensionCoef * alpha1 * (1. - alpha1)
        );
      
        Info << "mod Duret "  << endl;
    }

    sigmaEquil.max(0.0);

    sigma.primitiveFieldRef() = 
    (
       sigmaPrime.primitiveField() + sigmaMin.primitiveField()
    );

    sigma.primitiveFieldRef() = max(sigma.primitiveField() , 0.0);
    
    sigmaRatio.primitiveFieldRef() = max(
    (
        sigma.primitiveField() / 
        (sigmaEquil.primitiveField() + sigmaTiny.value()) 
    ), 
    1.);

    sigmaRatio.boundaryFieldRef() = sigmaPrime.boundaryField();

    fvScalarMatrix sigmaPrimeEqn
    (
        fvm::ddt(sigmaPrime)
      + fvm::div(phi,sigmaPrime)
      - fvm::laplacian
        (
            Dab + alphatab * turbulence->nut(), sigmaPrime,
            "laplacian(Dab,sigmaPrime)"
        )
     == (
            (cTurbElsaSigmaPrime * 
            (turbulence->k() / (elsaTinyNut + turbulence->nut())) * 
            (sigmaRatio)) * (sigmaEquil - sigmaMin)
            + fvm::Sp(-cTurbElsaSigmaPrime * 
            (turbulence->k() / (elsaTinyNut + turbulence->nut())) *
            (sigmaRatio), sigmaPrime)
        )
    );

    Info << "Solve elsaSigmaPrime ..." << endl;
    sigmaPrimeEqn.relax();
    sigmaPrimeEqn.solve();

    sigmaPrime.max(0.0); //?????

    sigma.primitiveFieldRef() = 
    (
        sigmaPrime.primitiveField() + sigmaMin.primitiveField()
    );

    D32=6.0*alpha1*(1.0-alpha1) / (sigma + elsaTiny);


    // ================ IRQ by gravity. done by JAP ========================

    // cell size
    volScalarField cbrtCellVol (pow(cellVol, 1.0/3.0));

    dimensionedScalar deltaN (1e-8/pow(average(cellVol),1.0/3.0));

    // Cell gradient of alpha
    volVectorField gradAlpha ( fvc::grad(alpha1) );

    // sigmaMin dense
    sigmaMinDense.primitiveFieldRef() = mag( gradAlpha );

    // Face unit interface normal
    volVectorField nHatfv ( gradAlpha / (mag(gradAlpha) + deltaN) ) ;
    
    // curavture
    curvature.primitiveFieldRef() = mixture.K();

    // turbulent kinetic energy
    tKE.primitiveFieldRef() = turbulence->k();

    // const dimensionedScalar& rhoConstant = ( (1/rho1) - (1/rho2) ) / (6*rho1);
    const dimensionedScalar rhoConstant = ((1.0/rho1)-(1.0/rho2)) / (6.0*rho1);

    // g*normal_vector_interface
    gNormalInt.primitiveFieldRef() = ( ( g/mag(g) ) & nHatfv );

    //IRQg without cell size nor constants
    magGradAlpha.primitiveFieldRef() = (pow(rho,2)*alpha1*(1.0-alpha1)*
                                       gNormalInt*mag(g)) /
                                       (tKE+tKETiny);

    // IRQ by grad of alpha
    IRQalpha.primitiveFieldRef() = mag( gradAlpha ) * cbrtCellVol;

    //IRQg
    IRQg.primitiveFieldRef() = -1*magGradAlpha*rhoConstant*cbrtCellVol;

    kappaDeltaX.primitiveFieldRef() = curvature*pow(mesh.V(),1.0/3.0);
    
    forAll(mesh.cells(),celli)
    {
      if (alpha1[celli] * (1-alpha1[celli]) < 1e-3 )
      {
        gNormalInt[celli] = 0;
        magGradAlpha[celli] = 0;
        IRQg[celli] = 1; //IRQg > 0.25 --> stabilization --> ICM
        IRQsigma[celli] = 1;
        cAlphaField[celli]=1;
      }
      else
      {
        IRQsigma[celli] = (sigmaMin[celli])/ (sigma[celli]);
        if (gNormalInt[celli] < 0)
            {
                IRQg[celli] = 0.0; //IRQ<g0.25 --> ELSA
            }
      }
    }

    forAll(mesh.cells(),celli)
    {
    //if ( IRQsigma[celli] < 0.5 || mag(kappaDeltaX[celli]) > 1 )
    
    if ( IRQsigma[celli] > 0.5 ) // --> ICM
      {
          cAlphaField[celli] = mixture.cAlpha(); // mixture.cAlpha==1 in general 
      }
    else // ELSA or ICM depend on IRQg
      {
          if (IRQg[celli] < 0.25 ) // ELSA
              
              cAlphaField[celli] = 0;

          else //>0.25 --> ICM
              cAlphaField[celli] = mixture.cAlpha();
      }

    }

    activateDiffusion=1;

    //=====to get variables at the faces.
    IRQgF = (fvc::interpolate(IRQg)) ;

    IRQalphaF = fvc::interpolate(IRQalpha) ;


    //forAll(mesh.faces(),faceI)
    IRQsigmaF = (fvc::interpolate(IRQsigma));

    kappaDeltaXF= (fvc::interpolate(kappaDeltaX));

    cAlphaFieldF = fvc::interpolate(cAlphaField);

    // loop for internal faces only. The internal faces are indexed first
    label nInternalFaces = mesh.nInternalFaces();

    for(label faceI=0; faceI < nInternalFaces; faceI++)
    {

      //if (IRQsigmaF[faceI] < 0.5 || mag(kappaDeltaXF[faceI]) > 1  ) 
      if (IRQsigmaF[faceI] > 0.5 ) 
        
        {
          cAlphaFieldF[faceI] = mixture.cAlpha();
        }

      else
        {
          if (IRQgF[faceI]<0.25 )
            {
              cAlphaFieldF[faceI] = 0;
            }
          else
              cAlphaFieldF[faceI] = mixture.cAlpha();
        }

    }

    //=======end==============================================

    IRQk.primitiveFieldRef() = 1/(mag(kappaDeltaX)+IRQkTiny);

}
