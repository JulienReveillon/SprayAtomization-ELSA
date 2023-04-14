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
    sigmaMin= min(sigmaMin, sigmaEquil * 0.001);

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

}
