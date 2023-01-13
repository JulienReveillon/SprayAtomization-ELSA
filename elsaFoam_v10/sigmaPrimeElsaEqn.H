{
     alpha1.max(0.0);

     alpha1.min(1.0);

     Info << "compute elsaSigmaMin ..." << endl;

     elsaSigmaMin.primitiveFieldRef() = cstElsaSigmaMin * 
                                        sqrt(alpha1 * (1.0 - alpha1)) * 
                                        pow(mesh.V(), -1.0/3.0);

     Info << "compute elsaSigmaEquil ..." << endl;

     if ( webEquilMod.value() < 0.5)
     {
        // WeberEquil Lebas 2009 IJMF
        elsaSigmaEquil.primitiveFieldRef() = rho * turbulence->k() / 
                                             elsaWeberEquil / 
                                             elsaSurfaceTensionCoef * 
                                             alpha1 * (1.0 - alpha1);
        Info << "mod Lebas" << endl;
     }
     else
     {
        // Modified WeberEquil Benjamin Duret et al IJMF 2013
        elsaSigmaEquil.primitiveFieldRef() = 2.0 * (rho1 + rho2) *
                                             turbulence->k() / elsaWeberEquil /
                                             elsaSurfaceTensionCoef * 
                                             alpha1 * (1.0 - alpha1);
        Info << "mod Duret "  << endl;
     }

     elsaSigmaEquil.max(0.0);


     elsaSigmaMin= min(elsaSigmaMin, elsaSigmaEquil * 0.001);

     Info << "compute Elsa Sigma ..." << endl;

     elsaSigma.primitiveFieldRef() = elsaSigmaPrime.primitiveField() + 
                                     elsaSigmaMin.primitiveField();

     elsaSigma.primitiveFieldRef() = max(elsaSigma.primitiveField(), 0.0);


     Info << "compute Elsa Ratio ..." << endl;

     elsaSigmaRatio.primitiveFieldRef() = (max(
                                      (elsaSigma / (elsaSigmaEquil+sigmaTiny)),
                                      1.0));
     elsaSigmaRatio.boundaryFieldRef()=elsaSigmaPrime.boundaryField();

     Info << "... done" << endl;

     fvScalarMatrix elsaSigmaPrimeEqn
     (
        fvm::ddt(elsaSigmaPrime)
      + fvm::div(phi,elsaSigmaPrime)
      - fvm::laplacian
        (
            Dab + alphatab*turbulence->nut(), elsaSigmaPrime,
            "laplacian(Dab,elsaSigmaPrime)"
        )
      == (
            cTurbElsaSigmaPrime * (turbulence->k() / 
            (elsaTinyNut + turbulence->nut())) * (elsaSigmaRatio)
         ) * 
         (
            elsaSigmaEquil - elsaSigmaMin
         ) +
         fvm::Sp(
            -cTurbElsaSigmaPrime * (turbulence->k() / 
            (elsaTinyNut + turbulence->nut()))*(elsaSigmaRatio),elsaSigmaPrime
            )
     );

     Info << "Solve elsaSigmaPrime ..." << endl;
     elsaSigmaPrimeEqn.relax();
     elsaSigmaPrimeEqn.solve();

     elsaSigmaPrime.max(0.0);

     elsaSigma = elsaSigmaPrime + elsaSigmaMin;
     Info << "...done" <<endl;
     elsaD32=6.0*alpha1*(1.0-alpha1)/(elsaSigmaMin+elsaSigmaPrime+elsaTiny);

}
