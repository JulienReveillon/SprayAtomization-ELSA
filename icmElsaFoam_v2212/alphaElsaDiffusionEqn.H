if (activateDiffusion==1){

      Info << "Diffusion term activated" << endl;

      volScalarField dTurb(Dab + alphatab*turbulence->nut());
      surfaceScalarField dTurbSurfaceField(fvc::interpolate(dTurb));


      fvScalarMatrix alpha1Eqn
	(
		  fvm::ddt(alpha1)
		- fvc::ddt(alpha1)
		- fvm::laplacian
      (
		      dTurbSurfaceField*(1-cAlphaFieldF),
		      alpha1,
		      "laplacian(Dab,alpha1)"
		  )
	 );

      alpha1Eqn.relax();
      alpha1Eqn.solve();

      alpha2 = 1.0 - alpha1;

      rhoPhi += alpha1Eqn.flux()*(rho1 - rho2);
}

rho = alpha1*rho1 + alpha2*rho2;
