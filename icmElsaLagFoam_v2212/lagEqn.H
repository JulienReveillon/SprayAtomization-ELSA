scalar NewDiam,NewVolume,NewSurface,NewNp;
scalar diamMax = 0;


forAll(mesh.cells(),celli)
{
    coeffAlpha[celli] = 0;
    particleCell[celli] = 0;
    volumeLag[celli] = 0;
    surfaceLag[celli] = 0;
    nbParticle[celli] = 0;
    nbVolume[celli] = 0;
    nbSurface[celli] = 0;
    Nbgoutte[celli] = 0;
}

forAll(alpha1,cellI)
{
  Ulag[cellI]=vector(0,0,0);
  USigmaGaz[cellI]=vector(0,0,0);
}

//Reconstruction des champs Volume,Surface,Vitesse via les particules
forAllIter(basicKinematicCollidingCloud,parcels2,pIter)
{
  basicKinematicParcel p = pIter();
  volumeLag[p.cell()] += PI/6*(p.d()*p.d()*p.d()*p.nParticle());
  surfaceLag[p.cell()] += p.areaS()*p.nParticle();
  nbParticle[p.cell()]++;
}

//Maj du diam et nParticle de particule afin d'avoir concordance avec la partie Eulerienne
forAllIter(basicKinematicCollidingCloud,parcels2,pIter)
{
  basicKinematicParcel p = pIter();
  scalar cellI  = p.cell();
  if (p.d() >diamMax) {diamMax = p.d();}
  if (alpha1[cellI]>0.011||alpha1[cellI]<0.000001)
    {
      parcels2.deleteParticle(pIter());
      tagDiffusion[cellI]=3;
    }
  else
  {
    tagDiffusion[cellI]=2;
    // Info << "2.1"<< endl;
    NewVolume = alpha1[cellI]*(mesh.cellVolumes()[cellI]/volumeLag[cellI])*(p.volume()*p.nParticle()) ;
    // Info << NewVolume << endl;
    NewSurface = ((sigma[cellI]+0.1)*mesh.cellVolumes()[cellI]/surfaceLag[cellI])*(p.areaS()*p.nParticle());
    //Info << NewSurface << endl;
    NewDiam = 6* NewVolume/NewSurface ;
    //Info << NewDiam << endl;
    NewNp = Foam::pow(NewSurface,3.)/(pi*Foam::pow(6*NewVolume,2.));
    //Info << NewNp << endl;
    //Info << " Diam : " << NewDiam << "         Nparticle : " << NewNp << endl;
    //ratio = alpha1[cellI]*mesh.cellVolumes()[cellI]/volumeLag[cellI];
    //NewDiam = p.d()*Foam::pow(ratio,1./3.);
    if (NewDiam > 1e-5) {
      Info <<"Warning PARTICLE TROP GROSSE :"<< NewDiam << "   ratio  : " << NewDiam/p.d() << " " <<  NewNp <<  endl;
      Info << "Alpha :  "<< alpha1[cellI] << "   Sigma :" << sigma[cellI] << " NBParticule :" << nbParticle[cellI] << endl;}

    pIter().d()=NewDiam;
    pIter().nParticle()=NewNp;
  }
}

Info <<"Diametre max particule :  "<< diamMax << endl;

forAll(mesh.cells(),celli)
{
    volumeLag[celli] = 0;
    surfaceLag[celli] = 0;
    particleCell[celli] = 0;
}

//Nouveau calcul du volume/Surface pour vérifier que notre code marche
// Calcul de la de la vitesse lagrangienne
forAllIter(basicKinematicCollidingCloud,parcels2,pIter)
{
  basicKinematicParcel p = pIter();
  volumeLag[p.cell()] += (PI*p.d()*p.d()*p.d()*p.nParticle()/(6.*mesh.cellVolumes()[p.cell()]));
  surfaceLag[p.cell()] += p.areaS()*p.nParticle()/mesh.cellVolumes()[p.cell()];
  Ulag[p.cell()] += p.U()*PI/6*(p.d()*p.d()*p.d())*p.nParticle();
  USigmaGaz[p.cell()] += p.U()*PI*(p.d()*p.d())*p.nParticle();
  particleCell[p.cell()]+=p.nParticle();
  nbVolume[p.cell()]+=PI/6*(p.d()*p.d()*p.d())*p.nParticle();
  nbSurface[p.cell()]+=PI*(p.d()*p.d())*p.nParticle();
}

//Calcule de Ulag :

forAll(alpha1,cellI)
{
  if (nbParticle[cellI]!=0 && particleCell[cellI]!= 0 )
    {
      errorVolume[cellI] = ( (volumeLag[cellI]-alpha1[cellI]) / alpha1[cellI] )*100;
      //Info << "errorVolume : " <<  errorVolume[cellI] << endl;

      errorSurface[cellI] = ( (surfaceLag[cellI]-sigma[cellI]) / sigma[cellI] )*100;
      //Info << "errorSurface : " <<  errorSurface[cellI] << endl;

      coeffAlpha[cellI]= max(1 - 1./Foam::sqrt(nbParticle[cellI]) , 0);

      Ulag[cellI]=Ulag[cellI]/nbVolume[cellI];
      Ugaz[cellI]=(U[cellI]-Ulag[cellI]);
      USigmaGaz[cellI]=U[cellI]-(USigmaGaz[cellI]/nbSurface[cellI]);
    }
}

coeffAlphaF = fvc::interpolate(coeffAlpha);

phiSigmaLag = fvc::interpolate(USigmaGaz*coeffAlpha) & mesh.Sf();
phiAlphaLag = fvc::interpolate(Ugaz*coeffAlpha*(1-cAlphaField)) & mesh.Sf();
