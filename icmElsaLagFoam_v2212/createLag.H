scalar nbParticleParcel = 10 ;
scalar PI = Foam::constant::mathematical::pi;
sigmaKineticEnergy = Foam::sqrt(2.0 * (turbulence->k()) / 3.0);

forAll (mesh.cells(),cellI)
{
// these two limits for alpha were selected to handle from one-way up to four-way coupling
//if (alpha1[cellI] > 0.0000001 && alpha1[cellI]<0.01 && compTerm[cellI] == 0 && nbParticle[cellI]==0 )

//==========fx proposal to include only fastest particles====
//scalar var( mag(rnd.GaussNormal()) );

//if (var > 2)
//{
//if (alpha1[cellI] > 0.000001 && alpha1[cellI]<0.01 && compTerm[cellI] == 0 && nbParticle[cellI]==0 )
 //   if (alpha1[cellI] > 0.0001 && alpha1[cellI] < 0.01 &&
 //       tagDiffusion[cellI] == 1 && nbParticle[cellI] == 0 )

    if (alpha1[cellI] > 0.0001 && alpha1[cellI] < 0.01 &&
        tagDiffusion[cellI] == 1 && nbParticle[cellI] == 0 )
    {
  // if (alpha1[cellI] > 0.00001 && alpha1[cellI]<0.011 && nbParticle[cellI]==0 )
      // creation of one parcel with the volume of the liquid within the cell
        scalar NewVolume = alpha1[cellI] * mesh.cellVolumes()[cellI];
        scalar NewSurface = (sigma[cellI] + 0.1) * mesh.cellVolumes()[cellI];

        scalar Diam = 6* NewVolume / NewSurface ;
        nbParticleParcel = Foam::pow(NewSurface,3.)/(PI*Foam::pow(6*NewVolume,2.));

        particleCell[cellI]+=nbParticleParcel;
        nbParticle[cellI]+=1;

        basicKinematicCollidingParcel* ptr = new basicKinematicCollidingParcel
        (
            mesh,
            mesh.C()[cellI],
            cellI
        );

        ptr->d() = Diam;
        //===============================
        // adding a fluctuating velocity to the mean particle velocity based on RANS
        //- Random angle between 0 - 2*PI
        const doubleScalar theta = rnd.sample01<scalar>()*2.0*PI;

        //- Random number between [-1,1]
        const doubleScalar u = 2*rnd.sample01<scalar>() - 1;

        doubleScalar a = Foam::sqrt(1 - Foam::pow(u,2));

        vector dir(a*Foam::cos(theta), a*Foam::sin(theta), u);

        ptr->U() = UaveLiq[cellI] + sigmaKineticEnergy[cellI]*mag(rnd.GaussNormal<scalar>())*dir;

        //========================================
        ptr->rho()=rho1.value();

        ptr->nParticle()=nbParticleParcel;

        parcels2.addParticle(ptr);
        //===============================

        //	delete ptr;
    }
}
