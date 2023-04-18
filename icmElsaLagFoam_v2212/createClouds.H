Info<< "\nConstructing cloud, create cloud routine start" << endl;


// Read kinematic viscosity
dimensionedScalar nu2("nu", mixture.nuModel2().viscosityProperties());

// create dynamic viscoisty
volScalarField mu
    (
        IOobject
        (
            "mu",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
	rho2*nu2
    );


Info<< "Constructing kinematicCloud " << endl;
basicKinematicCollidingCloud parcels2
(
"particle",
rho,
U,
mu,
g
);

Info<< "Construction de particleCell" << endl;

volScalarField particleCell
(
 IOobject
 (
  "particleCell",
  runTime.timeName(),
   mesh,
  IOobject::NO_READ,
  IOobject::AUTO_WRITE
  ),
 alpha1
 );

volScalarField volumeLag
(
 IOobject
 (
  "volumeLag",
  runTime.timeName(),
   mesh,
  IOobject::NO_READ,
  IOobject::AUTO_WRITE
  ),
 alpha1
 );

volScalarField surfaceLag
(
 IOobject
 (
  "surfaceLag",
  runTime.timeName(),
   mesh,
  IOobject::NO_READ,
  IOobject::AUTO_WRITE
  ),
 alpha1
 );

volScalarField nbParticle
(
 IOobject
 (
  "nbParticle",
  runTime.timeName(),
  mesh,
  IOobject::NO_READ,
  IOobject::AUTO_WRITE
  ),
 mesh,
 dimensionedScalar("nbParticle",dimless,0),
 "zeroGradient"
);

volVectorField Ulag
(
 IOobject
 (
  "Ulag",
  runTime.timeName(),
   mesh,
  IOobject::NO_READ,
  IOobject::AUTO_WRITE
  ),
  U
 );

volVectorField Ugaz
(
 IOobject
 (
  "Ugaz",
  runTime.timeName(),
   mesh,
  IOobject::NO_READ,
  IOobject::AUTO_WRITE
  ),
  U
 );

volVectorField USigmaGaz
(
 IOobject
 (
  "USigmaGaz",
  runTime.timeName(),
   mesh,
  IOobject::NO_READ,
  IOobject::AUTO_WRITE
  ),
  U
 );

//------------------------------
//Initialize the Lagrangian velocity

forAll (alpha1,cellI)
{
  Ulag[cellI]=vector(0,0,0);
  Ugaz[cellI]=vector(0,0,0);
  USigmaGaz[cellI]=vector(0,0,0);
}

//----------------------------------------------
volScalarField errorVolume
(
 IOobject
 (
  "errorVolume",
  runTime.timeName(),
   mesh,
  IOobject::NO_READ,
  IOobject::AUTO_WRITE
  ),
 alpha1
 );


volScalarField errorSurface
(
 IOobject
 (
  "errorSurface",
  runTime.timeName(),
   mesh,
  IOobject::NO_READ,
  IOobject::AUTO_WRITE
  ),
 alpha1
 );

volScalarField nbVolume
(
 IOobject
 (
  "NpVolume",
  runTime.timeName(),
   mesh,
  IOobject::NO_READ,
  IOobject::AUTO_WRITE
  ),
 alpha1
 );

volScalarField nbSurface
(
 IOobject
 (
  "nbSurface",
  runTime.timeName(),
   mesh,
  IOobject::NO_READ,
  IOobject::AUTO_WRITE
  ),
 alpha1
 );



//================
 // given seed to initialize the random class
Random rnd(1);
//================
