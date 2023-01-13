    volScalarField elsaRho
    (
        IOobject
        (
            "elsaRho",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        rho
    );

    Info<< "Reading field elsaSigmaPrime\n" << endl;
    volScalarField elsaSigmaPrime
    (
        IOobject
        (
            "elsaSigmaPrime",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh
    );

    volScalarField elsaSigmaMin
    (
        IOobject
        (
            "elsaSigmaMin",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        elsaSigmaPrime
    );

    volScalarField elsaSigma
    (
        IOobject
        (
            "elsaSigma",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        elsaSigmaPrime
    );

    volScalarField elsaSigmaEquil
    (
        IOobject
        (
            "elsaSigmaEquil",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        elsaSigmaPrime
    );

    volScalarField elsaD32
    (
        IOobject
        (
            "elsaSD32",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
       dimensionSet(0,1,0,0,0,0,0),
       "zeroGradient"
    );

    volScalarField elsaSigmaRatio
    (
        IOobject
        (
            "elsaSigmaRatio",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh,
       dimensionSet(0,0,0,0,0,0,0),
       "zeroGradient"
    );

    IOdictionary elsaProperties
    (
        IOobject
        (
            "elsaProperties",
            runTime.constant(),
            mesh,
            IOobject::MUST_READ_IF_MODIFIED,
            IOobject::NO_WRITE
        )
    );

    volScalarField cellVol
    (
        IOobject
        (
            "cellVol",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
	mesh,
       dimensionSet(0,3,0,0,0,0,0),
       "zeroGradient"
    );

    dimensionedScalar cstElsaSigmaMin(elsaProperties.lookup("cstElsaSigmaMin"));

    dimensionedScalar elsaSurfaceTensionCoef(elsaProperties.lookup("elsaSurfaceTensionCoef"));
    dimensionedScalar webEquilMod(elsaProperties.lookup("webEquilMod"));

    dimensionedScalar elsaWeberEquil(elsaProperties.lookup("elsaWeberEquil"));

    dimensionedScalar cTurbElsaSigmaPrime(elsaProperties.lookup("cTurbElsaSigmaPrime"));

    dimensionedScalar elsaTiny("elsaTiny",dimensionSet(0,-1,0,0,0,0,0),1.e-1);
    dimensionedScalar elsaTinyNut("elsaTinyNut",dimensionSet(0,2,-1,0,0,0,0),1.e-8);
    dimensionedScalar kTiny("kTiny",dimensionSet(0,2,-2,0,0,0,0),1.e-1);
    dimensionedScalar sigmaTiny("sigmaTiny",dimensionSet(0,-1,0,0,0,0,0),1.e-1);
