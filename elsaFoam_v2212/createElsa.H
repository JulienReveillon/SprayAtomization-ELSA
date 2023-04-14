Info<< "Reading field sigmaPrime\n" << endl;
volScalarField sigmaPrime
(
	IOobject
	(
		"sigmaPrime",
		runTime.timeName(),
  		mesh,
  		IOobject::MUST_READ,
 		IOobject::AUTO_WRITE
 	),
 	mesh
);

volScalarField sigma
(
	IOobject
	(
		"sigma",
		runTime.timeName(),
  		mesh,
  		IOobject::NO_READ,
 		IOobject::AUTO_WRITE
 	),
 	sigmaPrime
);

volScalarField sigmaMin
(
	IOobject
	(
  		"sigmaMin",
  		runTime.timeName(),
  		mesh,
  		IOobject::NO_READ,
  		IOobject::AUTO_WRITE
	),
	sigmaPrime
);

volScalarField sigmaEquil
(
	IOobject
	(
		"sigmaEquil",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	),
	sigmaPrime
 );

volScalarField D32
(
    IOobject
    (
        "D32",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionSet(0,1,0,0,0,0,0),
    "zeroGradient"
);

volScalarField sigmaRatio
(
	IOobject
	(
		"sigmaRatio",
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


dimensionedScalar cstElsaSigmaMin
(
    "cstElsaSigmaMin",
    dimensionSet(0,0,0,0,0,0,0),
    elsaProperties.lookup("cstElsaSigmaMin")
);

dimensionedScalar elsaSurfaceTensionCoef
(
    "elsaSurfaceTensionCoef",
    dimensionSet(1,0,-2,0,0,0,0),
    elsaProperties.lookup("elsaSurfaceTensionCoef")
);

dimensionedScalar webEquilMod
(
    "webEquilMod",
    dimensionSet(0,0,0,0,0,0,0),
    elsaProperties.lookup("webEquilMod")
);

dimensionedScalar elsaWeberEquil
(
    "elsaWeberEquil",
    dimensionSet(0,0,0,0,0,0,0),
    elsaProperties.lookup("elsaWeberEquil") 
);

dimensionedScalar cTurbElsaSigmaPrime
(
    "cTurbElsaSigmaPrime",
    dimensionSet(0,0,0,0,0,0,0),
    elsaProperties.lookup("cTurbElsaSigmaPrime") );

dimensionedScalar elsaTinyNut
(
    "elsaTinyNut",dimensionSet(0,2,-1,0,0,0,0),1.e-8
);

dimensionedScalar elsaTiny
(
    "elsaTiny",dimensionSet(0,-1,0,0,0,0,0),1.e-1
);

dimensionedScalar sigmaTiny
(
    "sigmaTiny",dimensionSet(0,-1,0,0,0,0,0),1.e-1
);


dimensionedScalar kTiny
(
    "kTiny",dimensionSet(0,2,-2,0,0,0,0),1.e-1
);
