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
  		"sigmaMinDilute",
  		runTime.timeName(),
  		mesh,
  		IOobject::NO_READ,
  		IOobject::AUTO_WRITE
	),
	sigmaPrime
);

volScalarField sigmaMinDense
(
	IOobject
	(
  		"sigmaMinDense",
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


volScalarField IRQalpha
(
	IOobject
	(
		"IRQalpha",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	),
	mesh,
	dimensionedScalar("IRQalpha", dimless, 0), //is set to 0 at time 0
	"zeroGradient"
);

volScalarField IRQsigma
(
	IOobject
	(
		"IRQsigma",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	),
	mesh,
	dimensionedScalar("IRQsigma",dimless,0),
	"zeroGradient"
);

surfaceScalarField IRQsigmaF
(
    "IRQsigmaF",
    (
        fvc::interpolate(IRQsigma)
    )
); 

IRQsigmaF = 0.0;
IRQsigma = mixture.cAlpha();

volScalarField cAlphaField                                                      
(                                                                               
    IOobject                                                                    
    (                                                                           
        "cAlpha",                                                               
        runTime.timeName(),                                                     
        mesh,                                                                   
        IOobject::NO_READ,                                                      
        IOobject::AUTO_WRITE                                                    
    ),                                                                          
    mesh,                                                                       
    dimensionedScalar("cAlphaField",dimless,0),                                 
    "zeroGradient"                                                              
); 

surfaceScalarField cAlphaFieldF
(
	"cAlphaFieldF",
	(
        fvc::interpolate(cAlphaField)
    )
);


surfaceScalarField IRQalphaF
(
    "IRQalphaF",
    (
        fvc::interpolate(IRQalpha)
    )
);

volScalarField curvature
(
	IOobject
	(
		"curvature",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	),
	mesh,
	dimensionedScalar("curvature",dimless/dimLength,0),
	"zeroGradient"
);

surfaceScalarField curvatureF
(
	"curvatureF",
	(
        fvc::interpolate(curvature)
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
        dimensionedScalar("zero",dimVolume,0.0)
);
cellVol.primitiveFieldRef() = mesh.V();

volScalarField gNormalInt
(
	IOobject
	(
		"gNormalInt",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	),
	mesh,
	dimensionSet(0,0,0,0,0,0,0),
	"zeroGradient"
);

volScalarField magGradAlpha
(
	IOobject
	(
		"magGradAlpha",
		runTime.timeName(),
		mesh,
		IOobject::NO_READ,
		IOobject::AUTO_WRITE
	),
	mesh,
	dimensionSet(2,-7,0,0,0,0,0),
	"zeroGradient"
);

volScalarField kappaDeltaX                                                      
(                                                                               
    IOobject                                                                    
    (                                                                           
        "kappaDeltaX",                                                          
        runTime.timeName(),                                                     
        mesh,                                                                   
        IOobject::NO_READ,                                                      
        IOobject::NO_WRITE                                                      
    ),                                                                          
    mesh,                                                                       
    dimensionedScalar("kappaDeltaX",dimless,0),                                 
    "zeroGradient"                                                              
);  

surfaceScalarField kappaDeltaXF                                                 
(                                                                               
    "kappaDeltaXF",                                                             
    (
        fvc::interpolate(kappaDeltaX)
    )                                             
);

volScalarField IRQg
(
    IOobject
    (
         "IRQg",
         runTime.timeName(),
         mesh,
         IOobject::NO_READ,
         IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionSet(0,0,0,0,0,0,0),
    "zeroGradient"
);

surfaceScalarField IRQgF
(
    "IRQgField",
    (
        fvc::interpolate(IRQg)
    )
);

volScalarField tKE
(
    IOobject
    (
         "tKE",
         runTime.timeName(),
         mesh,
         IOobject::NO_READ,
         IOobject::AUTO_WRITE
    ),
    turbulence->k()
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

volScalarField IRQk
(
    IOobject
    (
         "IRQk",
         runTime.timeName(),
         mesh,
         IOobject::NO_READ,
         IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("0", dimless, 0)
);

volVectorField  UaveLiq                                                      
(                                                                          
    IOobject                                                               
    (                                                                      
        "UaveLiq",                                                         
        runTime.timeName(),                                                
        mesh,                                                              
        IOobject::NO_READ,                                                 
        IOobject::AUTO_WRITE                                               
    ),                                                                     
    mesh,                                                                  
    dimensionedVector("UaveLiq",dimensionSet(0,1,-1,0,0),pTraits<vector>::zero)
);  

volScalarField coeffAlpha                                                      
(                                                                              
    IOobject                                                               
    (                                                                      
        "cLag",                                                               
         runTime.timeName(),                                                   
         mesh,                                                                 
         IOobject::NO_READ,                                                    
         IOobject::AUTO_WRITE                                                  
    ),                                                                     
    mesh,                                                                  
    dimensionedScalar("coeffAlpha",dimless,0),                              
    "zeroGradient"                                                          
);  


surfaceScalarField coeffAlphaF
(
    "coeffAlphaF",
    (fvc::interpolate(coeffAlpha)) 
);

surfaceScalarField phiSigmaLag
 (
  "phiSigmaLag",
  (phi*coeffAlphaF)
  );

surfaceScalarField phiAlphaLag
(
    "phiAlphaLag",
    (phi*coeffAlphaF)
);

volScalarField Nbgoutte
(
    IOobject
    (
        "Nbgoutte",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    mesh,
    dimensionedScalar("Nbgoutte",dimless,0),
    "zeroGradient"
);



volScalarField tagDiffusion
 (
  IOobject
  (
   "tagDiffusion",
   runTime.timeName(),
   mesh,
   IOobject::NO_READ,
   IOobject::AUTO_WRITE
   ),
   mesh,
  dimensionedScalar("tagDiffusion",dimless,0),
  "zeroGradient"
 );

label activateDiffusion = 0;

dimensionSet ds(0,0,0,0,0,0,0); //dimensionset for Lagrangian                   
                                                                                
volScalarField sigmaKineticEnergy                                               
    (                                                                           
        IOobject                                                                
        (                                                                       
            "sigmaKineticEnergy",                                               
            runTime.timeName(),                                                 
            mesh,                                                               
            IOobject::NO_READ,                                                  
            IOobject::NO_WRITE                                                  
        ),                                                                      
        mesh,                                                                   
       dimensionSet(0,1,-1,0,0,0,0),                                            
       "zeroGradient"                                                           
    );

dimensionedScalar Dab
(
    "Dab", 
    dimensionSet(0,2,-1,0,0,0,0),  
    mixture
);

dimensionedScalar alphatab
(
    "alphatab",
    dimless,
    mixture
);

dimensionedScalar cstElsaSigmaMin
(
    "cstElsaSigmaMin",
    dimensionSet(0,0,0,0,0,0,0),
    elsaProperties
);

dimensionedScalar elsaSurfaceTensionCoef
(
    "elsaSurfaceTensionCoef",
    dimensionSet(1,0,-2,0,0,0,0),
    elsaProperties
);

dimensionedScalar webEquilMod
(
    "webEquilMod",
    dimensionSet(0,0,0,0,0,0,0),
    elsaProperties
);

dimensionedScalar elsaWeberEquil
(
    "elsaWeberEquil",
    dimensionSet(0,0,0,0,0,0,0),
    elsaProperties 
);

dimensionedScalar cTurbElsaSigmaPrime
(
    "cTurbElsaSigmaPrime",
    dimensionSet(0,0,0,0,0,0,0),
    elsaProperties
);

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

dimensionedScalar IRQkTiny
(
    "IRQkTiny",dimensionSet(0,0,0,0,0,0,0),1.e-8
);

dimensionedScalar tKETiny
(
    "tKETiny",dimensionSet(0,2,-2,0,0,0,0),1.e-8
);
