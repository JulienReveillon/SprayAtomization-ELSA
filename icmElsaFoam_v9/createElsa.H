Info<< "Reading field elsaSigmaPrime\n" << endl;

// It must read sigmaPrime from time directory. Later a transport equation
// is solved for sigmaPrime
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

// sigmaMin is computed from alpha field and thus must not be read
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

// sigma is the sum of sigmaMin and sigmaPrime
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


// Is not re read the Sigma equil called also sigma* in literature
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

// Post-processing Sauter Mean Diameter Field
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


// Sigma/Sigma*
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

// IRQ alpha measure where the interface is
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
	dimensionedScalar("IRQalpha", dimless, 0),
	"zeroGradient"
);

// IRQsigma = sigmaMin/sigma
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

IRQsigma = 1;

// IRQ Sigma value interpolated to cell faces
surfaceScalarField IRQsigmaF
(
	"IRQsigmaF",
	(fvc::interpolate(IRQsigma))
);
IRQsigmaF = 1.0;

// IRQ alpha value interpolated to cell faces
surfaceScalarField IRQalphaF                                                    
(                                                                               
    "IRQalphaF",                                                                   
    (fvc::interpolate(IRQalpha))                                                   
);

// cAlpha value over the domain
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

// cAlpha values interpolated to cell faces
surfaceScalarField cAlphaFieldF
(
    "cAlphaFieldF",
    (fvc::interpolate(cAlphaField))
);

// IRQkappa inverse 
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

// IRQkappa inverse values interpolated to cell faces
surfaceScalarField kappaDeltaXF
(
	"kappaDeltaXF",
	(fvc::interpolate(kappaDeltaX))
);

// Cell volume field
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

//Dot product between gravity normal vector and surface normal vector 
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

// magnitude of gradiente of alpha G
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

// IRQg = magGradAlpha * deltaX
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

// IRQ g values interpolated to cell faces
surfaceScalarField IRQgF                                                    
 (                                                                              
   "IRQgF",                                                                 
   (fvc::interpolate(IRQg))                                                     
 );                                                                             

// turbulent kinetic energy 
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

// elsa properties dictionary in constant directory
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

// curvature calculated as the divergence of the surface normal vector
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
    dimensionedScalar("0", dimless/dimLength, 0)//,
);
 
// IRQ k = 1/kappa/deltaX
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
    dimensionedScalar("0", dimless, 0)//,
);

label activateDiffusion = 0;

dimensionedScalar Dab(mixture.lookup("Dab"));

dimensionedScalar alphatab(mixture.lookup("alphatab"));

dimensionedScalar cstElsaSigmaMin(elsaProperties.lookup("cstElsaSigmaMin"));

dimensionedScalar elsaSurfaceTensionCoef                                        
("elsaSurfaceTensionCoef",dimensionSet(1,0,-2,0,0,0,0),mixture.lookup("sigma"));

dimensionedScalar nDimension(mixture.lookup("nDimension"));

dimensionedScalar webEquilMod(elsaProperties.lookup("webEquilMod"));

dimensionedScalar elsaWeberEquil(elsaProperties.lookup("elsaWeberEquil"));

dimensionedScalar cTurbElsaSigmaPrime(elsaProperties.lookup("cTurbElsaSigmaPrime"));

dimensionedScalar elsaTinyNut("elsaTinyNut",dimensionSet(0,2,-1,0,0,0,0),1.e-8);

dimensionedScalar elsaTiny("elsaTiny",dimensionSet(0,-1,0,0,0,0,0),1.e-1);

dimensionedScalar sigmaTiny("sigmaTiny",dimensionSet(0,-1,0,0,0,0,0),1.e-1);

dimensionedScalar IRQkTiny("IRQkTiny",dimensionSet(0,0,0,0,0,0,0),SMALL);

dimensionedScalar tKETiny("tKETiny",dimensionSet(0,2,-2,0,0,0,0),1.e-8);
