Info << "Reading numericalParameters\n" << endl;
IOdictionary numericalParameters
(
    IOobject
    (
	"numericalParameters",
	runTime.constant(),
	mesh,
	IOobject::MUST_READ,
	IOobject::NO_WRITE
    )
);

scalar tol(numericalParameters.lookupOrDefault<scalar>("tol", 1e-8));
scalar gamma(numericalParameters.lookupOrDefault<scalar>("gamma", 1));
scalar c1(numericalParameters.lookupOrDefault<scalar>("c1", 1e-4));
scalar c2(numericalParameters.lookupOrDefault<scalar>("c2", 0.8));


Info<< "Reading field p\n" << endl;
volScalarField p
(
    IOobject
    (
        "p",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "createPhi.H"


label pRefCell = 0;
scalar pRefValue = 0.0;
setRefCell(p, simple.dict(), pRefCell, pRefValue);
mesh.setFluxRequired(p.name());


Info<< "Reading field pa\n" << endl;
volScalarField pa
(
    IOobject
    (
        "pa",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field Ua\n" << endl;
volVectorField Ua
(
    IOobject
    (
        "Ua",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

Info<< "Reading field Ud\n" << endl;
volVectorField Ud 
(
    IOobject
    (
        "Ud",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

#include "createPhia.H"


label paRefCell = 0;
scalar paRefValue = 0.0;
setRefCell
(
    pa,
    simple.dict(),
    paRefCell,
    paRefValue
);
mesh.setFluxRequired(pa.name());


singlePhaseTransportModel laminarTransport(U, phi);

autoPtr<incompressible::turbulenceModel> turbulence
(
    incompressible::turbulenceModel::New(U, phi, laminarTransport)
);


const labelList& inletCells = mesh.boundary()["inlet"].faceCells();
//const labelList& outletCells = mesh.boundary()["outlet"].faceCells();

volScalarField alpha
(
    IOobject
    (
        "alpha",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh //lambda*max(Ua & U, zeroSensitivity)
);

scalarField volField = mesh.V();
zeroCells(alpha, inletCells);
//zeroCells(alpha, outletCells);

#include "createFvOptions.H"

volScalarField alphaMax
(
    IOobject
    (
	"alphaMax",
	runTime.constant(),
	mesh,
	IOobject::MUST_READ,
	IOobject::NO_WRITE
    ),
    mesh
);


volScalarField alphaMin
(
    IOobject
    (
	"alphaMin",
	runTime.constant(),
	mesh,
	IOobject::MUST_READ,
	IOobject::NO_WRITE
    ),
    mesh
);



volScalarField alphak
(
    IOobject
    (
	"alphak",
	runTime.timeName(),
	mesh,
	IOobject::NO_READ,
	IOobject::NO_WRITE
    ),
    alpha
);


volScalarField sensitivity
(
    IOobject
    (
	"sensitivity",
	runTime.timeName(),
	mesh,
	IOobject::NO_READ,
	IOobject::AUTO_WRITE
    ),
    U & Ua
);



