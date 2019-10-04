
//gradAlpha = fvc::grad(alpha);

J =0 ;//  costLambda.value()*gSum(volField * (pow(mag(alpha),2)) );
// integrate in the boundary patches

scalar mass = gSum(volField*rho);
scalar maxVolFraction = 0.5;

J = pow(mass/maxVolFraction,10); // small value if the constraint is met and large value otherwise


forAll(mesh.boundary(), patchi)
{
    // only over non-empty patches
    if (!p.boundaryField()[patchi].empty())
    {
        J -= gSum(
            (p.boundaryField()[patchi] + 0.5*pow(mag(U.boundaryField()[patchi]),2) )
            *( U.boundaryField()[patchi] & mesh.Sf().boundaryField()[patchi])
        );
    }
}