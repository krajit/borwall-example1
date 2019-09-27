
label patchi = mesh.boundaryMesh().findPatchID("outlet");
J = 0.5*gSum(mesh.magSf().boundaryField()[patchi] * Foam::pow(mag(U.boundaryField()[patchi] - Ud.boundaryField()[patchi]) , 2));
