// Solve the energy equation

fvScalarMatrix TEqn
(
    fvm::ddt(T)
    + fvm::div(phi, T)
    == fvm::laplacian(a, T)
);
TEqn.solve();
