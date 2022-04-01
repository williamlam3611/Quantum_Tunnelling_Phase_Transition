# Quantum-Well-Resonance-in-SrTiO3-Heterostructure

This is a rewrite of an [old project](https://github.com/jwmmorley/MPhys-Project-Fortran).

## Libhqt
Models electron tunneling in Heterostructure by applying quantum wells. Can determine the Number density and band structure per electron band and/or energy state of the heterostructure.

Use the classic install method:
```
make
make install
```

Can be built to enable [openmpi](https://www.open-mpi.org/), else it uses [openmp](https://www.openmp.org/) [TODO].

## Examples
 - Determine band structure and number density of the dxy, dxz and dyz bands of a Strontium Titanait(SrTiO3) heterostructure with the aim of producing quantum resonance[TODO]
 - Apply Poisson-schrodinger equation to determine an impirical formula for fields applied at the surface of SrTiO3[TODO]
