```@meta
CurrentModule = Nmrc
```

# Nmrc

Documentation for [Nmrc](https://github.com/peakfind/Nmrc.jl).

## Usage

Nmrc.jl can be installed by opening a Julia REPL and typing:
```julia
]add https://github.com/peakfind/Nmrc.jl.git
```

## Methods Implemented in Nmrc.jl

- Use PML and FEM to estimate the exceptional values (quasi momentum).
- FEM with transparent boundary condition (DtN map).

## TBD

1. High order quadrature rules on facets may change the sparse pattern of the vector $\Theta^n$.

```@index
```

```@autodocs
Modules = [Nmrc]
```
