# Nmrc

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://peakfind.github.io/Nmrc.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://peakfind.github.io/Nmrc.jl/dev/)
[![Build Status](https://github.com/peakfind/Nmrc.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/peakfind/Nmrc.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/peakfind/Nmrc.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/peakfind/Nmrc.jl)

`Nmrc.jl` collects code snippets for numerical solution of the new radiation condition. All Finite element method codes are implemented by using [`Ferrite.jl`](https://github.com/Ferrite-FEM/Ferrite.jl).

## Structure
- Use PML and FEM to estimate the exceptional values (quasi momentum).
- Use FEM with transparent boundary condition (DtN map).

## Questions
1. High order quadrature rules on facets may change the sparse pattern of the vector $\Theta^n$.