# Finite Element Method with Perfectly Matched Layer

## Perfectly matched layer (PML)

Perfectly Matched Layer method was proposed in 1994. Then Chew and Weedon showed that 
the PML model can be obtained by using complex change of variables.

The complex coordinate approach is based on analytic continuation of wave equations (Maxwell's 
equation, Helmholtz equation) into complex spatial coordinates where the fields are exponentially 
decaying. The origin solution of wave equation depends on ``x \in \mathbb{R}``. Furthermore, 
they are analytic functions of ``x``. That means that we can analytically continue it, i.e., 
evaluating the solution at complex values.

## Finite element discretization