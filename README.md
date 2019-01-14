# CDCS

`CDCS.jl` is an interface to the **[CDCS](https://github.com/oxfordcontrol/CDCS)**
solver. It exports the `cdcs` function that is a thin wrapper on top of the
`cdcs` MATLAB function and use it to define the `CDCS.Optimizer` object that
implements the solver-independent `MathOptInterface` API.

To use it with JuMP, simply do
```julia
using JuMP
using CDCS
model = Model(with_optimizer(CDCS.Optimizer))
```
To suppress output, do
```julia
model = Model(with_optimizer(CDCS.Optimizer, fid=0))
```

## Installation

You can install `CDCS.jl` through the Julia package manager:
```julia
] add https://github.com/blegat/CDCS.jl.git
```
but you first need to make sure that you satisfy the requirements of the
[MATLAB.jl](https://github.com/JuliaInterop/MATLAB.jl) Julia package and that
the CDCS software is installed in your
[MATLAB™](http://www.mathworks.com/products/matlab/) installation.
