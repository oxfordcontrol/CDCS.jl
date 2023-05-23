# CDCS.jl

[CDCS.jl](https://github.com/oxfordcontrol/CDCS.jl) is a wrapper for the
[CDCS](https://github.com/oxfordcontrol/CDCS) solver. 

The wrapper has two components:

 * an exported `cdcs` function that is a thin wrapper on top of the `cdcs`
   MATLAB function
 * an interface to [MathOptInterface](https://github.com/jump-dev/MathOptInterface.jl)

## License

CDCS.jl is licensed under the [MIT license](https://github.com/oxfordcontrol/CDCS.jl/blob/master/LICENSE.md).

The underlying solver [oxfordcontrol/CDCS](https://github.com/oxfordcontrol/CDCS)
is licensed under the [LGPL v3 license](https://github.com/oxfordcontrol/CDCS/blob/master/LICENCE.txt).

In addition, CDCS requires an installation of MATLAB, which is a closed-source
commercial product for which you must [obtain a license](https://www.mathworks.com/products/matlab.html).

## Installation

First, make sure that you satisfy the requirements of the
[MATLAB.jl](https://github.com/JuliaInterop/MATLAB.jl) Julia package and that
the CDCS software is installed in your [MATLAB™](http://www.mathworks.com/products/matlab/)
installation.

Then, install `CDCS.jl` using `Pkg.add`:
```julia
import Pkg
Pkg.add("CDCS")
```

## Use with JuMP

To use CDCS with [JuMP](https://github.com/jump-dev/JuMP.jl), do:
```julia
using JuMP, CDCS
model = Model(CDCS.Optimizer)
set_attribute(model, "verbose", 0)
```

## Troubleshooting

If you get the error:
```
Undefined function or variable 'cdcs'.

Error using save
Variable 'jx_cdcs_arg_out_1' not found.

Linear Programming example: Error During Test at /home/blegat/.julia/dev/CDCS/test/lp.jl:5
  Got exception outside of a @test
  MATLAB.MEngineError("failed to get variable jx_cdcs_arg_out_1 from MATLAB session")
  Stacktrace:
    [1] get_mvariable(session::MATLAB.MSession, name::Symbol)
      @ MATLAB ~/.julia/packages/MATLAB/SVjnA/src/engine.jl:164
    [2] mxcall(::MATLAB.MSession, ::Symbol, ::Int64, ::Matrix{Float64}, ::Vararg{Any})
      @ MATLAB ~/.julia/packages/MATLAB/SVjnA/src/engine.jl:297
    [3] mxcall
      @ ~/.julia/packages/MATLAB/SVjnA/src/engine.jl:317 [inlined]
    [4] cdcs(A::Matrix{Float64}, b::Vector{Float64}, c::Vector{Float64}, K::CDCS.Cone; kws::Base.Pairs{Symbol, Int64, Tuple{Symbol}, NamedTuple{(:verbose,), Tuple{Int64}}})
```
The error means that we could not find the `cdcs` function with one output
argument using the MATLAB C API.

This most likely means that you did not add CDCS to the MATLAB's path (that is,
the `toolbox/local/pathdef.m` file).

If modifying `toolbox/local/pathdef.m` does not work, the following should work
where `/path/to/CDCS/` is the directory where the `CDCS` folder is located:
```julia
julia> using MATLAB

julia> cd("/path/to/CDCS/") do
           mat"cdcsInstall"
       end

julia> mat"savepath"
```
