# CDCS

`CDCS.jl` is an interface to the **[CDCS](https://github.com/oxfordcontrol/CDCS)**
solver. It exports the `cdcs` function that is a thin wrapper on top of the
`cdcs` MATLAB function and use it to define the `CDCS.Optimizer` object that
implements the solver-independent
[MathOptInterface](https://github.com/JuliaOpt/MathOptInterface.jl) API.

To use it with [JuMP](https://github.com/JuliaOpt/JuMP.jl), simply do
```julia
using JuMP
using CDCS
model = Model(optimizer_with_attributes(CDCS.Optimizer))
```
To suppress output, do either
```julia
set_silent(model)
```
or
```julia
model = Model(optimizer_with_attributes(CDCS.Optimizer, verbose=0))
```

## Installation

You can install `CDCS.jl` through the
[Julia package manager](https://docs.julialang.org/en/v1/stdlib/Pkg/index.html):
```julia
] add CDCS
```
but you first need to make sure that you satisfy the requirements of the
[MATLAB.jl](https://github.com/JuliaInterop/MATLAB.jl) Julia package and that
the CDCS software is installed in your
[MATLAB™](http://www.mathworks.com/products/matlab/) installation.

### Troubleshooting

#### CDCS not in PATH

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
The error means that we try to find the `cdcs` function with 1 output argument using the MATLAB C API but it wasn't found.
This most likely means that you did not add CDCS to the MATLAB's path (i.e. the `toolbox/local/pathdef.m` file).

If modifying `toolbox/local/pathdef.m` does not work, the following should work where `/path/to/CDCS/` is the directory where the `CDCS` folder is located:
```julia
julia> using MATLAB

julia> cd("/path/to/CDCS/") do
           mat"cdcsInstall"
       end
```
This should make `CDCS.jl` work for the Julia session in which this is run.
Alternatively, run
```julia
julia> mat"savepath"
```
to make `CDCS.jl` work for future Julia sessions.
