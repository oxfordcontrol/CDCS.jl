module CDCS

using SparseArrays
using MATLAB

export cdcs

# The fields should be integer values but the type should be `Float64` to work
# with CDCS
mutable struct Cone
    f::Float64 # number of free primal variables / linear dual equality constraints
    l::Float64 # length of LP cone
    q::Vector{Float64} # list of second-order cone dimension
    s::Vector{Float64} # list of semidefinite constraints
end
Cone(f::Real, l::Real) = Cone(f, l, Float64[], Float64[])

to_vec(x::Vector) = x
to_vec(x::Float64) = [x]
function to_vec(x::AbstractMatrix)
    @assert isempty(x) || isone(size(x, 2))
    return vec(x)
end

# Solve the primal/dual pair
# min c'x,      max b'y
# s.t. Ax = b,   c - A'x ∈ K
#       x ∈ K
function sedumi(A::Union{Matrix{Float64}, SparseMatrixCSC{Float64}},
                b::Vector{Float64}, c::Vector{Float64},
                K::Cone = Cone(0, size(A, 2));
                kws...)
    pars = Dict{String, Any}(string(key) => value for (key, value) in kws)
    # There are 3 output arguments: x, y and info so we use `3` above
    x, y, info = mxcall(:cdcs, 3, A, b, c, K, pars)
    return to_vec(x), to_vec(y), info
end

include("MOI_wrapper.jl")

end # module
