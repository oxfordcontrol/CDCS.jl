using LinearAlgebra

using MathOptInterface
const MOI = MathOptInterface
const CI = MOI.ConstraintIndex
const VI = MOI.VariableIndex

const MOIU = MOI.Utilities


# CDCS solves the primal/dual pair
# min c'x,       max b'y
# s.t. Ax = b,   c - A'x ∈ K
#       x ∈ K
# where K is a product of Zeros, Nonnegatives, SecondOrderCone,
# RotatedSecondOrderCone and PositiveSemidefiniteConeTriangle

# This wrapper copies the MOI problem to the CDCS dual so the natively
# supported supported sets are `VectorAffineFunction`-in-`S` where `S` is one
# of the sets just listed above.

const SF = Union{MOI.VectorOfVariables, MOI.VectorAffineFunction{Float64}}
const SS = Union{MOI.Zeros, MOI.Nonnegatives, MOI.Nonpositives,
                 MOI.SecondOrderCone, MOI.PositiveSemidefiniteConeTriangle}

mutable struct Solution
    x::Vector{Float64}
    y::Vector{Float64}
    slack::Vector{Float64}
    objval::Float64
    info::Dict{String, Any}
end

# Used to build the data with allocate-load during `copy_to`.
# When `optimize!` is called, a the data is passed to CDCS
# using `sedumi` and the `ModelData` struct is discarded
mutable struct ModelData
    m::Int # Number of rows/constraints of CDCS dual/MOI primal
    n::Int # Number of cols/variables of CDCS primal/MOI dual
    I::Vector{Int} # List of rows of A'
    J::Vector{Int} # List of cols of A'
    V::Vector{Float64} # List of coefficients of A
    c::Vector{Float64} # objective of CDCS primal/MOI dual
    objconstant::Float64 # The objective is min c'x + objconstant
    b::Vector{Float64} # objective of CDCS dual/MOI primal
end

# This is tied to CDCS's internal representation
mutable struct ConeData
    K::Cone
    sum_q::Int # cached value of sum(q)
    sum_s2::Int # cached value of sum(s.^2)
    nrows::Dict{Int, Int} # The number of rows of each vector sets, this is used by `constrrows` to recover the number of rows used by a constraint when getting `ConstraintPrimal` or `ConstraintDual`
    function ConeData()
        new(Cone(0, 0, Float64[], Float64[], Float64[]),
            0, 0, 0, Dict{Int, Int}())
    end
end

mutable struct Optimizer <: MOI.AbstractOptimizer
    cone::ConeData
    maxsense::Bool
    data::Union{Nothing, ModelData} # only non-Nothing between MOI.copy_to and MOI.optimize!
    sol::Union{Nothing, Solution}
    options::Iterators.Pairs
    function Optimizer(; options...)
        new(ConeData(), false, nothing, nothing, options)
    end

end

MOI.get(::Optimizer, ::MOI.SolverName) = "CDCS"

function MOI.is_empty(optimizer::Optimizer)
    !optimizer.maxsense && optimizer.data === nothing
end
function MOI.empty!(optimizer::Optimizer)
    optimizer.maxsense = false
    optimizer.data = nothing # It should already be nothing except if an error is thrown inside copy_to
    optimizer.sol = nothing
end

MOIU.supports_allocate_load(::Optimizer, copy_names::Bool) = !copy_names

function MOI.supports(::Optimizer,
                      ::Union{MOI.ObjectiveSense,
                              MOI.ObjectiveFunction{MOI.SingleVariable},
                              MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}})
    return true
end

MOI.supports_constraint(::Optimizer, ::Type{<:SF}, ::Type{<:SS}) = true

function MOI.copy_to(dest::Optimizer, src::MOI.ModelLike; kws...)
    return MOIU.automatic_copy_to(dest, src; kws...)
end

const LPCones = Union{MOI.Nonnegatives, MOI.Nonpositives}

# Computes cone dimensions
function constroffset(cone::ConeData,
                      ci::CI{<:MOI.AbstractFunction, MOI.Zeros})
    return ci.value
end
#_allocate_constraint: Allocate indices for the constraint `f`-in-`s`
# using information in `cone` and then update `cone`
function _allocate_constraint(cone::ConeData, f, s::MOI.Zeros)
    ci = Int(cone.K.f)
    cone.K.f += MOI.dimension(s)
    return ci
end
function constroffset(cone::ConeData,
                      ci::CI{<:MOI.AbstractFunction, <:LPCones})
    return Int(cone.K.f) + ci.value
end
function _allocate_constraint(cone::ConeData, f, s::LPCones)
    ci = cone.K.l
    cone.K.l += MOI.dimension(s)
    return ci
end
function constroffset(cone::ConeData,
                      ci::CI{<:MOI.AbstractVectorFunction,
                             <:MOI.SecondOrderCone})
    return Int(cone.K.f) + Int(cone.K.l) + ci.value
end
function _allocate_constraint(cone::ConeData, f, s::MOI.SecondOrderCone)
    ci = cone.sum_q
    push!(cone.K.q, s.dimension)
    cone.sum_q += s.dimension
    return ci
end
function constroffset(cone::ConeData,
                      ci::CI{<:MOI.AbstractFunction,
                             <:MOI.PositiveSemidefiniteConeTriangle})
    return Int(cone.K.f) + Int(cone.K.l) + cone.sum_q + ci.value
end
function _allocate_constraint(cone::ConeData, f,
                              s::MOI.PositiveSemidefiniteConeTriangle)
    ci = cone.sum_s2
    push!(cone.K.s, s.side_dimension)
    cone.sum_s2 += s.side_dimension^2
    return ci
end
function constroffset(optimizer::Optimizer, ci::CI)
    return constroffset(optimizer.cone, ci::CI)
end
function MOIU.allocate_constraint(optimizer::Optimizer, f::F, s::S) where {F <: MOI.AbstractFunction, S <: MOI.AbstractSet}
    return CI{F, S}(_allocate_constraint(optimizer.cone, f, s))
end

# Vectorized length for matrix dimension n
sympackedlen(n) = div(n*(n+1), 2)
# Matrix dimension for vectorized length n
sympackeddim(n) = div(isqrt(1+8n) - 1, 2)
sqrdim(n) = isqrt(n)
trimap(i::Integer, j::Integer) = i < j ? trimap(j, i) : div((i-1)*i, 2) + j
sqrmap(i::Integer, j::Integer, n::Integer) = i < j ? sqrmap(j, i, n) : i + (j-1) * n
function _copyU(x, n, mapfrom, mapto)
    y = zeros(eltype(x), mapto(n, n))
    for i in 1:n, j in 1:i
        y[mapto(i, j)] = x[mapfrom(i, j)]
    end
    return y
end
squareUtosympackedU(x, n=sqrdim(length(x))) = _copyU(x, n, (i, j) -> sqrmap(i, j, n), trimap)
sympackedUtosquareU(x, n=sympackeddim(length(x))) = _copyU(x, n, trimap, (i, j) -> sqrmap(i, j, n))

function sympackedUtosquareUidx(x::AbstractVector{<:Integer}, n)
    y = similar(x)
    map = squareUtosympackedU(1:n^2, n)
    for i in eachindex(y)
        y[i] = map[x[i]]
    end
    return y
end

# Scale coefficients depending on rows index on symmetric packed upper triangular form
# coef: List of coefficients
# minus: if true, multiply the result by -1
# rev: if true, we unscale instead (e.g. divide by √2 instead of multiply for PSD cone)
# rows: List of row indices
# d: dimension of set
function _scalecoef(coef, minus::Bool,
                    ::Type{<:MOI.AbstractSet},
                    rev, args...)
    return minus ? -coef : coef
end
function _scalecoef(coef, minus::Bool,
                    ::Type{MOI.Nonpositives},
                    rev, args...)
    return minus ? coef : -coef
end
function _scalecoef(coef::AbstractVector, minus::Bool,
                    ::Type{MOI.PositiveSemidefiniteConeTriangle},
                    rev::Bool,
                    rows::AbstractVector, d::Integer)
    scaling = minus ? -1 : 1
    scaling2 = rev ? scaling / 2 : scaling * 2
    output = copy(coef)
    diagidx = BitSet()
    for i in 1:d
        push!(diagidx, trimap(i, i))
    end
    for i in 1:length(output)
        if rows[i] in diagidx
            output[i] *= scaling
        else
            output[i] *= scaling2
        end
    end
    return output
end
# Unscale the coefficients in `coef` with respective rows in `rows` for a set `s` and multiply by `-1` if `minus` is `true`.
scalecoef(coef, minus, s, args...) = _scalecoef(coef, minus, typeof(s), false, args...)
function scalecoef(coef, minus, s::MOI.PositiveSemidefiniteConeTriangle,
                   rows)
    return _scalecoef(coef, minus, typeof(s), false, rows, s.side_dimension)
end
# Unscale the coefficients in `coef` for a set of type `S`
unscalecoef(coef, S::Type{<:MOI.AbstractSet}) = _scalecoef(coef, false, S, true)
# Unscale the coefficients of `coef` in symmetric packed upper triangular form
function unscalecoef(coef, S::Type{MOI.PositiveSemidefiniteConeTriangle})
    len = length(coef)
    return _scalecoef(coef, false, S, true, 1:len, sympackeddim(len))
end

output_index(t::MOI.VectorAffineTerm) = t.output_index
variable_index_value(t::MOI.ScalarAffineTerm) = t.variable_index.value
variable_index_value(t::MOI.VectorAffineTerm) = variable_index_value(t.scalar_term)
coefficient(t::MOI.ScalarAffineTerm) = t.coefficient
coefficient(t::MOI.VectorAffineTerm) = coefficient(t.scalar_term)
# constrrows: Recover the number of rows used by each constraint.
# When, the set is available, simply use MOI.dimension
constrrows(s::MOI.AbstractVectorSet) = 1:MOI.dimension(s)
constrrows(s::MOI.PositiveSemidefiniteConeTriangle) = 1:(s.side_dimension^2)
# When only the index is available, use the `optimizer.ncone.nrows` field
constrrows(optimizer::Optimizer, ci::CI{<:MOI.AbstractVectorFunction, <:MOI.AbstractVectorSet}) = 1:optimizer.cone.nrows[constroffset(optimizer, ci)]
MOIU.load_constraint(optimizer::Optimizer, ci, f::MOI.VectorOfVariables, s) = MOIU.load_constraint(optimizer, ci, MOI.VectorAffineFunction{Float64}(f), s)
orderval(val, s) = val
function orderval(val, s::MOI.PositiveSemidefiniteConeTriangle)
    return sympackedUtosquareU(val, s.side_dimension)
end
orderidx(idx, s) = idx
function orderidx(idx, s::MOI.PositiveSemidefiniteConeTriangle)
    sympackedUtosquareUidx(idx, s.side_dimension)
end
function MOIU.load_constraint(optimizer::Optimizer, ci, f::MOI.VectorAffineFunction, s::MOI.AbstractVectorSet)
    A = sparse(output_index.(f.terms), variable_index_value.(f.terms), coefficient.(f.terms))
    # sparse combines duplicates with + but does not remove zeros created so we call dropzeros!
    dropzeros!(A)
    I, J, V = findnz(A)
    offset = constroffset(optimizer, ci)
    rows = constrrows(s)
    optimizer.cone.nrows[offset] = length(rows)
    i = offset .+ rows
    # The SCS format is b - Ax ∈ cone
    # so minus=false for b and minus=true for A
    optimizer.data.c[i] = orderval(scalecoef(f.constants, false, s,
                                             1:MOI.dimension(s)),
                                   s)
    append!(optimizer.data.I, offset .+ orderidx(I, s))
    append!(optimizer.data.J, J)
    append!(optimizer.data.V, scalecoef(V, true, s, I))
end

function MOIU.allocate_variables(optimizer::Optimizer, nvars::Integer)
    optimizer.cone = ConeData()
    optimizer.sol = nothing
    return VI.(1:nvars)
end

function MOIU.load_variables(optimizer::Optimizer, nvars::Integer)
    cone = optimizer.cone
    m = Int(cone.K.f) + Int(cone.K.l) + cone.sum_q + cone.sum_s2
    I = Int[]
    J = Int[]
    V = Float64[]
    c = zeros(m)
    b = zeros(nvars)
    optimizer.data = ModelData(m, nvars, I, J, V, c, 0., b)
end

function MOIU.allocate(optimizer::Optimizer, ::MOI.ObjectiveSense, sense::MOI.OptimizationSense)
    optimizer.maxsense = sense == MOI.MAX_SENSE
end
function MOIU.allocate(::Optimizer, ::MOI.ObjectiveFunction,
                       ::MOI.Union{MOI.SingleVariable,
                                   MOI.ScalarAffineFunction{Float64}})
end

function MOIU.load(::Optimizer, ::MOI.ObjectiveSense, ::MOI.OptimizationSense)
end
function MOIU.load(optimizer::Optimizer, ::MOI.ObjectiveFunction,
                   f::MOI.SingleVariable)
    MOIU.load(optimizer,
              MOI.ObjectiveFunction{MOI.ScalarAffineFunction{Float64}}(),
              MOI.ScalarAffineFunction{Float64}(f))
end
function MOIU.load(optimizer::Optimizer, ::MOI.ObjectiveFunction,
                   f::MOI.ScalarAffineFunction)
    c0 = Vector(sparsevec(variable_index_value.(f.terms), coefficient.(f.terms),
                          optimizer.data.n))
    optimizer.data.objconstant = f.constant
    optimizer.data.b = optimizer.maxsense ? c0 : -c0
    return nothing
end

function MOI.optimize!(optimizer::Optimizer)
    cone = optimizer.cone
    m = optimizer.data.m
    n = optimizer.data.n
    At = sparse(optimizer.data.I, optimizer.data.J, optimizer.data.V)
    c = optimizer.data.c
    objconstant = optimizer.data.objconstant
    b = optimizer.data.b
    optimizer.data = nothing # Allows GC to free optimizer.data before At is loaded to CDCS

    x, y, info = sedumi(At, b, c, optimizer.cone.K; optimizer.options...)

    objval = (optimizer.maxsense ? 1 : -1) * dot(b, y) + objconstant
    optimizer.sol = Solution(x, y, c - At' * y, objval, info)
end

function MOI.get(optimizer::Optimizer, ::MOI.SolveTime)
    return optimizer.sol.info["cpusec"]
end

# Implements getter for result value and statuses
# SeDuMI returns one of the following values (based on SeDuMi_Guide_11 by Pólik):
# feasratio:  1.0 problem with complementary solution
#            -1.0 strongly infeasible problem
#             between -1.0 and 1.0 nasty problem
# pinf = 1.0 : y is infeasibility certificate => CDCS primal/MOI dual is infeasible
# dinf = 1.0 : x is infeasibility certificate => CDCS dual/MOI primal is infeasible
# pinf = 0.0 = dinf : x and y are near feasible
# numerr: 0 desired accuracy (specified by pars.eps) is achieved
#         1 reduced accuracy (specified by pars.bigeps) is achieved
#         2 failure due to numerical problems

function MOI.get(optimizer::Optimizer, ::MOI.TerminationStatus)
    if optimizer.sol isa Nothing
        return MOI.OPTIMIZE_NOT_CALLED
    end
    pinf      = optimizer.sol.info["pinf"]
    dinf      = optimizer.sol.info["dinf"]
    numerr    = optimizer.sol.info["numerr"]
    if numerr == 2
        return MOI.NUMERICAL_ERROR
    end
    @assert iszero(numerr) || isone(numerr)
    accurate = iszero(numerr)
    if isone(pinf)
        if accurate
            return MOI.DUAL_INFEASIBLE
        else
            return MOI.ALMOST_DUAL_INFEASIBLE
        end
    end
    if isone(dinf)
        if accurate
            return MOI.INFEASIBLE
        else
            return MOI.ALMOST_INFEASIBLE
        end
    end
    @assert iszero(pinf) && iszero(dinf)
    # TODO when do we return SLOW_PROGRESS ?
    #      Maybe we should use feasratio
    if accurate
        return MOI.OPTIMAL
    else
        return MOI.ALMOST_OPTIMAL
    end
end

MOI.get(optimizer::Optimizer, ::MOI.ObjectiveValue) = optimizer.sol.objval

function MOI.get(optimizer::Optimizer,
                 attr::Union{MOI.PrimalStatus, MOI.DualStatus})
    if optimizer.sol isa Nothing
        return MOI.NO_SOLUTION
    end
    pinf      = optimizer.sol.info["pinf"]
    dinf      = optimizer.sol.info["dinf"]
    numerr    = optimizer.sol.info["numerr"]
    if numerr == 2
        return MOI.UNKNOWN_RESULT_STATUS
    end
    @assert iszero(numerr) || isone(numerr)
    accurate = iszero(numerr)
    if isone(attr isa MOI.PrimalStatus ? pinf : dinf)
        if accurate
            return MOI.INFEASIBILITY_CERTIFICATE
        else
            return MOI.NEARLY_INFEASIBILITY_CERTIFICATE
        end
    end
    if isone(attr isa MOI.PrimalStatus ? dinf : pinf)
        return MOI.INFEASIBLE_POINT
    end
    @assert iszero(pinf) && iszero(dinf)
    if accurate
        return MOI.FEASIBLE_POINT
    else
        return MOI.NEARLY_FEASIBLE_POINT
    end
end
function MOI.get(optimizer::Optimizer, ::MOI.VariablePrimal, vi::VI)
    optimizer.sol.y[vi.value]
end
MOI.get(optimizer::Optimizer, a::MOI.VariablePrimal, vi::Vector{VI}) = MOI.get.(optimizer, a, vi)
reorderval(val, s) = val
function reorderval(val, ::Type{MOI.PositiveSemidefiniteConeTriangle})
    return squareUtosympackedU(val)
end
function MOI.get(optimizer::Optimizer, ::MOI.ConstraintPrimal,
                 ci::CI{<:MOI.AbstractFunction, S}) where S <: MOI.AbstractSet
    offset = constroffset(optimizer, ci)
    rows = constrrows(optimizer, ci)
    sqr = optimizer.sol.slack[offset .+ rows]
    tri = reorderval(sqr, S)
    return unscalecoef(tri, S)
end

function MOI.get(optimizer::Optimizer, ::MOI.ConstraintDual, ci::CI{<:MOI.AbstractFunction, S}) where S <: MOI.AbstractSet
    offset = constroffset(optimizer, ci)
    rows = constrrows(optimizer, ci)
    sqr = optimizer.sol.x[offset .+ rows]
    tri = reorderval(sqr, S)
    if S == MOI.PositiveSemidefiniteConeTriangle
        n = sqrdim(length(rows))
        for i in 1:n, j in 1:(i-1)
            # Add lower diagonal dual. It should be equal to upper diagonal dual
            # but `unscalecoef` will divide by 2 so it will do the mean
            tri[trimap(i, j)] += sqr[i + (j-1) * n]
        end
    end
    return unscalecoef(tri, S)
end

MOI.get(optimizer::Optimizer, ::MOI.ResultCount) = 1
