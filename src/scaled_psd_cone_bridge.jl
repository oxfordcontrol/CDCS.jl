struct ScaledPSDCone <: MOI.AbstractVectorSet
    side_dimension::Int
end

Base.copy(x::ScaledPSDCone) = ScaledPSDCone(x.side_dimension)

MOI.side_dimension(x::ScaledPSDCone) = x.side_dimension

function MOI.dimension(x::ScaledPSDCone)
    return x.side_dimension^2
end

function MOI.Utilities.set_with_dimension(::Type{ScaledPSDCone}, dim)
    return ScaledPSDCone(isqrt(dim))
end

struct ScaledPSDConeBridge{T,G} <: MOI.Bridges.Constraint.SetMapBridge{
    T,
    ScaledPSDCone,
    MOI.PositiveSemidefiniteConeTriangle,
    MOI.VectorAffineFunction{T},
    G,
}
    constraint::MOI.ConstraintIndex{MOI.VectorAffineFunction{T},ScaledPSDCone}
end

function MOI.Bridges.Constraint.concrete_bridge_type(
    ::Type{ScaledPSDConeBridge{T}},
    ::Type{G},
    ::Type{MOI.PositiveSemidefiniteConeTriangle},
) where {T,G<:Union{MOI.VectorOfVariables,MOI.VectorAffineFunction{T}}}
    return ScaledPSDConeBridge{T,G}
end

function MOI.Bridges.map_set(
    ::Type{<:ScaledPSDConeBridge},
    set::MOI.PositiveSemidefiniteConeTriangle,
)
    return ScaledPSDCone(set.side_dimension)
end

function MOI.Bridges.inverse_map_set(
    ::Type{<:ScaledPSDConeBridge},
    set::ScaledPSDCone,
)
    return MOI.PositiveSemidefiniteConeTriangle(set.side_dimension)
end

# Contrarily to SeDuMi, CDCS does not work if the A_i are not symmetric
# we move half of off-diagonal (i, j) coefficients to (j, i)

# Map ConstraintFunction from MOI -> CDCS
function MOI.Bridges.map_function(
    BT::Type{<:ScaledPSDConeBridge{T}},
    func::MOI.VectorOfVariables,
) where {T}
    new_f = MOI.Utilities.operate(*, Float64, 1.0, func)
    return MOI.Bridges.map_function(BT, new_f)
end
function MOI.Bridges.map_function(
    ::Type{<:ScaledPSDConeBridge},
    f::MOI.VectorAffineFunction,
)
    n = MOI.output_dimension(f)
    d = MOI.Utilities.side_dimension_for_vectorized_dimension(n)
    constants = triangle_to_square(f.constants, d)
    terms = copy(f.terms)
    triangle_to_square_indices!(terms, d)
    return MOI.VectorAffineFunction(terms, constants)
end

# Used to map the ConstraintPrimal from CDCS -> MOI
# No need to unscale (i, j) because half was moved to (j, i)
function MOI.Bridges.inverse_map_function(::Type{<:ScaledPSDConeBridge}, square)
    return square_to_triangle(square)
end

# Used to map the ConstraintDual from CDCS -> MOI
function MOI.Bridges.adjoint_map_function(::Type{<:ScaledPSDConeBridge}, square)
    triangle = square_to_triangle(square)
    n = isqrt(length(square))
    return triangle
end

function square_map(i::Integer, j::Integer, n::Integer)
    return i + (j - 1) * n
end

function copy_upper_triangle(x, n, map_from, map_to)
    y = zeros(eltype(x), map_to(n, n))
    for i in 1:n, j in 1:i
        y[map_to(i, j)] = x[map_from(i, j)]
    end
    return y
end
function square_to_triangle(x, n = isqrt(length(x)))
    return copy_upper_triangle(
        x,
        n,
        (i, j) -> square_map(i, j, n),
        MOI.Utilities.trimap,
    )
end
function triangle_to_square(x, n)
    y = zeros(eltype(x), square_map(n, n, n))
    k = 0
    for j in 1:n, i in 1:j
        k += 1
        y[square_map(i, j, n)] = y[square_map(j, i, n)] = x[k]
    end
    return y
end

function triangle_to_square_indices!(x::Vector{<:MOI.VectorAffineTerm}, n)
    map = Vector{Tuple{Int,Int}}(undef, MOI.Utilities.trimap(n, n))
    for j in 1:n, i in 1:j
        map[MOI.Utilities.trimap(i, j)] = (i, j)
    end
    for k in eachindex(x)
        i, j = map[x[k].output_index]
        t = x[k].scalar_term
        x[k] = MOI.VectorAffineTerm(square_map(i, j, n), t)
        if i != j
            push!(x, MOI.VectorAffineTerm(square_map(j, i, n), t))
        end
    end
end
