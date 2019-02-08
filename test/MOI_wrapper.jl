using Test

using MathOptInterface
const MOI = MathOptInterface
const MOIT = MOI.Test
const MOIU = MOI.Utilities
const MOIB = MOI.Bridges

# Iterations:
# linear5 : > 1000, < 2000
# linear9 : > 3000, < 4000
# linear15: > 20000, I don't know if ever converges so we exclude it
import CDCS
optimizer = CDCS.Optimizer(verbose=0, maxIter=4000)

@testset "SolverName" begin
    @test MOI.get(optimizer, MOI.SolverName()) == "CDCS"
end

@testset "supports_allocate_load" begin
    @test MOIU.supports_allocate_load(optimizer, false)
    @test !MOIU.supports_allocate_load(optimizer, true)
end

MOIU.@model(ModelData, (), (),
            (MOI.Zeros, MOI.Nonnegatives, MOI.SecondOrderCone,
             MOI.RotatedSecondOrderCone, MOI.PositiveSemidefiniteConeTriangle),
            (), (), (), (MOI.VectorOfVariables,), (MOI.VectorAffineFunction,))

# UniversalFallback is needed for starting values, even if they are ignored by CDCS
const cache = MOIU.UniversalFallback(ModelData{Float64}())
const cached = MOIU.CachingOptimizer(cache, optimizer)

# Essential bridges that are needed for all tests
const bridged = MOIB.Vectorize{Float64}(MOIB.NonposToNonneg{Float64}(cached))

config = MOIT.TestConfig(atol=3e-2, rtol=3e-2)

@testset "Unit" begin
    MOIT.unittest(MOIB.SplitInterval{Float64}(bridged),
                  config,
                  [# Need to investigate...
                   "solve_with_lowerbound", "solve_affine_deletion_edge_cases", "solve_blank_obj",
                   # Quadratic functions are not supported
                   "solve_qcp_edge_cases", "solve_qp_edge_cases",
                   # Integer and ZeroOne sets are not supported
                   "solve_integer_edge_cases", "solve_objbound_edge_cases"])
end

@testset "Continuous linear problems" begin
    MOIT.contlineartest(MOIB.SplitInterval{Float64}(bridged),
                        config, ["linear12", "linear15"])
end

@testset "Continuous conic problems" begin
    MOIT.contconictest(MOIB.SquarePSD{Float64}(MOIB.RootDet{Float64}(MOIB.GeoMean{Float64}(MOIB.RSOC{Float64}(bridged)))),
                       config, [# See https://github.com/JuliaOpt/MathOptInterface.jl/pull/632,
                                "rotatedsoc1v",
                                "rotatedsoc2", "rootdets", "exp", "logdet"])
end
