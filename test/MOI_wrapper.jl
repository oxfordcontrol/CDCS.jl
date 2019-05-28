using Test

using MathOptInterface
const MOI = MathOptInterface
const MOIT = MOI.Test
const MOIU = MOI.Utilities
const MOIB = MOI.Bridges

# Iterations:
# linear5 : > 1000, < 2000
# linear9 : > 3000, < 4000
# linear15: > 20000, Don't know if ever converges so we exclude it
import CDCS
optimizer = CDCS.Optimizer(maxIter=4000)
MOI.set(optimizer, MOI.Silent(), true)

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

const bridged = MOIB.full_bridge_optimizer(cached, Float64)

config = MOIT.TestConfig(atol=3e-2, rtol=3e-2)

@testset "Unit" begin
    MOIT.unittest(bridged, config,
                  [# Need to investigate...
                   "solve_with_lowerbound", "solve_affine_deletion_edge_cases", "solve_blank_obj",
                   # Quadratic functions are not supported
                   "solve_qcp_edge_cases", "solve_qp_edge_cases",
                   # Integer and ZeroOne sets are not supported
                   "solve_integer_edge_cases", "solve_objbound_edge_cases"])
end

@testset "Continuous linear problems" begin
    MOIT.contlineartest(bridged, config, ["linear12", "linear15"])
end

@testset "Continuous conic problems" begin
    MOIT.contconictest(bridged, config,
                       [# rotatedsoc2: Returns Inf and -Inf instead of infeasibility certificate
                        "rotatedsoc2", "rootdets", "exp", "logdet"])
end
