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
const OPTIMIZER_CONSTRUCTOR = MOI.OptimizerWithAttributes(
    CDCS.Optimizer,
    MOI.Silent() => true,
    "maxIter" => 4000,
)
const OPTIMIZER = MOI.instantiate(OPTIMIZER_CONSTRUCTOR)

@testset "SolverName" begin
    @test MOI.get(OPTIMIZER, MOI.SolverName()) == "CDCS"
end

@testset "supports_allocate_load" begin
    @test MOIU.supports_allocate_load(OPTIMIZER, false)
    @test !MOIU.supports_allocate_load(OPTIMIZER, true)
end

const BRIDGED =
    MOI.instantiate(OPTIMIZER_CONSTRUCTOR, with_bridge_type = Float64)
const CONFIG = MOIT.TestConfig(atol = 3e-2, rtol = 3e-2)

@testset "Unit" begin
    MOIT.unittest(
        BRIDGED,
        CONFIG,
        [
            # `NumberOfThreads` not supported.
            "number_threads",
            # `TimeLimitSec` not supported.
            "time_limit_sec",
            # Need to investigate...
            "solve_with_lowerbound",
            "solve_affine_deletion_edge_cases",
            "solve_blank_obj",
            "solve_single_variable_dual_min",
            "solve_single_variable_dual_max",
            # Need https://github.com/JuliaOpt/MathOptInterface.jl/issues/529
            "solve_qp_edge_cases",
            # Error using cdcs_hsde.preprocess (line 14)
            # No variables in your problem?
            "solve_unbounded_model",
            # Integer and ZeroOne sets are not supported
            "solve_integer_edge_cases",
            "solve_objbound_edge_cases",
            "solve_zero_one_with_bounds_1",
            "solve_zero_one_with_bounds_2",
            "solve_zero_one_with_bounds_3",
        ],
    )
end

@testset "Continuous linear problems" begin
    MOIT.contlineartest(BRIDGED, CONFIG, [
        # Need to investigate...
        "linear12",
        "linear15",
    ])
end

@testset "Continuous conic problems" begin
    MOIT.contconictest(
        BRIDGED,
        CONFIG,
        [
            # ITERATION_LIMIT
            "normone1v",
            "normone1f",
            # rotatedsoc2: Returns Inf and -Inf instead of infeasibility certificate
            "rotatedsoc2",
            # Need to investigate...
            "psdt3",
            "psds3",
            # Unsupported cones
            "pow",
            "dualpow",
            "rootdets",
            "exp",
            "dualexp",
            "logdet",
            "normspec",
            "normnuc",
            "relentr",
        ],
    )
end
