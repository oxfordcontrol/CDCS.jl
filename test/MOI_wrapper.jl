module TestCDCS

using Test
using MathOptInterface
import CDCS

const MOI = MathOptInterface

function runtests()
    for name in names(@__MODULE__; all = true)
        if startswith("$(name)", "test_")
            @testset "$(name)" begin
                getfield(@__MODULE__, name)()
            end
        end
    end
    return
end

function test_solver_name()
    @test MOI.get(CDCS.Optimizer(), MOI.SolverName()) == "CDCS"
end

function test_options()
    optimizer = CDCS.Optimizer()
    MOI.set(optimizer, MOI.RawOptimizerAttribute("printlevel"), 1)
    @test MOI.get(optimizer, MOI.RawOptimizerAttribute("printlevel")) == 1
end

function test_runtests()
    model = MOI.Utilities.CachingOptimizer(
        MOI.Utilities.UniversalFallback(MOI.Utilities.Model{Float64}()),
        MOI.instantiate(CDCS.Optimizer, with_bridge_type = Float64),
    )
    # `Variable.ZerosBridge` makes dual needed by some tests fail.
    MOI.Bridges.remove_bridge(
        model.optimizer,
        MathOptInterface.Bridges.Variable.ZerosBridge{Float64},
    )
    MOI.set(model, MOI.Silent(), true)
    MOI.set(model, MOI.RawOptimizerAttribute("maxIter"), 4000)
    MOI.Test.runtests(
        model,
        MOI.Test.Config(
            rtol = 3e-2,
            atol = 3e-2,
            exclude = Any[
                MOI.ConstraintBasisStatus,
                MOI.VariableBasisStatus,
                MOI.ObjectiveBound,
                MOI.SolverVersion,
            ],
        ),
        exclude = String[
            # `ITERATION_LIMIT`
            "test_conic_NormOneCone_VectorAffineFunction",
            "test_solve_VariableIndex_ConstraintDual_MAX_SENSE",
            "test_solve_VariableIndex_ConstraintDual_MIN_SENSE",
            "test_conic_NormOneCone_VectorOfVariables",
            "test_infeasible_MAX_SENSE",
            "test_linear_VectorAffineFunction_empty_row",
            # TODO CDCS just returns infinite and NaN values
            # See https://github.com/jump-dev/MathOptInterface.jl/issues/1759
            "test_conic_RotatedSecondOrderCone_INFEASIBLE",
            "test_linear_INFEASIBLE_2",
            "test_solve_DualStatus_INFEASIBILITY_CERTIFICATE_Interval_lower",
            "test_solve_DualStatus_INFEASIBILITY_CERTIFICATE_Interval_upper",
            "test_infeasible_MIN_SENSE",
            "test_infeasible_MIN_SENSE_offset",
            "test_infeasible_affine_MAX_SENSE",
            "test_infeasible_affine_MAX_SENSE_offset",
            "test_infeasible_affine_MIN_SENSE",
            "test_infeasible_affine_MIN_SENSE_offset",
            # No variables in your problem?
            "test_attribute_SolveTimeSec",
            "test_solve_TerminationStatus_DUAL_INFEASIBLE",
            "test_objective_ObjectiveFunction_blank",
            "test_attribute_RawStatusString",
        ],
    )
    return
end

end  # module

TestCDCS.runtests()
