module Gillespie

using LinearAlgebra
using Random
using StatsBase
using Plots
using Polynomials
using ProgressMeter

"""
     verify_working()
     
Verifies whether the library has been imported correctly. When called, prints "The package has been imported correctly" and returns true.

# Returns
- true
"""
function verify_working()
    println("The package has been imported correctly")
    return true
end

"""
    find_nearest(a,x)

Returns the nearest index in array a, assumed sorted, to the value x.

# Arguments
- `a`: array-like sorted object
- `x`: number

# Returns
- `n`: index in `a` of the closest element to `x`.
"""
function find_nearest(a,x)
    length(a) > 0 || return 0:-1
    r = searchsorted(a,x)
    length(r) > 0 && return r
    last(r) < 1 && return searchsorted(a,a[first(r)])
    first(r) > length(a) && return searchsorted(a,a[last(r)])
    x-a[last(r)] < a[first(r)]-x && return searchsorted(a,a[last(r)])
    x-a[last(r)] > a[first(r)]-x && return searchsorted(a,a[first(r)])
    return first(searchsorted(a,a[last(r)])):last(searchsorted(a,a[first(r)]))
end

function gillespie(
    H::Matrix{ComplexF64},
    M_l::Vector{Matrix{ComplexF64}},
    ψ0::Vector{ComplexF64},
    t_final::Float64,
    dt::Float64,
    number_trajectories::Int64,
    verbose::Bool=false)

    t_range = 0.:dt:t_final

    # Constructs the overall jump operator.
    J = zero(M_l[1])
    for M in M_l
        J += M' * M
    end
    # Effective (non-Hermitian) Hamiltonian.
    He = H - 1im/2. * J

    # Constructs the no-jump evolution operators for all the relevant times.
    V = Matrix{ComplexF64}[] # List of the no-jump evolution operators.
    Qs = Matrix{ComplexF64}[] # List of the non-state-dependent part of the waiting time distribution.
    for t in t_range
        ev_op = exp(-1im * He * t)
        push!(V, ev_op)
        nsd_wtd = ev_op' * J * ev_op
        push!(Qs, nsd_wtd)
    end

    # Prints the matrix norm for the latest Qs.
    error = norm(last(Qs))
    println("-> Truncation error given by norm of latest Qs matrix: " * string(error))
        
    # List for the results.
    trajectories_results = Array{Dict{String, Any}}[]

    # Cycle over the trajectories.
    @showprogress 1 "Gillespie evolution..." for trajectory in 1:number_trajectories
        
        # Initial state.
        ψ = ψ0
        # Absolute time.
        τ = 0
        
        results = Dict{String, Any}[]
        dict_initial = Dict("AbsTime" => 0,
            "TimeSinceLast" => 0,
            "JumpChannel" => nothing,
            "ψAfter" => ψ0)
        push!(results, dict_initial)
        
        while τ < t_final
            dict_jump = Dict()
            
            # Compute the waiting time distribution, exploiting the pre-computed part.
            Ps = Float64[]
            for Q in Qs
                wtd = real(ψ' * Q * ψ)
                push!(Ps, wtd)
            end
            
            # Sample from the waiting time distribution.
            n_T = sample(1:length(t_range), Weights(Ps))
                    
            # Increase the absolute time.
            τ += t_range[n_T]
            merge!(dict_jump, Dict("AbsTime" => τ, "TimeSinceLast" => t_range[n_T]))
            
            # Update the state.
            ψ = V[n_T] * ψ
            # Chooses where to jump.
            weights = Float64[]
            for M in M_l
                weight = real(ψ' * M' * M * ψ)
                push!(weights, weight)
            end
            n_jump = sample(1:length(M_l), Weights(weights))
            merge!(dict_jump, Dict("JumpChannel" => n_jump))
            # Update the state after the jump.
            ψ = M_l[n_jump] * ψ
            norm_state = norm(ψ)
            # Renormalize the state.
            ψ = ψ / norm_state
            merge!(dict_jump, Dict("ψAfter" => ψ))
            
            if verbose
                println(string(dict_jump))
            end
            
            push!(results, dict_jump)
        end
        
        push!(trajectories_results, results)
    end

    return trajectories_results, V, t_range
end

function state_at_time_on_trajectory(
    t_range::StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, Int64},
    relevant_times::StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, Int64},
    V::Vector{Matrix{ComplexF64}},
    trajectory_data::Vector{Dict{String, Any}})

    # Creates an array of states.
    v_states = Vector{ComplexF64}[]
    
    # Creates an array of jump times.
    jump_times = [trajectory_data[i]["AbsTime"] for i in eachindex(trajectory_data)]
    # Creates an array of states after the jumps.
    ψ_after_jumps = [trajectory_data[i]["ψAfter"] for i in eachindex(trajectory_data)]
    
    # Cycles over the jumps times.
    for n_jump in 1:length(jump_times)-1
        next_jump_time = jump_times[n_jump + 1]
        # Determines the set of relevant times between this jump and the following one.
        relevant_times_in_interval = [t for t in relevant_times if jump_times[n_jump] <= t < next_jump_time]
        # Cycles over the relevant times.
        for t_abs in relevant_times_in_interval
            ψ = ψ_after_jumps[n_jump]
            n_t = find_nearest(t_range, t_abs - jump_times[n_jump])[1]
            norm = sqrt(ψ' * V[n_t]' * V[n_t] * ψ)
            ψ = V[n_t] * ψ
            ψ = ψ / norm
            push!(v_states, ψ)
        end
    end

    # Now computes the state for all times after the latest jump.
    last_jump_absolute_time = last(jump_times)
    relevant_times_after_last_jump = [t for t in relevant_times if t >= last_jump_absolute_time]
    for t_abs in relevant_times_after_last_jump
        ψ = last(ψ_after_jumps)
        n_t = find_nearest(t_range, t_abs - last_jump_absolute_time)[1]
        norm = sqrt(ψ' * V[n_t]' * V[n_t] * ψ)
        ψ = V[n_t] * ψ
        ψ = ψ / norm
        push!(v_states, ψ)
    end

    return v_states
end

function expectation_at_time_on_trajectory(
    t_range::StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, Int64},
    relevant_times::StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, Int64},
    V::Vector{Matrix{ComplexF64}},
    trajectory_data::Vector{Dict{String, Any}},
    E_l::Vector{Matrix{ComplexF64}})

    # Creates an array of expectation values for each operator.
    expectations_v = []
    for E in E_l
        v = zeros(length(relevant_times))
        push!(expectations_v, v)
    end
    
    # Creates an array of jump times.
    jump_times = [trajectory_data[i]["AbsTime"] for i in eachindex(trajectory_data)]
    # Creates an array of states after the jumps.
    ψ_after_jumps = [trajectory_data[i]["ψAfter"] for i in eachindex(trajectory_data)]
    
    # Cycles over the jumps times.
    for n_jump in 1:length(jump_times)-1
        next_jump_time = jump_times[n_jump + 1]
        # Determines the set of relevant times between this jump and the following one.
        relevant_times_in_interval = [t for t in relevant_times if jump_times[n_jump] <= t < next_jump_time]
        # Cycles over the relevant times.
        for t_abs in relevant_times_in_interval
            ψ = ψ_after_jumps[n_jump]
            n_t = find_nearest(t_range, t_abs - jump_times[n_jump])[1]
            norm = sqrt(ψ' * V[n_t]' * V[n_t] * ψ)
            ψ = V[n_t] * ψ
            ψ = ψ / norm
            # Cycles over the operators to compute the expectation values.
            for n_E in eachindex(E_l)
                exp_val = ψ' * E_l[n_E] * ψ
                expectations_v[n_E][n_t] = exp_val
            end
        end
    end

    # Now computes the state for all times after the latest jump.
    last_jump_absolute_time = last(jump_times)
    relevant_times_after_last_jump = [t for t in relevant_times if t >= last_jump_absolute_time]
    for t_abs in relevant_times_after_last_jump
        ψ = last(ψ_after_jumps)
        n_t = find_nearest(t_range, t_abs - last_jump_absolute_time)[1]
        norm = sqrt(ψ' * V[n_t]' * V[n_t] * ψ)
        ψ = V[n_t] * ψ
        ψ = ψ / norm
        # Cycles over the operators to compute the expectation values.
        for n_E in eachindex(E_l)
            exp_val = ψ' * E_l[n_E] * ψ
            expectations_v[n_E][n_t] = exp_val
        end
    end

    return expectations_v
end

function compute_states_at_times(
    H::Matrix{ComplexF64},
    M_l::Vector{Matrix{ComplexF64}},
    ψ0::Vector{ComplexF64},
    t_final::Float64,
    dt::Float64,
    number_trajectories::Int64,
    verbose::Bool=false)

    trajectories_results, V, t_range = gillespie(H, M_l, ψ0, t_final, dt, number_trajectories, verbose)
    println()

    results = Vector{Vector{ComplexF64}}[]

    @showprogress 1 "Filling in the gaps..." for n_trajectory in eachindex(trajectories_results)
        v_states = state_at_time_on_trajectory(t_range, t_range, V, trajectories_results[n_trajectory])
        push!(results, v_states)
    end

    return results
end

function compute_expectation_values_at_times(
    H::Matrix{ComplexF64},
    M_l::Vector{Matrix{ComplexF64}},
    E_l::Vector{Matrix{ComplexF64}},
    ψ0::Vector{ComplexF64},
    t_final::Float64,
    dt::Float64,
    number_trajectories::Int64,
    verbose::Bool=false)

    # The operators of which the expectation value has to be computed have to be passed as the list E_l.

    # Gillespie evolution.
    trajectories_results, V, t_range = gillespie(H, M_l, ψ0, t_final, dt, number_trajectories, verbose)

    results = []

    # Holes filling and computation of expectation values.
    @showprogress 1 "Filling in the gaps..." for n_trajectory in eachindex(trajectories_results)
        v_expectations = expectation_at_times_on_trajectory(t_range, t_range, V, trajectories_results[n_trajectory], E_l)
        push!(results, v_expectations)
    end

    return results
end

end # module
