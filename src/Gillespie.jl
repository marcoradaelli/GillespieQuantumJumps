module Gillespie

using LinearAlgebra
using Random
using StatsBase
using Plots
using Polynomials
using ProgressMeter

"""
     verify_working() 
     
Verifies whether the library has been imported correctly. When called, prints "The package has been imported correctly" and a version number and returns true.

# Returns
- true
"""
function verify_working()
    println("The package has been imported correctly, version 0.2b")
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

"""
    gillespie(
        H::Matrix{ComplexF64},
        M_l::Vector{Matrix{ComplexF64}},
        ψ0::Vector{ComplexF64},
        t_final::Float64,
        dt::Float64,
        number_trajectories::Int64,
        verbose::Bool=false)

Simulates the jumps, according to the Gillespie algorithm, for the given dynamics.

# Arguments
- `H`: Hamiltonian matrix
- `M_l`: list of jump operators
- `ψ0`: initial state of the system
- `t_final`: final time of the evolution
- `dt`: time increment considered
- `number_trajectories`: number of trajectories of the simulation
- `verbose`: if true, gives more output. Verbose=true for large simulations can fill up completely the text buffer and cause a crash

# Returns
- `trajectories_results`: list of dictionaries with each jump channel, times and states after jumps
- `V`: list of pre-computed no-jump non-Hermitian evolution operators
- `t_range`: range of times at which the V operators are computed
"""
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

"""
    state_at_time_on_trajectory(
        t_range::StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, Int64},
        relevant_times::StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, Int64},
        V::Vector{Matrix{ComplexF64}},
        trajectory_data::Vector{Dict{String, Any}})

Taking as input the output of the `gillespie` function, fills the gaps between the jumps.

# Arguments
- `t_range`: range of times at which the V operators are computed
- `relevant_times`: list of times at which the state has to be computed (can also not coincide with `t_range`)
- `V`: list of non-Hermitian evolution operators, computed at times specified in `t_range`
- `trajectory_data`: a list of dictionaries in the form output by the `gillespie` function

# Returns 
- `v_states`: vector of quantum pure states (in vector form) at each of the times requested in `relevant_times`
"""
function state_at_time_on_trajectory(
    t_range::StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, Int64},
    relevant_times::StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, Int64},
    V::Vector{Matrix{ComplexF64}},
    trajectory_data::Vector{Dict{String, Any}},)

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

"""
    expectation_at_time_on_trajectory(
        t_range::StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, Int64},
        relevant_times::StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, Int64},
        V::Vector{Matrix{ComplexF64}},
        trajectory_data::Vector{Dict{String, Any}},
        E_l::Vector{Matrix{ComplexF64}})

Taking as input the output of the `gillespie` function, computes the expectation values of operators along the trajectory.
THIS FUNCTION HAS A RELEVANT BUG, STILL UNDER TEST

# Arguments
- `t_range`: range of times at which the V operators are computed
- `relevant_times`: list of times at which the state has to be computed (can also not coincide with `t_range`)
- `V`: list of non-Hermitian evolution operators, computed at times specified in `t_range`
- `trajectory_data`: a list of dictionaries in the form output by the `gillespie` function
- `E_l`: list of Hermitian measurement operators of which the expectation value has to be computed

# Returns
- `expectations_v`: list of lists of expectation values for all operators and times, such that `expectations_v[n_E][n_t]` is the expectation value for the operator indexed as n_E at time indexed by n_t.

"""
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
                expectations_v[n_E][n_t] = exp_val.re
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
            expectations_v[n_E][n_t] = exp_val.re
        end
    end

    return expectations_v
end

"""
    compute_states_at_times(
        H::Matrix{ComplexF64},
        M_l::Vector{Matrix{ComplexF64}},
        ψ0::Vector{ComplexF64},
        t_final::Float64,
        dt::Float64,
        number_trajectories::Int64,
        verbose::Bool=false,
        compute_V_each_step::Bool=false)

Function for external access, computes the states at the specified times (using both the `gillespie` and the `state_at_time_on_trajectory` when appropriate).

# Arguments
- `H`: system Hamiltonian
- `M_l`: list of jump operators
- `ψ0`: initial (pure) state of the system
- `t_final`: final time of the evolution
- `dt`: time step for the evolution
- `number_trajectories`: number of trajectories to be considered
- `verbose`: if true, gives more output. Verbose=true for large simulations can fill up completely the text buffer and cause a crash
- `compute_V_each_step`; if true, does not re-use pre-computed values for the no-jump evolution operator at all steps, but computes the operator directly.

# Returns
- `results`: list of lists of states along all the trajectories for all the considered times
"""
function compute_states_at_times(
    H::Matrix{ComplexF64},
    M_l::Vector{Matrix{ComplexF64}},
    ψ0::Vector{ComplexF64},
    t_final::Float64,
    dt::Float64,
    number_trajectories::Int64,
    verbose::Bool=false,
    compute_V_each_step=false)

    trajectories_results, V, t_range = gillespie(H, M_l, ψ0, t_final, dt, number_trajectories, verbose)
    println()

    results = Vector{Vector{ComplexF64}}[]

    # Constructs the overall jump operator.
    J = zero(M_l[1])
    for M in M_l
        J += M' * M
    end
    # Effective (non-Hermitian) Hamiltonian.
    He = H - 1im/2. * J

    if compute_V_each_step
        @showprogress 1 "Filling in the gaps..." for n_trajectory in eachindex(trajectories_results)
            v_states = state_at_time_on_trajectory_recomputing_V(t_range, t_range, trajectories_results[n_trajectory], He)
            push!(results, v_states)
        end
    else
        @showprogress 1 "Filling in the gaps..." for n_trajectory in eachindex(trajectories_results)
            v_states = state_at_time_on_trajectory(t_range, t_range, V, trajectories_results[n_trajectory])
            push!(results, v_states)
        end
    end
    

    return results
end

"""
    compute_expectation_values_at_times(
        H::Matrix{ComplexF64},
        M_l::Vector{Matrix{ComplexF64}},
        E_l::Vector{Matrix{ComplexF64}},
        ψ0::Vector{ComplexF64},
        t_final::Float64,
        dt::Float64,
        number_trajectories::Int64,
        verbose::Bool=false)

Function for external access, computes the expectation values of all the required operators at the specified times.
THIS FUNCTION HAS A RELEVANT BUG, STILL UNDER TEST
    
# Arguments
- `H`: system Hamiltonian
- `M_l`: list of jump operators
- `ψ0`: initial (pure) state of the system
- `t_final`: final time of the evolution
- `dt`: time step for the evolution
- `number_trajectories`: number of trajectories to be considered
- `verbose`: if true, gives more output. Verbose=true for large simulations can fill up completely the text buffer and cause a crash

# Returns
- `results`: list of lists of expectation values along all the trajectories for all the considered times for all the operators, such that `results[n_traj][n_E][n_t]` is the expectation value for the operator numbered by `n_E` on the trajectory `n_traj` at time `n_t`

"""
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
        v_expectations = expectation_at_time_on_trajectory(t_range, t_range, V, trajectories_results[n_trajectory], E_l)
        push!(results, v_expectations)
    end

    return results
end


function state_at_time_on_trajectory_recomputing_V(
    t_range::StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, Int64},
    relevant_times::StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, Int64},
    trajectory_data::Vector{Dict{String, Any}},
    He::Matrix{ComplexF64})

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
            Δt = t_abs - jump_times[n_jump]
            ev_op = exp(- 1im * He * Δt)
            norm = sqrt(ψ' * ev_op' * ev_op * ψ)
            ψ = ev_op * ψ
            ψ = ψ / norm
            push!(v_states, ψ)
        end
    end

    # Now computes the state for all times after the latest jump.
    last_jump_absolute_time = last(jump_times)
    relevant_times_after_last_jump = [t for t in relevant_times if t >= last_jump_absolute_time]
    for t_abs in relevant_times_after_last_jump
        ψ = last(ψ_after_jumps)
        Δt = t_abs - last_jump_absolute_time
        ev_op = exp(- 1im * He * Δt)
        norm = sqrt(ψ' * ev_op' * ev_op * ψ)
        ψ = ev_op * ψ
        ψ = ψ / norm
        push!(v_states, ψ)
    end

    return v_states
end


"""
    vectorize(x)

Vectorizes an operator.

# Arguments
- `x`: operator to be vectorized

# Returns
- the vectorized form of `x`
"""
function vectorize(x)
    return vec(x)
end

"""
    unvectorize(x)

Translates a vectorized operator into the matrix form.

# Arguments
- `x`; operator in vectorized form

# Returns
- the operator in matrix form
"""
function unvectorize(x)
    # Assumes the input object to be a reshapable to a square matrix.
    d = trunc(Int, sqrt(length(x)))
    return reshape(x, (d,d))
end

"""
    gillespie_partial_monitoring(
        H::Matrix{ComplexF64},
        M_l::Vector{Matrix{ComplexF64}},
        S_l::Vector{Matrix{ComplexF64}},
        ρ0::Matrix{ComplexF64},
        t_final::Float64,
        dt::Float64,
        number_trajectories::Int64,
        verbose::Bool=false)

Gillespie-based simulation with partial monitoring, and allowing for mixed initial states.

# Arguments
- `H`: Hamiltonian matrix
- `M_l`: list of monitored jump operators
- `S_l`: list of un-monitored jump operators
- `ρ0`: initial state of the system
- `t_final`: final time of the evolution
- `dt`: time increment considered
- `number_trajectories`: number of trajectories of the simulation
- `verbose`: if true, gives more output. Verbose=true for large simulations can fill up completely the text buffer and cause a crash

# Returns
- `trajectories_results`: list of dictionaries with each jump channel, times and states after jumps
- `V`: list of pre-computed no-jump non-Hermitian evolution operators
- `t_range`: range of times at which the V operators are computed
"""
function gillespie_partial_monitoring(
    H::Matrix{ComplexF64},
    M_l::Vector{Matrix{ComplexF64}},
    S_l::Vector{Matrix{ComplexF64}},
    ρ0::Matrix{ComplexF64},
    t_final::Float64,
    dt::Float64,
    number_trajectories::Int64,
    verbose::Bool=false)

    # Range of times.
    t_range = 0.:dt:t_final

    # M_l is the list of jumps corresponding to monitored channels.
    # S_l is the list of jumps corresponding to un-monitored channels.

    # Creates the appropriate identity matrix.
    d = size(H)[1]
    ide = Matrix{Float64}(I, d, d)

    # Creates the J operator.
    J = zero(ide)
    # Cycle over the monitored operators.
    for M in M_l
        J += M' * M
    end

    # Vectorized version of J.
    vect_J = vectorize(J)

    # Vectorized form of L_0.
    # Hamiltonian part.
    vect_L0 = -1im * kron(ide, H) + 1im * kron(transpose(H), ide)
    # Cycle over the un-monitored operators.
    for S in S_l
        vect_L0 += kron(conj.(S), S) - 0.5 * kron(ide, S' * S) - 0.5 * (transpose(S' * S), ide)
    end
    # Cycle over the monitored operators.
    for M in M_l
        vect_L0 += - 0.5 * kron(ide, M' * M) - 0.5 * kron(transpose(M' * M), ide)
    end
    
    # Vectorized form of L_0^\dagger.
    # Hamiltonian part.
    vect_L0_dagger = 1im * kron(ide, H) - 1im * kron(transpose(H), ide)
    # Cycle over the un-monitored operators.
    for S in S_l
        vect_L0_dagger += kron(transpose(S), S') - 0.5 * kron(transpose(S' * S), ide) - 0.5 * kron(ide, S' * S)
    end
    # Cycle over the monitored operators.
    for M in M_l
        vect_L0_dagger += - 0.5 * kron(transpose(M' * M), ide) - 0.5 * kron(ide, M' * M)
    end

    # Pre-computation stage.
    # Creates the list of the no-jump evolution operators and the non-state dependent part of the waiting time distribution.
    V = Matrix{ComplexF64}[] 
    Qs = Matrix{ComplexF64}[]
    for t in t_range
        ev_op = exp(vect_L0 * t)
        push!(V, ev_op)
        nsd_wtd = unvectorize(exp(vect_L0_dagger * t) * vect_J)
        push!(Qs, nsd_wtd)
    end

    # TODO: Some way of quantifying the error (like the norm of the latest Qs in the normal version)

    # Vector for the results of the computation.
    trajectories_results = Array{Dict{String, Any}}[]

    # Cycle over the trajectories.
    for trajectory in 1:number_trajectories
        # Initial state.
        ρ = ρ0
        # Absolute time.
        τ = 0

        # Creates the array of results for the single trajectory, and pushes the initial state as a fictitious fist jump.
        results = Dict{String, Any}[]
        dict_initial = Dict("AbsTime" => 0,
            "TimeSinceLast" => 0,
            "JumpChannel" => nothing,
            "ρAfter" => ρ0)
        push!(results, dict_initial)

        while τ < t_final 
            dict_jump = Dict()

            # Compute the waiting time distribution, exploiting the pre-computed part.
            Ps = Float64[]
            for Q in Qs
                wtd =  real(tr(Q * ρ))
                push!(Ps, wtd)
            end

            # Sample from the waiting time distribution.
            n_T = sample(1:length(t_range), Weights(Ps))

            # Increase the absolute time.
            τ += t_range[n_T]
            merge!(dict_jump, Dict("AbsTime" => τ, "TimeSinceLast" => t_range[n_T]))

            # Update the state.
            vect_ρ = V[n_T] * vectorize(ρ)
            ρ = unvectorize(vect_ρ)
            # Chooses where to jump.
            weights = Float64[]
            for M in M_l
                weight = real(tr(M' * M * ρ))
                push!(weights, weight)
            end
            n_jump = sample(1:length(M_l), Weights(weights))
            merge!(dict_jump, Dict("JumpChannel" => n_jump))
            # Update the state after the jump.
            ρ = M_l[n_jump] * ρ * (M_l[n_jump])'
            norm_state = real(tr(ρ))
            # Renormalize the state.
            ρ = ρ / norm_state
            merge!(dict_jump, Dict("ρAfter" => ρ))

            if verbose
                println(string(dict_jump))
            end

            push!(results, dict_jump)
        end

        push!(trajectories_results, results)
    end

    return trajectories_results, V, t_range
end

end # module