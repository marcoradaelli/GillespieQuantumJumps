# GillespieQuantumJumps v0.2

## Important!
At the moment, the repository is set as private. In the Julia documentation, it is not very clear whether it is possible to download from a private repository (see for instance [here](https://discourse.julialang.org/t/more-problems-trying-to-add-packages-from-private-repos/69059)). The best way is to temporarely turn the repo into a public one, follow the next steps, and then reset it to private.

## Installation instructions
Open a Julia terminal session, then press `]` to enter the `Packages` REPL. Insert the command:
```
add "https://github.com/marcoradaelli/GillespieQuantumJumps.git"
```
and the package will be installed. For further instructions on how to manage Julia packages, see [here](https://docs.julialang.org/en/v1/stdlib/Pkg/).

# Usage
Once installed, the package has to be imported in each project with the command:
```
using Gillespie
```

## Pure states, complete monitoring 
### Jumps only
```
Gillespie.gillespie(H, M_l, ψ0, t_final, dt, number_trajectories, verbose)
```
Returns a tuple `(trajectories_results, V, t_range)`, where:
* `trajectories_results` is a vector of vectors. Each exterior vector corresponds to a trajectory, and each internal vector corresponds to a specific jump on such trajectory. For instance, `trajectories_results[3][4]` is the fourth jump on the third trajectory. For each jump, a dictionary is recorded, with keys `AbsTime` (absolute time since the beginning of the evolution), `TimeSinceLast` (time since previous jump), `JumpChannel` (number of the jump channel), `ψAfter` (state after the jump). **Notice.** The initial point of the evolution is always recorded as a first jump, with `nothing` as `JumpChannel`. 
* `V` is a vector of no-jump evolution operators, computed at all times specified in `t_range`.
* `t_range` is a vector of times.


#### Arguments 
* `H`: Hermitian Hamiltonian of the system
* `M_l`: list of jump operators
* `ψ0`: initial state
* `t_final`: final time of the evolution
* `dt`: time step
* `number_trajectories`: number of trajectories to be computed
* `verbose`: if `true`, increases the amount of output to be printed out (for large calculations, can fill up completely the text buffer)


### States at all times
```
Gillespie.compute_states_at_times(H, M_l, ψ0, t_final, dt, number_trajectories, verbose)
```
Returns a vector of vectors `result`, such that `result[4][6]` is the state on the trajectory 4 at time `t_range[6]`, where `t_range = 0:dt:t_final`.

#### Arguments 
* `H`: Hermitian Hamiltonian of the system
* `M_l`: list of jump operators
* `ψ0`: initial state
* `t_final`: final time of the evolution
* `dt`: time step
* `number_trajectories`: number of trajectories to be computed
* `verbose`: if `true`, increases the amount of output to be printed out (for large calculations, can fill up completely the text buffer)

## Mixed states and/or partial monitoring
### Jumps only
```
Gillespie.gillespie_partial_monitoring(H, M_l, S_l, ρ0, t_final, dt, number_trajectories, verbose)
```
Returns a tuple `(trajectories_results, V, t_range)`, where:
* `trajectories_results` is a vector of vectors. Each exterior vector corresponds to a trajectory, and each internal vector corresponds to a specific jump on such trajectory. For instance, `trajectories_results[3][4]` is the fourth jump on the third trajectory. For each jump, a dictionary is recorded, with keys `AbsTime` (absolute time since the beginning of the evolution), `TimeSinceLast` (time since previous jump), `JumpChannel` (number of the jump channel), `ρAfter` (state after the jump). **Notice.** The initial point of the evolution is always recorded as a first jump, with `nothing` as `JumpChannel`. 
* `V` is a vector of no-jump evolution operators, computed at all times specified in `t_range`, already in vectorized form.
* `t_range` is a vector of times.


#### Arguments 
* `H`: Hermitian Hamiltonian of the system
* `M_l`: list of monitored jump operators
* `S_l`: list of un-monitored jump operators
* `ρ0`: initial state
* `t_final`: final time of the evolution
* `dt`: time step
* `number_trajectories`: number of trajectories to be computed
* `verbose`: if `true`, increases the amount of output to be printed out (for large calculations, can fill up completely the text buffer)

### States at all times
```
Gillespie.compute_states_at_times_mixed(H, M_l, S_l, ρ0, t_final, dt, number_trajectories, verbose)
```
Returns a vector of vectors `result`, such that `result[4][6]` is the state on the trajectory 4 at time `t_range[6]`, where `t_range = 0:dt:t_final`.

#### Arguments 
* `H`: Hermitian Hamiltonian of the system
* `M_l`: list of monitored jump operators
* `S_l`: list of un-monitored jump operators
* `ρ0`: initial state
* `t_final`: final time of the evolution
* `dt`: time step
* `number_trajectories`: number of trajectories to be computed
* `verbose`: if `true`, increases the amount of output to be printed out (for large calculations, can fill up completely the text buffer)
