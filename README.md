# GillespieQuantumJumps v0.2

## Important!
At the moment, the repository is set as private. In the Julia documentation, it is not very clear whether it is possible to download from a private repository (see for instance [here](https://discourse.julialang.org/t/more-problems-trying-to-add-packages-from-private-repos/69059)). The best way is to temporarely turn the repo into a public one, follow the next steps, and then reset it to private.

## Installation instructions
Open a Julia terminal session, then press `]` to enter the `Packages` REPL. Insert the command:
```
add "https://github.com/marcoradaelli/GillespieQuantumJumps.git"
```
and the package will be installed. For further instructions on how to manage Julia packages, see [here](https://docs.julialang.org/en/v1/stdlib/Pkg/).

## Usage
Once installed, the package has to be imported in each project with the command:
```
using Gillespie
```

### Pure states, complete monitoring 
```
Gillespie.gillespie(H, M_l, ψ0, t_final, dt, number_trajectories, false)
```
Returns a tuple `(trajectories_results, V, t_range)`, where:
* `trajectories_results` is a vector of vectors. Each exterior vector corresponds to a trajectory, and each internal vector corresponds to a specific jump on such trajectory. For instance, `trajectories_results[3][4]` is the fourth jump on the third trajectory. For each jump, a dictionary is recorder, with keys `AbsTime` (absolute time since the beginning), `TimeSinceLast` (time since previous jump), `JumpChannel` (number of the jump channel), `ψAfter` (state after the jump). 
