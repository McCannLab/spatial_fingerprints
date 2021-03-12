# Julia code part for simulation using ML

- `helpers.jl`: core functions,
- `scripts.jl`: call core function to do simulation,
- `install.packages.jl` install packages needed,
- `main.jl` for increasing number of bio-tracers,
- `main2.jl` was used for noise and training set size.

Plot were done in R.

## simulations

See `jobs` 

### nbio

Basically used [Compute Canada](https://www.computecanada.ca/) services,

```sh
julia -p 1 main.jl 0 b c d e f g
```

For simulation nbio :

c = 20

### ndistr

```sh
julia -p 1 main2.jl 2 b c d e f g
```

### noise

```sh
julia -p 1 main2.jl 3 b c d e f g
```