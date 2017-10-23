## GreenNMFk - submodule of NMF ##

### 1. Usage ###
Run `GreenNMFk.test()` to test the integrity of the current release.  
A simple example is present in `examples/` and more will be added shortly.  

### 2. Issues ###
At this point, it appears that the solver produces solutions approximate to the
Matlab version of the source. With the exception of the following two issues
(and the lack of a few addtl. examples), `GreenNMFk` should be in a state more or less
ready to ship.

1. Clustering (from lines 71-77) needs to be verified and bug tested. Currently
will crash at `NMFk.finalize()`.  I believe this is issue actually starts with
improperly formatted inputs to `NMFk.clustersolutions_old` - but the current
structure of the inputs (particularily `[vec(col_sources)']`)  are the only
format that will get through the function without crashing.

2. A true random sampling of `x_init` from `lb` to `ub` will often generate
inputs that do not converge (at least, in a resonable timescale) in the solver.
`lsqnonlin`  seems to be able to converge for generally random inputs.
In cases other than this, `lsqnonlin`  has comparable speeds to `Mads.minimize`.
