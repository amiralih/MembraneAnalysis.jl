# Test data

- `s_1.pdb`, `t_1.xtc`: Martini 2 bilayer of DOPC with 10 mol% DOPE (11 frames).
- `*.dat`: reference results computed with MembraneAnalysis.jl v0.3.2 (`fluctuation_spectrum`
  with grid spacing 15 Å, `q_max = 0.08`, curvature spectra between 0.01 and 0.21 1/Å in 10
  bins, TCB with `z_cutoff = 20`). `test/test_pipeline.jl` checks that v0.4 reproduces their
  means.
