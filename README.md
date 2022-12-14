# LTEI_TA

[![Build Status](https://travis-ci.com/sbadredd/LTEI_TA.jl.svg?branch=main)](https://travis-ci.com/sbadredd/LTEI_TA.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/sbadredd/LTEI_TA.jl?svg=true)](https://ci.appveyor.com/project/sbadredd/LTEI_TA-jl)
[![Coverage](https://codecov.io/gh/sbadredd/LTEI_TA.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/sbadredd/LTEI_TA.jl)
[![Coverage](https://coveralls.io/repos/github/sbadredd/LTEI_TA.jl/badge.svg?branch=main)](https://coveralls.io/github/sbadredd/LTEI_TA.jl?branch=main)

## Description 

The long-range Coulomb potential is a significant factor in molecular simulations, and its accurate evaluation is essential for accurate results. The present study proposes an approximation method for the numerical evaluation of the long-range Coulomb potential and the approximation of the resulting high-dimensional Two-Electron Integrals (TEI) tensor with long-range interactions. This method is based on the tensorized structure of the compressed two-electron integrals, which are obtained through two-dimensional Chebyshev interpolation combined with Gaussian quadrature. This is a prototype implementation in the Julia programming language.

## Dependencies
LTEI_TA uses a MATLAB function to find the optimal number of Chebyshev interpolation points

## DATA
Data are extracted from:<br />
-https://www.basissetexchange.org/ <br />
-https://quantumpackage.github.io/qp2/

## Licence
LTEI_TA is under  MIT licence.

## Credits
Funding from the European Research Council (ERC) under the European Union’s Horizon 2020 research
and innovation program (grant agreement No 810367).
