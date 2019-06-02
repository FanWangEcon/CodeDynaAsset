This is a work-in-progress [website](https://fanwangecon.github.io/CodeDynaAsset/) of collections of code for solving several infinite horizon exogenous incomplete dynamic asset models, produced by [Fan Wang](https://fanwangecon.github.io/). Materials gathered from various [papers](https://fanwangecon.github.io/research) on financial access with fixed costs, discrete and continuous asset choice grids and other features.

Generally, looped, vectorized, and optimized-vectorized implementations of the same solution algorithm with tabular, graphical and profiling results are shown. Looped codes are shown for clarity, vectorized codes are shown for speed. Codes are designed to not require special hardware or explicit parallelization. Codes tested on Windows 10 with [Matlab 2019a](https://www.mathworks.com/company/newsroom/mathworks-announces-release-2019a-of-matlab-and-simulink.html) for replicability. Please contact [FanWangEcon](https://fanwangecon.github.io/) for problems.

All functions are written with default parameters and are directly callable: 1, clone the project, see [here](docs/gitsetup.md); 2, add to path; 3, click run. There are three types of files:

1. **m**: matlab m file
2. **publish html**: html files generated by matlab publish, includes code mark-ups and tabular and graphical outputs from benchmark simulation
3. **profile**: html files generated by profiling the m file with timing results


# 1. The One Asset One Shock Problem (BZ)
<!-- https://fanwangecon.github.io/CodeDynaAsset/m_az -->

The *bz* problem: standard model with an asset and one shock, exogenous incomplete borrowing and savings, wage shocks follow AR1.

## 1.1 Main Optimization Solution Files (BZ)

Parameters can be adjusted [here](https://fanwangecon.github.io/CodeDynaAsset/m_az/paramfunc/html/ffs_az_set_default_param.html), for the benchmark simulation:

- savings only
- **750** grid points for asset states/choices
- **15** grid points for the AR1 shock

Using three algorithm that provide identical solutions:

1. *bz* model [looped solution](https://fanwangecon.github.io/CodeDynaAsset/m_az/solve/html/ff_az_vf.html): [**m**](https://github.com/FanWangEcon/CodeDynaAsset/blob/master/m_az/solve/ff_az_vf.m) \| [**publish html**](https://fanwangecon.github.io/CodeDynaAsset/m_az/solve/html/ff_az_vf.html) \| [**profile**](https://fanwangecon.github.io/CodeDynaAsset/m_az/solve/profile/ff_az_vf_default_p3/file0.html)
    * speed: **9765.7** seconds
    * loops: 1 for VFI, 1 for shocks, 1 for asset state, 1 for asset choice, 1 for future shocks
2. *bz* model [vectorized solution](https://fanwangecon.github.io/CodeDynaAsset/m_az/solve/html/ff_az_vf_vec.html): [**m**](https://github.com/FanWangEcon/CodeDynaAsset/blob/master/m_az/solve/ff_az_vf_vec.m) \| [**publish html**](https://fanwangecon.github.io/CodeDynaAsset/m_az/solve/html/ff_az_vf_vec.html) \| [**profile**](https://fanwangecon.github.io/CodeDynaAsset/m_az/solve/profile/ff_az_vf_vec_default_p3/file0.html)    
    * speed: **34.3** seconds
    * loops: 1 for VFI, 1 for shocks, vectorize remaining 3 loops
3. *bz* model [optimized vectorized solution](https://fanwangecon.github.io/CodeDynaAsset/m_az/solve/html/ff_az_vf_vecsv.html): [**m**](https://github.com/FanWangEcon/CodeDynaAsset/blob/master/m_az/solve/ff_az_vf_vecsv.m) \| [**publish html**](https://fanwangecon.github.io/CodeDynaAsset/m_az/solve/html/ff_az_vf_vecsv.html) \| [**profile**](https://fanwangecon.github.io/CodeDynaAsset/m_az/solve/profile/ff_az_vf_vecsv_default_p3/file0.html)    
    * speed: **1.2** seconds
    * loops: 1 for VFI, 1 for shocks, vectorize remaining 3 loops, reuse u(c)
    * reuse u(c) in cells, several speed improvements described [here](https://fanwangecon.github.io/M4Econ/)

## 1.2 Asset Distributions (BZ)

Solving for the asset distribution.

## 1.3 Solution Support Files (BZ)

**Parameters and Function Definitions**:
1. *bz* model [set default parameters](https://fanwangecon.github.io/CodeDynaAsset/m_az/paramfunc/html/ffs_az_set_default_param.html): [**m**](https://github.com/FanWangEcon/CodeDynaAsset/blob/master/m_az/paramfunc/ffs_az_set_default_param.m) \| [**publish html**](https://fanwangecon.github.io/CodeDynaAsset/m_az/paramfunc/html/ffs_az_set_default_param.html)
    * param_map: container map for carrying parameters across functions
    * support_map: container map for carrying programming instructions etc across functions
2. *bz* model [set functions](https://fanwangecon.github.io/CodeDynaAsset/m_az/paramfunc/html/ffs_az_set_functions.html): [**m**](https://github.com/FanWangEcon/CodeDynaAsset/blob/master/m_az/paramfunc/ffs_az_set_functions.m) \| [**publish html**](https://fanwangecon.github.io/CodeDynaAsset/m_az/paramfunc/html/ffs_az_set_functions.html)
    * functions: centrally define functions as function handles
3. *bz* model [generate states, choices, and shocks grids](https://fanwangecon.github.io/CodeDynaAsset/m_az/paramfunc/html/ffs_az_get_funcgrid.html): [**m**](https://github.com/FanWangEcon/CodeDynaAsset/blob/master/m_az/paramfunc/ffs_az_get_funcgrid.m) \| [**publish html**](https://fanwangecon.github.io/CodeDynaAsset/m_az/paramfunc/html/ffs_az_get_funcgrid.html)
    * func_map: container map containing all function handles
    * armt_map: container map containing matrixes for states and choices.

**Output Analysis**:
1.  *bz* model [solution results processing](https://fanwangecon.github.io/CodeDynaAsset/m_az/solvepost/html/ff_az_vf_post.html): [**m**](https://github.com/FanWangEcon/CodeDynaAsset/blob/master/m_az/solvepost/ff_az_vf_post.m) \| [**publish html**](https://fanwangecon.github.io/CodeDynaAsset/m_az/solvepost/html/ff_az_vf_post.html)
    * table: value and policy function by states and shocks
    * table: iteration convergence and percentage policy function change by shock
    * graph: value + policy functions with levels, logged levels and percentages
    * mat: store all workspace matrixes, arrays, scalar values to matrixes
2. *bz* model [solution results graphing](https://fanwangecon.github.io/CodeDynaAsset/m_az/solvepost/html/ff_az_vf_post_graph.html): [**m**](https://github.com/FanWangEcon/CodeDynaAsset/blob/master/m_az/solvepost/ff_az_vf_post_graph.m) \| [**publish html**](https://fanwangecon.github.io/CodeDynaAsset/m_az/solvepost/html/ff_az_vf_post_graph.html)
    * graph: value function by asset and shock
    * graph: consumption and asset choice levels
    * graph: consumption and asset logged levels
    * graph: consumption and asset as percentages of coh and assets


## 1.4 Cash-on-hand and One Shock (COH-Z)

The problem could be posed slightly differently with the asset state variable as cash-on-hand. The codes are basically identical, and are shown here in this [folder](https://fanwangecon.github.io/CodeDynaAsset/m_oz). Speeds and outcomes are the same. This is useful for thinking about finding asset distributions.

**Main Optimization Solution Files**:
- *coh-z* model [looped solution](https://fanwangecon.github.io/CodeDynaAsset/m_oz/solve/html/ff_oz_vf.html), [vectorized solution](https://fanwangecon.github.io/CodeDynaAsset/m_oz/solve/html/ff_oz_vf_vec.html), [optimized vectorized solution](https://fanwangecon.github.io/CodeDynaAsset/m_oz/solve/html/ff_oz_vf_vecsv.html)

**Parameters and Function Definitions**:
- *coh-z* model [set default parameters](https://fanwangecon.github.io/CodeDynaAsset/m_oz/paramfunc/html/ffs_oz_set_default_param.html), [set functions](https://fanwangecon.github.io/CodeDynaAsset/m_oz/paramfunc/html/ffs_oz_set_functions.html), [generate states, choices, and shocks grids](https://fanwangecon.github.io/CodeDynaAsset/m_oz/paramfunc/html/ffs_oz_get_funcgrid.html)

**Output Analysis**:
- *coh-z* model [solution results processing](https://fanwangecon.github.io/CodeDynaAsset/m_oz/solvepost/html/ff_oz_vf_post.html), [solution results graphing](https://fanwangecon.github.io/CodeDynaAsset/m_oz/solvepost/html/ff_oz_vf_post_graph.html)


# 2. The Risky + Safe Asset Problem

Two endogenous assets, one safe one risky. Risky asset could be stock with constant return to scale, or physical capital investment with depreciation and decreasing return to scale. Note that the utility function is CRRA, however, households do not have constant share of risky investment for any wealth (cash-on-hand) levels when risky asset has decreasing return to scale and when shock is highly persistent.

There are more analytical ways of solving the basic version of this problem. Here we stick to using this grid based solution algorithm which allows for flexibly solving non-differentiable and non-continuous problems. The grid based solution algorithm now with 2 endogenous choices and states requires exponentially more computation time than the *bz* model. Here I provide three sets of solution algorithms at increasing speeds:

- In **2.1**, solve the problem with the two asset choice explicitly
- In **2.2**, solve the problem in two stages
- In **2.3**, two stage solution with interpolation

## 2.1 Explicitly Solve for Safe and Risky Assets (BKZ)

The *bkz* problem. Parameters can be adjusted [here](https://fanwangecon.github.io/CodeDynaAsset/m_akz/paramfunc/html/ffs_akz_set_default_param.html), for the benchmark simulation:

- savings problem with alternative safe and risky assets
- **45** aggregate savings grid points, **1035** combinations of safe and risky asset choices. Choice grids shown [here](https://fanwangecon.github.io/CodeDynaAsset/m_akz/paramfunc/html/ffs_akz_get_funcgrid.html) at the end of the file. Note that the benchmark parameters for the *bz* model has *750* grid points of choices/states, the benchmark problem here is larger and hence takes more time.
- **15** grid points for the AR1 shock

1. *bkz* model [looped solution](https://fanwangecon.github.io/CodeDynaAsset/m_akz/solve/html/ff_akz_vf.html): [**m**](https://github.com/FanWangEcon/CodeDynaAsset/blob/master/m_akz/solve/ff_akz_vf.m) \|
    * speed: **19891.3** seconds
    * loops: 1 for VFI, 1 for shocks, 1 for coh(b,k,z), 1 for (b',k') choices, 1 for future shocks
2. *bkz* model [vectorized solution](https://fanwangecon.github.io/CodeDynaAsset/m_akz/solve/html/ff_akz_vf_vec.html): [**m**](https://github.com/FanWangEcon/CodeDynaAsset/blob/master/m_akz/solve/ff_akz_vf_vec.m) \| [**publish html**](https://fanwangecon.github.io/CodeDynaAsset/m_akz/solve/html/ff_akz_vf_vec.html) \| [**profile**](https://fanwangecon.github.io/CodeDynaAsset/m_akz/solve/profile/ff_akz_vf_vec_default_p3/file0.html)    
    * speed: **71.5** seconds
    * loops: 1 for VFI, 1 for shocks, vectorize remaining
3. *bkz* model [optimized vectorized solution](https://fanwangecon.github.io/CodeDynaAsset/m_akz/solve/html/ff_akz_vf_vecsv.html): [**m**](https://github.com/FanWangEcon/CodeDynaAsset/blob/master/m_akz/solve/ff_akz_vf_vecsv.m) \| [**publish html**](https://fanwangecon.github.io/CodeDynaAsset/m_akz/solve/html/ff_akz_vf_vecsv.html) \| [**profile**](https://fanwangecon.github.io/CodeDynaAsset/m_akz/solve/profile/ff_akz_vf_vecsv_default_p3/file0.html)    
    * speed: **2.1** seconds
    * loops: 1 for VFI, 1 for shocks, vectorize remaining
    * reuse u(c) in cells, several speed improvements described [here](https://fanwangecon.github.io/M4Econ/)

## 2.2 Two-stage Solution, First Aggregate Savings (WKZ)

The *wkz* problem, w=k'+b'. Takes significantly less time than *2.1*, produces identical results. Parameters can be adjusted [here](https://fanwangecon.github.io/CodeDynaAsset/m_akz/paramfunc/html/ffs_akz_set_default_param.html), for the benchmark simulation, same as *2.1*:

- savings problem with alternative safe and risky assets
- **45** aggregate savings grid points, **1035** combinations of safe and risky asset choices.
- **15** grid points for the AR1 shock

1. *wkz* model [looped solution](https://fanwangecon.github.io/CodeDynaAsset/m_akz/solve/html/ff_akz_vf.html): [**m**](https://github.com/FanWangEcon/CodeDynaAsset/blob/master/m_akz/solve/ff_akz_vf.m) \|
    * speed: **720.7** seconds (*72* times faster than *2.1*)
    * Step One solve k*(w,z); Step Two solve w*(z,coh(b,k,z)) given k*(w,z)
    * loops: 1 for VFI, 1 for shocks, 1 for coh(b,k,z), 1 for w(z)=k'+b'
2. *wkz* model [vectorized solution](https://fanwangecon.github.io/CodeDynaAsset/m_akz/solve/html/ff_akz_vf_vec.html): [**m**](https://github.com/FanWangEcon/CodeDynaAsset/blob/master/m_akz/solve/ff_akz_vf_vec.m) \| [**publish html**](https://fanwangecon.github.io/CodeDynaAsset/m_akz/solve/html/ff_akz_vf_vec.html) \| [**profile**](https://fanwangecon.github.io/CodeDynaAsset/m_akz/solve/profile/ff_akz_vf_vec_default_p3/file0.html)    
    * speed: **4.1** seconds (*17* times faster than *2.1*)
    * Step One solve k*(w,z); Step Two solve w*(z,coh(b,k,z)) given k*(w,z)
    * loops: 1 for VFI, 1 for shocks, vectorize remaining
3. *wkz* model [optimized vectorized solution](https://fanwangecon.github.io/CodeDynaAsset/m_akz/solve/html/ff_akz_vf_vecsv.html): [**m**](https://github.com/FanWangEcon/CodeDynaAsset/blob/master/m_akz/solve/ff_akz_vf_vecsv.m) \| [**publish html**](https://fanwangecon.github.io/CodeDynaAsset/m_akz/solve/html/ff_akz_vf_vecsv.html) \| [**profile**](https://fanwangecon.github.io/CodeDynaAsset/m_akz/solve/profile/ff_akz_vf_vecsv_default_p3/file0.html)    
    * speed: **0.8** seconds (*2.6* times faster than *2.1*)
    * Step One solve k*(w,z); Step Two solve w*(z,coh(b,k,z)) given k*(w,z)
    * loops: 1 for VFI, 1 for shocks, vectorize remaining
    * store u(c) in cells, update when k*(w,z) changes, several speed improvements described [here](https://fanwangecon.github.io/M4Econ/)

## 2.3 Two-stage Solution, First Aggregate Savings (WKZ)


## 2.4 Asset Distributions (BKZ)

Solving for the asset distribution.

## 2.5 Solution Support Files (Shared)

All solution algorithms share the same support files.

**Parameters and Function Definitions**:
1. *bkz* model [set default parameters](https://fanwangecon.github.io/CodeDynaAsset/m_akz/paramfunc/html/ffs_akz_set_default_param.html): [**m**](https://github.com/FanWangEcon/CodeDynaAsset/blob/master/m_akz/paramfunc/ffs_akz_set_default_param.m) \| [**publish html**](https://fanwangecon.github.io/CodeDynaAsset/m_akz/paramfunc/html/ffs_akz_set_default_param.html)
    * param_map: container map for carrying parameters across functions
    * support_map: container map for carrying programming instructions etc across functions
2. *bkz* model [set functions](https://fanwangecon.github.io/CodeDynaAsset/m_akz/paramfunc/html/ffs_akz_set_functions.html): [**m**](https://github.com/FanWangEcon/CodeDynaAsset/blob/master/m_akz/paramfunc/ffs_akz_set_functions.m) \| [**publish html**](https://fanwangecon.github.io/CodeDynaAsset/m_akz/paramfunc/html/ffs_akz_set_functions.html)
    * functions: centrally define functions as function handles
3. *bkz* model [generate states, choices, and shocks grids](https://fanwangecon.github.io/CodeDynaAsset/m_akz/paramfunc/html/ffs_akz_get_funcgrid.html): [**m**](https://github.com/FanWangEcon/CodeDynaAsset/blob/master/m_akz/paramfunc/ffs_akz_get_funcgrid.m) \| [**publish html**](https://fanwangecon.github.io/CodeDynaAsset/m_akz/paramfunc/html/ffs_akz_get_funcgrid.html)
    * func_map: container map containing function handles
    * armt_map: container map containing matrixes for states and choices.

**Output Analysis**:
1.  *bkz* model [solution results processing](https://fanwangecon.github.io/CodeDynaAsset/m_akz/solvepost/html/ff_akz_vf_post.html): [**m**](https://github.com/FanWangEcon/CodeDynaAsset/blob/master/m_akz/solvepost/ff_akz_vf_post.m) \| [**publish html**](https://fanwangecon.github.io/CodeDynaAsset/m_akz/solvepost/html/ff_akz_vf_post.html)
    * table: value and policy function by states and shocks
    * table: convergence and percentage policy function change by shock
    * graph: value + policy functions with levels, logged levels and percentages
    * mat: store all workspace matrixes, arrays, scalar values to matrixes
2. *bkz* model [solution results graphing](https://fanwangecon.github.io/CodeDynaAsset/m_akz/solvepost/html/ff_akz_vf_post_graph.html): [**m**](https://github.com/FanWangEcon/CodeDynaAsset/blob/master/m_akz/solvepost/ff_akz_vf_post_graph.m) \| [**publish html**](https://fanwangecon.github.io/CodeDynaAsset/m_akz/solvepost/html/ff_akz_vf_post_graph.html)
    * graph: value function by asset and shock
    * graph: consumption and asset choice levels
    * graph: consumption and asset logged levels
    * graph: consumption and asset as percentages of coh and assets
