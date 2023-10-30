# Scattering Results {#scattering}


## Results {#results}

After compiling the example program, we can run it with

```bash
./scattering parameter.prm
```

With the [parameter file](@ref parameter) given below, we expect the following
output:

```bash
Sapphire::Start example scattering with parameter file "parameter.prm" on 1 processor(s) [2023/10/30 13:17:39]
Sapphire:VFP::Run VFP equation solver.          [13:17:39]
Sapphire:VFP::Time step    1 at t = 0.00000     [13:17:39]
Sapphire:VFP::Time step    2 at t = 0.100000    [13:17:39]
Sapphire:VFP::Time step    3 at t = 0.200000    [13:17:39]
Sapphire:VFP::Time step    4 at t = 0.300000    [13:17:39]
Sapphire:VFP::Time step    5 at t = 0.400000    [13:17:39]
Sapphire:VFP::Time step    6 at t = 0.500000    [13:17:39]
Sapphire:VFP::Time step    7 at t = 0.600000    [13:17:39]
Sapphire:VFP::Time step    8 at t = 0.700000    [13:17:39]
Sapphire:VFP::Time step    9 at t = 0.800000    [13:17:39]
Sapphire:VFP::Time step   10 at t = 0.900000    [13:17:39]
Sapphire:VFP::Time step   11 at t = 1.00000     [13:17:39]
Sapphire:VFP::Time step   12 at t = 1.10000     [13:17:39]
Sapphire:VFP::Time step   13 at t = 1.20000     [13:17:39]
Sapphire:VFP::Time step   14 at t = 1.30000     [13:17:39]
Sapphire:VFP::Time step   15 at t = 1.40000     [13:17:39]
Sapphire:VFP::Time step   16 at t = 1.50000     [13:17:39]
Sapphire:VFP::Time step   17 at t = 1.60000     [13:17:39]
Sapphire:VFP::Time step   18 at t = 1.70000     [13:17:39]
Sapphire:VFP::Time step   19 at t = 1.80000     [13:17:39]
Sapphire:VFP::Time step   20 at t = 1.90000     [13:17:39]
Sapphire:VFP::The simulation ended.             [13:17:39]

+------------------------------------------+------------------+------------+------------------+
| Total wallclock time elapsed             |    0.1548s     0 |    0.1548s |    0.1548s     0 |
|                                          |                  |                               |
| Section                      | no. calls |   min time  rank |   avg time |   max time  rank |
+------------------------------------------+------------------+------------+------------------+
| DG matrix                    |        21 |   0.07757s     0 |   0.07757s |   0.07757s     0 |
| FE system                    |         1 |   0.01056s     0 |   0.01056s |   0.01056s     0 |
| Grid setup                   |         1 |  0.001882s     0 |  0.001882s |  0.001882s     0 |
| Mass matrix                  |         1 |  0.002406s     0 |  0.002406s |  0.002406s     0 |
| Output                       |        21 |    0.0196s     0 |    0.0196s |    0.0196s     0 |
| Project f onto the FEM space |         1 |  0.004495s     0 |  0.004495s |  0.004495s     0 |
| Setup                        |         1 |   0.02981s     0 |   0.02981s |   0.02981s     0 |
| Theta method                 |        20 |    0.1114s     0 |    0.1114s |    0.1114s     0 |
+------------------------------------------+------------------+------------+------------------+
Sapphire:VFP::Compute the global error
Sapphire::L2 error: 0.000654138
```

The results of the run will be saved in the 'results/scattering` folder.


### Parameter file {#parameter}

We run the program with the following parameter file:

\include scattering/parameter.prm


## The plain program {#plainCode}


### config.h {#plainConfig}

\include scattering/config.h


### scattering.cpp {#plainMain}

\include scattering/scattering.cpp

