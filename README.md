# MMM1D
Steady-state, non-isothermal, two-phase simulation of low-temperature PEM fuel cells in 1D

<img src="https://romanvetter.ethz.ch/software/mmm1d/image1.png" alt="Screenshot 1" height="200"> <img src="https://romanvetter.ethz.ch/software/mmm1d/image2.png" alt="Screenshot 2" height="200"> <img src="https://romanvetter.ethz.ch/software/mmm1d/image3.png" alt="Screenshot 3" height="200">

## Original publication
**Free open reference implementation of a two-phase PEM fuel cell model**  
Roman Vetter, Juergen O. Schumacher  
*Computer Physics Communications* 234, 222-234 (2019)  
<https://doi.org/10.1016/j.cpc.2018.07.023>

## License
See [LICENSE](LICENSE).

## Installation & system requirements
No installation is needed. Just open [MMM1D.m](MMM1D.m) in MATLAB R2018b or newer.

## Execution
To run the program, open it in MATLAB and execute the function either by clicking the "Run" button on the MABLAB GUI or by typing `MMM1D` in MATLAB's command window.

- The program generates three figures with plots for all fields, fluxes and the polarization curve.
- By default, the program iterates over all cell voltages in steps of 50 mV.
- The expected runtime on a modern PC is approximately two seconds for the default setup.
- By default, moderate error tolerances are used (relative: `1e-4`, absolute: `1e-6`). Decrease these values (e.g., relative: `1e-6`, absolute: `1e-10`) for high-precision results.

## Input & output
There is no input needed to run a simulation. To change material parameters or operating conditions, edit the source code directly.
Running the program by typing `[IU, SOL] = MMM1D;` in MATLAB's command window will return the following two output data structures:
- `IU`: a two-column matrix containing the list of computed current densities and corresponding cell voltages (i.e., the polarization curve) of the fuel cell
- `SOL`: a cell array of solution structures (one for each row in `IU`) containing all fields and fluxes

## Changes since original publication
- Better readable and generalizable code: named rather than indexed quantities and layers
- Dependency on the Curve Fitting Toolbox removed
- Redundancy in formulas reduced
- Code simplified and modernized (now requires MATLAB R2018b or newer)
- Various improvements (comments etc.)
