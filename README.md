This repository provides the complete code and analysis tools for the paper  

**Distinct Roles of AMPA, NMDA, and GABA Kinetics in Shaping Macroscopic Cortical Dynamics**  
Hongsheng Deng, Xinkun Zhang, Hongjie Bi, Xiyun Zhang  

The repository contains:  
- A C++ implementation of a large-scale, sparsely connected excitatory–inhibitory network of quadratic integrate-and-fire (QIF) neurons with separate AMPA, NMDA, and GABA synaptic kinetics.  
- Tools for batch submission on PBS clusters and power-spectrum analysis in Matlab.  
- Detailed instructions for reproducing the bifurcation and continuation analysis of the exact low-dimensional mean-field model using MatCont.

### Default network parameters (used in all figures unless otherwise stated)
- Excitatory neurons: N_E = 10,000  
- Inhibitory neurons: N_I = 2,500  
- Mean in-degree: K = 1,000  
- Synaptic weights:  
  J^ee_AMPA = 0.6,          J^ie_AMPA = 0.6 / 1.01  
  J^ee_NMDA = 0.0178,       J^ie_NMDA = 0.0178 / 1.01  
  J^ei_GABA = 0.0155,       J^ii_GABA = 0.0155 / 1.005  
- Degree heterogeneity: Δ_k^ee = 5, Δ_k^ii = 30  
- Membrane time constant: τ_m = 20 ms  

---

## 1. Spiking Network Simulation (net.cpp)

### Compilation
```bash
g++ net.cpp -o net -std=c++11
Running a single simulation
Bash./net -tau_ampa=VALUE -tau_nmda=VALUE -tau_gaba=VALUE -i=VALUE

-tau_ampa, -tau_nmda, -tau_gaba : synaptic decay time constants (ms)
-i : external drive current I_0 (μA)

Example (parameters used in Fig. 1 of the manuscript):
Bash./net -tau_ampa=4 -tau_nmda=90 -tau_gaba=5 -i=0.1
Output files
After each run the program generates three files:

data_e.csv – spike times of all excitatory neurons
data_i.csv – spike times of all inhibitory neurons
(format: first column = dimensionless time; remaining columns = indices of neurons that fired at that time)
data.csv – population-averaged quantities (13 columns):
t 2. rall 3. vall 4. re 5. ve 6. ri 7. vi
saee 9. snee 10. sgei 11. saie 12. snie 13. sgii


Units

Time inside the simulation is dimensionless: 1 time unit = 20 ms
To convert firing rates to physical units: frequency (Hz) = rate × 50

Batch submission on PBS clusters

Edit parameters.csv (one parameter set per line, values in physical units).
Submit all combinations:Bash./submit_jobs.sh output_directory_nameEach parameter set is submitted as an independent PBS job; results are stored in the specified directory.

Power-spectrum analysis
matlabScript/Frequency.m
Loads data.csv and computes the power spectral density (used to generate Fig. 5).

2. Bifurcation and Continuation Analysis with MatCont (Mean-Field Model)
The dimensionless mean-field equations are presented as Eq. (1)–(7) in the paper.
Entering the system in MatCont

Select → System → New
Define the seven state variables in this exact order:
re, ve, ri, vi, S_ampa, S_nmda, S_gaba
Equations 1–3: define symbolically
Equations 4–7: define numerically

Computing transient dynamics

Type → Initial Point → Point
Provide a small non-zero initial condition (e.g., re = 0.01, ve = -2)
Open a 2D plot window
Compute → Forward

One-parameter continuation of equilibria

To obtain a robust starting steady state: set τ_NMDA > 516 ms, integrate forward (Interval > 500), then select the Last point.
Switch to Type → Equilibrium
Choose the continuation parameter (e.g., tau_nmda)
Set the y-axis to re (most sensitive indicator of network activity)
Compute → Forward / Backward

Branch points and limit cycles

Double-click a Branch Point (BP) → set as new Initial Point → continue in both directions.
For Hopf points:
Limit cycles: Initializer → Limit cycle (init H LC), select Period as ordinate.
Two-parameter Hopf loci: Initializer → Hopf (init H H).


Extracting oscillation frequencies
Limit-cycle frequency
matlabh = findobj(gca,'Type','line');
x = get(h,'XData');     % parameter value
T = get(h,'YData');     % period (dimensionless)
f_Hz = 50 ./ T;         % frequency in Hz
Linear frequency near fixed points (imaginary part of eigenvalues)
matlabh = findobj(gca,'Type','line');
x = get(h,'XData');
imag_part = get(h,'YData');           % rad / dimensionless time
f_Hz = imag_part / (2*pi) * 50;       % Hz

Validation and Reproducibility
Spiking network simulations closely match the mean-field predictions (see Fig. 1). All phase diagrams (Fig. 2), one-parameter bifurcation diagrams (Fig. 3), frequency analyses (Fig. 4), and power spectra from finite-size networks (Fig. 5) were produced using the exact procedures described above.
For additional plotting scripts or questions, please open an issue or contact the corresponding author:
hongsheng deng (derandeng@stu.jnu.edu.cn)
