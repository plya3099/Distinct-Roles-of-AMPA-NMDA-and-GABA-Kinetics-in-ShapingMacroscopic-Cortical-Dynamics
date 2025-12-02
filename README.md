
# README.md

This repository contains the complete code and analysis tools for the paper  

**Distinct Roles of AMPA, NMDA, and GABA Kinetics in Shaping Macroscopic Cortical Dynamics**  

The repository includes:  
- A C++ spiking network simulator with separate AMPA, NMDA, and GABA synaptic kinetics  
- PBS batch submission scripts  
- Matlab power-spectrum analysis tool  
- Detailed instructions for MatCont bifurcation analysis of the exact mean-field model  

### Default network parameters (used in all figures unless stated otherwise)
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
```

### Running a single simulation
```bash
./net -tau_ampa=VALUE -tau_nmda=VALUE -tau_gaba=VALUE -i=VALUE
```
- `-tau_ampa`, `-tau_nmda`, `-tau_gaba`: decay time constants (ms)  
- `-i`: external drive current I₀ (μA)  

**Example** (parameters used in Fig. 1 of the paper):
```bash
./net -tau_ampa=4 -tau_nmda=90 -tau_gaba=5 -i=0.1
```

### Output files
- `data_e.csv` — spike times of excitatory neurons  
- `data_i.csv` — spike times of inhibitory neurons  
  (each row: first column = dimensionless time; remaining columns = indices of neurons firing at that time)  
- `data.csv` — population-averaged variables (13 columns):  
  1. t 2. rall 3. vall 4. re 5. ve 6. ri 7. vi  
  8. saee 9. snee 10. sgei 11. saie 12. snie 13. sgii  

### Units
- Simulation time is dimensionless: **1 time unit = 20 ms**  
- To obtain physical firing rates: **frequency (Hz) = rate × 50**

### Batch runs on PBS clusters
1. Edit `parameters.csv` (one parameter combination per line, values in physical units)  
2. Submit all jobs:
```bash
./submit_jobs.sh output_directory_name
```
Each combination is submitted as an independent PBS job; results are saved in the specified directory.

### Power-spectrum analysis
```matlab
Script/Frequency.m
```
Loads `data.csv` and computes the power spectral density (used for Fig. 5).

---

## 2. Bifurcation Analysis with MatCont (Mean-Field Model)

The dimensionless mean-field equations are given as Eq. (1)–(7) in the manuscript.

### Entering the system in MatCont
- `Select → System → New`  
- Define the seven state variables **in this exact order**:  
  `re, ve, ri, vi, S_ampa, S_nmda, S_gaba`  
- Equations 1–3: define symbolically  
- Equations 4–7: define numerically  

### Computing transient dynamics
1. `Type → Initial Point → Point`  
2. Set small non-zero initial conditions (e.g., re = 0.01, ve = -2)  
3. Open a 2D plot window  
4. `Compute → Forward`

### One-parameter continuation of equilibria
1. For a robust starting steady state: set τ_NMDA > 516 ms, integrate forward (Interval > 500), then take the Last point  
2. Switch to `Type → Equilibrium`  
3. Choose continuation parameter (e.g., `tau_nmda`)  
4. Set y-axis to `re`  
5. `Compute → Forward / Backward`

### Branch points and limit cycles
- Double-click a Branch Point (BP) → set as new Initial Point → continue both directions  
- For Hopf points:  
  - Limit cycles: `Initializer → Limit cycle (init H LC)`, select `Period` as ordinate  
  - Two-parameter Hopf curves: `Initializer → Hopf (init H H)`

### Extracting oscillation frequencies

#### Limit-cycle frequency
```matlab
h = findobj(gca,'Type','line');
x = get(h,'XData');     % parameter value
T = get(h,'YData');     % period (dimensionless)
f_Hz = 50 ./ T;         % frequency in Hz
```

#### Linear frequency near fixed points
```matlab
h = findobj(gca,'Type','line');
x = get(h,'XData');
imag_part = get(h,'YData');
f_Hz = abs(imag_part) / (2*pi) * 50;   % Hz
```

---

### Validation
Spiking network simulations accurately reproduce the mean-field predictions (Fig. 1). All figures in the paper (phase diagrams, bifurcation diagrams, frequency analyses, and power spectra) were generated using the procedures described above.

For questions or additional scripts, please open an issue or contact the corresponding author:  
**HongSheng deng** (derandeng@stu.jnu.edu.cn)


