# README.md

This repository contains the code and analysis tools accompanying the paper  

**Distinct Roles of AMPA, NMDA, and GABA Kinetics in Shaping Macroscopic Cortical Dynamics**  
Hongsheng Deng, Xinkun Zhang, Hongjie Bi, Xiyun Zhang  

The repository includes:  
- A C++ implementation of a large-scale spiking network of quadratic integrate-and-fire (QIF) neurons with separate AMPA, NMDA, and GABA synaptic currents.  
- Scripts for batch submission on PBS clusters.  
- Matlab tools for power-spectrum analysis.  
- Instructions for reproducing the bifurcation and continuation analysis of the exact low-dimensional mean-field (neural mass) model using MatCont.

### Default network parameters (used throughout the manuscript)
- Excitatory population size: N_E = 10,000  
- Inhibitory population size: N_I = 2,500  
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
