# Electromagnetic and Thermo-Hydraulic Design of a 15 T HTS CICC Solenoid

![Politecnico di Torino](https://img.shields.io/badge/University-Politecnico_di_Torino-blue.svg)
![Field](https://img.shields.io/badge/Field-Superconducting_Magnet_Engineering-green.svg)
![Technology](https://img.shields.io/badge/Technology-REBCO_HTS_CICC-red.svg)

---

## Overview
This repository contains the computational models, technical documentation, experimental data analysis, and calculation scripts developed for the design of a **high-field HTS solenoid** based on **REBCO coated conductors** arranged in a **Cable-in-Conduit Conductor (CICC)** configuration.

The project focuses on the electromagnetic and thermo-hydraulic design of a solenoid composed of stacked **double pancakes**, with particular attention to the central double pancake located on the equatorial plane. The target operating conditions are a magnetic field of approximately **15 T** and a nominal transport current of **50 kA**, reached through a linear current ramp of **200 A/s**.

The analysis combines superconducting-tape characterization, magnetic-field calculations, nonlinear electrical modeling, AC-loss estimation, and transient helium cooling simulations in order to demonstrate the feasibility of the proposed conductor and winding configuration while maintaining a minimum current-sharing temperature margin of **5 K**.

## Key Objectives & Analysis Features
- **REBCO Tape Characterization:** Uses HERMES coated-conductor data from Faraday Factory to determine the dependence of critical current and power-law exponent on temperature and magnetic field, supported by laboratory I-V measurements.
- **Solenoid Geometry Optimization:** Determines the number of double pancakes and radial turns required to reproduce the target magnetic field while respecting the imposed winding dimensions.
- **Three-Dimensional Magnetic-Field Evaluation:** Computes external, mutual, and self-field contributions along the central double pancake and evaluates the local field experienced by the superconducting tape stacks.
- **Critical Current & Current-Sharing Analysis:** Determines the spatial distribution of critical current and current-sharing temperature along the conductor, including a conservative assessment of supplier variability.
- **Nonlinear DAE Electrical Model:** Models the six equivalent REBCO tape stacks through a differential-algebraic system including inductive coupling and nonlinear superconducting resistance, verifying natural current equipartition during the ramp.
- **CICC Cross-Section Design:** Defines a six-slot copper former containing stacks of 80 REBCO tapes, surrounded by a stainless-steel jacket and cooled through a central supercritical-helium channel.
- **AC Loss Evaluation:** Calculates magnetization hysteresis and inter-stack coupling losses during the 0–50 kA current ramp and provides the distributed heat source required by the thermal model.
- **Transient Thermo-Hydraulic Analysis:** Implements a reduced GANDALF-type model coupling the superconducting stacks, copper former, steel jacket, and forced helium flow.
- **Thermal-Margin Verification:** Determines the minimum helium mass flow required to maintain a current-sharing margin above 5 K throughout the complete current ramp and post-ramp transient.
- **Operating-Window Assessment:** Evaluates the trade-off between helium inlet temperature and mass flow, identifying the feasible operating frontier of the cooling system.

## Methodology & Modeling Approach
1. **Superconducting Material Database:** Critical-current lift factors are obtained from manufacturer data for HERMES REBCO tapes, while a separate experimental dataset is used to characterize the nonlinear power-law exponent \(n(T,B)\). Laboratory I-V measurements provide an independent qualitative validation of the adopted analysis procedure.

2. **Electromagnetic Design:** The solenoid geometry is constructed as an axial stack of double pancakes. The final configuration consists of **29 double pancakes with 7 radial turns each**, producing a peak conductor field of approximately **15.07 T** at the nominal current of 50 kA.

3. **Current Distribution Model:** The central double pancake is represented using six twisted equivalent strands, each corresponding to a stack of 80 REBCO tapes. A nonlinear DAE model follows the current distribution during the 250 s ramp and predicts a maximum deviation from ideal current equipartition of only about **0.65%**.

4. **AC-Loss Model:** Magnetization hysteresis losses are evaluated using a Bean critical-state approximation, while inter-stack coupling losses account for the twist pitch, copper magnetoresistance, and effective transverse electrical conductivity of the CICC former. The resulting ramp-averaged AC power is approximately **179.5 W**, corresponding to about **44.9 kJ** deposited in the central double pancake.

5. **Thermo-Hydraulic Model:** The conductor is represented by four coupled one-dimensional thermal regions: REBCO stacks, copper former, stainless-steel jacket, and helium coolant. Axial conduction, inter-material heat transfer, forced convection, and helium advection are solved during both the current ramp and the subsequent thermal relaxation.

6. **Cooling-System Verification:** With helium entering at **4.5 K**, the minimum total mass flow required to maintain the prescribed 5 K thermal margin is approximately **4.0 g/s**, significantly below the available 15 g/s limit. The minimum transient margin is approximately **5.00 K**, with a maximum superconducting-stack temperature close to **16 K**.

7. **Operating Frontier:** Increasing the helium mass flow allows progressively higher inlet temperatures. At the maximum available flow of **15 g/s**, an inlet temperature of approximately **11.2 K** can still satisfy the required thermal margin.

## Repository Structure
```text
.
├── src/                # Electromagnetic, DAE, AC-loss, and thermo-hydraulic models
├── data/               # REBCO supplier data, material properties, and experimental data
├── docs/               # Technical report and project documentation
├── results/            # Magnetic-field, current-sharing, AC-loss, and thermal results
└── README.md           # Project summary and documentation
