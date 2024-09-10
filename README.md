# Signal Reconstruction in Diffusion-Based Molecular Communication

**Overview**

This repository contains MATLAB codes developed to simulate the reconstruction of molecular signals in diffusion-based molecular communication (MC) systems. Additionally, the simulation codes for the optimal receiver nanomachine design are provided. These codes are based on the model proposed in the following paper:

**Citation:**

Atakan, B., & Gulec, F. (2019). "Signal reconstruction in diffusion-based molecular communication." *Transactions on Emerging Telecommunications Technologies*, 30(12), e3699. doi: [10.1002/ett.3699](https://doi.org/10.1002/ett.3699)

**Background**

Molecular communication is a prominent nanoscale communication paradigm where nanomachines (NMs) form nanonetworks by transmitting information using molecules. In this model:

- A Receiver Nanomachine (RN) senses and reconstructs the molecular signal by measuring the molecule concentration.
- The molecular signal around the RN is treated as a Gaussian random process, rather than the deterministic models often used in the literature.
- The reconstructed signal is derived as a doubly stochastic Poisson process, with signal distortion introduced as a key performance parameter.
- This distortion is used to derive RN design parameters such as the RN radius and sampling period.

The MATLAB codes simulate the random walk and reconstruction process, validating the derived distortion function and optimizing RN design parameters to minimize signal distortion.

**Code Overview**

The code file names correspond to the figure names in the paper for easier reference.

- **Random Walk Simulation:** Simulates the diffusion process where molecules propagate in the reception volume of the RN.
- **Signal Reconstruction:** Reconstructs the signal based on molecule concentrations sensed by the RN and calculates signal distortion.
- **Receiver Design:** Uses derived formulas to find the optimal RN design parameters that minimize signal distortion.

**Results**

The simulation results demonstrate that:

- The original signal can be satisfactorily reconstructed with low distortion.
- There is a trade-off among RN design parameters, allowing for joint optimization to achieve the desired level of signal distortion.

**How to Cite**

If you use or build upon these codes in your research, please ensure that you cite the paper:

Atakan, B., & Gulec, F. (2019). "Signal reconstruction in diffusion-based molecular communication." *Transactions on Emerging Telecommunications Technologies*, 30(12), e3699. doi: [10.1002/ett.3699](https://doi.org/10.1002/ett.3699)
