# MATLAB Version

This folder contains the MATLAB implementation of **causal preconditioning filters** for multichannel active noise control (ANC).

---

## Requirements
- MATLAB R2020b or newer (for `tiledlayout`, `exportgraphics`)
- Control System Toolbox (for `idare`)

---

## Folder Structure
```text
Matlab_version/
├── src/             <-- core algorithms
│   ├── precond_obtain_filter.m
│   ├── mul_specFact.m
├── example/         <-- demonstration script and figure outputs
│   ├── example.m
│   ├── figures/     <-- generated automatically after running
└── README.md        <-- this file
```

## Running the Example
1. Open MATLAB.
2. Navigate to `matlab/example`.
3. Run:
   ```matlab
   example
4. The script will:
    - Compute preconditioning filters.
    - Generate magnitude/phase and impulse response figures.
    - Save all figures automatically to matlab/example/figures.
5. In the example code, the secondary path was set to the same setup as in the paper. You can compare the all-pass component of the secondary path with the analytical solution: only the 6th coefficients of the diagonal elements in G_{all pass} impulse response matrix are 1, the other elements should be close to 0:
[Impulse response plot of the decomposed all-pass filter for secondary path matrix](Matlab_version/examples/figures/G_all_FreqResp.png)

### License
This repository is released under the [MIT License](LICENSE).

### Citation
If you use this code in your research, please cite:

```text
@article{WANG2025110950,
title = {Causal preconditioning filters design for real-time multichannel active noise control},
author = {Yiming Wang and Yongjie Zhuang and Yangfan Liu},
journal = {Applied Acoustics},
volume = {240},
pages = {110950},
year = {2025},
issn = {0003-682X},
doi = {https://doi.org/10.1016/j.apacoust.2025.110950},
url = {https://www.sciencedirect.com/science/article/pii/S0003682X25004220},
}
```

