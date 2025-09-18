# Causal Preconditioning Filters Design for Real-Time Multichannel ANC

This repository contains the code accompanying the paper:

> **Yiming Wang, Yongjie Zhuang, and Yangfan Liu**,  
> "Causal Preconditioning Filters Design for Real-Time Multichannel Active Noise Control",  
> *Applied Acoustics*, 2025.
>
> DOI: https://doi.org/10.1016/j.apacoust.2025.110950
> 
> Here is a free preprint copy:
> [https://yongjie-zhuang.com/files/preprint_2025_preconditioning_ANC.pdf](https://yongjie-zhuang.com/files/preprint_2025_preconditioning_ANC.pdf)

The code implements causal preconditioning filters for multichannel active noise control (ANC),  
including both **MATLAB** and **Python** versions. Both versions are expected to produce the same results  
(with minor numerical precision differences possible).

---

## Folder Structure
```text
.
├── LICENSE
├── README.md                <-- this file
├── Matlab_version/          <-- MATLAB implementation
│   ├── src/                 <-- main algorithms
│   ├── example/             <-- example usage and figure outputs
│   └── README.md            <-- MATLAB-specific instructions
└── Python_version/          <-- Python implementation
    ├── src/                 <-- main algorithms
    ├── example/             <-- example usage and figure outputs
    └── README.md            <-- Python-specific instructions
```

---

## Quick Start

### MATLAB Version
1. Navigate to the `Matlab_version/example` folder.
2. Open `example.m` in MATLAB and directly run the `example.m` file.
Figures will be generated and saved in example/figures.
(See Matlab/README.md for more details.)

### Python Version
1. Navigate to the `Python_version/example` folder.
2. Open `example.py` and directly run the `example.py` file.
Figures will be generated and saved in example/figures.
(See python/README.md for details.)

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
