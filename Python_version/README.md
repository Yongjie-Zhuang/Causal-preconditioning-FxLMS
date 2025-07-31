# Python Implementation: Causal Preconditioning Filters for Real-Time Multichannel ANC

This folder contains the **Python implementation** of the algorithms described in:

> **Yiming Wang, Yongjie Zhuang, and Yangfan Liu**,  
> *Causal Preconditioning Filters Design for Real-Time Multichannel Active Noise Control*,  
> Applied Acoustics, 2025.  
> [Free access link (valid until September 16, 2025)](https://authors.elsevier.com/a/1lVwp,5Mxwwgy)

The Python version reproduces the same functionality as the MATLAB version, using only open-source Python packages.

---

## Requirements

Python 3.10 or newer is recommended.  
Dependencies include:
- numpy
- scipy >= 1.12
- matplotlib

---
## Folder Structure

```text
Python_version/
├── README.md <-- this file
├── example/ <-- example script & output figures
│ ├── example.py <-- main usage example
│ └── figures/ <-- generated plots from example.py
└── src/ <-- core algorithm implementations
├── precondition.py <-- main preconditioning filter function
└── init.py
```

## Running the Example
1. Activate your Python virtual environment (optional but recommended). Install necessary dependencies.
2. Navigate to example folder:
    ```bash
    cd Python_version/example
3. Run:
   ```bash
   python example.py
4. The script will:
    - Compute preconditioning filters.
    - Generate magnitude/phase and impulse response figures.
    - Save all figures automatically to Python_version/example/figures.

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

