# Observational-Astro-Analysis: Measuring the Habitability of Exoplanets Detected Using Transit Photometry

## Project Overview

This repository documents an observational astrophysics project focused on applying key data analysis algorithms to **Kepler-4** flux data to detect orbiting exoplanets and assess their potential habitability.

The core objective was to perform **transit photometry** to identify exoplanets, determine their orbital characteristics (like period and radius), and analyse their physical properties (like mass-period relations and density) to conclude their likelihood of being habitable.

## Key Findings

The analysis successfully identified and characterized four distinct exoplanets orbiting the Kepler-4 star system:

| Planet | Orbital Period (days) | Estimated Radius ($R_{\oplus}$) | Consistent Mass-Period Relation |
| :---: | :---: | :---: | :---: |
| **F** | $13.1754 \pm 0.0002$ | $2.92 \pm 0.14$ | Rocky Planet |
| **E** | $21.7760 \pm 0.0004$ | $5.11 \pm 0.22$ | Super-Earth-like |
| **L** | $31.784 \pm 0.0016$ | $3.83 \pm 0.17$ | Super-Earth-like |
| **S** | $41.023 \pm 0.0021$ | $4.05 \pm 0.18$ | Warm-Neptune |

**Habitability Conclusion:** Based on an assessment across the **habitability zone**, **exoplanet density**, and **escape velocity**, it was concluded that **none of the detected Kepler-4 exoplanets are currently habitable.**

---

## Methodology & Algorithms

The exoplanet detection and characterization process involved the application of several key algorithms and techniques:

1.  **Data Preprocessing:** **Data Normalization** and **Savitzky-Golay Filtering** were used to clean the raw flux data, remove noise, and prepare it for analysis.
2.  **Period Determination:** **Lomb-Scargle Periodograms** were implemented to search for periodic signals indicative of planetary transits, thus determining the orbital periods.
3.  **Lightcurve Analysis:** **Curve Fitting** techniques were applied to the phased lightcurves to model the transit shape and estimate key parameters, such as the planetary radius.
4.  **Habitability Assessment:** Habitability was determined by calculating the location of the planets relative to the star's **Habitable Zone (HZ)**, estimating **planetary density**, and computing the planet's **escape velocity**.

---

##  Repository Structure

The main files within this repository are structured as follows:

| File/Folder | Type | Description |
| :--- | :--- | :--- |
| `code/` | Directory | Contains supporting Python scripts (`mytools.py`, `mytools2.py`) for custom functions and data processing. |
| `Coursework_B_25406.ipynb` | Jupyter Notebook | The main project file containing all the analysis steps, code, visualizations, and commentary. |
| `Exoplanet_Report.pdf` | PDF Document | The complete formal written report detailing the background, methodology, results, and discussion. **(Recommended starting point)** |
| `nasa.csv` | Data File | The primary dataset used for the analysis (likely containing the Kepler-4 flux measurements). |
| `image*.jpg` | Image Files | Supporting images and figures used in the analysis or report. |

---

## Requirements and Installation

This project is written primarily in Python and utilizes standard scientific libraries.

### Prerequisites

* Python 3

### Installation

1.  **Clone the repository:**
    ```bash
    git clone [https://github.com/FunmiLS/Observational-Astro-Analysis.git](https://github.com/FunmiLS/Observational-Astro-Analysis.git)
    cd Observational-Astro-Analysis
    ```
2.  **Install dependencies:**
    You will likely need common scientific packages such as `numpy`, `pandas`, `matplotlib`, and `astropy` (or `scipy` for Lomb-Scargle/Savitzky-Golay).

    ```bash
    pip install numpy pandas matplotlib scipy astropy jupyter
    ```
3.  **Launch the notebook:**
    ```bash
    jupyter notebook Coursework_B_25406.ipynb
    ```

## License

This project is licensed under the **MIT License** - see the included `LICENSE` file for details (if not present, consider adding one).

---

