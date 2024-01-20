# CervixRadiomic# Project Title: Cervix Signature Generation from Radiomic Biomarkers

## Overview

This GitHub repository contains code for generating a cervix signature based on radiomic biomarkers. The analysis utilizes pyradiomics features from two datasets: TCGIA and Norway. Feature selection is performed using mrmrE, and the final model is constructed using 100 combinations of 2 radiomics features. The resulting signature is then validated on an internal PMH cohort.

## Table of Contents

1. [Introduction](#introduction)
2. [Datasets](#datasets)
3. [Code Structure](#code-structure)
4. [Usage](#usage)
5. [Dependencies](#dependencies)
6. [Results](#results)
7. [Validation](#validation)
8. [License](#license)
9. [Acknowledgements](#acknowledgements)

## Introduction

Radiomic biomarkers offer valuable insights into cervical cancer. This project aims to generate a cervix signature by leveraging pyradiomics features from the TCGIA and Norway datasets. Feature selection is crucial, and mrmrE is employed to identify the most relevant radiomics features. The final model combines 100 sets of 2 radiomics features to create a robust signature.

## Datasets

- **TCGIA Dataset:** The Cancer Genome Atlas's cervical cancer dataset.
- **Norway Dataset:** Dataset from Norway containing relevant radiomic information.
- **PMH Dataset:** Internal Dataset from PMH used to validate prognostic potential of radiomics signature

## Code Structure

The code is organized as follows:

- **train:** Contains data, notebook used to generate final signature. Multiple potential signatures were tested using CV of training set. A final model using an ensemble of 100 sets of 2 radiomics features was extracted as final signature
- **archive** Scripts for any extra analysis conducted.
- **validation:** Validates the generated signature on an internal PMH cohort.

## Usage

To use the code, follow these steps:

1. Clone the repository:

   ```bash
   git clone https://github.com/bhklab/CervixRadiomics
   ```

2. Navigate to the project directory:

   ```bash
   cd CervixRadiomics
   ```

3. Run a jupyter notebook server to conduct analysis.

   ```bash
   jupyter notebook
   ```

## Dependencies

Ensure you have the following dependencies installed:

- Python (>=3.6)
- pymrmre (1.0.7)

Install dependencies using:

```bash
pip install matplotlib pandas numpy seaborn lifelines scikit-learn scikit-image pymrmre
```

## Results

The generated cervix signature is located in the `results` directory. Detailed results and performance metrics are also provided.

## Validation

The signature is validated on an internal PMH cohort to assess its generalizability.

## License

This project is licensed under the [MIT License](LICENSE).

## Acknowledgements

Special thanks to Dr. Kathy Han & Dr. Benjamin Haibe Kans, Dr. Philip Ye, Jessica Weiss & all corresponding contributors and data providers for making this research possible.

Feel free to contribute, report issues, or suggest improvements!s