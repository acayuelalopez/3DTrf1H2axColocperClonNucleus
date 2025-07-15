# 3D Quantification of TRF1, γH2AX, and Their Colocalization

This repository provides a 3D image analysis pipeline for quantifying TRF1 and γH2AX signals and their spatial colocalization within nuclei, using confocal microscopy data and deep learning-based segmentation.

## Overview

The pipeline includes:

- **Segmentation**: 3D segmentation of DAPI, TRF1, and γH2AX using custom-trained Cellpose models.
- **Quantification**:
  - Number and intensity of TRF1 and γH2AX foci per nucleus.
  - Colocalization events based on centroid overlap within nuclear volumes.
- **Statistical Analysis**: Automated t-tests comparing wild-type and knockout conditions.

## Requirements

- Fiji (ImageJ) with Groovy scripting
- MCIB3D library
- Cellpose (3D model training and inference)
- Java 8+

## Usage

1. Segment DAPI, TRF1, and γH2AX channels using your trained Cellpose models.
2. Run `3DTrf1H2axColocperClonNucleus.groovy` in Fiji.
3. Results are saved in `output/csv/`, including per-nucleus quantification and statistical comparisons.

## Output Metrics

- Number of TRF1 and γH2AX foci per nucleus
- Mean, sum, and standard deviation of intensity
- Number of colocalized foci (TRF1–γH2AX)
- p-values from t-tests comparing WT vs KO



