# Spectral hardening tests from the appendix from "Quantifying Lyman-α emissions from reionization fronts" [(Wilson et al. 2024)](https://arxiv.org/abs/2406.14622

Bayu Wilson

Note that this is not documented well. It is mainly for Bayu's personal use and/or if someone is especially curious how the spectral hardening plots were created in my paper.

| **Script**    | **Description**                                                                                         |
|------------------|---------------------------------------------------------------------------------------------------------|
| `get_neutral_island_location.py`     | Reads ndot catalog. Then makes 3 plots. 1) 3D_scatterplot.png. Self-expanatory. 2) GammaHI heatmap using gas file with source positions overlayed. map_GammaHI_with_source_positions_.png. 3) A sightline poked through the most neutral region. A vertical line is placed to find the most neutral location along the skewer. See sightline.png |
| power_law_tests.py | Testing validity of power law spectrum for hardened spectra. See figure 5 in the appendix of paper I. |
| hardening_tests/plot_GammaHI_alpha.py| Shows degree of spectral hardening by residual HI absorption in the IGM at z=5.7. See figure 6 or paper I. | 
