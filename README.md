# AP3ESM

**AP3ESM** (AI-Powered, Performance-Portable Earth System Model) is a kilometer-scale coupled Earth system model. It integrates the GRIST atmosphere, the LICOM ocean, CICE sea ice, and the CESM coupler, and is built to run at global resolutions down to about 1 km on heterogeneous supercomputers.

## Overview

Kilometer-scale Earth system models have to integrate the atmosphere, ocean, land, and sea ice on tens to hundreds of billions of grid points, over periods from days to decades. AP3ESM does this inside the CESM 2.2 CPL7/MCT coupler: GRIST supplies the atmosphere, LICOM (the LASG/IAP climate-system ocean model) supplies the ocean, and CICE supplies sea ice. In the SC '25 experiments the fully coupled 1 km configuration reached **0.54 SYPD** on Sunway OceanLight. The standalone 1 km atmosphere reached 0.85 SYPD. A 3v2 forecast (about 3 km atmosphere, about 2 km ocean) reproduced Super Typhoon Doksuri (2023) and the associated extreme rainfall over China.

This repository is the model source.

## Model components

| Component | Code | Role in AP3ESM |
| --- | --- | --- |
| Coupler | `cime/`, CPL7 + MCT | Owns the clocks, component concurrency, and inter-component communication. GRIST and LICOM connect through the MCT init / run / finalize and import / export interfaces. |
| Atmosphere | `components/grist/` | GRIST (Global–Regional Integrated Forecast System). Dynamical core on an unstructured hexagonal mesh. Physics packages live in `src/atmosphere/physics` (`amipc_physics`, `amipw_physics`, `simple_physics`). On Sunway, MPI plus OpenMP target directives (SWGOMP) offload loops onto the computing processing elements. |
| Ocean | `components/licom/` | LICOM3 on a tripolar grid with 80 levels. Hot loops in this tree run on the Sunway Athread computing processing elements (`SW_ATHREAD`). Each ocean grid has a header `source/def-undef.h-<grid>`, copied to `def-undef.h` at build time. |
| Sea ice | `components/cice/` | The CESM CICE component. The published experiments use CICE4 physics, selected in CIME as `CICE%CICE4`. |
| Land | Stub land (`SLND`) on the coupler side | In the published configuration GRIST exchanges directly with the land surface, so the coupler does not own a land component. The fully coupled compset therefore uses `SLND`. |

The AI physics suite reported in the paper is separate from the sources in this tree. A one-dimensional convolutional network along the vertical column (five ResUnits, about 11 layers, about \(5\times10^5\) parameters) predicts tendencies from wind, temperature, specific humidity, and pressure. A 7-layer residual MLP diagnoses surface downward shortwave and longwave fluxes. It was trained on 5 km GRIST fields from 80 days, 20 from each season. Physics sources in this repository are the conventional GRIST packages listed above.

## Layout

```
ap3esm/
├── cime/                  CIME driver, scripts, machines, grids
├── cime_config/           compsets for the coupled model
├── components/
│   ├── grist/             atmosphere (GRIST dynamical core and physics)
│   ├── licom/             ocean (LICOM3; source/ is the model code)
│   ├── cice/              sea ice
│   ├── clm/  cam/  pop/   CESM 2.2 components kept for the CIME build
│   └── ...
├── manage_externals/      CESM external checkout tool
├── Externals.cfg          upstream CESM component pins
└── LICENSE.txt            CESM 2.2 and bundled-component licenses
```

`components/licom/source` is a git submodule (branch `canuto2010`). The other CESM components are pinned in `Externals.cfg` near the CESM 2.2 tags (CIME `cime5.8.32`, CTSM `release-cesm2.2.01`, CICE `cice5_20200430`) and are already present in this tree.

## License

The coupler and the components brought in with CESM 2.2 follow the UCAR/NCAR terms in [`LICENSE.txt`](LICENSE.txt), together with the third-party terms listed there (CIME, CISM, FATES, ESMF, and others). CAM and WRF code inside the GRIST physics packages follows the license of the original host model; GRIST modifications are identified in the source files. Please cite the SC '25 paper below when publishing results from this model, and retain the CESM copyright notice.

## Acknowledgements

AP3ESM builds on the kilometer-scale GRIST atmosphere, the kilometer-scale LICOM ocean model, and the CESM 2.2 coupler ported to Sunway. The performance experiments used Sunway OceanLight and ORISE. CESM is developed by NCAR and its partners.

## Reference

> Kai Xu, Maoxue Yu, Yuhu Chen, Jie Gao, et al. **Kilometer-Scale AI-Powered and Performance-Portable Earth System Model (AP3ESM) to Achieve Year-Scale Simulation Speed on Heterogeneous Supercomputers.** *Proceedings of the International Conference for High Performance Computing, Networking, Storage and Analysis (SC '25)*, 2025.  
> https://doi.org/10.1145/3712285.3771788

```bibtex
@inproceedings{xu2025ap3esm,
  author    = {Xu, Kai and Yu, Maoxue and Chen, Yuhu and Gao, Jie and others},
  title     = {Kilometer-Scale {AI}-Powered and Performance-Portable {Earth} System Model ({AP3ESM}) to Achieve Year-Scale Simulation Speed on Heterogeneous Supercomputers},
  booktitle = {Proceedings of the International Conference for High Performance Computing, Networking, Storage and Analysis},
  series    = {SC '25},
  year      = {2025},
  doi       = {10.1145/3712285.3771788},
  url       = {https://dl.acm.org/doi/10.1145/3712285.3771788}
}
```
