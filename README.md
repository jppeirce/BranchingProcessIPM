# IPM with Branching Processes and Allee Effects

Code accompanying the paper:

> *Combining branching processes and Allee effects into an integral projection model to assess invasion risk*

**Authors:** James P. Peirce[1,2],
Cynthia S. Tam[3],
Joseph J. Parkos III[4],
Grace L. Loppnow[5],
Mark W. Fritts[6],
Alison A. Coulter[7],
Richard A. Erickson[8]

[1] Mathematics & Statistics Department, University of Wisconsin - La Crosse, La Crosse, WI, USA

[2] River Studies Center, University of Wisconsin - La Crosse, La Crosse, WI, USA

[3] U.S.~Geological Survey, Ecosystems Mission Area, Reston, VA, USA

[4] Kaskaskia Biological Station, Illinois Natural History
  Survey, Prairie Research Institute, University of Illinois at
  Urbana-Champaign, Sullivan, IL, USA

[5] Minnesota Department of Natural Resources, St. Paul, MN, USA

[6] U.S. Fish and Wildlife Service, La Crosse Fish and Wildlife
  Conservation Office, Onalaska, WI, USA

[7] Department of Natural Resource Management, South Dakota
  State, Brookings, SD, USA

[8] U.S. Geological Survey, Upper Midwest Environmental Sciences
  Center, La Crosse, WI, USA

---

## Abstract

Invasive species present a major threat to both ecological systems and economy. Their successful establishment is governed, in part, by demographic stochasticity where random individual-level variations can significantly influence invasion dynamics. A critical factor in this process is the Allee effect, which describes reduced population growth at low organism densities, potentially serving as a barrier to species persistence. Additionally, size-dependent growth and survival rates further shape population dynamics, leading to variable outcomes based on the distribution of individual sizes within the population. To address these complexities, we developed a dynamic population-level model that utilizes an Integral Projection Model (IPM) for the expected distribution of sizes. This framework incorporates an Allee effect for survival and recruitment dynamics, along with a branching process to account for demographic stochasticity. We applied this framework to Silver carp (*Hypophthalmichthys molitrix*), a species that has established populations in the Upper Mississippi River system and poses a significant risk of expanding into the Great Lakes. Our analysis compared the temporal dynamics and population-size distributions of large adult fish and smaller subadult fish under hypothetical introduction scenarios. In both simulations, the initial population size of Silver carp relative to the Allee threshold was a critical factor in determining the probability of fish establishment. These findings highlight the model’s utility for invasion assessments, offering valuable insights for predicting and managing the spread of invasive species.

## Model overview

This model tracks an invasive fish population over time by combining three modeling frameworks:

**Integral Projection Model (IPM):** The population is represented as a continuous distribution over body length. Each time step, individuals survive, grow (following von Bertalanffy growth), and reproduce according to length-dependent vital rates.

**Branching process:** Rather than advancing the population deterministically via matrix multiplication, each individual fish is treated as an independent random variable at each time step. Each fish draws one of four outcomes:
- Dies (probability depends on survival and reproductive outcome)
- Survives and grows
- Spawns *k* young-of-year (YOY) recruits and dies
- Survives, grows, and spawns *k* YOY recruits

This individual-level stochasticity is especially important at low abundances, where demographic noise can drive extinction even when average conditions favor growth.

**Allee effect:** Reproduction is modulated by a biomass-based Allee threshold *C*. The proportion of mature fish that reproduces (*A*) depends on how much biomass is "available" to fill the gap between current survivor biomass and the carrying-capacity trajectory. When total biomass is below *C*, reproduction is suppressed, creating a strong Allee effect that can cause small founding populations to fail.

### Main functions

**`ipm_bp_allee(...)`** — runs a single stochastic trajectory forward in time. Returns a `(n_len × n_time+1)` array of population density across the length mesh.

**`stoch_wrapper(...)`** — runs `ipm_bp_allee` for `n_iter` replicate trajectories and returns a long-format `DataFrame` with total population size per year per iteration, suitable for plotting extinction probability or stochastic envelopes.

### Parameters

| Parameter | Description |
|---|---|
| `gamma` | Intrinsic growth rate (biomass-based) |
| `C` | Allee threshold biomass |
| `K` | Carrying capacity biomass |
| `k` | Number of YOY recruits produced per reproductive event |
| `min_len`, `max_len`, `n_len` | Length mesh bounds and resolution |
| `n_time` | Number of annual time steps |
| `yoy_mean`, `yoy_sd` | Mean and SD of YOY length at recruitment |
| `init_pop` | Dictionary of initial fish counts and lengths |

### Key vital rate submodels

| Submodel | Description |
|---|---|
| `length_weight` | Log-10 length–weight relationship |
| `prob_survival` | Annual survival as a function of body weight |
| `node` (von Bertalanffy) | Stochastic growth kernel — projects a fish from length *z* to length *z'* |
| `mature` | Logistic maturity-at-length |

---

## Recreating the environment

Python 3.9 is required (the `.venv` was created with the system Python 3.9).

```bash
# From the Code/ directory
python3.9 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

The notebook also uses `plotnine` for figures, which is not in `requirements.txt` because it is only needed for visualization:

```bash
pip install plotnine
```

To launch the notebook:

```bash
jupyter notebook IPM_Allee.ipynb
```

---

## File overview

| File | Description |
|---|---|
| `IPMbpA.py` | Model implementation (all classes and functions) |
| `IPM_Allee.ipynb` | Notebook reproducing figures from the paper |
| `requirements.txt` | Pinned package versions for the core environment |
| `figures/` | Output figures |
