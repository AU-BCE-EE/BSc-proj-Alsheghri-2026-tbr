# BSc-proj-Alsheghri-2026-tbr

Osman Alsheghri's trickle-bed reactor BSc project.

The repository contains a numerical model for CO2 absorption in NaOH solution. The model tracks pH evolution in both the packed column and reservoir.

- Main model file: `mod_co2_main.py` 
- Speciation model file: `speciation_model.py`
- Enhancmenet factor file: `enhancement_mod.py`


# Maintainer
Osman Alsheghri and Sasha Hafner (contact information here: <https://au.dk/sasha.hafner@bce.au.dk>).

# Contents

-   `laboratory`: raw experimental data and Python scripts for model validation.
-   `mods`: python package containing the main model code, speciation model code, enhancement factor and a simplified version of the main model.
-   `report_results`: Python scripts used to generate report plots, and other report results
-   `shared`: subdirectory for sharing example code and more. Contents, names, and organization flexible.
-   `tests`: model tests used during model development.
-   `verification`: verification and debugging files.

# Running the model

To run the model for a given set of input, create a python script in either the repository root or a folder.(for example: example.py)

First import the model in the script

```{python}
import mods.mod_co2_main as md
```

Then call the main model function (tfmod)

```{python}
result = md.tfmod(inputs)
# extract the different results from the output dictionary 
# e.g.
gas_concentration = results['gas_conc']
pH = results['pH_profile]
```

Then run the script from the repository root in a terminal:
``` bash
python -m example
```