=====
PET-FBA
=====

PET-FBA (pH-, Enzyme allocation-, and Thermodynamics-constrained Flux Balance Analysis) is a computational framework that integrates catalytic enzyme protein availability and reaction themodynamic feasibility, conveniently incorporating the pH effect on them. It provides a more comprehensive mechanistically grounded approach to modeling metabolic fluxes.

This repository host the Python implementation, released as etfba, which allows applying these constraints either individually or jointly within a model. This flexibility enables users to solve various optimization problems, including:

- FBA: Traditional flux balance analysis.
- EFBA: FBA with enzyme protein allocation constraints.
- TFBA: FBA with thermodynamic constraints.
- ETFBA: FBA with both enzyme protein allocation and thermodynamic constraints.

Variability analysis is supported, allowing users to evaluate the potential ranges of metabolic fluxes, enzyme protein costs, and reaction Gibbs energy changes while ensuring that the objective function remains within a specified level of optimality. This analysis includes:

- FVA: Flux variability analysis.
- EFVA: Enzyme-constrained FVA.
- TFVA: Thermodynamic FVA.
- EVA: Enzyme variability analysis.
- TEVA: Thermodynamic EVA.
- TVA: Thermodynamic variability analysis.
- ETVA: Enzyme-constrained TVA.

For further details, please refer to our `documentation <https://etfba.readthedocs.io/en/latest/index.html>`__. For example of applying pH-dependent flux analysis, see the `script <https://github.com/Chaowu88/etfba/tree/main/scripts>`__.

Installation
============

The package has been tested with Python versions 3.8, 3.9, 3.10 and 3.11. It can be installed using *pip* from PyPI:

.. code-block:: python

  python -m pip install --upgrade pip
  pip install etfba

Alternatively, you can install it from source (assuming `git <https://git-scm.com/>`__ is installed):

.. code-block:: python

  git clone https://github.com/Chaowu88/etfba.git /path/to/etfba
  pip install /path/to/etfba

Note: It is recommended to install within a `virtual environment <https://docs.python.org/3.8/tutorial/venv.html>`__.

Solver installation
===================

The package uses the modeling language `Pyomo <https://www.pyomo.org/>`__ to formulate linear programming (LP) and mixed integer linear programming (MILP) problems. You can install the freely available solver GLPK via conda:

.. code-block:: python

  conda install -c conda-forge glpk

For larger and more complex models, such as genome-scale models with thermodynamic constraints, it is highly recommended to use the commercial optimizer `Gurobi <https://www.gurobi.com/>`__ with its Python interface installed:

.. code-block:: python

  conda install -c gurobi gurobi

Example Usage
=============

You can build a model from scratch or convert it from a `COBRA <https://cobrapy.readthedocs.io/en/latest/io.html>`__ model (refer to `here <https://etfba.readthedocs.io/en/latest/building_model.html>`__ for more details). Below is an example of estimating the flux distribution constrained by enzyme protein allocation and thermodynamics:

.. code-block:: python

  from etfba import Model

  model = Model.load('/path/to/model.bin')
  
  res = model.optimize(
      'etfba',
      objective=objective,      # typically the growth rate
      flux_bound=flux_bound,    # bounds for metabolic fluxes 
      conc_bound=conc_bound,    # bounds for metabolite concentrations
      preset_flux=preset_flux,  # preset values for specific metabolic fluxes
      preset_conc=preset_conc,  # preset values for specific metabolite concentrations
      ex_thermo_cons=ex_rxns,   # reactions excluded from thermodynamic constraint
      inc_enz_cons=eff_rxns,    # reactions included in enzyme protein constraint
      enz_prot_lb=enz_ub,       # upper bound on enzyme protein allocation
      parsimonious=True         # to obtain parsimonious flux distributions
  ).solve(solver='gurobi')

  opt_growth_rate = res.opt_objective
  opt_metabolic_fluxes = res.opt_fluxes

To estimate the variability of fluxes:

.. code-block:: python

  res = model.evaluate_variability(
      'etfva',
      objective=objective,
      obj_value=obj_value,   # optimal objective value obtained by "optimize"
      gamma=gamma,           # fraction of the optimum objective to achieve
      flux_bound=flux_bound,
      conc_bound=conc_bound,
      preset_flux=preset_flux,
      preset_conc=preset_conc,
      ex_thermo_cons=ex_rxns,
      inc_enz_cons=eff_rxns,
      enz_prot_lb=enz_ub
  ).solve(solver='gurobi', n_jobs=100)

  metabolic_flux_ranges = res.flux_ranges

For more detailed information, please refer to the complete `documentation <https://etfba.readthedocs.io/en/latest/index.html>`__.

Citation
========

`Chao Wu, Jeffrey N. Law, Onyeka Onyenemezu, Jetendra K. Roy, Peter C. St. John, Robert L. Jernigan, Yannick J. Bomble, Laura Jarboe. PET-FBA: A lightweight enzyme allocation and thermodynamics-constrained flux analysis approach to explore Escherichia coli metabolic adaptation to intracellular acidification. Metabolic Engineering. 2025. <https://doi.org/10.1016/j.ymben.2025.12.003>`__

