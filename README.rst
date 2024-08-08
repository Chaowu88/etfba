ETFBA
=======================================

ETFBA is a Python package designed for performing **e**\ nzyme protein allocation and **t**\ hermodynamics constraint-based **f**\ lux **b**\ alance **a**\ nalysis. It extends traditional flux balance analysis (FBA) by incorporating constraints based on the availability of catalytic enzyme proteins and the second law of thermodynamics, offering a more comprehensive approach to modeling metabolic fluxes.

ETFBA allows for the application of enzyme protein allocation and thermodynamic constraints either individually or jointly within the model. This flexibility enables users to solve various optimization problems, including:

- FBA: Traditional flux balance analysis.
- EFBA: FBA with enzyme protein allocation constraints.
- TFBA: FBA with thermodynamic constraints.
- ETFBA: FBA with both enzyme protein allocation and thermodynamic constraints.

Variability analysis allows you to evaluate the potential ranges of metabolic fluxes, enzyme protein costs, and reaction Gibbs energy changes while ensuring that the objective function remains within a specified level of optimality. This analysis includes:

- FVA: Flux variability analysis.
- EFVA: Enzyme-constrained FVA.
- TFVA: Thermodynamic FVA.
- EVA: Enzyme variability analysis.
- TEVA: Thermodynamic EVA.
- TVA: Thermodynamic variability analysis.
- ETVA: Enzyme-constrained TVA.

For further details, refer to our `documentation <https://etfba.readthedocs.io/en/latest/index.html>`__.

Installation
============

ETFBA has been tested with Python versions 3.8, 3.9, 3.10 and 3.11. It can be installed using *pip* from PyPI:

.. code-block:: python

  python -m pip install --upgrade pip
  pip install etfba

Alternatively, you can install it from source (assuming `git <https://git-scm.com/>`__ is installed):

.. code-block:: python

  git clone https://github.com/Chaowu88/etfba.git /path/to/etfba
  pip install /path/to/etfba

Note: It is recommended to install ETFBA within a `virtual environment <https://docs.python.org/3.8/tutorial/venv.html>`__.

Solver installation
===================

ETFBA uses the modeling language `Pyomo <https://www.pyomo.org/>`__ to formulate linear programming (LP) and mixed integer linear programming (MILP) problems. You can install the freely available solver GLPK via conda:

.. code-block:: python

  conda install -c conda-forge glpk

For larger models, such as genome scale models, it is highly recommended to use the commercial optimizer `Gurobi <https://www.gurobi.com/>`__ and install the Python support:

.. code-block:: python

  conda install -c gurobi gurobi

Example Usage
=============

You can build an ETFBA model from scratch or convert it from a `COBRA <https://cobrapy.readthedocs.io/en/latest/io.html>`__ model (refer to `here <https://etfba.readthedocs.io/en/latest/building_model.html>`__ for more details). Below is an example of estimating the flux distribution constrained by enzyme protein allocation and thermodynamics:

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

For more detailed information, please refer to the complete `documentation <https://etfba.readthedocs.io/en/latest/index.html>`__.



