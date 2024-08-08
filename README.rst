ETFBA
=======================================

ETFBA is a Python package designed for performing **e**\ nzyme protein allocation and **t**\ hermodynamics constraint-based **f**\ lux **b**\ alance **a**\ nalysis. It extends traditional flux balance analysis (FBA) by incorporating constraints based on the availability of catalytic enzyme proteins and the second law of thermodynamics, offering a more comprehensive approach to modeling metabolic fluxes.

ETFBA allows for the application of enzyme protein allocation and thermodynamic constraints either individually or jointly within the model. This flexibility enables users to solve various optimization problems, including:

- FBA: Traditional flux balance analysis.
- EFBA: FBA with enzyme protein allocation constraints activated.
- TFBA: FBA with thermodynamic constraints activated.
- ETFBA: FBA with both enzyme protein allocation and thermodynamic constraints activated.

The constraints of enzyme protein allocation and thermodynamics can be used seperately or jointly in the model, thus yielding various optimization problems, such as: FBA (traditional flux balance analysis), EFBA (enzyme protein allocation constraint activated), TFBA (thermodynamic constraints activated) and ETFBA (both type of constraints activated). Variability analysis can further be performed to assess the possible range of metabolic fluxes (FVA, EFVA, TFVA and ETFVA), enzyme protein costs (EVA and TEVA) and reaction Gibbs energy change (TVA and ETVA) while still achieving a certain level of optimality in the objective function. For more information, refer to our `documentation <https://etfba.readthedocs.io/en/latest/index.html>`__.

Installation
============

ETFBA was tested in Python 3.8, 3.9, 3.10 and 3.11. It can be installed using *pip* from PyPI:

.. code-block:: python

  python -m pip install --upgrade pip
  pip install etfba

or from source (assuming you have `git <https://git-scm.com/>`__ installed):

.. code-block:: python

  git clone https://github.com/Chaowu88/etfba.git /path/to/etfba
  pip install /path/to/etfba

Installation within an `virtual environment <https://docs.python.org/3.8/tutorial/venv.html>`__ is recommendated.

Solver installation
===================

ETFBA uses the modeling language `Pyomo <https://www.pyomo.org/>`__ to formulate linear programming (LP) and mixed integer linear programming (MILP) problems. Freely available solver glpk can be installed by:

.. code-block:: python

  conda install -c conda-forge glpk

For larger models, such as genome scale models, it is highly recommended to use the commercial optimizer `Gurobi <https://www.gurobi.com/>`__ and install the Python support:

.. code-block:: python

  conda install -c gurobi gurobi

Example Usage
=============

The ETFBA model can be built from stratch or translated from a `COBRA <https://cobrapy.readthedocs.io/en/latest/io.html>`__ model as examplified `here <https://etfba.readthedocs.io/en/latest/building_model.html>`__. A typicall usage to estimate the flux distribution constrained by enzyme protein allocation and thermodynamics is as below:

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

And estimate the variability of fluxes:

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

For more information, please refer to the `documentation <https://etfba.readthedocs.io/en/latest/index.html>`__.



