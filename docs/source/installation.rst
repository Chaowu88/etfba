Installation
============

Using PIP
---------

ETFBA is compatible with Python versions 3.8 through 3.11 and can be easily installed using *pip* from PyPI. To get started, first upgrade pip using the following command:

.. code-block:: python

  python -m pip install --upgrade pip

Next, install ETFBA with the command:

.. code-block:: python

  pip install etfba  

Alternatively, you can install FreeFlux from the source code by cloning the GitHub repository using the following command (assuming you have `git <https://git-scm.com/>`__ installed):

.. code-block:: python

  git clone https://github.com/Chaowu88/etfba.git /path/to/etfba

Then, install ETFBA using *pip*:

.. code-block:: python

  pip install /path/to/etfba
  
.. Note::
  It's recommended to install ETFBA within a virtual environment to avoid conflicts with other Python packages. Refer to these `instructions <https://docs.python.org/3.8/tutorial/venv.html>`_ on creating virtual environments or `here <https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#creating-an-environment-with-commands>`_ for Conda environments.

Solver Installation
-------------------
 
ETFBA uses the modeling language `Pyomo <http://www.pyomo.org/>`__ to formulate linear programming (LP) and mixed integer linear programming (MILP) problems. For small-sized models, the freely available solver glpk is capable of handling the workload, which can be installed by:

.. code-block:: python
  
  conda install -c conda-forge glpk  

For larger models, such as genome scale models, it is highly recommended to use the commercial optimizer `Gurobi <https://www.gurobi.com/>`_ and install the Python support:

.. code-block:: python

  conda install -c gurobi gurobi
  

