  ----------------- Summary ------------------------------------
  HACTLib is a library of functions and classes designed to
  allow for the creation of a heterogenous agent continuous
  time model.

  ----------------- Basic Components ---------------------------
  - The Params class is used to store model parameters.
  - The Income class is used to store the income process.
  - The Grid class is used to store the asset grids.
  - The Preferences class is used to store utility functions.
  - The Model class initializes and solves the model.
  - Additional classes are provided to solve the HJB and KFE
    and to compute or simulate MPCs.

  ----------------- Conventions --------------------------------
  Frequently, we repeat the use of certain variables across this
  library without redefinition. They are as follows:

    z : An extra state variable used for preference and returns
    heterogeneity.

    nb : Number of points on the liquid asset grid for the HJB.

    na : Number of points on the illiquid asset grid.

    nz : Number of points on the z-grid.

    n_states : Total number of states, i.e. nb*na*nz*ny, or for
    the KFE, this is nb_KFE*na_KFE*nz*ny.

    ny : Number of points on the income grid.

    nb_KFE : Number of points on the liquid asset grid for the KFE.

    na_KFE : Number of points on the illiquid asset grid for the KFE.

    grids_HJB : A Grid object for the HJB.

    grids_KFE : A Grid object for the KFE.

    ytrans : The square income transition matrix, of shape (ny, ny).
    