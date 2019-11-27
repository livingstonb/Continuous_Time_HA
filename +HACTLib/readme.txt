

  ----------------- Conventions --------------------------------
  Frequently, we repeat the use of certain variables across this
  toolbox without redefinition. They are as follows:

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
    