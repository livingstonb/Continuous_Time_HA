# Continuous Time HA

Code for a continuous time heterogenous agent model.

## Repository structure

This repository consists of two main parts:

* A library with classes and functions to construct a heterogenous agent continuous time model (HACTLib).

* Code that implements the library. HACTLib components can be used independently, so the user may implement the library with his/her own code.

## SuiteSparse
We use SuiteSparse, a package written by Timothy A. Davis.

-- <cite> Direct Methods for Sparse Linear Systems, Timothy A. Davis, SIAM, Philadelphia, Sept. 2006. </cite>

# HACTLib

## Summary
HACTLib is a library of functions and classes designed to allow for the creation of a heterogenous agent continuous time model. The *Model* class utilizes most of the library objects and the source code for *Model* can be used as a guideline for how the different model objects can be used together. Details omitted in this document can for most classes can be found in comments written in the source code, which can often be viewed directly in MATLAB with the `help` command. The main components of the library are the following:

* The Params class, used to store model parameters.

* The Income class, used to store the income process.

* The Grid class, used to store the asset grids.

* The Preferences class, used to store utility functions.

* Additional classes provided to solve the HJB and/or KFE and to compute or simulate MPCs.

## Conventions
Frequently, we repeat the use of certain variables across this library without redefinition. They are as follows:

* z : An extra state variable used for preference and returns heterogeneity.

* nb : Number of points on the liquid asset grid for the HJB.

* na : Number of points on the illiquid asset grid.

* nz : Number of points on the z-grid.

* n_states : Total number of states, i.e. nb\*na\*nz\*ny, or for the KFE, this is nb\_KFE\*na\_KFE\*nz\*ny.

* ny : Number of points on the income grid.

* nb\_KFE : Number of points on the liquid asset grid for the KFE.

* na\_KFE : Number of points on the illiquid asset grid for the KFE.

* grids\_HJB : A Grid object for the HJB.

* grids\_KFE : A Grid object for the KFE.

* ytrans : The square income transition matrix, of shape (ny, ny).

## Params

### Instantiation

The Params object stores the model parameters and contains various methods related to these parameters. Depending on use case, various parameters will be unused and can be ignored. Parameters take their default values from *ParamsDefaults.m* and can be overriden in either or both of the following two ways:

* Instantiate the Params object with a structure with fields matching the desired parameters to be overriden. E.g. to override the default number of points on the illiquid grid, one could pass the class constructor a structure defined by `p = struct('na', 125)`, when creating the Params object. Any number of parameters can be overriden this way.

* Pass name, value pairs to the class constructor, e.g. `new_params = Params('na', 125)`.

Note that some parameters set with the Params object are intended to be used to instantiate the properties of other objects, such as the Grid object. See the appropriate sections of this readme for more information on those parameters.

Note that some of the properties of the Params object are structures containing multiple parameters, which I use to pass these parameters to other objects, such as HJBSolver. These structures should not be set manually; instead, the user should override the appropriate parameters, and the Params object will update these structures accordingly. E.g. to override the default value for the delta passed to HJBSolver (the default value for which can be seen in *HJBSolver.m*), one should override the HJB\_delta property when instantiating the Params object.

### Parameters that *must* be set by the user

The user must set the directory parameters. These parameters are *out_dir*, *temp_dir*, and *direc*, which should be set to the paths of the output directory, a temp directory, and the path of the directory containing the +HACTLib folder.

### Frequency intepretation of parameters

Parameters related to the passage of time should be set to their equivalent quarterly values. E.g. if the user wants the probability of death within a given year to be 0.02, the *deathrate* parameter should be set to 0.005.

### Model features to be turned on or off

Some parameters are expected to take true/false values, which may activate or disactivate certain model features; some of these only apply if the *main.m* script is used. Some of these parameters are as follows:

* *OneAsset*

	When true, the illiquid asset is effectively turned off and the standard one-asset model is recovered.

* *Bequests*

	When true, household liquid and illiquid assets are maintained upon death. When false, both assets are set to zero upon death.

* *perfectannuities*

	When true, the *Bequests* parameter is set to false automatically and returns on both the liquid and illiquid asset, where they are used, are augmented by the death rate.

* *SDU*

	When true, stochastic-differential utility is used. The *invies* and *riskaver* parameters should both be set by the user if SDU is used.

* *DealWithSpecialCase*

	This feature has gone out of date and should be left as false at the moment.

* *NoRisk*

	When true, the model is solved both with and without income risk. In some cases, turning off income risk leads convergence to fail, so set this to false if it's causing problems and you don't need the model solved for the case of no income risk.

* *ComputeMPCS* and *ComputeMPCS_illiquid*

	When true, MPCs out of liquid and/or illiquid assets are computed.

* *SimulateMPCS*

	When true, MPCs out of the liquid asset are simulated. Note that some model features may not be supported for MPC simulations; simulation code was written primarily to test the results from direct computation.

* *ComputeMPCS_news* and *SimulateMPCS_news*

	When true, compute or simulate MPCs out of news of a future transfer.

### Endogenous labor supply

This feature should be more-or-less complete, but it has not been tested.

### Discount rate heterogeneity

Discount rate heterogeneity is achieved by selecting an array of constants to be added to the model parameter *rho*. The grid actually used in the code will be the parameter rhos = *rho* + *rho_grid*. The parameter rhos should not be overriden, it will be automatically generated based on the selected values for *rho* and *rho_grid*. The user can override the *rho_grid* parameter to be a vector such as [-0.01, 0, 0.01]. Then the grid of discount factors will be [*rho* - 0.01, *rho*, *rho* + 0.01]. Note that when iterating over the discount factor, each element of the grid will therefore shift by the change in rho at each iteration, but the spacing will remain constaint.

### Risk aversion heterogeneity

A grid can be used for the risk aversion coefficient by setting the *riskaver* parameter equal to a row vector, e.g. [0.5, 1, 1.5]. Unless using stochastic-differential utility, the *invies* parameter will be automatically set equal to the *riskaver* parameter, and should be ignored by the user.


## Grid

The *Grid* object stores grid parameters, constructs the asset grid parameters, and contains the asset grids as properties.

### Instantiation

The best way to create the *Grid* object, is to pass the constructor a *Params* instance with grid parameters overriden as desired. The constructor accepts two more arguments: The number of income grid points, and the grid type. The grid type should be passed as either `'HJB'` or `'KFE'`, and this determines which grid parameters are used from the *Params* object, either the HJB parameters or the KFE parameters, in case the user wants to use different grids to solve the HJB and KFE.

### Borrowing constraints

The illiquid asset grid is always constructed with a hard constraint at zero illiquid assets. The liquid asset grid can contain negative values, in which case the number of points below zero will be determined by the number of total grid points chosen less the number of positive grid points chosen, both of which can be overridden by the user. The number of negative grid points should not be set manually. Curved grids can be used, in which case more points can be concentrated near the constraints. If borrowing is allowed in the liquid asset, the borrowing constraint is set by the parameter *bmin* and the soft constraint (presumably zero) is set by the parameter *b_soft_constraint*.

### Methods to create the grid variables

Once the *Grid* object is constructed, the grids will be constructed by calling the *autoconstruct* method, e.g.

`new_grid = Grid(params, ny, 'KFE').auto_construct()`

If a custom grid is desired, the *autoconstruct* method must be skipped and grids can be constructed by passing your grid vectors to the *create_agrids* and *create_bgrids* methods and then calling the *generate_variables* method.

Grid spacings used for finite-difference approximations as well as variables used for trapezoidal integration are created within the above methods.

### Grid curvature

The *Grid* object allows for curved grids. For a grid between zero and a positive value, two linearly-spaced grids in the interval [0,1] are constructed, and each is raised to some power to induce curvature. A linear combination of these vectors is then used as the grid, after rescaling to the desired grid min and max. For the first term I use a very low weight, e.g. 0.01, and very low curvature, e.g. a power of around 1.1. This first term prevents too many grid points from clustering right next to the constraint when high curvature is used on the second term, but an ordinary power-spaced grid can be used by setting the weight of the first term to zero. Curvature parameters are set in the *Params* object, and end with *\_gcurv*, and the exponent used in grid creation will be the inverse of these values. The weight on the first term is set by parameters ending with *\_term1\_weight*, and the weight on the second term is always one.

## AdjustmentCost

A class used for computing adjustment costs, with related methods.
The model can accomodate a single, liquid asset, or a liquid asset in combination with an illiquid asset. In the latter case, the household pays a cost for transferring funds in and out of the illiquid account.

The form of the adjustment cost function can be viewed in the *AdjustmentCost* class. There are two functional forms in the code, which are equivalent, but we use only the second.

## HJBSolver, KFESolver, MPCs, and MPCsNews

Each of these classes have their own set of options, which can be set manually or by passing a *Params* object into the class constructor. These classes have methods that solve various parts of the model and return the result(s).

## Statistics

Used for computing moments and other statistics associated with the stationary distribution. Most properties of this class are structures containing three fields:

1. *value*

	The value taken by the given statistic.

2. *label*

	The description of the given statistic.

3. *indicator*

	An indicator of which asset the statistic applies to. Takes the value of zero if the statistic is relevant for both one- and two-asset models, one if the statistic is only relevant for the one-asset model, and two if the statistic is only relevant for the two-asset model. This is used by the StatsTable object which may selectively present statistics based on the type of model used.

## StatsTable and BaseTable

The *StatsTable* class, which inherits from *BaseTable*, is used to produce results tables based on values provided by an instance of the *Statistics* class.

## Calibrator

Used to calibrate the model to certain targets set by the user by calling a Matlab solver to iterate over the values of one or more parameters.

## Rate of return risk

Rate of return risk can be enabled by setting the *sigma\_r* parameter to a positive value. Set *retrisk\_KFE* to true to include returns risk in the KFE, or false to exclude it from the KFE, but still include it in the HJB. If the two-asset model is used, returns risk is applied to the illiquid asset only. Otherwise, returns risk is applied to the liquid asset.