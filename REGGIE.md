
## Flexi Regressioncheck 

## How to execute

## Which examples?

| **No.** |            **Example**            | **When** |     **CMAKE-CONFIG**    |      **Feature**               |         **Execution**                       |           **Comparing**          |
|:-------:|:---------------------------------:|:--------:|:-----------------------:|:------------------------------:|:-------------------------------------------:|:--------------------------------:|
|    1    |             run_basic             |  checkin |       maxwell,RK4       |      DG-Operator               |        nProcs=1,2,5,8                       |              L2,Linf             |
|    2    |    feature_poisson_powerdensity   |  nightly | poisson, Crank-Nicolson | Implicit, CalcTimeAvg          |  DoRefMapping=T/F, nProcs=2                 |       Final TimeAvg, h5diff      |
|    3    |   feature_tracking_DSMC_ANSA_box  |  nightly |      maxwell, DSMC      |        Tracking                | DoRefMapping=T,F, nProcs=1,2                | PartInt, PartPos in Bounding Box |
|    4    |    feature_tracking_DSMC_curved   |  nightly |       maxwell,DSMC      |        Tracking                |  DoRefMapping=T, nProcs=1,2                 |              PartInt             |
|    5    |    feature_tracking_DSMC_mortar   |  nightly |       maxwell,DSMC      |  Tracking with Mortar          | DoRefMapping=T,F, nProcs=1,2                | PartInt, PartPos in Bounding Box |
|    6    |   feature_tracking_DSMC_periodic  |  nightly |       maxwell,DSMC      |        Tracking                | DoRefMapping=T,F, nProcs=1,2                | PartInt, PartPos in Bounding Box |
|    7    | feature_tracking_DSMC_sphere_soft |  nightly |       maxwell,DSMC      |        Tracking                |                                             |                                  |
|    8    | feature_poisson_powerdensity      |  nightly |       poisson,CN        |      CalcTimeAvg               | DoRefMapping=1,2, nProcs=2, CN implicit     |  TimeAvg                         |
|    9    |   feature_emission_gyrotron       |  nightly |       maxwell,RK4       | Part-Inflow,TimeDep            | N=1,3,6,9,10, nProcs=1,2,10,25, gyro-circle |  LineIntegration of nPartIn      |
|   10    |      feature_TWT_recordpoints     |  nightly |       maxwell,RK4       | RPs, ExactFlux                 | nProcs=1,4, RPs, interior TE-Inflow         |  RP_State, RP_Daata              |
|   11    |    feature_PIC_HDG_plasma_wave    |  nightly |       poisson,RK4,CN    | Poisson-PIC,Shape-Function-1D  | nProcs=2, Imex for CN                       |  W_el LineIntegration over 2Per  |
|   12    |  feature_PIC_maxwell_plasma_wave  |  weekly  | maxwell,RK4,ImplicitO4  | Maxwell-PIC,SF1D, FastPeriodic | nProcs=2, IMEX for ImplicitO4               |  W_el LineIntegration over 2Per  |
