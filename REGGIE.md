## Flexi Regressioncheck 

## How to execute

python reggie.py --help

## Checks/Examples

| **No.** |            **Check**              | **When** |     **CMAKE-CONFIG**    |      **Examples**              |      **Feature**               |         **Execution**                       |           **Comparing**          |
|:-------:|:----------------------------------|:--------:|:-----------------------:|:------------------------------:|:------------------------------:|:-------------------------------------------:|:--------------------------------:|
|    1    | run_basic (flexi)                 | checkin  | default                 | freestream_2D                  |  DG-Operator                   |  MPI=1,2                                    | L2                               |
|         |                                   |          |                         | freestream_3D                  |  DG-Operator                   |  MPI=1,2                                    | L2                               |
|    2    | convtest (flexi)                  | nighlty  | FLEXI_2D=ON             | h_2D                           |  h-convergece                  |  single                                     | L2                               |
|         |                                   |          |                         | h_3D                           |  h-convergece                  |  single                                     | L2                               |
|         |                                   |          | FLEXI_FV=ON             | h_3D_FV                        |  h-convergece                  |  single                                     | L2                               |
|         |                                   |          |                         | h_3D_mortar                    |  h-convergece                  |  single                                     | L2                               |
|         |                                   |          | FLEXI_PARABOLIC=OFF     | h_3D_parabolic_off             |  h-convergece                  |  single                                     | L2                               |
|         |                                   |          | FLEXI_2D=ON             | p_2D                           |  p-convergece                  |  single                                     | L2                               |
|         |                                   |          |                         | p_3D                           |  p-convergece                  |  single                                     | L2                               |
|         |                                   |          |                         | p_3D_mortar                    |  p-convergece                  |  single                                     | L2                               |
|         |                                   |          | FLEXI_PARABOLIC=OFF     | p_3D_parabolic_off             |  p-convergece                  |  single                                     | L2                               |
|         |                                   |          |                         |                                |                                |                                             |                                  |
|         |                                   |          |                         |                                |                                |                                             |                                  |
|         |                                   |          |                         |                                |                                |                                             |                                  |


## Analyze routines

|**analyze**        | **options**                          | **values**                                            |
|:-----------------:|:-------------------------------------|:------------------------------------------------------|
|L2 error           | analyze\_L2                          | 1e-5                                                  |
|h-convergence test | analyze\_Convtest\_h\_cells          | 1,2,4,8                                               |
|                   | analyze\_Convtest\_h\_tolerance      | 0.3                                                   |
|                   | analyze\_Convtest\_h\_rate           | 1.0                                                   |
|p-convergence test | analyze\_Convtest\_p\_rate           | 0.6                                                   |
|h5diff             | h5diff\_file                         | single-particle\_State\_000.00000000000000000.h5      |
|                   | h5diff\_reference\_file              | single-particle\_State\_000.00000000000000000.h5\_ref |
|                   | h5diff\_data\_set                    | DG\_Solution                                          |
|                   | h5diff\_tolerance\_value             | 1.0e-2                                                |
|                   | h5diff\_tolerance\_type              | relative                                              |
|data file line     | compare\_data\_file\_name            | Database.csv                                          |
|                   | compare\_data\_file\_reference       | Database.csv\_ref                                     |
|                   | compare\_data\_file\_tolerance       | 6e-2                                                  |
|                   | compare\_data\_file\_tolerance\_type | relative                                              |
