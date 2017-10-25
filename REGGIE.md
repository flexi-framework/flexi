## Flexi Regressioncheck 

## How to execute

## Which examples?

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
