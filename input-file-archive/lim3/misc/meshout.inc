$----------------------------------------------------------------------
$                                                                               
$                                                                               
$   MESH / CONTOUR PLOTS                                                        
$   ====================                                                        
$
$    These plots are only valid for LIM runs using 3D. Integrations are 
$    only over the defined 3D space.
$                                                                               
$    Ensure all titles begin with a reference from the following list :-        
$      MC         MR            ML            MP          MI                    
$    Clouds    PowerLoss   Line Radiation    PLRPs    Ionisation                
$                                                                               
$ I             : Plot Option 0-off  1-contour  2-mesh
$ P             : Axis / Plane of Integration (0-X  1-Y  2-P)
$                 ( 0 - Y vs P,  1 - X vs P,  2 - X vs Y)
$ X range       : (An over estimate is allowed, eg -99.0 to 99.0)               
$ Y range       : (An over estimate is allowed, eg -99.0 to 99.0)               
$ P range       : (An over estimate is allowed, eg -99.0 to 99.0)               
$ N             : Number of grid points for first dimension of mesh
$ M             : Number of grid points for second dimension of mesh
$ Z             : Ionisation state to be plotted (-1:prim neuts,                
$                 0:neutrals, 1 to nizs: normal, nizs+1: total line)            
$ F             : Yfold option   0:off   1:fold -Y results into +Y region    
$                                                                               
$                                                                               
$ Ref Title        I  P    X range       Y range       P range     N   M   Z  F 
$ -----------------------------------------------------------------------------
'MC Clouds     3D' 0  0  -99.0  99.0   -10.0  10.0   -0.50  0.50   40  40  1  0 
'MC Clouds     3D' 0  0  -99.0  99.0   -10.0  10.0   -0.10  0.10   50  50  2  1 
'MC Clouds     3D' 0  0  -0.10  0.10   -0.40  0.40   -0.40  0.40   40  40  1  0 
'MR Power Loss 3D' 0  0  -0.10  0.10   -0.40  0.40   -0.40  0.40   40  40  1  0 
'ML Line Rad.  3D' 0  0  -0.10  0.10   -0.40  0.40   -0.40  0.40   40  40  1  0 
'MP PLRPs      3D' 0  0  -0.10  0.10   -0.40  0.40   -0.40  0.40   40  40  1  0 
'MI Ionisation 3D' 0  0  -0.10  0.10   -0.40  0.40   -0.40  0.40   40  40  1  0 
$ -----------------------------------------------------------------------------
