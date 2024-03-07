$                                                                               
$                                                                               
$                                                                               
$   3 DIMENSIONAL PLOTS                                                        
$   ===================                                                         
$                                                                               
$    Ensure all titles begin with a reference from the following list :-        
$      3C         3R            3L            3P          3I                    
$    Clouds    PowerLoss   Line Radiation    PLRPs    Ionisation                
$                                                                               
$ X range       : (An over estimate is allowed, eg -99.0 to 99.0)               
$ Y range       : (An over estimate is allowed, eg -99.0 to 99.0)               
$ npts option   : 0:no plot >0 (eg 100): no. of gridpoints in each dim.         
$ State option  : Ionisation state to be plotted (-1:prim neuts,                
$                 0:neutrals, 1 to nizs: normal, nizs+1: total line)            
$ P plane option: 99:integrate over P  -4 to 4: Which P plane required          
$ Yfold option  : 0:off   1:fold -Y results into +Y region                      
$ Smooth        : Apply a flat 3 X 3 convolution to the surface
$                                                                               
$                                                                               
$ Ref   Title             X range      Y range   npts State Plane Yfold Smoth        
$ ---------------------------------------------------------------------------         
'3C     Clouds      3D' -0.1   0.1    -6.0   6.0    0    1    99     0    0
'3C     Clouds      3D' -0.1   0.1    -6.0   6.0    0    2    99     0    0      
'3C     Clouds      3D' -0.1   0.1    -0.4   0.4    0    1    99     0    0      
'3C     Clouds      3D' -0.1   0.1    -0.4   0.4    0    2    99     0    0      
'3R   Power Loss    3D' -0.02  0.08   -0.05  0.15   0   99    99     0    0      
'3L  Line Radiation 3D' -0.02  0.08   -0.1   0.1    0   99    99     0    0      
'3P     PLRPs       3D' -0.02  0.08   -0.1   0.1    0   99    99     0    0      
'3I   Ionisation    3D' -0.02  0.08   -0.05  0.15   0    3    99     0    0      
$----------------------------------------------------------------------------
