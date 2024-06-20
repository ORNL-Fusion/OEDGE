$                                                                          
$                                                                          
$                                                                          
$   2 DIMENSIONAL PLOTS                                                    
$   ===================                                                    
$                                                                          
$ Fill in plotting range and option flags for the graphs required          
$ Graphs can be in any order and repeated if required                          
$ Ensure all titles begin with a reference from the following list :-      
$                                                                          
$ 2CY 2CX 2CP 2CT   2DX   2NX   2WY   2RY 2RX 2RP  2LY 2LX 2LP             
$     Clouds        Deps  NetE  Wall    PowerLoss    Line Rad                 
$                                                                          
$ 2PY 2PX 2PP 2PT   2IY 2IX 2IP   2ZY 2ZX   2TY 2TX   2YY 2YX              
$     PLRPs          Ionisation    Zeffs    Tempture  Y steps              
$                                                                          
$ Plotting range: (An over estimate is allowed, eg -99.0 to 99.0)          
$ Integration   : Range to sum over for plots. Overestimates are 
$               : allowed. Also, a range of 0.0 to 0.0 is taken to be 
$               : the full scale. Any Y integration greater than 2 * CL 
$               : is assumed to be the full scale only.
$ Plot option   : 0:no plot    1:plot required     2:normalised plot          
$ Smooth option : Z:smooth states Z upwards (eg 0:all 5:some 99:none)         
$ Maxiz option  : 0:plot neutrals only  >0:max state to plot                  
$                -1:plot primary neutrals only  99:plot all states            
$ P plane option: 99:integrate over P  -4 to 4: Which P plane required        
$ Yfold option  : 0:off   1:fold -Y results into +Y region, & vice versa      
$ Alliz option  : 0:off   1:include "total of all ionisations" line           
$ View option   : 0:off   1:Xview position   2:Calculate from Sv,Tv,etc       
$                                                                             
$                                                                             
$ Ref   Title    Plotting  Integration  Plot Smooth Maxiz Plane Yfld Allz Vu  
$ -------------------------------------------------------------------------- 
'2CY Clouds Y=' -0.25 0.25    0.0  0.0    0     1     8     99    0   0   0   
'2CY Clouds Y=' -2.5  2.5     0.0  0.0    0     1     8     99    0   0   0  
'2CY Clouds Y=' -99.0  99.0   0.0  0.0    2    99     6     99    0   0   0  
'2CY Clouds Y=' -1.0   1.0    0.0  0.0    0    99     6     99    0   0   0
'2CY Clouds Y=' -99.0  99.0   0.0  0.0    0    99     6     -2    0   0   0  
'2CY Clouds Y=' -99.0  99.0   0.0  0.0    0    99     6     -1    0   0   0
'2CY Clouds Y=' -99.0  99.0   0.0  0.0    0    99     6      0    0   1   0
'2CY Clouds Y=' -99.0  99.0   0.0  0.0    0    99     6      1    0   0   0
'2CY Clouds Y=' -99.0  99.0   0.0  0.0    0    99     6      2    0   0   0
'2CY Clouds Y=' -99.0  99.0   0.0  0.0    0    99     6      4    0   0   0
'2CY Clouds Y=' -99.0  99.0   0.0  0.0    0    99     6      2    0   0   0
'2CY Clouds Y=' -99.0  99.0   0.0  0.0    0    99     6      2    0   0   0
'2CX Clouds X=' -0.02  0.08   0.0  0.0    0    99     8     99    0   0   0
'2CX Clouds X=' -99.0  99.0   0.0  0.0    2    99     6     99    0   0   0
'2CX Clouds X=' -0.1   0.4    0.0  0.0    0    99     4     99    0   0   0
'2CX Clouds X=' -99.0  99.0   0.0  0.0    0    99     6      2    1   0   0
'2CX Clouds X=' -99.0  99.0   0.0  0.0    0    99     6      0    1   0   0
'2CP Clouds P=' -99.0  99.0   0.0  0.0    0    99    28     99    0   0   0
'2CT Clouds T='  0.0   .0005  0.0  0.0    0    99    28     99    0   0   0
'2CT Clouds T='  0.0   .0001  0.0  0.0    0    99    28     99    0   1   0
'2PY PLRPs  Y=' -0.10  0.4    0.0  0.0    0    99     4     99    0   0   2
'2PY PLRPs  Y=' -0.15  0.15   0.0  0.0    0    99     2     0     0   0   0
'2PY PLRPs  Y=' -99.0  99.0   0.0  0.0    0    99     6     99    1   0   0
'2PX PLRPs  X=' -0.1   0.15   0.0  0.0    0    99     4     99    0   0   0
'2PX PLRPs  X=' -0.1   0.4    0.0  0.0    0    99     2     99    0   0   2
'2PP PLRPs  P=' -99.0  99.0   0.0  0.0    0    99     2     99    0   0   0
'2ZY Z eff  Y=' -0.4   0.4    0.0  0.0    0    99     6     99    0   0   0
'2ZY Z eff  Y=' -0.1   0.4    0.0  0.0    0    99     6     99    0   0   2
'2ZY Z eff  Y=' -99.0  99.0   0.0  0.0    0    99    28     99    0   0   0
'2ZX Z eff  X=' -0.05  0.15   0.0  0.0    0    99     6     99    0   0   0
'2ZX Z eff  X=' -99.0  99.0   0.0  0.0    0    99     6     99    0   0   0
'2DX Deps   X=' -0.1   0.0    0.0  0.0    0    99    28     99    0   0   0
'2NX NetEro X=' -0.1   0.0    0.0  0.0    0    99    28     99    0   0   0
'2NX NetEro X=' -0.1   0.0    0.0  0.0    0    99    28     99    0   0   0
'2YX Av dY  X=' -0.05  0.15   0.0  0.0    0    99    28     99    0   0   0
'2WY WalDep Y='  0.0   0.2    0.0  0.0    0    99    28     99    0   0   0
'2WY WalDep Y='  0.0   99.0   0.0  0.0    0    99    28     99    0   0   0
'2RY PowLos Y=' -0.4   0.4    0.0  0.0    0     3     6     99    0   1   0
'2RY PowLos Y=' -0.1   0.4    0.0  0.0    0     3     6     99    0   1   2
'2RY PowLos Y=' -99.0  99.0   0.0  0.0    0    99    28     99    0   0   0
'2RX PowLos X=' -0.05  0.15   0.0  0.0    0    99     6     99    0   1   0
'2RX PowLos X=' -99.0  99.0   0.0  0.0    0    99     6     99    0   1   0
'2RP PowLos P=' -99.0  99.0   0.0  0.0    0    99    28     99    0   0   0
'2LY LinRad Y=' -0.4   0.4    0.0  0.0    0     3     6     99    0   1   0
'2LY LinRad Y=' -0.1   0.4    0.0  0.0    0     3     6     99    0   1   2
'2LY LinRad Y=' -99.0  99.0   0.0  0.0    0    99    28     99    0   0   0
'2LX LinRad X=' -0.05  0.15   0.0  0.0    0    99     6     99    0   1   0
'2LX LinRad X=' -99.0  99.0   0.0  0.0    0    99     6     99    0   1   0
'2LP LinRad P=' -99.0  99.0   0.0  0.0    0    99    28     99    0   0   0
'2IY Ionis  Y=' -0.4   0.4    0.0  0.0    0     0     3     99    0   0   0
'2IY Ionis  Y=' -0.1   0.4    0.0  0.0    0     0     3     99    0   0   2
'2IX Ionis  X=' -0.05  0.15   0.0  0.0    0     0     3     99    0   0   0
'2IX Ionis  X=' -99.0  99.0   0.0  0.0    0     0     4     99    0   0   0
'2IP Ionis  P=' -99.0  99.0   0.0  0.0    0    99    28     99    0   0   0
'2TY CloudT Y=' -99.0  99.0   0.0  0.0    0     0    28     99    0   0   0
'2TX CloudT X=' -99.0  99.0   0.0  0.0    0    99    28     99    0   0   0
'2SY CloudT Y=' -0.5   0.5    0.0  0.0    0    99     3     99    0   0   0
'2SX CloudT X=' -0.1   0.15   0.0  0.0    0    99     3     99    0   0   0
'2YY Av dY  Y=' -99.0  99.0   0.0  0.0    0     0    28     99    0   0   0
'2YX Av dY  X=' -0.05  0.15   0.0  0.0    0    99    28     99    0   0   0
'2CY Clouds Y=' -0.1   0.1    0.0  0.0    0    99     3     -1    0   0   0
'2CY Clouds Y=' -0.05  0.15   0.0  0.0    0    99    28      0    0   0   0
'2CY Clouds Y=' -0.1   0.1    0.0  0.0    0    99    28      1    0   0   0
'2CY Clouds Y=' -0.1   0.1    0.0  0.0    0    99    28      2    0   0   0
$                                                                          
$                                                                          
