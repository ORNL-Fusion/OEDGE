$                                                                               
$                                                                               
$                                                                               
$   2 DIMENSIONAL PLOTS - R,THETA COORDINATE SYSTEM 
$   ===============================================                                                         
$                                                                               
$ Fill in plotting range and option flags for the graphs required               
$ Graphs can be in any order and repeated if required                          
$ Ensure all titles begin with a reference from the following list :-           
$                                                                               
$  RCA RCR RCS     RRA RRR RRS  RLA RLR RLS                  
$     Clouds        PowerLoss    Line Rad                 
$                                                                               
$ RPA RPR RPS      RIA RIR RIS   RZA RZR RZS                  
$     PLRPs         Ionisation    Zeffs                      
$                                                                               
$ Plotting range: (An over estimate is allowed, eg -99.0 to 99.0)               
$               : for R plots [0,Rwall] - R 
$               : for Theta plots [0,360] or [-180,180] - A 
$               : for Psi (S) plots [-pi,pi] - S 
$ Integration   : The range of values over which the
$    Range      : integration is to occur. The same value for upper
$               : and lower bounds will plot only one line of values.
$               : for R plots     : [0,360] or [-180,180] - R
$               : for Theta plots : [0,Rwall] - A
$               : for Psi(S) plots: ignored - S
$               : overestimates are allowed
$ Plot option   : 0:no plot    1:plot required     2:normalised plot
$ Smooth option : Z:smooth states Z upwards (eg 0:all 5:some 99:none)           
$ Maxiz option  : 0:plot neutrals only  >0:max state to plot                    
$                -1:plot primary neutrals only  99:plot all states              
$ P plane option: 99:integrate over P  -4 to 4: Which P plane required          
$                 (NOT CURRENTLY SUPPORTED FOR R,THETA PLOTS)
$ Yfold option  : 0:off   1:fold -Y results into +Y region, & vice versa        
$ Alliz option  : 0:off   1:include "total of all ionisations" line             
$ View option   : 0:off   1:Xview position   2:Calculate from Sv,Tv,etc
$                 (NOT CURRENTLY SUPPORTED FOR R,THETA PLOTS)         
$                                                                               
$                                                                               
$ Ref   Title      Plotting  Integration  Plot Sm. Maxiz Plane Yfld Allz Vu
$ -------------------------------------------------------------------------
'RCA Clouds A='  -45.0  45.0   0.0  99.0    0   1     8    99    0   0   0
'RCA Clouds A='  -0.0   180.0  0.0  99.0    0   1     8    99    0   0   0
'RCA Clouds A='  -90.0  90.0   0.0  99.0    0  99     6    99    0   0   0
'RCA Clouds A='  -0.1   0.4    0.0  99.0    0  99     6    99    0   0   2
'RCR Clouds R='  -0.02  0.08   0.0 360.0    0  99     8    99    0   0   0
'RCR Clouds R='  -99.0  99.0   0.0 360.0    0  99     6    99    0   0   0
'RCR Clouds R='   0.1   0.8    0.0 360.0    0  99     4    99    0   0   0
'RCS Clouds S='  -99.0  99.0   0.0   0.0    0  99     6     2    1   0   0
'RCS Clouds S='  -99.0  99.0   0.0   0.0    0  99     6     0    1   0   0
'RPA PLRPs  A='  -0.10  0.4    0.0  99.0    0  99     4    99    0   0   2
'RPA PLRPs  A='  -0.15  0.15   0.0  99.0    0  99     2    0     0   0   0
'RPA PLRPs  A='  -99.0  99.0   0.0  99.0    0  99     6    99    1   0   0
'RPR PLRPs  R='  -0.1   0.15   0.0 360.0    0  99     4    99    0   0   0
'RPR PLRPs  R='  -0.1   0.4    0.0 360.0    0  99     2    99    0   0   2
'RPS PLRPs  S='  -0.1   0.15   0.0   0.0    0  99     4    99    0   0   0
'RPS PLRPs  S='  -0.1   0.4    0.0   0.0    0  99     2    99    0   0   2
'RZA Z effs A='  -0.4   0.4    0.0  99.0    0  99     6    99    0   0   0
'RZA Z effs A='  -0.1   0.4    0.0  99.0    0  99     6    99    0   0   2
'RZA Z effs A='  -99.0  99.0   0.0  99.0    0  99    28    99    0   0   0
'RZR Z effs R='  -0.05  0.15   0.0 360.0    0  99     6    99    0   0   0
'RZR Z effs R='  -99.0  99.0   0.0 360.0    0  99     6    99    0   0   0
'RZS Z effs R='  -0.05  0.15   0.0   0.0    0  99     6    99    0   0   0
'RZS Z effs R='  -99.0  99.0   0.0   0.0    0  99     6    99    0   0   0
'RRA PowLos A='  -0.4   0.4    0.0  99.0    0   3     6    99    0   1   0
'RRA PowLos A='  -0.1   0.4    0.0  99.0    0   3     6    99    0   1   2
'RRA PowLos A='  -99.0  99.0   0.0  99.0    0  99    28    99    0   0   0
'RRR PowLos R='  -0.05  0.15   0.0 360.0    0  99     6    99    0   1   0
'RRR PowLos R='  -99.0  99.0   0.0 360.0    0  99     6    99    0   1   0
'RRS PowLos S='  -99.0  99.0   0.0   0.0    0  99    28    99    0   0   0
'RLA LinRad A='  -0.4   0.4    0.0  99.0    0   3     6    99    0   1   0
'RLA LinRad A='  -0.1   0.4    0.0  99.0    0   3     6    99    0   1   2
'RLA LinRad A='  -99.0  99.0   0.0  99.0    0  99    28    99    0   0   0
'RLR LinRad R='  -0.05  0.15   0.0 360.0    0  99     6    99    0   1   0
'RLR LinRad R='  -99.0  99.0   0.0 360.0    0  99     6    99    0   1   0
'RLS LinRad S='  -99.0  99.0   0.0   0.0    0  99    28    99    0   0   0
'RIA Ionis. A='  -0.4   0.4    0.0  99.0    0   0     3    99    0   0   0
'RIA Ionis. A='  -0.1   0.4    0.0  99.0    0   0     3    99    0   0   2
'RIR Ionis. R='  -0.05  0.15   0.0 360.0    0   0     3    99    0   0   0
'RIR Ionis. R='  -99.0  99.0   0.0 360.0    0   0     4    99    0   0   0
'RIS Ionis. S='  -99.0  99.0   0.0   0.0    0  99    28    99    0   0   0