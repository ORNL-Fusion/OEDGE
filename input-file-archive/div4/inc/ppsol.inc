$
$ The following lines specify the parameters for the linearly
$ interpolated specified BG SOL option. This is purely empirical.
$                    S-value      Function Value
$ The form is:         0.0            F0  
$                       S1            F1
$                       S2            F2
$ For S>S2 F=F2 
$
'TN     Te S1 - First  S -value = ctes1 * SMAX       '    0.0
'TN     Te F1 - First  Te-value = ctef1 * te0        '    1.0
'TN     Te S2 - Second S -value = ctes2 * SMAX       '    0.0
'TN     Te F2 - Second Te-value = ctef2 * te0        '    1.0
$
'TN     Ti S1 - First  S -value = ctis1 * SMAX       '    0.0
'TN     Ti F1 - First  Te-value = ctif1 * ti0        '    1.0
'TN     Ti S2 - Second S -value = ctis2 * SMAX       '    0.0
'TN     Ti F2 - Second Te-value = ctif2 * ti0        '    1.0
$
'TN     Nb S1 - First  S -value = cnes1 * SMAX       '    0.0
'TN     Nb F1 - First  Te-value = cnef1 * ne0        '    1.0
'TN     Nb S2 - Second S -value = cnes2 * SMAX       '    0.0
'TN     Nb F2 - Second Te-value = cnef2 * ne0        '    1.0
$
'TN     vb S1 - First  S -value = cvbs1 * SMAX       '    0.0
'TN     vb F1 - First  Te-value = cvbf1 * ti0        '    1.0
'TN     vb S2 - Second S -value = cvbs2 * SMAX       '    0.0
'TN     vb F2 - Second Te-value = cvbf2 * ti0        '    1.0
