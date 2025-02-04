#!/bin/ksh
echo rundiv V 6.27 : DIV Input = $1 : OUT Input = $2 : Grid = $3 : BG Data = $4 : CFD Data = $5 : DIV Data = $6
if [ $# -lt 3 ] 
then 
   echo "usage: $0 <DIV input file> <OUT input file> <geometry file name> <fluid plasma filename extention - optional> <CFD solution - optional> <DIVIMP solution - optional>"
   echo "Use the word none as the argument for uneeded options" 
   exit 1
fi
#
# Set Run-time environment variables
#
# These should be all that require changing for most installations
#
#-------------------------------------
# 
export CASENAME=$1
export XLFRTEOPTS=namelist=old
export ADASCENT=/afs/ipp/home/a/adas/adas
#
VER=sun
RUNROOT=$HOME/divimp
EXEROOT=$HOME/divimp
DATAROOT=$HOME/divimp
#
#  Loctions of the various DIVIMP files 
#
#  DIVRUNDIR - Directory where DIVIMP will run - must have read/write access
#  DIVDATDIR - Directory containing the input data files
#  DIVRESDIR - Directory where the results from the case will be stored
#
#  DIVMAINDIR- Main DIVIMP directory containing the source code tree 
#  DIVEXEDIR - Directory containing the DIVIMP executable
#  DIVOUTDIR - Directory containing the OUT excutable
#  PINEXEDIR - Directory containing the PIN executable
#  EIREXEDIR - Directory containing the EIRENE executable 
#
#  EQUDIR    - Directory containg the equilibrium, background plasma and pump files
#
#  DIVEXE    - DIVIMP executable
#  OUTEXE    - OUT executable
#  PINEXE    - PIN executable
#  EIREXE    - EIRENE executable
#
export DIVRUNDIR=/tmp
export DIVDATDIR=$RUNROOT/data
export DIVRESDIR=$RUNROOT/results
#
export DIVMAINDIR=$EXEROOT
export DIVEXEDIR=$EXEROOT/bin/@sys
export OUTEXEDIR=$EXEROOT/bin/@sys
export PINEXEDIR=$EXEROOT/pin6
export PIN5EXEDIR=$EXEROOT/pin5
export EIREXEDIR=$EXEROOT/eirene
export EIREXE99DIR=$EXEROOT/bin/@sys
export EQUDIR=$DATAROOT/shots
#
export DIVEXE=$DIVEXEDIR/div6O$VER
export OUTEXE=$OUTEXEDIR/out6O$VER
export PINEXE=$PINEXEDIR/pin6O
export PIN5EXE=$PIN5EXEDIR/pin5O
export EIREXE=$EIREXEDIR/Sources/eirene
export EIREXE99=$EIREXE99DIR/eirene
#
#------------------------------------- 
#
#
# Change to execution directory
#
cd $DIVRUNDIR
mkdir $1
cd $1
#
# Link in the EIRENE databases
#
#
  ln -s $EIREXE99DIR/a+m/HYDHEL HYDHEL
  ln -s $EIREXE99DIR/a+m/SPUTER SPUTER
  ln -s $EIREXE99DIR/a+m/AMJUEL.TEX AMJUEL
  ln -s $EIREXE99DIR/a+m/METHANE METHANE
#
#  Unit  4 is the equilibrium input file
#  Unit 13 is the xy grid file (if available)
#  Unit 14 is the DIVIMP input file - assigned to unit 14 for NIMBUS purposes
#                                   - now connected in the runpin script
#
if [[ -f  $EQUDIR/$3 ]] then
  ln -s $EQUDIR/$3 $DIVRUNDIR/$1/fort.4
else
  echo "Grid file " $EQUDIR/$3 " does not exist."
  echo "Abnormal Script Termination - Script Exiting"
  cd ..
  rm -rf $1
  exit 1
fi
#
# Info file for Eirene option related print outs
#
if [[ -f  $EQUDIR/info.dat ]] then 
  ln -s $EQUDIR/info.dat $DIVRUNDIR/$1/info.dat
fi  
#
# JET hybrid wall data file
#
if [[ -f  $EQUDIR/hybrid.dat ]] then 
  ln -s $EQUDIR/hybrid.dat $DIVRUNDIR/$1/fort.28
#else
#  echo "JET hybrid wall data file not found."
fi  
#
#  CFD OSM files
#
if [[ -f  $DIVRESDIR/$5.cfd.gz ]] then 
  gunzip $DIVRESDIR/$5.cfd.gz
fi
#
if [[ -f  $DIVRESDIR/$5.cfd ]] then 
  echo "Linking ..." $DIVRESDIR/$5.cfd
  ln -s $DIVRESDIR/$5.cfd $DIVRUNDIR/$1/fort.74
fi  
#
#  Other files that are not usually used.
#
# cp $EQUDIR/$3.grd $DIVRUNDIR/$1/fort.13
# cp $EQUDIR/$3.dat $DIVRUNDIR/$1/fort.12
# cp $EQUDIR/$3.prn $DIVRUNDIR/$1/fort.10
#
#  divimp_plasma.dat is a background plasma solution - if required.   
#
if [[ -f  $DIVRESDIR/$6.bgp ]] then 
  echo "Linking ..." $DIVRESDIR/$6.bgp
  ln -s $DIVRESDIR/$6.bgp $DIVRUNDIR/$1/divimp_plasma.dat
fi  
#
#  divimp_aux_data.dat contains various data from a previous DIVIMP run
#
if [[ -f  $DIVRESDIR/$6.auxdata ]] then 
  echo "Linking ..." $DIVRESDIR/$6.auxdata
  ln -s $DIVRESDIR/$6.auxdata $DIVRUNDIR/$1/divimp_aux_data.dat
fi  
#
#  A fluid plasma solution is linked to fort.11 - if one is available.
#
if [[ -n $4 ]] then     
  if [[ -f $EQUDIR/$4 ]] then 
    ln -s $EQUDIR/$4 $DIVRUNDIR/$1/fort.11
  else
    if [[ -f $EQUDIR/$3.$4 ]] then 
      ln -s $EQUDIR/$3.$4 $DIVRUNDIR/$1/fort.11
    else
      if [[ -f $EQUDIR/$3.g80 ]] then 
         ln -s $EQUDIR/$3.g80 $DIVRUNDIR/$1/fort.11
      fi 
    fi 
  fi
else
  if [[ -f $EQUDIR/$3.g80 ]] then 
    ln -s $EQUDIR/$3.g80 $DIVRUNDIR/$1/fort.11
  fi 
fi
#
#  Unit 12 is an auxiliary file for the background plasma solution - if required.   
#
if [[ -n $4 ]] then     
  if [[ -f $EQUDIR/$4.aux ]] then 
    ln -s $EQUDIR/$4.aux $DIVRUNDIR/$1/fort.11
  else
    if [[ -f $EQUDIR/$3.$4.aux ]] then 
      ln -s $EQUDIR/$3.$4.aux $DIVRUNDIR/$1/fort.12
    else
      if [[ -f $EQUDIR/$3.aux ]] then 
         ln -s $EQUDIR/$3.aux $DIVRUNDIR/$1/fort.12
      fi 
    fi
  fi 
else
  if [[ -f $EQUDIR/$3.aux ]] then 
    ln -s $EQUDIR/$3.aux $DIVRUNDIR/$1/fort.12
  fi 
fi
#
#  Unit 13 contains experimental data for the case - if available. 
#
if [[ -n $4 ]] then     
  if [[ -f $EQUDIR/$4.experiment ]] then 
    ln -s $EQUDIR/$4.experiment $DIVRUNDIR/$1/fort.11
  else
    if [[ -f $EQUDIR/$3.$4.experiment ]] then 
      ln -s $EQUDIR/$3.$4.experiment $DIVRUNDIR/$1/fort.13
    else
      if [[ -f $EQUDIR/$3.experiment ]] then 
         ln -s $EQUDIR/$3.experiment $DIVRUNDIR/$1/fort.13
      fi 
    fi
  fi
else
  if [[ -f $EQUDIR/$3.experiment ]] then 
    ln -s $EQUDIR/$3.experiment $DIVRUNDIR/$1/fort.13
  fi 
fi
#
# If an external flux file is available for the shot it will 
# be linked to the "ext_flux.dat" file.
#
#
if [[ -n $4 ]] then     
  if [[ -f $EQUDIR/$4.extflx ]] then 
    ln -s $EQUDIR/$4.extflx $DIVRUNDIR/$1/ext_flux.dat
  else
    if [[ -f $EQUDIR/$3.$4.extflx ]] then 
      ln -s $EQUDIR/$3.$4.extflx $DIVRUNDIR/$1/ext_flux.dat
    else
      if [[ -f $EQUDIR/$3.extflx ]] then 
         ln -s $EQUDIR/$3.extflx $DIVRUNDIR/$1/ext_flux.dat
      fi 
    fi
  fi 
else
  if [[ -f $EQUDIR/$3.extflx ]] then 
    ln -s $EQUDIR/$3.extflx $DIVRUNDIR/$1/ext_flux.dat
  fi 
fi
#
# Link to the toroidal camera inversion data (SL)
#
if [[ -f $EQUDIR/$3.camera ]] then
   ln -s $EQUDIR/$3.camera $DIVRUNDIR/$1/camera.dat
fi
#
# These are the scripts that will run the apprpriate hydrogenic
# neutral code depending on the type of grid. PIN/NIMBUS or EIRENE.
#
#
cp $DIVMAINDIR/runpin rpindiv
#cp $DIVMAINDIR/runpin5 rpindiv5
#cp $DIVMAINDIR/runeire reirediv
cp $DIVMAINDIR/runeire99 reire99
#
# Execute DIVIMP
#
#  Unit  7 is the text output from DIVIMP
#  Unit  9 is an echo of the DIVIMP input
#  Unit 17 is the passing file from DIVIMP to PIN
#  Unit 19 is the temporary file used to record PIN iteration information
#  Unit 21 is the print file from the CALCSOL option
#  Unit 22 is the HTML version of the .dat file
#  Unit 23 is a temporary file used to store HTML as the code runs
#  Unit 25 is a diagnostic print file of grid information
#  Unit 26 is used by OUT to print the numerics for plots
#  Unit 27 is another grid print file
#  Unit 41 is the TRAN file output for JET post-processors
#  Unit 50 is a diagnostic file from the Theta routines and Steve's code
#  Unit 88 is a SOL22 (and other stuff) diagnostic file
#
#  Link DIVIMP input data file - so name can be extracted. 
#
if [[ -f  $DIVDATDIR/$1.d6i ]] then
  ln -s $DIVDATDIR/$1.d6i fort.5
else
  echo "DIVIMP Input file " $DIVDATDIR/$1.d6i " does not exist."
  echo "Abnormal Script Termination - Script Exiting"
  cd ..
  rm -rf $1
  exit 1
fi
#
# Copy EIRENE input file template:
#
if [[ -f $EIREXE99DIR/$1.dat ]] then 
  ln -s $EIREXE99DIR/$1.dat $DIVRUNDIR/$1/fort.80
else
  ln -s $EIREXE99DIR/eirene99.dat $DIVRUNDIR/$1/fort.80
fi
#
#  Run DIVIMP - Use optimized version.
#
echo Starting DIVIMP ...  
#
$DIVEXE < fort.5 > divout1
#
echo DIVIMP Finished.
#
if [[ -n $PINCASE ]] then
  echo "PIN/NIMBUS case: " $PINCASE
fi
#
#  Save SOL22 plot file 
#
if [[ -f POSTSCPT.LIS ]] then 
   mv POSTSCPT.LIS $DIVRESDIR/$1.sol22.ps 
fi
#
#  Copy or move results.
#  .lim Debugging file
#
if [[ -f divout1 ]] then 
   mv divout1 $DIVRESDIR/$1.lim
fi
#
#  .dat case output file
#
if [[ -f fort.7 ]] then 
   mv fort.7  $DIVRESDIR/$1.dat
fi
#
#  .html case output file
#
if [[ -f fort.22 ]] then 
   mv fort.22  $DIVRESDIR/$1.html
fi
#
#  .inp input echo file + other information
#
if [[ -f fort.9 ]] then 
   mv fort.9 $DIVRESDIR/$1.inp
fi
#
#  PIN input plasma file
#
if [[ -f fort.17 ]] then 
   mv fort.17  $DIVRESDIR/$1.pin
fi
#
#  Additional information from SOL option 22.
#
if [[ -f fort.21 ]] then 
   mv fort.21  $DIVRESDIR/$1.sol22
fi
#
#  Print-out from PIN
#
if [[ -f fort.24 ]] then 
  mv fort.24  $DIVRESDIR/$1.pinprn
fi
#
# TRAN file for JET postprocessors
#
if [[ -f fort.41 ]] then 
   mv fort.41  $DIVRESDIR/$1.tran
fi
#
#  Diagnostic information from THETA module (among others)
#
if [[ -f fort.50 ]] then 
   mv fort.50  $DIVRESDIR/$1.theta
fi
#
#  Print out for Kevin
#
if [[ -f fort.56 ]] then 
   mv fort.56  $DIVRESDIR/$1.probe
fi
#
#  Additional information from SOL option 23.
#
if [[ -f fort.71 ]] then 
   mv fort.71  $DIVRESDIR/$1.sol23
fi
#
#  Excel print information from sol23
#
if [[ -f fort.73 ]] then 
   mv fort.73  $DIVRESDIR/$1.exl23
fi
#
#  CFD solution file
#
if [[ -f fort.75 ]] then 
   mv fort.75   $DIVRESDIR/$1.cfd
fi
#
#  Diagnostic files (SL)
#
if [[ -f fort.88 ]] then 
   mv fort.88  $DIVRESDIR/$1.src
fi
#
if [[ -f fort.87 ]] then 
   mv fort.87  $DIVRESDIR/$1.g3 
fi
#
#  DIVIMP format background plasma file - if writing one
#  has been requested. 
#
if [[ -f divimp_plasma.out ]] then 
  mv divimp_plasma.out $DIVRESDIR/$1.bgp
fi
#
#  Additional DIVIMP data file - output in a standardized tagged format. 
#
if [[ -f divimp_aux_data.out ]] then 
  mv divimp_aux_data.out  $DIVRESDIR/$1.auxdata
fi
#
# Core file in case DIVIMP crashed
#
if [[ -f core ]] then 
  mv core  $DIVRESDIR/$1.core
fi
#
# mv fort.25  $DIVRESDIR/$1.lst.grd
# mv fort.27  $DIVRESDIR/$1.cell.grd
#
#  Convert the DIVIMP output data file to postscript. 
#
a2ps $DIVRESDIR/$1.dat --output=$DIVRESDIR/$1.ps 
#
# a2ps -nn $DIVRESDIR/$1.sol > $DIVRESDIR/$1.sol.ps 
# mv POSTSCPT.LIS $DIVRESDIR/$1.sol.psg
#
#  Execute OUT
#
echo Starting OUTput processing ...  
#
$OUTEXE < $DIVDATDIR/$2.d6o > outout1 
#
echo OUTput processing complete.
#
#  Copy or move results of OUT run.
#
#  .out - OUT debugging file 
#
if [[ -f outout1 ]] then 
   mv outout1 $DIVRESDIR/$1.out
fi
#  .raw complete results file
#
if [[ -f fort.8 ]] then 
   mv fort.8 $DIVRESDIR/$1.raw
fi  
#
# AUG file
#
if [[ -f fort.59 ]] then
   mv fort.59 $DIVRESDIR/$1.AUGdiv
fi
#
# Supplimentary .raw files (SL)
#
if [[ -f fort.89 ]] then 
   mv fort.89  $DIVRESDIR/$1.raw.src
fi
if [[ -f fort.94 ]] then 
   mv fort.94  $DIVRESDIR/$1.raw.pla
fi
if [[ -f fort.95 ]] then 
   mv fort.95  $DIVRESDIR/$1.raw.geo
fi
if [[ -f fort.79 ]] then 
   mv fort.79  $DIVRESDIR/$1.raw.rel
fi
#
#  .dag OUT output file
#
if [[ -f fort.7 ]] then 
   mv fort.7  $DIVRESDIR/$1.dag
fi
#
#  .psg - Postscript plots 
#
if [[ -f POSTSCPT.LIS ]] then 
  ps2pdf POSTSCPT.LIS $DIVRESDIR/$1.pdf
  mv POSTSCPT.LIS $DIVRESDIR/$1.psg 
fi
#
#  .ing - Echo of input to graphing routines
#
if [[ -f fort.9 ]] then 
   mv fort.9 $DIVRESDIR/$1.ing
fi
#
#  Special printed listing of some plot information.
#
if [[ -f fort.26 ]] then
   mv fort.26 $DIVRESDIR/$1.grp
fi 
#
#  Special printed listing of some plot information.
#
if [[ -f fort.49 ]] then
   mv fort.49 $DIVRESDIR/$1.plt
fi
#
# EIRENE related output:
#
if [[ -f fort.81 ]] then
   mv fort.81 $DIVRESDIR/$1.eirdat
fi
if [[ -f eirtrac.dat ]] then
   mv eirtrac.dat $DIVRESDIR/$1.eirtrc
fi
#if [[ -f fort.80 ]] then
#   mv fort.80 $DIVRESDIR/$1.eirtmp
#fi
if [[ -f fort.52 ]] then
   mv fort.52 $DIVRESDIR/$1.eirgeo
fi
if [[ -f fort.50 ]] then
   mv fort.50 $DIVRESDIR/$1.debug
fi
if [[ -f fort.85 ]] then
   mv fort.85 $DIVRESDIR/$1.g1
fi
if [[ -f addsur.dat ]] then
   mv addsur.dat $DIVRESDIR/$1.eirsur
fi
#
# Move image files if any ...
#
if [[ -f $1_image01.jpg ]] then 
  mv *image*.jpg $DIVRESDIR
fi
#
# cp fort.13 $EQUDIR/$3.grd 
# a2ps -nn -p $DIVRESDIR/$1.plt > $DIVRESDIR/$1.plt.ps
#
#
# Clean-up
#
rm *.dat
rm fort.*
rm -f rpindiv
rm -f rpindiv5
rm -f reirediv
rm -f reire99
#
rm -f HYDHEL
rm -f SPUTER
rm -f AMJUEL
rm -f METHANE
#
# Remove camera file if linked
#
if [[ -f camera.dat ]] then 
   rm -f camera.dat
fi
#
# Remove info.dat file if linked
#
if [[ -f info.dat ]] then 
  rm -f info.dat
fi 
#
cd $DIVRUNDIR
rmdir $1
cd $DIVRESDIR
#
# Compress the larger output files - check for file existence
#                                    before compressing.
#
if [[ -f $1.raw ]] then
   compress -f -v $1.raw
fi
# (SL)
if [[ -f $1.raw.src ]] then
   compress -f -v $1.raw.src
fi
if [[ -f $1.raw.pla ]] then
   compress -f -v $1.raw.pla
fi
if [[ -f $1.raw.geo ]] then
   compress -f -v $1.raw.geo
fi
if [[ -f $1.raw.rel ]] then
   compress -f -v $1.raw.rel
fi
#
if [[ -f $1.lim ]] then
   compress -f -v $1.lim
fi
#
if [[ -f $1.sol22 ]] then
   compress -f -v $1.sol22
fi
#
if [[ -f $1.sol22 ]] then
   compress -f -v $1.sol23
fi
#
#if [[ -f $1.sol22 ]] then
#   compress -f -v $1.sol22
#fi
#
if [[ -f $1.out ]] then
   compress -f -v $1.out
fi
#
if [[ -f $1.pinout ]] then
   compress -f -v $1.pinout
fi 
#
if [[ -f $1.pinraw ]] then
   compress -f -v $1.pinraw
fi 
#
if [[ -f $1.pinnim ]] then
   compress -f -v $1.pinnim
fi 
#
if [[ -f $1.pinmc ]] then
   compress -f -v $1.pinmc
fi 
#
#  Print
#
# lp $1.ps
# lp $1.sol.ps
# lp $1.plt.ps
# lp $1.psg
# lp $1.sol.ps
#
#  Return to base directory.
#
cd $DIVRUNDIR





