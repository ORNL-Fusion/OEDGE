#!/bin/ksh
echo "runout V 6.27 : DIV Input =" $1 ": OUT Input =" $2 ": <Optional Grid = $3>"
if [ $# -lt 2 ]
then 
   echo "usage: $0 <DIV case name> <OUT input file> <Optional - Grid file name>"
   exit 1
fi
#
# Set Run-time environment variables
#
#-------------------------------------
#
export CASENAME=$1
export XLFRTEOPTS=namelist=old
export ADASCENT=/afs/ipp/home/a/adas/adas
#
VER=sun
RUNROOT=/tmp
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
export DIVRUNDIR=$RUNROOT
export DIVDATDIR=$RUNROOT/data
export DIVRESDIR=$RUNROOT/results
#
export DIVMAINDIR=$EXEROOT
export OUTEXEDIR=$EXEROOT/bin/@sys
export EQUDIR=$DATAROOT/shots
#
export OUTEXE=$OUTEXEDIR/out6O$VER
#
#
#-------------------------------------
#
#  Change to Execution directory
#
cd $DIVRUNDIR
mkdir $1
cd $1
#
#  Execute OUT
#
#  Unit  4 is the equilibrium input file
#  Unit  8 is the binary output file from DIVIMP
#  Unit  9 is an echo of the input to OUT
#  Unit 13 is the xy grid
#  Unit 26 is a print file containing the data from all plots
#  Unit 49 is a summary of all plots
#
#  Copy and uncompress the results of the DIVIMP run 
#
cp $DIVRESDIR/$1.raw.Z $DIVRUNDIR/$1/$1.raw.Z
uncompress -f -v $1.raw.Z
mv $1.raw fort.8
#
# Supplimentary raw files (SL)
# Some consolidation of these files needs to
# take place, but it has not been done yet (Sep 6, 2000)
#
if [[ -f $DIVRESDIR/$1.raw.src.Z ]] then
   cp $DIVRESDIR/$1.raw.src.Z $DIVRUNDIR/$1
   uncompress -f -v $1.raw.src.Z
   mv $1.raw.src fort.89
fi
if [[ -f $DIVRESDIR/$1.raw.pla.Z ]] then
   cp $DIVRESDIR/$1.raw.pla.Z $DIVRUNDIR/$1
   uncompress -f -v $1.raw.pla
   mv $1.raw.pla fort.94
fi
if [[ -f $DIVRESDIR/$1.raw.geo.Z ]] then
   cp $DIVRESDIR/$1.raw.geo.Z $DIVRUNDIR/$1
   uncompress -f -v $1.raw.geo
   mv $1.raw.geo fort.95
fi
if [[ -f $DIVRESDIR/$1.raw.rel.Z ]] then
   cp $DIVRESDIR/$1.raw.rel.Z $DIVRUNDIR/$1
   uncompress -f -v $1.raw.rel.Z
   mv $1.raw.rel fort.79
fi
if [[ -f $DIVRESDIR/$1.eirtrc ]] then
   cp $DIVRESDIR/$1.eirtrc eirtrac.dat
fi
#
# Info file for output
#
if [[ -f  $EQUDIR/info.dat ]] then 
  ln -s $EQUDIR/info.dat $DIVRUNDIR/$1/info.dat
fi 
#
#  Connect an experimental data input file to unit 13.
#  If one has been specified. 
#
#  Unit 13 contains experimental data for the case - if available. 
#
if [[ -n $4 ]] then     
  if [[ -f $EQUDIR/$4.experiment ]] then 
    ln -s $EQUDIR/$4.experiment $DIVRUNDIR/$1/fort.13
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
# Link to the toroidal camera inversion data (SL)
#
if [[ -f $EQUDIR/$3.camera ]] then
   ln -s $EQUDIR/$3.camera $DIVRUNDIR/$1/camera.dat
fi
#
#if [[ -n $3 ]] then 
#   if [[ -f $EQUDIR/$3.experiment ]] then 
#      ln -s $EQUDIR/$3.experiment $DIVRUNDIR/$1/fort.13
#   fi 
#fi
#
#
#  Execute OUT
#
$OUTEXE < $DIVDATDIR/$2.d6o > outouta
#
#  Copy or move results.
#
#
#  .outa - OUT debugging file 
#
if [[ -f outouta ]] then 
   mv outouta $DIVRESDIR/$1.outa
fi
#
#  .psga - Postscript plots 
#
if [[ -f POSTSCPT.LIS ]] then 
  ps2pdf POSTSCPT.LIS $DIVRESDIR/$1-a.pdf
  mv POSTSCPT.LIS $DIVRESDIR/$1.psga 
fi
#
#  .daga OUT output file
#
if [[ -f fort.7 ]] then 
   mv fort.7  $DIVRESDIR/$1.daga
fi
#
#  Echo of input
#
if [[ -f fort.9 ]] then
   mv fort.9  $DIVRESDIR/$1.inga
fi
#
#  Special print file 
#
if [[ -f fort.26 ]] then
   mv fort.26 $DIVRESDIR/$1.grpa
fi
#
#  Special print file 
#
if [[ -f fort.49 ]] then
   mv fort.49 $DIVRESDIR/$1.plta
fi
#
#  Print out for Kevin
#
if [[ -f fort.56 ]] then 
   mv fort.56  $DIVRESDIR/$1.probe
fi
#
# AUG file
#
if [[ -f fort.59 ]] then
   mv fort.59 $DIVRESDIR/$1.AUGdiv
fi
#
# Core file in case OUT crashed
#
if [[ -f core ]] then 
  mv core  $DIVRESDIR/$1.core
fi
#
# Move image files if any ...
#
if [[ -f $1_image01.jpg ]] then 
  mv *image*.jpg $DIVRESDIR
fi
#
#  Clean-up
#
rm fort.*
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
# Remove eirtrac.dat if linked
#
if [[ -f eirtrac.dat ]] then
   rm -f eirtrac.dat
fi
#
cd $DIVRUNDIR
rmdir $1
cd $DIVRESDIR
#
#  Compress and generate any remaining print files
#
#compress $1.outa
if [[ -f $1.plta ]] then
   a2ps -R $1.plta --output=$1.plta.ps
fi
#
#  Print
#
#lp $1.plt.ps 
#lp $1.psga
#
#  Return to starting directory. 
#
cd $DIVRUNDIR



