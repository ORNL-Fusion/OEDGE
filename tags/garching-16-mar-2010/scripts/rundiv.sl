#!/bin/tcsh 
#
#
# ======================================================================
# PROCESS INPUT FLAGS
# ======================================================================
#
# ----------------------------------------------------------------------
# Remote execution:
#
  setenv MACHINE "local"
  setenv USER    "none"
  if      ( "$1" == "-r" ) then
    shift 
    setenv MACHINE $1
    shift
    setenv USER $1
    shift
  endif
#
# ----------------------------------------------------------------------
# Retrieve data from remote execution:
#
  set FETCH = no
  if      ( "$1" == "-f" ) then
    shift 
    setenv MACHINE $1
    shift
    setenv USER $1
    shift
    set FETCH = yes
  endif
#
# ----------------------------------------------------------------------
# Decide which version of EIRENE to run, and if photon transport is 
# included:
# 
  set EIRTEMP = eirene99
  set EIRFILE = eirene99.dat
  set EIR13 = 
  set EIRPHOTONS = 
  if      ( "$1" == "-04" ) then
    shift 
    set EIRTEMP = eirene04
    set EIRFILE = eirene04.dat
    if      ( "$1" == "-photons" ) then
      shift 
      set EIRPHOTONS = yes
      set EIRFILE = eirene04.dat.photons
    endif
    if      ( "$1" == "-photons.cont" ) then
      shift 
      set EIRPHOTONS = yes
      set EIRFILE = eirene04.dat.photons.cont
    endif
    if      ( "$1" == "-photons.bgk" ) then
      shift 
      set EIRPHOTONS = yes
      set EIRFILE = eirene04.dat.photons.bgk
    endif
    if      ( "$1" == "-13" ) then
      shift 
      set EIR13 = $1
      shift 
    endif
  endif
#
# ----------------------------------------------------------------------
# Identify fluid code data to be loaded:
# 
  set FLUID = none
  if ( "$1" == "-solps" ) then
    shift 
    set FLUID = solps
    set FLUIDCASE = nodrifts
  endif
#
# ----------------------------------------------------------------------
# Just run OUT, the post-processing and graphics package:
# 
  set OUTONLY = no
  if ( "$1" == "-o" ) then
    shift 
    set OUTONLY = yes
  endif
#
# ----------------------------------------------------------------------
# Just run OUT, the post-processing and graphics package:
# 
  set CLEAN = yes
  if ( "$1" == "-nc" ) then
    shift 
    set CLEAN = no
  endif
#
# ----------------------------------------------------------------------
# Case parameter summary to the screen:
#
                             echo 
                             echo  "  DIVIMP case         : $1      "
                             echo  "  OUT input file      : $2      "
                             echo  "  Geometry file       : $3      "
  if ( $4 != "")             echo  "  Saved plasma case   : $4      "
  if ( $OUTONLY == "no")     echo  "  EIRENE version      : $EIRTEMP"
  if ( $MACHINE != "local" ) echo  "  Remote execution on : $MACHINE"
  if ( $MACHINE != "local" ) echo  "  User                : $USER"
                             echo 
#
# ----------------------------------------------------------------------
# Validate input (weak):
#
  if ( $1 == "" || $2 == "" || $3 == "" )  then
    echo "Parmeter error"
    exit
  endif
#
#
# ======================================================================
# SET SHELL VARIABLES
# ======================================================================
#
# DIVIMP root directory:
  setenv DIVHOME ~/divimp

# Set case variables:
  setenv CASENAME $1
  setenv GRIDNAME $3
  setenv SAVENAME $4
  if ( $OUTONLY == "yes" ) setenv SAVENAME $CASENAME

  setenv ADASCENT $DIVHOME/adas/adas
  setenv RESDIR   $DIVHOME/results
  setenv EIRDIR   $DIVHOME/$EIRTEMP

# Location of graphics library used by EIRENE:
  setenv GLI_HOME $DIVHOME/libsrc

# Set local directories:
  set RUNDIR = $DIVHOME
  set BATDIR = $DIVHOME/scripts
  set EXEDIR = $DIVHOME/cases/$CASENAME
  set EQUDIR = $DIVHOME/shots
  set DATDIR = $DIVHOME/data
  set OUTDIR = $DIVHOME/out6
#
#
# ======================================================================
# COPY SETUP FILES
# ======================================================================
#
# ----------------------------------------------------------------------
# Set DIVIMP/OUT executables:
#
  set OUTEXE = out6

  if ( $OUTONLY == "yes" ) then
#   Don't worry about which ... bad... 
    set DIVDIR = $DIVHOME/div6
    set DIVEXE = div6
    set DIVEXT = d6i

  else if ( (-e $DATDIR/$CASENAME.d5i) && (-e $DATDIR/$CASENAME.d6i) ) then

    echo Both .d5i and .d6i cases exist - not sure what code to run
    exit

  else if ( -e $DATDIR/$CASENAME.d6i ) then

#   Execute the main version of DIVIMP, version 6 currently:
    echo "  Executing DIVIMP version 6"
    echo 
    set DIVEXE = div6
    set DIVDIR = $DIVHOME/div6
    set DIVEXT = d6i

  else if ( -e $DATDIR/$CASENAME.d5i ) then

#   Execute thesis version (almost obsolete):
    echo "  Executing DIVIMP version 5 (thesis)"
    echo 
    set DIVEXE = div5   
    set DIVDIR = $DIVHOME/div5
    set DIVEXT = d5i
    if ( $EIRTEMP != 'eirene04' ) setenv EIRDIR $DIVHOME/$EIRTEMP.big

  else

    echo ERROR: $CASENAME not found

  endif
#
# ----------------------------------------------------------------------
# Create and enter execution directory:
#
  cd $DIVHOME/cases
  if ( -e $CASENAME ) rm -rf $CASENAME
  mkdir $CASENAME
  cd    $CASENAME
#
# ----------------------------------------------------------------------
# Copy DIVIMP executable and input file:
#
  cp  $DIVDIR/$DIVEXE           $EXEDIR
  cp  $DATDIR/$CASENAME.$DIVEXT fort.5
  cpc $DATDIR/$2.d6o            outinput.dat

  if ( $EIR13 != "" ) then
    echo "  Getting .eir13 and .eir for $EIR13"
    cp ~/divimp/results/$EIR13.eir13 eirene.13
    cp ~/divimp/results/$EIR13.eir eirene.transfer
  endif

#
# ----------------------------------------------------------------------
# Copy data files from a previous case, if requested:
#
  if ( $SAVENAME != "" ) then

    if ( $OUTONLY == "yes" ) then
      cpc $RESDIR/$SAVENAME.raw.zip .
      unzipc $SAVENAME.raw.zip
      mvc $SAVENAME.raw fort.8

      cpc $RESDIR/$CASENAME.dump      dump.dat
      cpc $RESDIR/$CASENAME.eirtrc    eirtrac.dat
      cpc $RESDIR/$CASENAME.triangles triangles.dat
      cpc $RESDIR/$CASENAME.eir       eirene.transfer
    endif

#   Supplimental .raw files:
    cpc $RESDIR/$SAVENAME.raw.pla.zip .
    cpc $RESDIR/$SAVENAME.raw.src.zip .
    cpc $RESDIR/$SAVENAME.raw.geo.zip .
    cpc $RESDIR/$SAVENAME.raw.vac.zip .
    cpc $RESDIR/$SAVENAME.raw.jum.zip .
    cpc $RESDIR/$CASENAME.raw.tri.zip .

    unzipc $SAVENAME.raw.pla.zip
    unzipc $SAVENAME.raw.src.zip
    unzipc $SAVENAME.raw.geo.zip
    unzipc $SAVENAME.raw.vac.zip
    unzipc $SAVENAME.raw.jum.zip
    unzipc $CASENAME.raw.tri.zip

    mvc $SAVENAME.raw.pla plasma.dat   
    mvc $SAVENAME.raw.src source.dat   
    mvc $SAVENAME.raw.geo geomty.dat  
    mvc $SAVENAME.raw.vac vac-grid.dat  
    mvc $SAVENAME.raw.jum jummap.dat  
    mvc $CASENAME.raw.tri triangles.raw

#   Decide which DIVIMP format background plasma file to copy:
    if      ( -e $EQUDIR/$SAVENAME.bgp ) then 

      cp $EQUDIR/$SAVENAME.bgp divimp_plasma.dat

    else if ( -e $RESDIR/$SAVENAME.bgp ) then 

      cp $RESDIR/$SAVENAME.bgp divimp_plasma.dat

    else if ( -e $EQUDIR/$3.$SAVENAME ) then 

      cp $EQUDIR/$3.$SAVENAME fort.11
#    Use $SAVENAME = g80 for this one
#    else if ( -e $EQUDIR/$3.g80) then
#      cp $EQUDIR/$3.g80 fort.11

    else if ( $SAVENAME != "none" ) then

      echo "Warning: background plasma file $SAVENAME.bgp not found"

    endif

#   Unit 12 is an auxiliary file for the background plasma solution, if required:
    cpc $EQUDIR/$3.$SAVENAME.aux fort.12

#    (These are for continuing a previous EIRENE run, not often used):
#    cpc $RESDIR/$SAVENAME.eird10.zip  . 
#    cpc $RESDIR/$SAVENAME.eird11.zip  .
#    unzipc $SAVENAME.eird11.zip
#    unzipc $SAVENAME.eird10.zip
#    mvc $SAVENAME.eird10 eirene.d10
#    mvc $SAVENAME.eird11 eirene.d11

  endif
#
# ----------------------------------------------------------------------
# Process optional data (most of this is obsolete):
#
  cpc $EQUDIR/$3.pdt fort.91
  cpc $EQUDIR/$3.dat.line match.dat
#  mvc $EQUDIR/casedata/probe.dat  fort.91
#  mvc $EQUDIR/casedata/halpha.dat fort.92
#  mvc $EQUDIR/casedata/camera.dat .
#  mvc $EQUDIR/casedata/pressure.dat .
#  mvc $EQUDIR/casedata/plasma.dat .
#  mvc $EQUDIR/casedata/source.dat .
#  mvc $EQUDIR/casedata/geomty.dat .
  cpc $EQUDIR/hybrid.dat fort.28
  cpc $EQUDIR/info.dat .
#
# ----------------------------------------------------------------------
# Copy fluid code data (from SOLPS, UEDGE, etc.) to run directory: 
#
  if ( $FLUID != "none" ) then
    echo $EQUDIR/$GRIDNAME.$FLUID.$FLUIDCASE
    if ( -e $EQUDIR/$GRIDNAME.$FLUID.$FLUIDCASE ) then
      cp $EQUDIR/$GRIDNAME.$FLUID.$FLUIDCASE/* . 
    else
      echo "  Fluid code $FLUID data directory for case $FLUIDCASE not found"
    endif
  endif
#
# ----------------------------------------------------------------------
# Copy EIRENE related files:
#
# Data for an EIRENE time dependent run (not often used):
  cpc $RESDIR/$CASENAME.eirtim .
#
# Custom EIRENE data file (not often used):
  if ( -e $DATDIR/$CASENAME.eirdat && $FETCH == "no" ) then
    echo "  Copying custom EIRENE input file from data directory"
    cp $DATDIR/$CASENAME.eirdat .
  endif
#
# Copy the EIRENE input file template, which DIVIMP currently uses
# to generate the actual input file.  This will hopefully become
# unnecessary for EIRENE04 in the semi-near future:   
  cpc $EIRDIR/$EIRFILE fort.80

  cpc $GLI_HOME/gksfont.dat .
  cpc $BATDIR/runeire   reire99
  cpc $BATDIR/runeire02 reire02
  cpc $BATDIR/runeire04 reire04

  if ( $EIRTEMP == 'eirene04' ) then

#    cp ~/divimp/eirene99/a+m/HYDHEL     .
#    cp ~/divimp/eirene99/a+m/SPUTER     .
#    cp ~/divimp/eirene99/a+m/AMJUEL.TEX AMJUEL
#    cp ~/divimp/eirene99/a+m/METHANE    .
#    cp ~/divimp/eirene99/a+m/H2VIBR.TEX H2VIBR
#    mkdir TRIM
#    cp ~/divimp/eirene99/TRIM/* TRIM

    cp $EIRDIR/Eirene_04/Database/Surfacedata/SPUTER SPUTER
    cp $EIRDIR/Eirene_04/Database/AMdata/amjuel.tex  AMJUEL
    cp $EIRDIR/Eirene_04/Database/AMdata/hydhel.tex  HYDHEL
    cp $EIRDIR/Eirene_04/Database/AMdata/h2vibr.tex  H2VIBR
    cp $EIRDIR/Eirene_04/Database/AMdata/methane.tex METHANE
    mkdir TRIM
    cp -r $EIRDIR/Eirene_04/Database/Surfacedata/TRIM/* TRIM

#   Executable:
    cp $EIRDIR/eirene $EXEDIR/eirene04

  else

    cp $EIRDIR/a+m/HYDHEL     .
    cp $EIRDIR/a+m/SPUTER     .
    cp $EIRDIR/a+m/AMJUEL.TEX AMJUEL
    cp $EIRDIR/a+m/METHANE    .
    cp $EIRDIR/a+m/H2VIBR.TEX H2VIBR
    mkdir TRIM
    cp $EIRDIR/TRIM/* TRIM
#   Executable:
    cp $EIRDIR/eirene $EXEDIR/eirene99

  endif

#
# ----------------------------------------------------------------------
# SOL23 plasma file (not currently supported, at least not properly):
#
  cpc $RESDIR/$5.cfd $EXEDIR/fort.74
#
# ----------------------------------------------------------------------
# Copy triangle grid related source files:
#
#  if ( $TRIANGLE_FILES == "none" ) then
#  else if !( -e $EQUDIR/"$TRIANGLE_FILES"_triageom.elemente ) then
#    echo "Requested triangle files not found"
#    exit
#  else
#    cpc $EQUDIR/"$TRIANGLE_FILES"_triageom.elemente  triageom.elemente
#    cpc $EQUDIR/"$TRIANGLE_FILES"_triageom.neighbor  triagepm.neighbor
#    cpc $EQUDIR/"$TRIANGLE_FILES"_triageom.npco_char triageom.npco_char
#    cpc $EQUDIR/"$TRIANGLE_FILES"_triageom.plasma    triageom.plasma
#  endif

# Copy TRIANGLE executable, the triangular mesh generator
# currently used in the setup of EIRENE04:
  cp $DIVHOME/triangle/triangle .
#
# ----------------------------------------------------------------------
# Copy magnetic grid files:
#
  if !( -e $EQUDIR/$GRIDNAME ) then
    echo "Error: magnetic grid file $GRIDNAME not found"
    exit
  else
    cp  $EQUDIR/$GRIDNAME fort.4
    cpc $EQUDIR/$GRIDNAME.experiment fort.13

    cpc $EQUDIR/$GRIDNAME.dat        grid.dat
    cpc $EQUDIR/$GRIDNAME.dat.1      grid.dat.1
    cpc $EQUDIR/$GRIDNAME.dat.2      grid.dat.2
    cpc $EQUDIR/$GRIDNAME.dat.3      grid.dat.3
#    cpc $EQUDIR/$GRIDNAME.facelift   facelift.dat
#    cpc $EQUDIR/$GRIDNAME.facelift.2 facelift.dat.2
#    cpc $EQUDIR/$GRIDNAME.facelift.3 facelift.dat.3
#    cpc $EQUDIR/$GRIDNAME.facelift.4 facelift.dat.4
  endif
#
#
# ======================================================================
# EXECUTE
# ======================================================================
#
# ----------------------------------------------------------------------
# OSM-EIRENE-DIVIMP:
#
  if ( $FETCH == "yes" ) then
#
#   --------------------------------------------------------------------
#   Download results from a case that was run remotely: 

    scp "$USER@$MACHINE":divimp/results/$CASENAME.results.zip .
    unzipc $CASENAME.results.zip
    rm -f  $CASENAME.results.zip
  endif

  if ( $FETCH == "no" && $MACHINE != "local" ) then
#
#   --------------------------------------------------------------------
#   Run remotely on linux boxes:

#   Determine if binaries are to be uploaded (whether or not they have
#   been rebuilt since the last remote execution).  Avoiding this can
#   save a lot of time/bandwidth:

    if ( $OUTONLY == "yes" ) then
      echo "  Uploading OUT binary"
      cp $DATDIR/$2.d6o out.input
      cp $RUNDIR/out6/$OUTEXE .
      rm -f $DIVEXE
      rm -f eirene99
      rm -f eirene04
    else
      set OUTONLY = "done"
      if ( -e $DIVDIR/new.$MACHINE ) then
        echo "  Uploading DIVIMP binary"
      else
        rm -rf $DIVEXE
      endif
      if ( -e $EIRDIR/new.$MACHINE ) then
        echo "  Uploading EIRENE binary"
      else
        rm -rf eirene99
        rm -rf eirene04
      endif
      if ( -e $OUTDIR/new.$MACHINE ) then
        echo "  Uploading OUT binary"
        cp $RUNDIR/out6/$OUTEXE .
      endif
    endif

#   Prepare and upload case files:
    zip -r -q $CASENAME.zip *
    cp  $BATDIR/rundiv-remote rundiv-$CASENAME
    scp rundiv-$CASENAME $CASENAME.zip $USER@$MACHINE':'
    rm -rf * 

#   Register that binaries have been uploaded, so that they are not uploaded again
#   until they are rebuilt:
    if ( $OUTONLY == "no" || $OUTONLY == "done" ) then
      if ( -e $DIVDIR/new.$MACHINE ) rm -rf $DIVDIR/new.$MACHINE 
      if ( -e $EIRDIR/new.$MACHINE ) rm -rf $EIRDIR/new.$MACHINE
      if ( -e $OUTDIR/new.$MACHINE ) rm -rf $OUTDIR/new.$MACHINE
    endif

    if ( $2 == "upload" || $2 == "upload-run" ) then

#     Upload the case and drop the ssh link, with the option of starting the
#     case with a remote call to the at command:
      echo "./rundiv-$CASENAME $MACHINE $CASENAME $DIVEXE $2 $OUTONLY" > command-$CASENAME.txt

      chmod u+x command-$CASENAME.txt
      scp command-$CASENAME.txt $USER@$MACHINE':'.
      if ( $2 == "upload-run" ) then
        ssh -l $USER -n $MACHINE "at -f ./command-$CASENAME.txt now"
      endif

#     Clean and stop this script:
      cd $DIVHOME/cases
      if ( -e $CASENAME ) then
        rm -r -f $CASENAME
      endif
      echo Case uploaded, stopping run script
      exit

    else

#     Execute in active shell:
      ssh -l $USER -n $MACHINE "~/rundiv-$CASENAME $MACHINE $CASENAME $DIVEXE $2 $OUTONLY"

#     Download results:
      scp    $USER@$MACHINE':'$CASENAME.results.zip .
      unzipc $CASENAME.results.zip 
      rm     *.zip

#     Clean remote computer:
      ssh -l $USER -n $MACHINE "rm ~/rundiv-$CASENAME"
      ssh -l $USER -n $MACHINE "rm ~/$CASENAME.results.zip"

    endif

  else

    if ( $OUTONLY == "no" && $FETCH == "no" ) then
#
#     --------------------------------------------------------------------
#     Execute DIVIMP locally:
#
#      Xpgdbg $DIVEXE < fort.5 > divout1
      $DIVEXE < fort.5 > divout1

#     Process this output now, since these file units conflict with OUT:
      mvc POSTSCPT.LIS old.LIS
      mvc fort.50 $CASENAME.slm
      mvc fort.7  $CASENAME.dat
      mvc fort.22 $CASENAME.html
      mvc fort.88 $CASENAME.src

# Clean this up, may not be necessary, get rid of the raw-'s and just use the renamed names...:
      mvc raw-plasma.dat   plasma.dat   
      mvc raw-sources.dat  source.dat   
      mvc raw-geometry.dat geomty.dat  
    endif
#
#   ----------------------------------------------------------------------
#   Execute OUT:
#
    if ( $2 == "none" ) then
      echo "OUT not executed (by request)"

    else if !( -e $DATDIR/$2.d6o ) then
      echo "OUT not executed (input file $2.d6o not found)"

    else if !( -e fort.8 ) then
      echo "OUT not executed (RAW data file not found)"

    else

      cp $DATDIR/$2.d6o fort.99

#     Decide to run big OUT (thesis) or little OUT:
      if ( (-e $DATDIR/$CASENAME.d5i) && (-e $DATDIR/$CASENAME.d6i) ) then

        echo Both .d5i and .d6i case exist - not sure what code to run
        exit

      else if ( -e $DATDIR/$CASENAME.d6i ) then

        echo
        echo "  Executing OUT version 6"
        echo 
        $RUNDIR/out6/$OUTEXE < fort.99 > outout1

      else if ( -e $DATDIR/$CASENAME.d5i ) then

        echo 
        echo "  Executing OUT version 6 - BIG"
        echo 
        $RUNDIR/out6.big/$OUTEXE < fort.99 > outout1

      endif

    endif

  endif
#
#
# ======================================================================
# PROCESS RESULTS
# ======================================================================
#
# ----------------------------------------------------------------------
# Rename DIVIMP/OUT files:
#

#  if ( -e POSTSCPT.LIS ) then
#    ps2pdf POSTSCPT.LIS $CASENAME.pdf
#  endif

  mvc fort.22 $CASENAME.html

#  mvc fort.99 $CASENAME.divdat
  mvc divout1 $CASENAME.lim
  mvc fort.17 $CASENAME.pin
  mvc fort.21 $CASENAME.sol
  mvc fort.62 $CASENAME.bgp
  mvc divimp_plasma.out $CASENAME.bgp

  mvc outout1      $CASENAME.out
  mvc POSTSCPT.LIS $CASENAME.ps
  mvc fort.49      $CASENAME.plt
  mvc outdata.dat  $CASENAME.dat.out

  mvc fort.8     $CASENAME.raw
  mvc plasma.dat $CASENAME.raw.pla
  mvc source.dat $CASENAME.raw.src
  mvc geomty.dat $CASENAME.raw.geo
  mvc triangles.raw $CASENAME.raw.tri
  mvc fort.79       $CASENAME.raw.rel
  mvc vac-grid.dat  $CASENAME.raw.vac
  mvc jummap.dat    $CASENAME.raw.jum

  mvc fort.90 $CASENAME.osm.log
  mvc fort.85 $CASENAME.g1
  mvc fort.86 $CASENAME.g2
  mvc fort.87 $CASENAME.g3
  mvc fort.65 $CASENAME.e1
  mvc fort.66 $CASENAME.e2
  mvc fort.67 $CASENAME.e3
  mvc fort.25 $CASENAME.sonnet

  mvc dump.dat $CASENAME.dump

  mvc fort.52 $CASENAME.eig

  mvc line.dat   $CASENAME.dat.line
  mvc strata.dat $CASENAME.dat.strata
  mvc bgk.dat    $CASENAME.dat.bgk
#  mvc acd.dat    $CASENAME.dat.acd
  mvc rec.dat    $CASENAME.dat.rec

  mvc fort.9     $CASENAME.inp

  mvc addsur.dat    $CASENAME.addsur

  mvc triangles.dat    $CASENAME.triangles
  mvc triangles.points $CASENAME.triangles.points
  mvc triangles.sides  $CASENAME.triangles.sides
  mvc triangles.map    $CASENAME.triangles.map
  mvc triangles.plasma $CASENAME.triangles.plasma
  mvc triangles.efield $CASENAME.triangles.efield

  mvc quadrangles.geometry $CASENAME.quadrangles.geometry
  mvc quadrangles.plasma   $CASENAME.quadrangles.plasma
  mvc quadrangles.efield   $CASENAME.quadrangles.efield

  mvc plasma-transfer.dat $CASENAME.rzplasma

  mvc celldata.dat $CASENAME.celldata

  mvc eirene.transfer $CASENAME.eir

  mv *image*.jpg $RESDIR

  if ($EIRPHOTONS == "yes") mvc eirene.13 $CASENAME.eir13
#
# ----------------------------------------------------------------------
# Rename some EIRENE files:
#
  mvc eirene.d10  $CASENAME.eird10
  mvc eirene.d11  $CASENAME.eird11
  mvc eirtrac.dat $CASENAME.eirtrc
#  rm *.eirtim
#
# ----------------------------------------------------------------------
# Rename SOL23 files:
#
  mvc fort.71 $CASENAME.sol23
  mvc fort.73 $CASENAME.exl23
  mvc fort.75 $CASENAME.cfd
#
# ----------------------------------------------------------------------
# Compress large files:
#
  zipc $CASENAME.raw
  zipc $CASENAME.raw.src
  zipc $CASENAME.raw.pla
  zipc $CASENAME.raw.geo
  zipc $CASENAME.raw.rel
  zipc $CASENAME.raw.vac
  zipc $CASENAME.raw.jum
  zipc $CASENAME.raw.tri
#  zipc $CASENAME.eir
#  zipc $CASENAME.pin
#  zipc $CASENAME.eird10
#  zipc $CASENAME.eird11
#
# ----------------------------------------------------------------------
# Move files for upload to the web:
#
#  mvc $CASENAME.html $RESDIR/transfer
#  mvc $CASENAME.pdf  $RESDIR/transfer
#
# ----------------------------------------------------------------------
#
# Move all files to results directory, over-writing any previous results:
#
  mv -f $CASENAME.* $RESDIR
#
# ----------------------------------------------------------------------
# Clean:
#
  if ( $CLEAN == "yes" ) then
    cd $DIVHOME/cases
    rm -r -f $CASENAME
  else
    echo "  Not cleaning execution directory"
  endif
#
#
#





