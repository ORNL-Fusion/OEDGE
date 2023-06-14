c slmod begin
c
c maxnfla added to gfsub3w and maptodiv calls
c nfla changed to nfla2
c
c slmod end
      subroutine eireprn
      use mod_params
      use mod_pindata
      use mod_cgeom
      implicit none
c     include 'params'
c     include 'pindata'
c     include 'cgeom'

      call prrmatdiv(pinion,maxnks,nks(irsep),nrs,24,'ionization')
      call prrmatdiv(pinatom,maxnks,nks(irsep),nrs,24,'neutrals')
      call prrmatdiv(pinalpha,maxnks,nks(irsep),nrs,24,'halpha')
      call prrmatdiv(pinqi,maxnks,nks(irsep),nrs,24,'pinqe')
      call prrmatdiv(pinqe,maxnks,nks(irsep),nrs,24,'pinqi')
      call prrmatdiv(pinmol,maxnks,nks(irsep),nrs,24,'molecules')

      return
      end

c     ------------------------------------------------------------------

      subroutine readeire_97

      use mod_params
      use mod_io_units
      use mod_pindata
      use mod_comtor
      use mod_cgeom
      use mod_slcom
      implicit none
c     include 'params'
c     include 'pindata'
c     include 'comtor'
c     include 'cgeom'
c slmod begin
c     INCLUDE 'slcom'
c slmod end
c
c     readeire: this subroutine is designed to read in the data
c     required from EIRENE output files in the specific data format
c     used by those files and to scale them to the units used in
c     DIVIMP
c
      integer ir,ik,topik
      real sum,sum1,sum2,sum1a,sum2a,sum3,sum3a
c
c     nfla2 = number of ion species in the Braams data file;
c     for the EIRENE file (which is in similar format, I guess
c     it has to be 1)

      integer maxix,maxiy,maxis
      parameter (maxix=100,maxiy=50,maxis=7)
      integer nfla2
      parameter (nfla2=1)
      real ndummy(0:maxix+1,0:maxiy+1,maxis)
      real tdummy(0:maxix+1,0:maxiy+1)
      real signv
      integer nx,ny,nxd,nyd
c
c     needed for temporary vars
c
      real pinbuf(maxnks,maxnrs)
c slmod begin
      INTEGER i1
      REAL    rdum1
c slmod end
c
c     calls to read routine
c
c     BEWARE: the EIRENE pin file only writes out the fields without
c             the dead cells, i.e the number of rings and knots is
c             the numbers of the used grid -2 !
c
c     NOTE: -4 is used here because -
c              1) There are two dead cells on each ring that are
c                 not used in the Eirene version of the B2 file.
c                 This accounst for -2
c              2) The numbers maxkpts and maxrings are set
c                 in DIVIMP for use in loops from 1 to maxrings or
c                 maxkpts. Whereas the reading routine gfsub3r uses
c                 loops running from 0 to nx +1 - and from 0 to
c                 ny +1 - thus to get the correct number of both
c                 knots and rings in the read routine for the Eirene
c                 data - another -2 must be taken from the indices.
c
c
      nx = maxkpts-4
      ny = maxrings-4
      nxd = maxix
      nyd = maxiy
c
c-----------------------------------------------------------------------
c
c     read in data from pin passing file
c
c-----------------------------------------------------------------------
c
c     zero arrays (should be done below)
c
      call rzero (pinvdist,3*14*maxnks*maxnrs)
c
      rewind(ipinout)
c
c-----------------------------------------------------------------------
c
c     read in the pin data and map them from B2 storing convention
c     into divimp format arrays. note also that eirene passes in
c     units per cell so conversion with cell volume is necessary
*     Krieger IPP/96: I think we made an error here:
*     1. ion source had wrong unit (see below, does not really matter
*                                   if flux tubes are renormalized)
*        from Ampere to particles/s -> divide by 1.602e-19
*                                      (done in EIRENE now)
*     2. as cell volumes, we have to use kvol2/kvols instead of
*        karea2/kareas as in our first version
*
*        this was corrected IPP/96 ...

c     Krieger, IPP/97

c     the structure of the pin data file written by EIRENE is now as
c     follows:
c     - H+ source                [particles/(cell x sec)]     PINION
c     - Impurity ion source      [particles/(cell x sec)]     PINIONZ
c     - Hydrogen momentum source [(kg x m)/(cell x sec^2)]    PINMP
c       This is a cell edge quantity and needs special treatment
c       for DIVIMP presumably.
c       Sign convention like for velocities.
c     - Impurity momentum source [(kg x m)/(cell x sec^2)]
c     - H neutral density        [particles/cell]             PINATOM
c     - Impurity neutral density [particles/cell]             PINZ0
c     - H2 density               [molecules/cell]             PINMOL
c     - H2+ density              [molecules/cell]            +PINMOL
c     - Avg. neutr. H energy     [eV]                         PINENA
c     - Avg. neutr. imp. energy  [ev]                         PINENZ
c     - Avg. H2 energy           [eV]                         PINENM
c     - Avg. H2+ energy          [eV]  !!! (not yet added to DIVIMP) !!!
c     - Halpha emissivity        [photons/(cell x sec)]       PINALPHA
c       (5 blocks for different generation mechanisms)
c     - H rec. source            [particles/(cell x sec)]     PINREC
c     - H+ (+imp.) energy source [Watt/cell]                  PINQI
c     - Electron energy source   [Watt/cell]                  PINQE

c     Neutral particle velocity distribution information from nimbus
c     (pinvdist) has no EIRENE equivalent yet
c
c     NB: the EIRENE output comes already with the dead cells removed
c         at both targets. Therefore, if we have a target option with
c         dead cells included (like ctargopt=4), we need to add these!
c
c-----------------------------------------------------------------------
c
c     H+ source
c
      call rzero (pinion,maxnks*maxnrs)
      call gfsub3r(ipinout,nx,ny,nxd,nyd,nfla2,maxnfla,tdummy(0,0),0)
      call deadadd(tdummy,nx,ny,nxd,nyd)
      call maptodiv(cutring,cutpt1,cutpt2,nx+2,ny+2,nxd,nyd,
     >             nfla2,maxnfla,tdummy(0,0),pinion,maxnks,maxnrs,1.0,0)
c
c     normalize from "per cell" to "per volume"
c
      call eire_renorm(pinion,kvols,nrs,nks,irsep,sum,0)
c
      write(6,*) 'Total ionizations per second (PIN): ',sum
c
c-----------------------------------------------------------------------
c
c     Impurity ion source
c
      call rzero (pinionz,maxnks*maxnrs)
      call gfsub3r(ipinout,nx,ny,nxd,nyd,nfla2,maxnfla,tdummy(0,0),0)
      call deadadd(tdummy,nx,ny,nxd,nyd)
      call maptodiv(cutring,cutpt1,cutpt2,nx+2,ny+2,nxd,nyd,
     >            nfla2,maxnfla,tdummy(0,0),pinionz,maxnks,maxnrs,1.0,0)
c
c     normalize from "per cell" to "per volume"
c
      call eire_renorm(pinionz,kvols,nrs,nks,irsep,sum,0)
      write(6,*) 'Total imp. ionizations per second (PIN): ',sum
c
c-----------------------------------------------------------------------
c
c     Hydrogen momentum source
c
      call rzero (pinqi,maxnks*maxnrs)
      call gfsub3r(ipinout,nx,ny,nxd,nyd,nfla2,maxnfla,tdummy(0,0),0)
      call deadadd(tdummy,nx,ny,nxd,nyd)
      call maptodiv(cutring,cutpt1,cutpt2,nx+2,ny+2,nxd,nyd,
     >              nfla2,maxnfla,tdummy(0,0),pinmp,maxnks,maxnrs,1.0,0)
c
c     normalize from "per cell" to "per volume"
c
      call eire_renorm(pinmp,kvols,nrs,nks,irsep,sum,0)
c
c-----------------------------------------------------------------------
c
c     Impurity momentum source  (NOT USED YET)
c
      call rzero (pinbuf,maxnks*maxnrs)
      call gfsub3r(ipinout,nx,ny,nxd,nyd,nfla2,maxnfla,tdummy(0,0),0)
      call deadadd(tdummy,nx,ny,nxd,nyd)
      call maptodiv(cutring,cutpt1,cutpt2,nx+2,ny+2,nxd,nyd,
     >             nfla2,maxnfla,tdummy(0,0),pinbuf,maxnks,maxnrs,1.0,0)
c
c     normalize from "per cell" to "per volume"
c
      call eire_renorm(pinbuf,kvols,nrs,nks,irsep,sum,1)
c
c-----------------------------------------------------------------------
c
c     H neutral density
c
      call rzero (pinatom,maxnks*maxnrs)
      call gfsub3r(ipinout,nx,ny,nxd,nyd,nfla2,maxnfla,tdummy(0,0),0)
      call deadadd(tdummy,nx,ny,nxd,nyd)
      call maptodiv(cutring,cutpt1,cutpt2,nx+2,ny+2,nxd,nyd,
     >            nfla2,maxnfla,tdummy(0,0),pinatom,maxnks,maxnrs,1.0,0)
c
c     normalize from "per cell" to "per volume"
c
      call eire_renorm(pinatom,kvols,nrs,nks,irsep,sum,0)
c
c-----------------------------------------------------------------------
c
c     Impurity neutral density (optional for certain launch options)
c
      call rzero (pinz0,maxnks*maxnrs)
      call gfsub3r(ipinout,nx,ny,nxd,nyd,nfla2,maxnfla,tdummy(0,0),0)
      call deadadd(tdummy,nx,ny,nxd,nyd)
      call maptodiv(cutring,cutpt1,cutpt2,nx+2,ny+2,nxd,nyd,
     >              nfla2,maxnfla,tdummy(0,0),pinz0,maxnks,maxnrs,1.0,0)
c
c     normalize from "per cell" to "per volume"
c
      call eire_renorm(pinz0,kvols,nrs,nks,irsep,sum,0)
c
c-----------------------------------------------------------------------
c
c     H2 and H2+ molecular density, summed up in one array
c
c     first H2
c
      call rzero (pinmol,maxnks*maxnrs)
      call gfsub3r(ipinout,nx,ny,nxd,nyd,nfla2,maxnfla,tdummy(0,0),0)
      call deadadd(tdummy,nx,ny,nxd,nyd)
      call maptodiv(cutring,cutpt1,cutpt2,nx+2,ny+2,nxd,nyd,
     >             nfla2,maxnfla,tdummy(0,0),pinmol,maxnks,maxnrs,1.0,0)
c
c     then H2+
c
      call rzero (pinbuf,maxnks*maxnrs)
      call gfsub3r(ipinout,nx,ny,nxd,nyd,nfla2,maxnfla,tdummy(0,0),0)
      call deadadd(tdummy,nx,ny,nxd,nyd)
      call maptodiv(cutring,cutpt1,cutpt2,nx+2,ny+2,nxd,nyd,
     >             nfla2,maxnfla,tdummy(0,0),pinbuf,maxnks,maxnrs,1.0,0)
c     add up all contributions, then normalize ...
      call sum_contrib(pinmol,pinbuf,nrs,nks,2)
c
      call eire_renorm(pinmol,kvols,nrs,nks,irsep,sum,0)
c
c-----------------------------------------------------------------------
c
c     Average neutral H energy
c
      call rzero (pinena,maxnks*maxnrs)
      call gfsub3r(ipinout,nx,ny,nxd,nyd,nfla2,maxnfla,tdummy(0,0),0)
      call deadadd(tdummy,nx,ny,nxd,nyd)
      call maptodiv(cutring,cutpt1,cutpt2,nx+2,ny+2,nxd,nyd,
     >             nfla2,maxnfla,tdummy(0,0),pinena,maxnks,maxnrs,1.0,0)
c
c     normalize from "per cell" to "per volume"
c
      call eire_renorm(pinena,kvols,nrs,nks,irsep,sum,1)
c
c-----------------------------------------------------------------------
c
c     Average neutral imp. energy
c
      call rzero (pinenz,maxnks*maxnrs)
      call gfsub3r(ipinout,nx,ny,nxd,nyd,nfla2,maxnfla,tdummy(0,0),0)
      call deadadd(tdummy,nx,ny,nxd,nyd)
      call maptodiv(cutring,cutpt1,cutpt2,nx+2,ny+2,nxd,nyd,
     >             nfla2,maxnfla,tdummy(0,0),pinenz,maxnks,maxnrs,1.0,0)
c
c     normalize from "per cell" to "per volume"
c
      call eire_renorm(pinenz,kvols,nrs,nks,irsep,sum,1)
c
c-----------------------------------------------------------------------
c
c     average H2 molecule energy
c
      call rzero (pinenm,maxnks*maxnrs)
      call gfsub3r(ipinout,nx,ny,nxd,nyd,nfla2,maxnfla,tdummy(0,0),0)
      call deadadd(tdummy,nx,ny,nxd,nyd)
      call maptodiv(cutring,cutpt1,cutpt2,nx+2,ny+2,nxd,nyd,
     >             nfla2,maxnfla,tdummy(0,0),pinenm,maxnks,maxnrs,1.0,0)
c
c     normalize from "per cell" to "per volume"
c
      call eire_renorm(pinenm,kvols,nrs,nks,irsep,sum,1)
c
c-----------------------------------------------------------------------
c
c     average H2+ molecule energy      !!! (not yet added to DIVIMP) !!!
c
      call gfsub3r(ipinout,nx,ny,nxd,nyd,nfla2,maxnfla,tdummy(0,0),0)
      call deadadd(tdummy,nx,ny,nxd,nyd)
      call maptodiv(cutring,cutpt1,cutpt2,nx+2,ny+2,nxd,nyd,
     >             nfla2,maxnfla,tdummy(0,0),pinbuf,maxnks,maxnrs,1.0,0)
c
c     normalize from "per cell" to "per volume"
c
      call eire_renorm(pinbuf,kvols,nrs,nks,irsep,sum,1)
c
c-----------------------------------------------------------------------
c
c     Halpha emissivity (various contributions summed up)
c
c     halpha emission 1 (emitted by H up to ionization)
c
      call rzero (pinalpha,maxnks*maxnrs)
      call gfsub3r(ipinout,nx,ny,nxd,nyd,nfla2,maxnfla,tdummy(0,0),0)
      call deadadd(tdummy,nx,ny,nxd,nyd)
      call maptodiv(cutring,cutpt1,cutpt2,nx+2,ny+2,nxd,nyd,
     >           nfla2,maxnfla,tdummy(0,0),pinalpha,maxnks,maxnrs,1.0,0)
c
c     halpha emission 2 (emitted in H+ recombination)
c
      call rzero (pinbuf,maxnks*maxnrs)
      call gfsub3r(ipinout,nx,ny,nxd,nyd,nfla2,maxnfla,tdummy(0,0),0)
      call deadadd(tdummy,nx,ny,nxd,nyd)
      call maptodiv(cutring,cutpt1,cutpt2,nx+2,ny+2,nxd,nyd,
     >             nfla2,maxnfla,tdummy(0,0),pinbuf,maxnks,maxnrs,1.0,0)
c
c     add up contributions
c
      call sum_contrib(pinalpha,pinbuf,nrs,nks,2)
c
c     halpha emission 3 (emitted by H2 up to dissociation)
c
      call rzero (pinbuf,maxnks*maxnrs)
      call gfsub3r(ipinout,nx,ny,nxd,nyd,nfla2,maxnfla,tdummy(0,0),0)
      call deadadd(tdummy,nx,ny,nxd,nyd)
      call maptodiv(cutring,cutpt1,cutpt2,nx+2,ny+2,nxd,nyd,
     >             nfla2,maxnfla,tdummy(0,0),pinbuf,maxnks,maxnrs,1.0,0)
c
c     add up contributions
c
      call sum_contrib(pinalpha,pinbuf,nrs,nks,2)
c
c     halpha emission 4 (emitted by H2+ up to dissociation)
c
      call rzero (pinbuf,maxnks*maxnrs)
      call gfsub3r(ipinout,nx,ny,nxd,nyd,nfla2,maxnfla,tdummy(0,0),0)
      call deadadd(tdummy,nx,ny,nxd,nyd)
      call maptodiv(cutring,cutpt1,cutpt2,nx+2,ny+2,nxd,nyd,
     >             nfla2,maxnfla,tdummy(0,0),pinbuf,maxnks,maxnrs,1.0,0)
c
c     add up contributions
c
      call sum_contrib(pinalpha,pinbuf,nrs,nks,2)
c
c     halpha emission 5 (emitted in CX of H and H+)
c
      call rzero (pinbuf,maxnks*maxnrs)
      call gfsub3r(ipinout,nx,ny,nxd,nyd,nfla2,maxnfla,tdummy(0,0),0)
      call deadadd(tdummy,nx,ny,nxd,nyd)
      call maptodiv(cutring,cutpt1,cutpt2,nx+2,ny+2,nxd,nyd,
     >             nfla2,maxnfla,tdummy(0,0),pinbuf,maxnks,maxnrs,1.0,0)
c
c     add up contributions, then normalize ...
c
      call sum_contrib(pinalpha,pinbuf,nrs,nks,2)
c
      call eire_renorm(pinalpha,kvols,nrs,nks,irsep,sum,0)
c
c slmod begin
c
c-----------------------------------------------------------------------
c
c     Hgamma emissivity (various contributions summed up)
c
c     Hgamma emission 1 (emitted by H up to ionization)
c
      call rzero (pinline,maxnks*maxnrs*6*2)
      call gfsub3r(ipinout,nx,ny,nxd,nyd,nfla2,maxnfla,tdummy(0,0),0)
      call deadadd(tdummy,nx,ny,nxd,nyd)
      call maptodiv(cutring,cutpt1,cutpt2,nx+2,ny+2,nxd,nyd,
     >           nfla2,maxnfla,tdummy(0,0),pinline(1,1,1,H_BGAMMA),
     .           maxnks,maxnrs,1.0,0)
c
c     Hgamma emission 2 (emitted in H+ recombination)
c
      call rzero (pinbuf,maxnks*maxnrs)
      call gfsub3r(ipinout,nx,ny,nxd,nyd,nfla2,maxnfla,tdummy(0,0),0)
      call deadadd(tdummy,nx,ny,nxd,nyd)
      call maptodiv(cutring,cutpt1,cutpt2,nx+2,ny+2,nxd,nyd,
     >             nfla2,maxnfla,tdummy(0,0),pinbuf,maxnks,maxnrs,1.0,0)
c
c     add up contributions
c
      call sum_contrib(pinline(1,1,2,H_BGAMMA),pinbuf,nrs,nks,2)
c
c     Hgamma emission 3 (emitted by H2 up to dissociation)
c
      call rzero (pinbuf,maxnks*maxnrs)
      call gfsub3r(ipinout,nx,ny,nxd,nyd,nfla2,maxnfla,tdummy(0,0),0)
      call deadadd(tdummy,nx,ny,nxd,nyd)
      call maptodiv(cutring,cutpt1,cutpt2,nx+2,ny+2,nxd,nyd,
     >             nfla2,maxnfla,tdummy(0,0),pinbuf,maxnks,maxnrs,1.0,0)
c
c     add up contributions
c
      call sum_contrib(pinline(1,1,3,H_BGAMMA),pinbuf,nrs,nks,2)
c
c     Hgamma emission 4 (emitted by H2+ up to dissociation)
c
      call rzero (pinbuf,maxnks*maxnrs)
      call gfsub3r(ipinout,nx,ny,nxd,nyd,nfla2,maxnfla,tdummy(0,0),0)
      call deadadd(tdummy,nx,ny,nxd,nyd)
      call maptodiv(cutring,cutpt1,cutpt2,nx+2,ny+2,nxd,nyd,
     >             nfla2,maxnfla,tdummy(0,0),pinbuf,maxnks,maxnrs,1.0,0)
c
c     add up contributions
c
      call sum_contrib(pinline(1,1,4,H_BGAMMA),pinbuf,nrs,nks,2)
c
c     Hgamma emission 5 (emitted in CX of H and H+)
c
      call rzero (pinbuf,maxnks*maxnrs)
      call gfsub3r(ipinout,nx,ny,nxd,nyd,nfla2,maxnfla,tdummy(0,0),0)
      call deadadd(tdummy,nx,ny,nxd,nyd)
      call maptodiv(cutring,cutpt1,cutpt2,nx+2,ny+2,nxd,nyd,
     >             nfla2,maxnfla,tdummy(0,0),pinbuf,maxnks,maxnrs,1.0,0)
c
c     add up contributions, then normalize ...
c
      call sum_contrib(pinline(1,1,5,H_BGAMMA),pinbuf,nrs,nks,2)
c
      DO i1 = 1, 5
        call sum_contrib(pinline(1,1, 6,H_BGAMMA),
     .                   pinline(1,1,i1,H_BGAMMA),nrs,nks,2)
        call eire_renorm(pinline(1,1,i1,H_BGAMMA),kvols,nrs,nks,irsep,
     .                   sum,0)
      ENDDO
      call eire_renorm(pinline(1,1,6,H_BGAMMA),kvols,nrs,nks,irsep,
     .                 sum,0)
c
c slmod end
c-----------------------------------------------------------------------
c
c     H recombination source rate
c
      call rzero (pinrec,maxnks*maxnrs)
      call gfsub3r(ipinout,nx,ny,nxd,nyd,nfla2,maxnfla,tdummy(0,0),0)
      call deadadd(tdummy,nx,ny,nxd,nyd)
      call maptodiv(cutring,cutpt1,cutpt2,nx+2,ny+2,nxd,nyd,
     >             nfla2,maxnfla,tdummy(0,0),pinrec,maxnks,maxnrs,1.0,0)
c
c     normalize from "per cell" to "per volume"
c
      call eire_renorm(pinrec,kvols,nrs,nks,irsep,sum,0)
c
      write(6,*) 'Total recombinations per second (PIN): ',sum
c
c-----------------------------------------------------------------------
c
c     H+ (and impurities) energy source
c
      call rzero (pinqi,maxnks*maxnrs)
      call gfsub3r(ipinout,nx,ny,nxd,nyd,nfla2,maxnfla,tdummy(0,0),0)
      call deadadd(tdummy,nx,ny,nxd,nyd)
      call maptodiv(cutring,cutpt1,cutpt2,nx+2,ny+2,nxd,nyd,
     >              nfla2,maxnfla,tdummy(0,0),pinqi,maxnks,maxnrs,1.0,0)
c
c     normalize from "per cell" to "per volume"
c
      call eire_renorm(pinqi,kvols,nrs,nks,irsep,sum,0)
      write(6,*) 'Total ion energy source (PIN): ',sum
c
c-----------------------------------------------------------------------
c
c     Electron energy source
c
      call rzero (pinqe,maxnks*maxnrs)
      call gfsub3r(ipinout,nx,ny,nxd,nyd,nfla2,maxnfla,tdummy(0,0),0)
      call deadadd(tdummy,nx,ny,nxd,nyd)
      call maptodiv(cutring,cutpt1,cutpt2,nx+2,ny+2,nxd,nyd,
     >            nfla2,maxnfla,tdummy(0,0),pinqe,maxnks,maxnrs,1.0,0)
c
c     normalize from "per cell" to "per volume"
c
      call eire_renorm(pinqe,kvols,nrs,nks,irsep,sum,0)
      write(6,*) 'Total electron energy source (PIN): ',sum
c
c-----------------------------------------------------------------------
c
c     Other quantities ... not read in yet.
c
c
c     Neutral particle velocity distribution information from nimbus
c     (pinvdist)
c
c     (not read in from EIRENE yet)
c
c-----------------------------------------------------------------------
c
c     Close the files that have been used for input
c
      close (ipinout)
c
c     Calculate the DIVIMP values for the recombination.
c
c slmod begin
      call calc_divrec(rdum1)
c
c      call calc_divrec
c slmod end
c
c     Calculate the target fluxes and compare particle balances
c
      call targflux
c
c     Print out some of the Eirene data in DIVIMP format
c
      call eireprn
c
      return
      end
c
c
c
      subroutine deadadd(tdummy,nx,ny,nxd,nyd)
      implicit none
c
c     DEADADD:
c
c     this routine is used to add dead boundary rings to the EIRENE
c     output; the dead cells are set to the values in adjacent points
c
c
      integer ir,ik,nx,ny,nxd,nyd
      real tdummy(0:nxd+1,0:nyd+1)
c
c     Add knots to ring ends - set virtual cells to zero so that
c                              integrations and balances work
c
      do ir=0,ny+1
        do ik=nx+1,0,-1
           tdummy(ik+1,ir)=tdummy(ik,ir)
        enddo
c
c        tdummy(0   ,ir)=tdummy(1,ir)
c        tdummy(nx+3,ir)=tdummy(nx+2,ir)
c
        tdummy(0   ,ir)=0.0
        tdummy(nx+3,ir)=0.0
      enddo
c
c     Add virtual rings - first move data
c
      do ir=ny+1,0,-1
        do ik=0,nx+ 3
          tdummy(ik,ir+1)=tdummy(ik,ir)
        enddo
      enddo
c
c     Copy values from adjacent rings to virtual rings.
c     Actually - set to zero so that summations work correctly.
c
      do ik=0,nx+3
c
c         tdummy(ik,0   )=tdummy(ik,1)
c         tdummy(ik,ny+3)=tdummy(ik,ny+2)
c
         tdummy(ik,0   )=0.0
         tdummy(ik,ny+3)=0.0

      enddo
c
      return
      end
c
c
c
      subroutine wrteirene_97
      use mod_params
      use mod_cgeom
      use mod_comtor
      use mod_slcom
      implicit none
c     include 'params'
c     include 'cgeom'
c     include 'comtor'
c slmod begin
c     INCLUDE 'slcom'

      INTEGER ik,ir
c slmod end
c     the purpose of this routine is to bring first the arrays
c
c     bratio(ik,ir) : magnetic field ratio B_pol/B
c     knbs  (ik,ir) : density of background plasma
c     kvhs  (ik,ir) : parallel flow velocity background plasma
c     ktibs (ik,ir) : ion temperature of background plasma
c     ktebs (ik,ir) : electron temperature
c     rs    (ik,ir) : radius of bin center
c
c     to B2 standard form, i.e. to the DIVIMP arrays we add virtual
c     grid points at the end of the rings.
c     In the following call of b2wrpl, maptodiv converts the arrays
c     to the B2 ordering scheme; gfsub3w writes the arrays conforming
c     to the B2 file structure.
c
c
c     ALSO - the velocities and fluxes need to be specified as
c            CELL EDGE values - this means that extra values
c            for the ends - at the targets - must be passed to
c            the mapping routines.
c
c     the following variables are needed in addition:
c
c     from comtor:
c     cutring  = number of separatrix ring
c     cutpt1   = number of trap knots in first half
c     cutpt2   = number of trap knots in second half of trap rings
c     maxrings = total number of rings
c     maxkpts  = total number of knots on ring
c
c     from tauin1 -> raug:
c
c     nrs     = number of rings in divimp-array
c     nks(ir) = number of knots on ring ir
c     irsep   = innermost sol-ring (separatrix)
c
c     from tauin1 -> raug -> b2repl -> (gfsub3r & maptodiv):
c
c     bratio, knbs, kvhs, ktibs, ktebs
c slmod begin
      IF (eirgeom.EQ.1) THEN
        cutpt1   = ikto + 1
        cutpt2   = ikti + 1
        cutring  = irsep - 1
        maxkpts  = nks(irsep) + 2
        maxrings = irwall
      ENDIF
c
      write (6,*) 'Cut Values:',cutring,cutpt1,cutpt2,maxrings

c slmod end
      do 392 ir = irsep,nrs
         do 394 ik = nks(ir),1,-1
            bratio(ik+1,ir) = bratio(ik,ir)
            knbs  (ik+1,ir) = knbs  (ik,ir)
            kvhs  (ik+1,ir) = kvhs  (ik,ir)
            ktibs (ik+1,ir) = ktibs (ik,ir)
            ktebs (ik+1,ir) = ktebs (ik,ir)
            rs    (ik+1,ir) = rs    (ik,ir)
c
c           Also adjust kss and ksb since they are needed
c           to generate the cell boundary quantities correctly.
c
            kss (ik+1,ir)  = kss(ik,ir)
            ksb (ik+1,ir)  = ksb(ik,ir)
c
  394    continue
c
c        Setting values in the virtual points
c
c        The velocity is treated differently because it will
c        be converted to a CELL EDGE quantity by averaging the
c        contents of adjacent cells - as a result - it is
c        necessary to load the first and last cells with
c        values that will be used for the actual target velocity.
c
         bratio(1,ir) = bratio(2,ir)
         knbs  (1,ir) = knbs(2,ir)
c         kvhs  (1,ir) = kvds(idds(ir,2))*qtim*csolvf
         kvhs  (1,ir) = kvds(idds(ir,2))
         ktibs (1,ir) = ktibs(2,ir)
         ktebs (1,ir) = ktebs(2,ir)
         rs    (1,ir) = rs(2,ir)
c
         bratio(nks(ir)+2,ir) = bratio(nks(ir)+1,ir)
         knbs  (nks(ir)+2,ir) = knbs  (nks(ir)+1,ir)
c         kvhs  (nks(ir)+2,ir) = kvds(idds(ir,1))*qtim*csolvf
         kvhs  (nks(ir)+2,ir) = kvds(idds(ir,1))
         ktibs (nks(ir)+2,ir) = ktibs (nks(ir)+1,ir)
         ktebs (nks(ir)+2,ir) = ktebs (nks(ir)+1,ir)
         rs    (nks(ir)+2,ir) = rs (nks(ir)+1,ir)
c
c        Put values in kss and ksb - even if they will not be used.
c
         kss(1,ir)         = kss(2,ir)
         kss(nks(ir)+2,ir) = kss(nks(ir)+1,ir)
c
c        ksb is a 0-index array
c
         ksb(1,ir)         = ksb(0,ir)
         ksb(nks(ir)+2,ir) = ksb(nks(ir)+1,ir)
c
         nks(ir)=nks(ir)+2           ! because we have 2 more points!
c
  392 continue

*     added northopt to argument list; Krieger, IPP 6/95
      call b2wrpl_97(maxrings,maxkpts,cutring,cutpt1,cutpt2,
     >               qtim,csolvf,rizb,crmb,northopt,cioptf,cmachno)
c
c     Restore original values in these arrays
c
      do 492 ir = irsep,nrs
         do 494 ik = 1,nks(ir)-1
            bratio(ik,ir) = bratio(ik+1,ir)
            knbs  (ik,ir) = knbs  (ik+1,ir)
            kvhs  (ik,ir) = kvhs  (ik+1,ir)
            ktibs (ik,ir) = ktibs (ik+1,ir)
            ktebs (ik,ir) = ktebs (ik+1,ir)
            rs    (ik,ir) = rs    (ik+1,ir)
c
c           Readjust kss and ksb - ksb(0,...) doesn't need fixing
c                                  because it wasn't changed.
c
            kss (ik,ir) = kss(ik+1,ir)
            ksb (ik,ir) = ksb(ik+1,ir)
c
  494    continue
         nks(ir)=nks(ir)-2           ! because we deleted 2 points
  492 continue

      return
      end
c
c
c
      subroutine b2wrpl_97(mrings,mkpts,cutring,cutpt1,cutpt2,
     >                     qtim,csolvf,rizb,crmb,northopt,
     >                     cioptf,cmachno)
      use mod_params
      use mod_cgeom
      implicit none
c
c     include 'params'
c     include 'cgeom'
c
      integer mrings,mkpts,cutring,cutpt1,cutpt2,northopt
      integer cioptf
      real qtim,csolvf,rizb,crmb
      real cmachno(maxnrs,2)
c
c     B2WRPL: (B2 WRite PLasma)
c
c     this subroutine is a "mirror" of b2repl to write the arrays
c     from wrteirene to a B2 conformant background plasma file
c
c     this can of course be adapted to other B2 conformant files
c     with DIVIMP array quantities
c
c
      integer maxix,maxiy,maxis
      integer i,j,k
      integer ix,iy,is,id,ir,ik,nplasf,nx,ny,nxd,nyd
c     parameter nfla2 must be initialized before initialization of maxis
c     changed by Krieger IPP/97
      integer nfla2
      parameter (nfla2=2)
      parameter (maxix=100,maxiy=50,maxis=nfla2)
      real cs
      real ndummy(0:maxix+1,0:maxiy+1,maxis)
      real tdummy(0:maxix+1,0:maxiy+1)
      real uu(maxnks,maxnrs),fniy(maxnks,maxnrs)
c
      real sum1,sum2,sum
c
      parameter (nplasf=17)
c
c slmod begin
      CALL OutputGrid(87,'Writing EIRENE data file')
c slmod end
      rewind(nplasf)
c
      nx = mkpts-2
      ny = mrings-2
      nxd = maxix
      nyd = maxiy
c
c-----------------------------------------------------------------------
c
c     Plasma density knbs(ik,ir) (ni)
c     scaled by 1.0
c
      WRITE(nplasf,*) '[plasma density]'
      call maptob2(cutring,cutpt1,cutpt2,nx,ny,nxd,nyd,
     >                    nfla2,ndummy(0,0,1),knbs,maxnks,maxnrs,1.0,0)
      call gfsub3w(nplasf,nx,ny,nxd,nyd,nfla2,ndummy(0,0,1))

c
c-----------------------------------------------------------------------
c
c     poloidal velocity (uu)
c     scaled by qtim*csolvf to obtain m/s
c
      do i=1,maxnrs
         do j=1,maxnks
            uu(j,i)=kvhs(j,i)*bratio(j,i)

            write (6,*) 'Poloidal velocity:',j,i,kvhs(j,i),bratio(j,i),
     >                     uu(j,i)

         enddo
      enddo
c
      call maptob2(cutring,cutpt1,cutpt2,nx,ny,nxd,nyd,
     >                    nfla2,ndummy(0,0,1),uu,maxnks,maxnrs,
     >                    1.0,1)
      WRITE(nplasf,*) '[poloidal velocity]'
      call gfsub3w(nplasf,nx,ny,nxd,nyd,nfla2,ndummy(0,0,1))
c
c-----------------------------------------------------------------------
c
c     Radial velocity (vv)
c     set to zero
c
      call rzero(ndummy,(maxix+1)*(maxiy+1)*maxis)
      WRITE(nplasf,*) '[radial velocity]'
      call gfsub3w(nplasf,nx,ny,nxd,nyd,nfla2,ndummy(0,0,1))
c
c-----------------------------------------------------------------------
c
c     Electron temperature ktebs(ik,ir) (te)
c     scaled by 1.0/1.6e-19 (eV -> Joule)
c
      call maptob2(cutring,cutpt1,cutpt2,nx,ny,nxd,nyd,
     >             1,tdummy(0,0),ktebs,maxnks,maxnrs,1.0/1.6e-19,0)
      WRITE(nplasf,*) '[electron temperature]'
      call gfsub3w(nplasf,nx,ny,nxd,nyd,1,tdummy(0,0))
c
c
c-----------------------------------------------------------------------
c
c     ion temperature ktibs(ik,ir) (ti)
c     scaled by 1.0/1.6e-19 (eV -> Joule)
c
      call maptob2(cutring,cutpt1,cutpt2,nx,ny,nxd,nyd,
     >             1,tdummy(0,0),ktibs,maxnks,maxnrs,1.0/1.6e-19,0)
      WRITE(nplasf,*) '[ion temperature]'
      call gfsub3w(nplasf,nx,ny,nxd,nyd,1,tdummy(0,0))
c
c
c-----------------------------------------------------------------------
c
c     (pr)
c     set to zero
c
      call rzero(tdummy,(maxix+1)*(maxiy+1))
      WRITE(nplasf,*) '[pr - set to zero]'
      call gfsub3w(nplasf,nx,ny,nxd,nyd,1,tdummy(0,0))
c
c
c-----------------------------------------------------------------------
c
c     parallel velocity kvhs(ik,ir) (up)
c     scaled by qtim*csolvf to obtain m/s
c
      call maptob2(cutring,cutpt1,cutpt2,nx,ny,nxd,nyd,
     >                    nfla2,ndummy(0,0,1),kvhs,maxnks,maxnrs,
     >                    1.0,1)
      WRITE(nplasf,*) '[parallel velocity]'
      call gfsub3w(nplasf,nx,ny,nxd,nyd,nfla2,ndummy(0,0,1))
c
c
c-----------------------------------------------------------------------
c
c     B-field ratio bratio(ik,ir) (pit)
c     scaled by 1
c
      call maptob2(cutring,cutpt1,cutpt2,nx,ny,nxd,nyd,
     >                    nfla2,tdummy(0,0),bratio,maxnks,maxnrs,1.0,0)
      WRITE(nplasf,*) '[B-field ratio]'
      call gfsub3w(nplasf,nx,ny,nxd,nyd,1,tdummy(0,0))
c
c
c-----------------------------------------------------------------------
c
c     Particle current on target cells (fnix)
c     i.e. parallel flux * target segment length * cos(incl. angle) *
c                          2 * pi * major radius * B-field ratio
c
c     Naming inconsistency in use of fniy - one of (fnix,fniy)
c     represents the parallel flux/current direction - we believe
c     that it is fnix at the moment so the use of fniy below is
c     inconsistent - this conclusion may be backward :-)
c
      call rzero(ndummy,(maxix+1)*(maxiy+1)*maxis)
      call rzero(fniy,maxnrs*maxnks)
c
      sum1 = 0.0
      sum2 = 0.0
      sum = 0.0
c
c     These fluxes should only need to be defined on the actual
c     edges of the end cells - this would be index ik=1 and ik=nks(ir)-1
c     However, at the moment they are also being defined (as the same
c     values) for ik=2 and ik=nks(ir) - so that the value may be obtained
c     no matter how Eirene is interpreting  an array of cell edge data.
C     Note that the array is then mapped as a cell-centred one - since the
c     "Averaging" that is used to calculate the cell boundary data for
c     an entire ring are unneeded here.
c
c
      do 640 ir = irsep, nrs
c
c         First target
c
            ik = 2
            id = idds(ir,2)
c
c           Handle possibly supersonic flow from SOL 22
c
            if (cioptf.eq.22) then
                CS = 9.79E3 * SQRT (0.5*(KTEDS(ID)+KTIDS(ID))*
     >                (1.0+RIZB)/CRMB) * cmachno(ir,2)
            else
                cs = 9.79e3 * sqrt (0.5*(kteds(id)+ktids(id))*
     >                    (1.0+rizb)/crmb)
            endif
c
*           write (6,*) 'testfl-i:',ir,ik,id,dds(id),costet(id)
*           added code for northopt=2; Krieger 6/95
c
            if (northopt.eq.0) then
c             minus sign f. B2
              fniy(ik,ir) =-knds(id)*cs*dds(id)*
     >                      bratio(ik,ir)*2.*pi*rp(id)
            else if (northopt.eq.1.or.northopt.eq.2.or.northopt.eq.3)
     >                          then
c             minus sign f. B2
              fniy(ik,ir) =-knds(id)*cs*dds(id)*
     >                      bratio(ik,ir)*2.*pi*rp(id)*costet(id)
            endif
c
            fniy(ik-1,ir)=fniy(ik,ir)
c
c         Second target
c
            ik = nks(ir) -1
            id = idds(ir,1)
c
c           Handle possibly supersonic flow from SOL 22
c
            if (cioptf.eq.22) then
               CS = 9.79E3 * SQRT (0.5*(KTEDS(ID)+KTIDS(ID))*
     >                (1.0+RIZB)/CRMB) * cmachno(ir,1)
            else
               cs = 9.79e3 * sqrt (0.5*(kteds(id)+ktids(id))*
     >                    (1.0+rizb)/crmb)
            endif
c
*           write (6,*) 'testfl-o:',ir,ik,id,dds(id),costet(id)
*           added code for northopt=2; Krieger 6/95
c
            if (northopt.eq.0) then
              fniy(ik,ir) = knds(id)*cs*dds(id)*
     >                      bratio(ik,ir)*2.*pi*rp(id)
            else if (northopt.eq.1.or.northopt.eq.2.or.northopt.eq.3)
     >                     then
              fniy(ik,ir) = knds(id)*cs*dds(id)*
     >                      bratio(ik,ir)*2.*pi*rp(id)*costet(id)
            endif
c
            fniy(ik+1,ir)=fniy(ik,ir)
c
            sum1 = sum1 + fniy(2,ir)
            sum2 = sum2 + fniy(nks(ir)-1,ir)
            sum = sum - fniy(2,ir) +  fniy(nks(ir)-1,ir)

            write(6,*) 'Target Flux: ',ir,idds(ir,2),fniy(2,ir),
     >                 idds(ir,1),fniy(nks(ir)-1,ir),sum1,sum2,sum

c
  640 continue
c
c     Although fniy should be a cell-edge based quantity - it is being
c     mapped as a cell centered quantity because it is already being
c     calculated above as a cell-edge quantity and shouldn't need
c     any further mapping.
c
      call maptob2(cutring,cutpt1,cutpt2,nx,ny,nxd,nyd,
     >                    nfla2,ndummy(0,0,1),fniy,maxnks,maxnrs,1.0,0)
      WRITE(nplasf,*) '[particle current on targets]'
      call gfsub3w(nplasf,nx,ny,nxd,nyd,nfla2,ndummy(0,0,1))
c
c
c-----------------------------------------------------------------------
c
c     (fnix)
c     set to zero
c
      call rzero(ndummy,(maxix+1)*(maxiy+1)*maxis)
      WRITE(nplasf,*) '[fnix - set to zero]'
      call gfsub3w(nplasf,nx,ny,nxd,nyd,nfla2,ndummy(0,0,1))
c
c
c-----------------------------------------------------------------------
c
c     (feix)
c     set to zero

      call rzero(tdummy,(maxix+1)*(maxiy+1))
      WRITE(nplasf,*) '[feix - set to zero]'
      call gfsub3w(nplasf,nx,ny,nxd,nyd,1,tdummy(0,0))
c
c
c-----------------------------------------------------------------------
c
c     (feiy)
c     set to zero
c
      call rzero(tdummy,(maxix+1)*(maxiy+1))
      WRITE(nplasf,*) '[feiy - set to zero]'
      call gfsub3w(nplasf,nx,ny,nxd,nyd,1,tdummy(0,0))
c
c
c-----------------------------------------------------------------------
c
c     (feex)
c     set to zero
c
      call rzero(tdummy,(maxix+1)*(maxiy+1))
      WRITE(nplasf,*) '[feex - set to zero]'
      call gfsub3w(nplasf,nx,ny,nxd,nyd,1,tdummy(0,0))

c
c-----------------------------------------------------------------------
c
c     (feey)
c     set to zero
c
      call rzero(tdummy,(maxix+1)*(maxiy+1))
      WRITE(nplasf,*) '[feey - set to zero]'
      call gfsub3w(nplasf,nx,ny,nxd,nyd,1,tdummy(0,0))
c
c-----------------------------------------------------------------------
c
c     we need to close the file because the buffers must be flushed
c
      close(nplasf)
c
      return
      end
c
c
c
c
c SUBROUTINE MAPTOB2 MOVED TO EIRENE.D6A
c

c     ------------------------------------------------------------------


      subroutine gfsub3w(kard,nx,ny,ndimx,ndimy,ndims,dummy)

c     write array dummy(0:nx+1,0:ny+1,1:ndims)
c     in channel kard in format 5(e16.8)
c     pattern in file :
c     x(0,0,1) x(1,0,1) x(2,0,1) x(3,0,1) x(4,0,1)
c     x(5,0,1) x(6,0,1) x(7,0,1) x(8,0,1) x(9,0,1)
c                            .....
c     x(nx+1,0,1) x(nx+2,0,1)
c     x(0,1,1) x(1,1,1) x(2,1,1) x(3,1,1) x(4,1,1)
c     x(5,1,1) x(6,1,1) x(7,1,1) x(8,1,1) x(9,1,1)
c                            .....
c     x(nx+1,1,1) x(nx+2,1,1)
c                            .....
c                            .etc.

c     thus every block represents a ring ir

      implicit none

      integer kard,nx,ny,ndimx,ndimy,ndims,nd1,lim,is,ix,iy,iii
      real dummy(0:ndimx+1,0:ndimy+1,1:ndims)

      nd1 = nx+2
      lim = (nd1/5)*5 - 4
      do 110 is = 1,ndims
        do 110 iy = 0,ny+1
          do 100 ix = 1,lim,5
             write(kard,910) (dummy(ix-2+iii,iy,is),iii=1,5)
 100          continue
          if ((lim+4).eq.nd1) goto 110
          write(kard,910) (dummy(ix-1,iy,is),ix=lim+5,nd1)
 110  continue

      return
 910  format(5(e16.8))
      end
c
c
c
      subroutine eire_renorm(pinarray,renormf,nrs,nks,irsep,sum,renorm)
      use mod_params
      implicit none
c     include 'params'
      integer nrs,nks(nrs),irsep
      integer renorm
      real sum
      real pinarray(maxnks,maxnrs),renormf(maxnks,maxnrs)
c
c     EIRE_RENORM: Renormalizes the arrays of data read in from an
c                  Eirene transfer file.
c
c     If the argument renorm is set to 0 - then the renormalization is
c     performed - if renorm is set to 1 - then only the shift of
c     array elements is performed - the contents are not affected.
c
c

      integer ik,ir
c
      sum = 0.0
c
      do ir = 1,nrs
         do ik = 1,nks(ir)
c
c           Shift the array into the DIVIMP coordinate system -
c           the cell volume/renormalization factor should already
c           be in this system.
c
c slmod begin - bug...?
            IF (ir.GE.irsep) pinarray(ik,ir) = pinarray(ik+1,ir)
c
c            pinarray(ik,ir) = pinarray(ik+1,ir)
c slmod end
c
            if (renorm.eq.0) then
               if (renormf(ik,ir).ne.0.0) then
                  sum = sum + pinarray(ik,ir)
                  if (.not.(ir.lt.irsep.and.ik.eq.nks(ir))) then
                     pinarray(ik,ir) = pinarray(ik,ir) / renormf(ik,ir)
                  endif
               else
                  pinarray(ik,ir) = 0.0
               endif
            endif
         enddo
      enddo
c
      return
      end
c
c
c
      subroutine sum_contrib(pinarrsum,pinarradd,nrs,nks,offset)
      use mod_params
      implicit none
c     include 'params'
      integer nrs,nks(nrs),offset
c     corrected dimensions, Krieger, IPP/97
      real pinarrsum(maxnks,maxnrs),pinarradd(maxnks,maxnrs)
c
c     SUM_CONTRIB: Sums array components
c
c
      integer ik,ir
c
      do ir = 1,nrs
         do ik = 1,nks(ir) + offset
            pinarrsum(ik,ir) = pinarrsum(ik,ir)+pinarradd(ik,ir)
         enddo
      enddo
c
      return
      end

