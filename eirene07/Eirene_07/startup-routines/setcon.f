C  06.12.05:  added: au_to_cm2 (was previously defined in fpatha)
C  25.08.06:  use double precision constant for definition of PMASSA
C  22.12.06:  set up periodic table of elements
c
      SUBROUTINE SETCON
      USE PRECISION
      USE PARMMOD
      USE CCONA

      IMPLICIT NONE

      REAL(DP) :: AMUAKG

      CALL ALLOC_CCONA

C   CONSTANTS
C  ELECTRON CHARGE [EV]
      ELCHA=1.6022D-19
C  ENERGY CONVERSION
      EV_TO_J=ELCHA
      J_TO_EV=1.D0/EV_TO_J
      J_TO_ERG=1.0D7
      ERG_TO_J=1.D0/J_TO_ERG
      EV_TO_ERG=EV_TO_J*J_TO_ERG
      ERG_TO_EV=1.D0/EV_TO_ERG
      EVKEL=8.6173D-5
      EV2HZ  = 2.41798834D14
C  PLANCK CONSTANT, [ERG S]
      HPLANCK=6.6260755D-27
C  PLANCK CONSTANT, [EV S]
      HPLNK=HPLANCK*ERG_TO_EV
C  SPEED OF LIGHT, [CM/S]
      CLIGHT=2.99792458D+10
c h*c [eV*cm] FOR CONVERSION FROM ENERGY TO WAVELENGTH UNITS
      HPCL   = CLIGHT*HPLNK
C  BOHR MAGNETON, [EV/TESLA]
      MUB=5.788381804D-5
C  ATOMIC MASS UNIT [G]
      AMUA=1.6606D-24
C  PROTON, ELECTRON MASS
      PMASSA=1.0073D0
      PMASSE=5.448D-4
C  CONVERT CROSS SECTIONS FROM ATOMIC UNITS TO CM**2
      AU_TO_CM2=5.29177E-9**2
C  NUMERICAL PRECISION PARAMETERS
      EPS60=1.D-60
      EPS30=1.D-30
      EPS12=1.D-12
      EPS10=1.D-10
      EPS6=1.D-6
      EPS5=1.D-5
C
      PIA=4._DP*ATAN(1._DP)
      PI2A=2._DP*PIA
      PIHA=PIA/2._DP
      PIQU=PIA/4._DP
      PISQ=SQRT(PIA)
      PIAI=1._DP/PIA
C
      SQ2=SQRT(2._DP)
      SQ2I=1._DP/SQ2
      DEGRAD=PIA/180._DP
      RADDEG=180._DP/PIA
C  EIRENE UNITS CONVERSIONS
      AMUAKG=AMUA*1.D-3
C  V_THERMAL=CVEL2A*SQRT(T(EV)/RMASS(AMU)),  CVEL2A=0.98227E6
      CVEL2A=SQRT(1.D4*ELCHA/AMUAKG)
C  VEL(CM/S)=CVELAA*SQRT(E0(EV)/RMASS(AMU)), CVELAA=1.38912E6
      CVELAA=SQ2*CVEL2A

      CVELI2=1._DP/CVELAA/CVELAA
      EFACT=CVELI2*PMASSA
      EFCT23=EFACT*2._DP/3._DP
      HPLNK_BAR=HPLNK/PI2A
C
C  IONIZATION POTENTIAL OF NEUTRAL HYDOGEN ATOM
      EIONH=13.6_DP
C  IONIZATION POTENTIAL OF NEUTRAL HYDROGEN MOLECULE
      EIONH2=15.4_DP
C  IONIZATION POTENTIAL OF NEUTRAL HELIUM ATOM
      EIONHE=24.588_DP
C
C
C  SETUP PERIODIC TABLE OF ELEMENTS

      CALL SET_PTE_ELEMENT( 1,'Hydrogen     ','H ',1.0_DP,1.0_DP)
      CALL SET_PTE_ELEMENT( 2,'Helium       ','He',4.0_DP,2.0_DP)
      CALL SET_PTE_ELEMENT( 3,'Lithium      ','Li',6.9_DP,3.0_DP)
      CALL SET_PTE_ELEMENT( 4,'Beryllium    ','Be',9.0_DP,4.0_DP)
      CALL SET_PTE_ELEMENT( 5,'Bor          ','B ',10.8_DP,5.0_DP)
      CALL SET_PTE_ELEMENT( 6,'Carbon       ','C ',12.0_DP,6.0_DP)
      CALL SET_PTE_ELEMENT( 7,'Nitrogen     ','N ',14.0_DP,7.0_DP)
      CALL SET_PTE_ELEMENT( 8,'Oxygen       ','O ',16.0_DP,8.0_DP)
      CALL SET_PTE_ELEMENT( 9,'Fluorine     ','F ',19.0_DP,9.0_DP)
      CALL SET_PTE_ELEMENT(10,'Neon         ','Ne',20.2_DP,10.0_DP)
      CALL SET_PTE_ELEMENT(11,'Sodium       ','Na',23.0_DP,11.0_DP)
      CALL SET_PTE_ELEMENT(12,'Magnesium    ','Mg',24.2_DP,12.0_DP)
      CALL SET_PTE_ELEMENT(13,'Aluminium    ','Al',27.0_DP,13.0_DP)
      CALL SET_PTE_ELEMENT(14,'Silicium     ','Si',28.1_DP,14.0_DP)
      CALL SET_PTE_ELEMENT(15,'Phosphor     ','P ',31.0_DP,15.0_DP)
      CALL SET_PTE_ELEMENT(16,'Sulphur      ','S ',32.1_DP,16.0_DP)
      CALL SET_PTE_ELEMENT(17,'Chlorine     ','Cl',35.5_DP,17.0_DP)
      CALL SET_PTE_ELEMENT(18,'Argon        ','Ar',40.0_DP,18.0_DP)
      CALL SET_PTE_ELEMENT(19,'Potassium    ','K ',39.1_DP,19.0_DP)
      CALL SET_PTE_ELEMENT(20,'Calcium      ','Ca',40.1_DP,20.0_DP)
      CALL SET_PTE_ELEMENT(21,'Scandium     ','Sc',45.0_DP,21.0_DP)
      CALL SET_PTE_ELEMENT(22,'Titan        ','Ti',47.9_DP,22.0_DP)
      CALL SET_PTE_ELEMENT(23,'Vanadium     ','V ',50.9_DP,23.0_DP)
      CALL SET_PTE_ELEMENT(24,'Chromium     ','Cr',52.0_DP,24.0_DP)
      CALL SET_PTE_ELEMENT(25,'Manganese    ','Mn',54.9_DP,25.0_DP)
      CALL SET_PTE_ELEMENT(26,'Iron         ','Fe',55.9_DP,26.0_DP)
      CALL SET_PTE_ELEMENT(27,'Cobalt       ','Co',58.9_DP,27.0_DP)
      CALL SET_PTE_ELEMENT(28,'Nickel       ','Ni',58.7_DP,28.0_DP)
      CALL SET_PTE_ELEMENT(29,'Copper       ','Cu',63.6_DP,29.0_DP)
      CALL SET_PTE_ELEMENT(30,'Zinc         ','Zn',65.4_DP,30.0_DP)
      CALL SET_PTE_ELEMENT(31,'Gallium      ','Ga',69.7_DP,31.0_DP)
      CALL SET_PTE_ELEMENT(32,'Germanium    ','Ge',72.6_DP,32.0_DP)
      CALL SET_PTE_ELEMENT(33,'Arsenic      ','As',74.9_DP,33.0_DP)
      CALL SET_PTE_ELEMENT(34,'Selenium     ','Se',79.0_DP,34.0_DP)
      CALL SET_PTE_ELEMENT(35,'Bromine      ','Br',79.9_DP,35.0_DP)
      CALL SET_PTE_ELEMENT(36,'Krypton      ','Kr',83.8_DP,36.0_DP)
      CALL SET_PTE_ELEMENT(37,'Rubidium     ','Rb',85.5_DP,37.0_DP)
      CALL SET_PTE_ELEMENT(38,'Strontium    ','Sr',87.6_DP,38.0_DP)
      CALL SET_PTE_ELEMENT(39,'Yttrium      ','Y ',88.9_DP,39.0_DP)
      CALL SET_PTE_ELEMENT(40,'Zirconium    ','Zr',91.2_DP,40.0_DP)
      CALL SET_PTE_ELEMENT(41,'Niobium      ','Nb',92.9_DP,41.0_DP)
      CALL SET_PTE_ELEMENT(42,'Molybdenum   ','Mo',95.9_DP,42.0_DP)
      CALL SET_PTE_ELEMENT(43,'Technetium   ','Tc',97.9_DP,43.0_DP)
      CALL SET_PTE_ELEMENT(44,'Ruthenium    ','Ru',101.1_DP,44.0_DP)
      CALL SET_PTE_ELEMENT(45,'Rhodium      ','Rh',102.9_DP,45.0_DP)
      CALL SET_PTE_ELEMENT(46,'Palladium    ','Pd',106.4_DP,46.0_DP)
      CALL SET_PTE_ELEMENT(47,'Silver       ','Ag',108.9_DP,47.0_DP)
      CALL SET_PTE_ELEMENT(48,'Cadmium      ','Cd',112.4_DP,48.0_DP)
      CALL SET_PTE_ELEMENT(49,'Indium       ','In',114.8_DP,49.0_DP)
      CALL SET_PTE_ELEMENT(50,'Tin          ','Sn',118.7_DP,50.0_DP)
      CALL SET_PTE_ELEMENT(51,'Antimony     ','Sb',121.8_DP,51.0_DP)
      CALL SET_PTE_ELEMENT(52,'Tellurium    ','Te',127.6_DP,52.0_DP)
      CALL SET_PTE_ELEMENT(53,'Iodine       ','I ',126.9_DP,53.0_DP)
      CALL SET_PTE_ELEMENT(54,'Xenon        ','Xe',131.3_DP,54.0_DP)
      CALL SET_PTE_ELEMENT(55,'Caesium      ','Cs',132.9_DP,55.0_DP)
      CALL SET_PTE_ELEMENT(56,'Barium       ','Ba',137.3_DP,56.0_DP)
      CALL SET_PTE_ELEMENT(57,'Lanthanum    ','La',138.9_DP,57.0_DP)
      CALL SET_PTE_ELEMENT(58,'Cerium       ','Ce',140.1_DP,58.0_DP)
      CALL SET_PTE_ELEMENT(59,'Praseodymium ','Pr',140.9_DP,59.0_DP)
      CALL SET_PTE_ELEMENT(60,'Neodymium    ','Nd',144.2_DP,60.0_DP)
      CALL SET_PTE_ELEMENT(61,'Promethium   ','Pm',144.9_DP,61.0_DP)
      CALL SET_PTE_ELEMENT(62,'Samarium     ','Sm',150.4_DP,62.0_DP)
      CALL SET_PTE_ELEMENT(63,'Europium     ','Eu',152.0_DP,63.0_DP)
      CALL SET_PTE_ELEMENT(64,'Gadolinium   ','Gd',157.3_DP,64.0_DP)
      CALL SET_PTE_ELEMENT(65,'Terbium      ','Tb',158.9_DP,65.0_DP)
      CALL SET_PTE_ELEMENT(66,'Dysprosium   ','Dy',162.5_DP,66.0_DP)
      CALL SET_PTE_ELEMENT(67,'Holmium      ','Ho',164.9_DP,67.0_DP)
      CALL SET_PTE_ELEMENT(68,'Erbium       ','Er',167.3_DP,68.0_DP)
      CALL SET_PTE_ELEMENT(69,'Thulium      ','Tm',168.9_DP,69.0_DP)
      CALL SET_PTE_ELEMENT(70,'Ytterbium    ','Yb',173.0_DP,70.0_DP)
      CALL SET_PTE_ELEMENT(71,'Lutetium     ','Lu',175.0_DP,71.0_DP)
      CALL SET_PTE_ELEMENT(72,'Hafnium      ','Hf',178.5_DP,72.0_DP)
      CALL SET_PTE_ELEMENT(73,'Tantalum     ','Ta',181.0_DP,73.0_DP)
      CALL SET_PTE_ELEMENT(74,'Tungsten     ','W ',183.9_DP,74.0_DP)
      CALL SET_PTE_ELEMENT(75,'Rhenium      ','Re',186.2_DP,75.0_DP)
      CALL SET_PTE_ELEMENT(76,'Osmium       ','Os',190.2_DP,76.0_DP)
      CALL SET_PTE_ELEMENT(77,'Iridium      ','Ir',192.2_DP,77.0_DP)
      CALL SET_PTE_ELEMENT(78,'Platinum     ','Pt',195.1_DP,78.0_DP)
      CALL SET_PTE_ELEMENT(79,'Gold         ','Au',197.0_DP,79.0_DP)
      CALL SET_PTE_ELEMENT(80,'Mercury      ','Hg',200.6_DP,80.0_DP)
      CALL SET_PTE_ELEMENT(81,'Thallium     ','Tl',204.4_DP,81.0_DP)
      CALL SET_PTE_ELEMENT(82,'Lead         ','Pb',207.2_DP,82.0_DP)
      CALL SET_PTE_ELEMENT(83,'Bismuth      ','Bi',209.0_DP,83.0_DP)
      CALL SET_PTE_ELEMENT(84,'Polonium     ','Po',209.0_DP,84.0_DP)
      CALL SET_PTE_ELEMENT(85,'Astatine     ','At',210.0_DP,85.0_DP)
      CALL SET_PTE_ELEMENT(86,'Radon        ','Rn',222.0_DP,86.0_DP)
      CALL SET_PTE_ELEMENT(87,'Francium     ','Fr',223.0_DP,87.0_DP)
      CALL SET_PTE_ELEMENT(88,'Radium       ','Ra',226.0_DP,88.0_DP)
      CALL SET_PTE_ELEMENT(89,'Actinium     ','Ac',227.0_DP,89.0_DP)
      CALL SET_PTE_ELEMENT(90,'Thorium      ','Th',232.0_DP,90.0_DP)
      CALL SET_PTE_ELEMENT(91,'Protactinium ','Pa',231.0_DP,91.0_DP)
      CALL SET_PTE_ELEMENT(92,'Uranium      ','U ',238.0_DP,92.0_DP)
      CALL SET_PTE_ELEMENT(93,'Neptunium    ','Np',237.0_DP,93.0_DP)
      CALL SET_PTE_ELEMENT(94,'Plutonium    ','Pu',244.0_DP,94.0_DP)
      CALL SET_PTE_ELEMENT(95,'Americium    ','Am',243.0_DP,95.0_DP)
      CALL SET_PTE_ELEMENT(96,'Curium       ','Cm',247.0_DP,96.0_DP)
      CALL SET_PTE_ELEMENT(97,'Berkelium    ','Bk',247.0_DP,97.0_DP)
      CALL SET_PTE_ELEMENT(98,'Californium  ','Cf',251.0_DP,98.0_DP)
      CALL SET_PTE_ELEMENT(99,'Einsteinium  ','Es',252.0_DP,99.0_DP)
      CALL SET_PTE_ELEMENT(100,'Fermium      ','Fm',257.0_DP,100.0_DP)
      CALL SET_PTE_ELEMENT(101,'Mendelevium  ','Md',258.0_DP,101.0_DP)
      CALL SET_PTE_ELEMENT(102,'Nobelium     ','No',259.0_DP,102.0_DP)
      CALL SET_PTE_ELEMENT(103,'Lawrencium   ','Lr',262.0_DP,103.0_DP)
      CALL SET_PTE_ELEMENT(104,'Rutherfordium','Rf',261.0_DP,104.0_DP)
      CALL SET_PTE_ELEMENT(105,'Dubnium      ','Db',262.0_DP,105.0_DP)
      CALL SET_PTE_ELEMENT(106,'Seaborgium   ','Sg',266.0_DP,106.0_DP)
      CALL SET_PTE_ELEMENT(107,'Bohrium      ','Bh',264.0_DP,107.0_DP)
      CALL SET_PTE_ELEMENT(108,'Hassium      ','Hs',277.0_DP,108.0_DP)
      CALL SET_PTE_ELEMENT(109,'Meitnerium   ','Mt',268.0_DP,109.0_DP)
      CALL SET_PTE_ELEMENT(110,'Darmstadtium ','Ds',271.0_DP,110.0_DP)
      CALL SET_PTE_ELEMENT(111,'Roentgenium  ','Rg',272.0_DP,111.0_DP)

      RETURN
C
      END
