load chi_001.txt
load chi_002.txt
load chi_003.txt
load chi_004.txt
load chi_005.txt
load chi_006.txt
load chi_007.txt
load chi_008.txt
load chi_009.txt
load chi_010.txt
load chi_011.txt
load chi_012.txt
load chi_013.txt
load chi_014.txt
load chi_015.txt
load chi_016.txt


load R_001.txt
load R_002.txt
load R_003.txt
load R_004.txt
load R_005.txt
load R_006.txt
load R_007.txt
load R_008.txt
load R_009.txt
load R_010.txt
load R_011.txt
load R_012.txt
load R_013.txt
load R_014.txt
load R_015.txt
load R_016.txt

load Z_001.txt
load Z_002.txt
load Z_003.txt
load Z_004.txt
load Z_005.txt
load Z_006.txt
load Z_007.txt
load Z_008.txt
load Z_009.txt
load Z_010.txt
load Z_011.txt
load Z_012.txt
load Z_013.txt
load Z_014.txt
load Z_015.txt
load Z_016.txt

load ../Mesh/meshx_001.txt
load ../Mesh/meshx_002.txt
load ../Mesh/meshx_003.txt
load ../Mesh/meshx_004.txt
load ../Mesh/meshx_005.txt
load ../Mesh/meshx_006.txt
load ../Mesh/meshx_007.txt
load ../Mesh/meshx_008.txt
load ../Mesh/meshx_009.txt
load ../Mesh/meshx_010.txt
load ../Mesh/meshx_011.txt
load ../Mesh/meshx_012.txt
load ../Mesh/meshx_013.txt
load ../Mesh/meshx_014.txt
load ../Mesh/meshx_015.txt
load ../Mesh/meshx_016.txt

load ../Mesh/meshz_001.txt
load ../Mesh/meshz_002.txt
load ../Mesh/meshz_003.txt
load ../Mesh/meshz_004.txt
load ../Mesh/meshz_005.txt
load ../Mesh/meshz_006.txt
load ../Mesh/meshz_007.txt
load ../Mesh/meshz_008.txt
load ../Mesh/meshz_009.txt
load ../Mesh/meshz_010.txt
load ../Mesh/meshz_011.txt
load ../Mesh/meshz_012.txt
load ../Mesh/meshz_013.txt
load ../Mesh/meshz_014.txt
load ../Mesh/meshz_015.txt
load ../Mesh/meshz_016.txt

figure(1)
hold on
pcolor(R_001,Z_001,chi_001);
pcolor(R_002,Z_002,chi_002);
pcolor(R_003,Z_003,chi_003);
pcolor(R_004,Z_004,chi_004);
pcolor(R_005,Z_005,chi_005);
pcolor(R_006,Z_006,chi_006);
pcolor(R_007,Z_007,chi_007);
pcolor(R_008,Z_008,chi_008);
pcolor(R_009,Z_009,chi_009);
pcolor(R_010,Z_010,chi_010);
pcolor(R_011,Z_011,chi_011);
pcolor(R_012,Z_012,chi_012);
pcolor(R_013,Z_013,chi_013);
pcolor(R_014,Z_014,chi_014);
pcolor(R_015,Z_015,chi_015);
pcolor(R_016,Z_016,chi_016);
axis equal

figure(2)
hold on
pcolor(meshz_001,meshx_001,chi_001);
pcolor(meshz_002,meshx_002,chi_002);
pcolor(meshz_003,meshx_003,chi_003);
pcolor(meshz_004,meshx_004,chi_004);
pcolor(meshz_005,meshx_005,chi_005);
pcolor(meshz_006,meshx_006,chi_006);
pcolor(meshz_007,meshx_007,chi_007);
pcolor(meshz_008,meshx_008,chi_008);
pcolor(meshz_009,meshx_009,chi_009);
pcolor(meshz_010,meshx_010,chi_010);
pcolor(meshz_011,meshx_011,chi_011);
pcolor(meshz_012,meshx_012,chi_012);
pcolor(meshz_013,meshx_013,chi_013);
pcolor(meshz_014,meshx_014,chi_014);
pcolor(meshz_015,meshx_015,chi_015);
pcolor(meshz_016,meshx_016,chi_016);
axis equal


