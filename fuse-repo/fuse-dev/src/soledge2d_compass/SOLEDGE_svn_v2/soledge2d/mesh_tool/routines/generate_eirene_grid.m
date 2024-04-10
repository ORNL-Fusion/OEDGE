global eirene;
global eireneOK;
addpath('./eirene_grid');

eirene.npuff=0;

move_wall_inward;

Nsteps=13;
heir=waitbar(0,'Generating eirene knots','Name','Eirene grid generation');
generate_eirene_knots;
direct_or_indirect;
disp('knots generated in the plasma')

waitbar(1/Nsteps,heir,'Detecting wall intersections','Name','Eirene grid generation');
detect_cut;
detect_aligned_mesh;


waitbar(2/Nsteps,heir,'Removing knots in the wall','Name','Eirene grid generation');
remove_knots_in_solid;
waitbar(3/Nsteps,heir,'Generating knots on the wall','Name','Eirene grid generation');
generate_knots_on_wall;

finalization_remove;
disp('knots generated on the wall')

if(eirene.direct==0) %indirect triangles
    waitbar(4/Nsteps,heir,'Generating triangles in the plasma','Name','Eirene grid generation');
    generate_indirect_triangles_in_the_plasma;
    disp('triangles generated in the plasma')
    waitbar(5/Nsteps,heir,'Generating triangles on the wall','Name','Eirene grid generation');
    generate_indirect_triangles_on_wall;
    disp('triangles generated on the wall')
else %direct triangles
    waitbar(4/Nsteps,heir,'Generating triangles in the plasma','Name','Eirene grid generation');
    generate_direct_triangles_in_the_plasma;
    disp('triangles generated in the plasma')
    waitbar(5/Nsteps,heir,'Generating triangles on the wall','Name','Eirene grid generation');
    generate_direct_triangles_on_wall;
    disp('triangles generated on the wall')
end

waitbar(6/Nsteps,heir,'Finding triangle neighbors','Name','Eirene grid generation');
find_triangle_area;
find_neighbors;
waitbar(7/Nsteps,heir,'Determining triangle connexity','Name','Eirene grid generation');
find_connexity;
remove_non_connex_triangles;

eirene.ntriangles_=ntriangles_;
eirene.triangles_=triangles_;
waitbar(8/Nsteps,heir,'Writing eirene files','Name','Eirene grid generation');
write_npco_char_file;
write_elemente_file;
write_neighbors_file;

waitbar(9/Nsteps,heir,'Generating soledge-eirene interpolations','Name','Eirene grid generation');
generate_data_for_interpolation_pass1_2;
generate_data_for_interpolation_pass3_4;
check_interpolation_triangles;

waitbar(10/Nsteps,heir,'Finding wall triangle sequence','Name','Eirene grid generation');
identify_wall_triangles;
generate_interpolation_flux_on_wall;
find_sequence;

waitbar(11/Nsteps,heir,'Determining parallel wall flux interpolations','Name','Eirene grid generation');
assign_triangles_to_trans_mesh2;
compute_weights;
remapping_para;

waitbar(12/Nsteps,heir,'Determining perpendicular wall flux interpolations','Name','Eirene grid generation');
assign_triangles_to_trans_mesh2_perp;
compute_weights_perp;
remapping_perp;



eirene.ntriangles_=ntriangles_;
eirene.nknots_=nknots_;
eirene.R_=R_;
eirene.Z_=Z_;
eirene.knots_interp=knots_interp;
eirene.triangle=triangle;
eirene.triangles_=triangles_;
eirene.triR=R;
eirene.triZ=Z;


eirene.wall.nummaterial=1;
eirene.wall.material(1)='W';
eirene.wall.numpump=0;


waitbar(13/Nsteps,heir,'Writing triangle.h5 file','Name','Eirene grid generation');
write_knot_files;

disp('Eirene grid generated')
eireneOK=true;

delete(heir);