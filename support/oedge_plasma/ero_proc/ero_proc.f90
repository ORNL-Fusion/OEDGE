program ero_proc
use ero_part_process
implicit none


! read transformation data including file names
call read_transform_data

! process the data file from ERO into one usable in DIVIMP
call process_ero_part_data

end program ero_proc

