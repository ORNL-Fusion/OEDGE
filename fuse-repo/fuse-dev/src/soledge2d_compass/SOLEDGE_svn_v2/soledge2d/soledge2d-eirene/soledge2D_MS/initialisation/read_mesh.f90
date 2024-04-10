subroutine read_mesh()
  use all_variables, only : flags
  use Mlog_message
  implicit none
  call read_mesh_all()
  call allocate_masks()
  if(flags%is_Pen) then
     call read_chis()
  end if
  call allocate_geometry()
  if(flags%is_SLAB) then
     call read_geometry_slab()
  else
     call read_geometry()
  end if
  call MD_broadcast_mesh()
  call MD_broadcast_corners_mesh()
  call MD_broadcast_corners_B()
  call MD_broadcast_geometry()
  call compute_index()
  call write_log_message(msg_mesh_loaded)
end subroutine read_mesh

