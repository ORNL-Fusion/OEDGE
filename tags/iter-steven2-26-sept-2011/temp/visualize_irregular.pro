;+
;NAME: visualize_irregular
;
;PURPOSE: This code takes data from an irregular spaced grid and interpolates it onto a regular
;        grid for volume visualization with ivolume. Creates a IMAGE PLANE on the volume with a
;        corresponding image in the same iTool, so the IMAGE PLANE can be moved on the volume to
;        show different image slices of the volume
;
;
;CALLING SEQUENCE:
;      visualize_irregular
;
;INPUTS:    x    fltarr(N)    1D array of x coordinates of vertices of tetrahedral mesh
;           y    fltarr(N)    1D array of y coordinates of vertices of tetrahedral mesh
;           z    fltarr(N)    1D array of z coordinates of vertices of tetrahedral mesh
;           data fltarr(N)    1D array of data at each vertex of tetrahedral mesh
;
;
;OUTPUTS:
;    
;
;KEYWORDS:
;
;AUTHOR:
;	rmchurch  Aug 12, 2010
;-
pro visualize_irregular ; , x,y,z,data

tstart=systime(1)

restore,'tetrahedron_test.sav'
x=tet.x & y=tet.y & z=tet.z
data=tet.n_d

x=double(x) & y=double(y) & z=double(z)
data=double(data)

;for large datasets, take off places where not defined
inds=where(data gt 0.0)
x=x[inds] & y=y[inds] & z=z[inds] & data=data[inds]

;visualize only the main puff region
inds=where(abs(z) lt 0.2)
x=x[inds] & y=y[inds] & z=z[inds] & data=data[inds]

;create regular grid to plot on
Nx=400. & Ny=400. & Nz=400.
ngrid=[Nx,Ny,Nz]

;qhull, qgrid3 to place irregular spaced grid on a regular spaced grid.
;ignore points outside original mesh
;qhull takes sometime. to save time, can do once for a particular mesh,
;then restore the triangular connection array (tr). Independent of data.
qhull,x,y,z,tr,/delaunay
datag=qgrid3(x,y,z,data,tr,dimension=ngrid,missing=!values.f_nan)

datag=alog10(datag)

;transpose for correct orientation in volume
datag=transpose(datag,[0,2,1])

;opacity values
;this will give a good 3D visualization of the puff,
;but makes the 2D image hard to see
;n=2
;max=25.
;tmp=findgen(256)^n
;opacity=tmp/max(tmp)*max

;visualize 3D volume
ivolume,datag,rgb_table0=33,/auto_render,$;opacity_table0=opacity,$
    xtitle='x',ytitle='z',ztitle='y',view_grid=[2,1],/insert_colorbar;,$
;    volume_dim=[max(x)-min(x),max(z)-min(z),max(y)-min(y)],
;    volume_loc=[min(x),min(z),min(y)],
;the above commented keywords don't work in IDL 7.0 when there is an IMAGE PLANE (prevents it from stepping)
;if you have IDL 7.1 or higher, the bug is fixed; these will set the x, y, and z coordinates correctly
;(instead of just index coordinates)

id=itgetcurrent(tool=otool)

;take off opacity colorbar
colorid=otool->findIdentifiers('*COLORBAR_1')
result=otool->removeByIdentifier(colorid)

;name colorbar
colorid=otool->findIdentifiers('*COLORBAR',/annot)
ocolor=otool->getByIdentifier(colorid)
ocolor->getProperty,Title=otitle
otitle->setProperty,strings='LOG n_D'

;create an image plane to show slices
imgplnop_id=otool->findIdentifiers('*IMAGEPLANE*',/operations)
result=otool->doAction(imgplnop_id)
planeid=otool->findIdentifiers('*IMAGE PLANE',/visualizations)
;orient the image plane to the X-Y plane
result = otool->doSetProperty(planeId, 'orientation', 2)

;show the image in the view next to the volume. place a contour on it
selectedItems = otool->getSelectedItems(count=nSelectedItems)
data = selectedItems->getParameter('IMAGEPIXELS')
imageId = IDLitSys_CreateTool('Image Tool', /view_next, view_zoom=400/Nx,opacity=findgen(256),initial_data=data, $
                               visualization_type='IMAGE') 
contId = IDLitSys_CreateTool('Contour Tool', /over, initial_data=data, $
                               visualization_type='CONTOUR') 
otool->commitactions

print,'Time Elapsed: '+strtrim(systime(1)-tstart,2)+' sec'
end