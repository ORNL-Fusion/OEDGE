function varargout = domain_definition(varargin)
% DOMAIN_DEFINITION MATLAB code for domain_definition.fig
%      DOMAIN_DEFINITION, by itself, creates a new DOMAIN_DEFINITION or raises the existing
%      singleton*.
%
%      H = DOMAIN_DEFINITION returns the handle to a new DOMAIN_DEFINITION
%      or the handle to
%      the existing singleton*.
%
%      DOMAIN_DEFINITION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DOMAIN_DEFINITION.M with the given input arguments.
%
%      DOMAIN_DEFINITION('Property','Value',...) creates a new DOMAIN_DEFINITION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before domain_definition_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to domain_definition_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help domain_definition

% Last Modified by GUIDE v2.5 14-Dec-2015 15:09:09

% Begin initialization code - DO NOT EDIT

gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @domain_definition_OpeningFcn, ...
    'gui_OutputFcn',  @domain_definition_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before domain_definition is made visible.
function domain_definition_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to domain_definition (see VARARGIN)

% Choose default command line output for domain_definition
handles.output = hObject;
global X_points;
global equ_OK;
equ_OK=false;
global wall_OK;
wall_OK=false;
global psicoreOK;
psicoreOK=false;
global psioutOK;
psioutOK=false;
global show_ortho;
show_ortho=false;
global limOK;
limOK=false;
global frontiers;
frontiers.num=0;
global eireneOK
eireneOK=false;
global zone_elements;
global zones;
global megazone;
global Pmegazone;
global extrapol_val;
global sel_mode_d;
global resize_mode;
global raise_value;
raise_value=1.e-4;
resize_mode=false;
sel_mod_d=false;
extrapol_val=0.;
global newgraph;
newgraph=true;
clc
addpath('./perso_toolbox'); %if necessary
addpath('./routines'); %if necessary
addpath('./eirene_grid'); 
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes domain_definition wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = domain_definition_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global r2D;
global z2D;
global flux2D;
global Br2D;
global Bz2D;
global Bphi2D;
global r2Do;
global z2Do;
global flux2Do;
global Br2Do;
global Bz2Do;
global Bphi2Do;
global X_points;
global equ_OK;
global rred;
global zred;
global dp_dR;
global dp_dZ;
global resize_mode;
[baseName, folder]=uigetfile();
equ_file=fullfile(folder, baseName);
load(equ_file);
r=r2D;
z=z2D;
r2D=linspace(min(min(r2D)),max(max(r2D)),400);
z2D=linspace(min(min(z2D)),max(max(z2D)),400);
[r2D,z2D]=meshgrid(r2D,z2D);
flux2D=interp2(r,z,flux2D,r2D,z2D);
Br2D=interp2(r,z,Br2D,r2D,z2D);
Bz2D=interp2(r,z,Bz2D,r2D,z2D);
Bphi2D=interp2(r,z,Bphi2D,r2D,z2D);
X_points=find_Xpoints(r2D,z2D,flux2D);
for k=1:X_points.num
    X_points.sel(k)=false;
end
surround_Xpoints;
equ_OK=true;
set(handles.text3,'String',num2str(min(min(flux2D)),'%3.2e'));
set(handles.text4,'String',num2str(max(max(flux2D)),'%3.2e'));
resize_mode=true;
r2Do=r2D;
z2Do=z2D;
flux2Do=flux2D;
Br2Do=Br2D;
Bz2Do=Bz2D;
Bphi2Do=Bphi2D;
update_figure1(hObject, eventdata, handles);
update_tableX(hObject, eventdata, handles);


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Rwall;
global Zwall;
global wall_OK;
global resize_mode;
[baseName, folder]=uigetfile();
wall_file=fullfile(folder, baseName);
load(wall_file);
wall_OK=true;
resize_mode=true;
update_figure1(hObject, eventdata, handles);



function update_figure1(hObject, eventdata, handles)
global r2D;
global z2D;
global flux2D;
global Rwall;
global Zwall;
global equ_OK;
global wall_OK;
global manual_psi_val;
global psicore;
global psicoreOK;
global psiout;
global psioutOK;
global show_ortho;
global X_points;
global limOK;
global frontiers;
global zones;
global sel_mode_d;
global X_point_sel;
global cut_sel;
global resize_mode;
global newgraph;
axes(handles.axes1)
set( gcf, 'toolbar', 'none' )
set( gcf, 'toolbar', 'figure' )
Xlims=get(gca,'xlim');
Ylims=get(gca,'ylim');
cla reset
hold on
colorlist={'b','g','c','m','y','k'};
if(equ_OK)
    for k=1:X_points.num
        contour(r2D,z2D,flux2D,[X_points.psi(k) X_points.psi(k)],'r-','LineWidth',2)
        text(X_points.R(k),X_points.Z(k),num2str(X_points.index(k)),'FontSize',18)
    end
end
if(wall_OK)
    plot(Rwall,Zwall,'k-','LineWidth',2)
end
if(get(handles.checkbox1,'Value')==true)
    contour(r2D,z2D,flux2D,[manual_psi_val manual_psi_val]);
end
if(show_ortho)
    for k=1:X_points.num
        for n=1:4
            plot(X_points.cut(k).arc(n).R,X_points.cut(k).arc(n).Z,'b-','LineWidth',2)
            plot(X_points.cut(k).arc(n).psiR,X_points.cut(k).arc(n).psiZ,'r--')
        end
    end
end
if(get(handles.checkbox2,'Value'))
    for k=1:X_points.num
        %     text(X_points.branch(k).R(1),X_points.branch(k).Z(1),'1')
        %     text(X_points.branch(k).R(2),X_points.branch(k).Z(2),'2')
        %     text(X_points.branch(k).R(3),X_points.branch(k).Z(3),'3')
        %     text(X_points.branch(k).R(4),X_points.branch(k).Z(4),'4')
        text(X_points.cut(k).R(1),X_points.cut(k).Z(1),'A')
        text(X_points.cut(k).R(2),X_points.cut(k).Z(2),'B')
        text(X_points.cut(k).R(3),X_points.cut(k).Z(3),'C')
        text(X_points.cut(k).R(4),X_points.cut(k).Z(4),'D')
    end
end
for k=1:frontiers.num
    if(frontiers.lim(k).sel)
        plot(frontiers.lim(k).R,frontiers.lim(k).Z,'c.--','LineWidth',2)
    else
        plot(frontiers.lim(k).R,frontiers.lim(k).Z,'b.--','LineWidth',2)
    end
    plot(frontiers.lim(k).P1(:,1),frontiers.lim(k).P1(:,2),'mv')
end
if(limOK)
    for k=1:zones.num
        R=(zones.zone(k).south.R(floor(end/2))+zones.zone(k).north.R(floor(end/2)))*0.5;
        Z=(zones.zone(k).south.Z(floor(end/2))+zones.zone(k).north.Z(floor(end/2)))*0.5;
        plot(zones.zone(k).south.R,zones.zone(k).south.Z,'r-','LineWidth',2)
        plot(zones.zone(k).north.R,zones.zone(k).north.Z,'b--','LineWidth',2)
        plot(zones.zone(k).east.R,zones.zone(k).east.Z,'k-','LineWidth',2)
        plot(zones.zone(k).west.R,zones.zone(k).west.Z,'g--','LineWidth',2)
        text(R,Z,num2str(k))
    end
end
if(sel_mode_d)
    plot(X_points.cut(X_point_sel).arc(cut_sel).R,...
        X_points.cut(X_point_sel).arc(cut_sel).Z,'m-','Linewidth',3)
end
axis equal
if((~newgraph)&&(~resize_mode))
    xlim(Xlims)
    ylim(Ylims)
else
    resize_mode=false;
end
newgraph=false;



% --- Executes on slider movement.
function slider2_Callback(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global flux2D;
global manual_psi_val;
val=get(hObject,'Value');
manual_psi_val=min(min(flux2D))+(max(max(flux2D))-min(min(flux2D)))*val;
set(handles.edit2,'String',num2str(manual_psi_val));
update_figure1(hObject, eventdata, handles);
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function edit2_Callback(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global manual_psi_val;
manual_psi_val=str2num(get(hObject,'String'));
update_figure1(hObject, eventdata, handles);
% Hints: get(hObject,'String') returns contents of edit2 as text
%        str2double(get(hObject,'String')) returns contents of edit2 as a double


% --- Executes during object creation, after setting all properties.
function edit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
update_figure1(hObject, eventdata, handles);
% Hint: get(hObject,'Value') returns toggle state of checkbox1



function edit3_Callback(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit3 as text
%        str2double(get(hObject,'String')) returns contents of edit3 as a double


% --- Executes during object creation, after setting all properties.
function edit3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global psicore;
global psicoreOK;
psicore = str2num(get(handles.edit3,'String'));
set(handles.text5,'String',['psi core = ',num2str(psicore)]);
psicoreOK=true;
update_figure1(hObject, eventdata, handles);



function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit4 as text
%        str2double(get(hObject,'String')) returns contents of edit4 as a double


% --- Executes during object creation, after setting all properties.
function edit4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global psiout;
global psioutOK;
psiout = str2num(get(handles.edit4,'String'));
set(handles.text6,'String',['psi out = ',num2str(psiout)]);
psioutOK=true;
update_figure1(hObject, eventdata, handles);


% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global show_ortho;
global r2D;
global z2D;
global flux2D;
global X_points;
global psicore;
global psiout;
global psidiv1;
global rred;
global zred;
global dp_dR;
global dp_dZ;
d=get(handles.uitable1,'Data');
b=zeros(1,X_points.num);
for k=1:X_points.num
    b(k)=X_points.index(k);
end
for k=1:X_points.num
    X_points.cut(find(b==k)).psilim(1)=d(k,2);
    X_points.cut(find(b==k)).psilim(2)=psiout;
    X_points.cut(find(b==k)).psilim(3)=d(k,3);
    X_points.cut(find(b==k)).psilim(4)=psiout;
end
plot_ortho_lines;
getArcLim;
show_ortho=true;
update_figure1(hObject, eventdata, handles);



function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit5 as text
%        str2double(get(hObject,'String')) returns contents of edit5 as a double


% --- Executes during object creation, after setting all properties.
function edit5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function update_tableX(hObject, eventdata, handles)
global X_points;
d=zeros(X_points.num,3);
for k=1:X_points.num
    d(X_points.index(k),1)=X_points.psi(k);
    d(X_points.index(k),2)=X_points.cut(k).psilim(1);
    d(X_points.index(k),3)=X_points.cut(k).psilim(3);
end
set(handles.uitable1,'Data',d)

function update_tableLim(hObject, eventdata, handles)
global frontiers;
d=zeros(frontiers.num,2);
for k=1:frontiers.num
    d(k,1)=frontiers.lim(k).psimin;
    d(k,2)=frontiers.lim(k).psimax;
end
set(handles.uitable2,'Data',d)

% --- Executes on button press in checkbox2.
function checkbox2_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
update_figure1(hObject, eventdata, handles);
% Hint: get(hObject,'Value') returns toggle state of checkbox2


% --- Executes on button press in pushbutton8.
function pushbutton8_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global r2D;
global z2D;
global flux2D;
global X_points;
global limOK;
global frontiers;
global zone_elements;
global psicore;
global zones;
global megazone;
global Pmegazone;
find_zone_corners;
find_segment_theta;
find_cote_psi;
remove_double_cote;
find_zones;
limOK=true;
update_figure1(hObject, eventdata, handles);
analyze_topology;
hsub=zones_prop();
WinOnTop(hsub);
mesher();

% --- Executes on button press in pushbutton9.
function pushbutton9_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global r2D;
global z2D;
global flux2D;
global psiout;
hsub=select_pointlist(handles.figure1);
WinOnTop(hsub);
uiwait(handles.figure1);
update_frontiers();
update_figure1(hObject, eventdata, handles);
update_tableLim(hObject, eventdata, handles);
% [x,y]=gradient_psi(pt,psiout,r2D,z2D,flux2D);


function update_frontiers()
global frontiers;
global r2D;
global z2D;
global flux2D;
for k=1:frontiers.num
    R1=frontiers.lim(k).R(1);
    Z1=frontiers.lim(k).Z(1);
    frontiers.lim(k).psiA=interp2(r2D,z2D,flux2D,R1,Z1);
    Rend=frontiers.lim(k).R(end);
    Zend=frontiers.lim(k).Z(end);
    frontiers.lim(k).psiB=interp2(r2D,z2D,flux2D,Rend,Zend);
end


% --- Executes on button press in pushbutton10.
function pushbutton10_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes when selected cell(s) is changed in uitable2.
function uitable2_CellSelectionCallback(hObject, eventdata, handles)
% hObject    handle to uitable2 (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) currently selecteds
% handles    structure with handles and user data (see GUIDATA)
ind=eventdata.Indices;
global frontiers;
for k=1:frontiers.num
    if(size(find(k==ind(:,1)))>0)
        frontiers.lim(k).sel=true;
    else
        frontiers.lim(k).sel=false;
    end
end
update_figure1(hObject, eventdata, handles);


% --- Executes on button press in pushbutton11.
function pushbutton11_Callback(hObject, eventdata, handles)
global extrapol_val;
global r2D;
global z2D;
global flux2D;
global Br2D;
global Bz2D;
global Bphi2D;
global resize_mode;
global X_points;
global equ_OK;
extrapolate_psi_out;
resize_mode=true;
X_points=find_Xpoints(r2D,z2D,flux2D);
for k=1:X_points.num
    X_points.sel(k)=false;
end
surround_Xpoints;
equ_OK=true;
set(handles.text3,'String',num2str(min(min(flux2D)),'%3.2e'));
set(handles.text4,'String',num2str(max(max(flux2D)),'%3.2e'));
resize_mode=true;
update_tableX(hObject, eventdata, handles);
update_figure1(hObject, eventdata, handles);
% hObject    handle to pushbutton11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on slider movement.
function slider3_Callback(hObject, eventdata, handles)
% hObject    handle to slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global extrapol_val;
extrapol_val=get(hObject,'Value');
set(handles.text10,'String',num2str(get(hObject,'Value')))
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in pushbutton12.
function pushbutton12_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global X_points;
for k=1:X_points.num
    if(X_points.sel(k))
        indrem=X_points.index(k);
        for k1=k+1:X_points.num
            X_points.sel(k1-1)=X_points.sel(k1);
            X_points.R(k1-1)=X_points.R(k1);
            X_points.Z(k1-1)=X_points.Z(k1);
            X_points.psi(k1-1)=X_points.psi(k1);
            X_points.branch(k1-1)=X_points.branch(k1);
            X_points.cut(k1-1)=X_points.cut(k1);
            X_points.index(k1-1)=X_points.index(k1);
        end
        X_points.num=X_points.num-1;
        for k1=1:X_points.num
            if(X_points.index(k1)>=indrem)
                X_points.index(k1)=X_points.index(k1)-1;
            end
        end
    end
end
update_figure1(hObject, eventdata, handles);
update_tableX(hObject, eventdata, handles);


% --- Executes when selected cell(s) is changed in uitable1.
function uitable1_CellSelectionCallback(hObject, eventdata, handles)
% hObject    handle to uitable1 (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) currently selecteds
% handles    structure with handles and user data (see GUIDATA)
ind=eventdata.Indices;
global X_points;
for k=1:X_points.num
    if(size(find(X_points.index(k)==ind(:,1)))>0)
        X_points.sel(k)=true;
    else
        X_points.sel(k)=false;
    end
end


% --- Executes on button press in pushbutton13.
function pushbutton13_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton13 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global X_points;
global equ_OK;
global wall_OK;
global psicoreOK;
global psioutOK;
global show_ortho;
global limOK;
global frontiers;
global zone_elements;
global zones;
global megazone;
global extrapol_val;
global r2D;
global z2D;
global flux2D;
global Br2D;
global Bz2D;
global Bphi2D;
global Rwall;
global Zwall;
global psicore;
global psiout;
global Pmegazone;
global eirene;
global eireneOK;
mesh_file=uiputfile();
save(mesh_file,'zones','megazone','r2D','z2D','flux2D','Br2D','Bz2D','Bphi2D','Rwall','Zwall','psicore','psiout','X_points',...
    'equ_OK','wall_OK','psicoreOK','psioutOK','show_ortho','limOK','frontiers','zone_elements','extrapol_val','Pmegazone','eirene','eireneOK')

% --- Executes on button press in pushbutton14.
function pushbutton14_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton14 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global X_points;
global equ_OK;
global wall_OK;
global psicoreOK;
global psioutOK;
global show_ortho;
global limOK;
global frontiers;
global zone_elements;
global zones;
global megazone;
global extrapol_val;
global r2D;
global z2D;
global flux2D;
global Br2D;
global Bz2D;
global Bphi2D;
global Rwall;
global Zwall;
global psicore;
global psiout;
global Pmegazone;
global eireneOK;
global eirene;
[baseName, folder]=uigetfile();
mesh_file=fullfile(folder, baseName);
load(mesh_file);
update_figure1(hObject, eventdata, handles);
update_tableX(hObject, eventdata, handles);
update_tableLim(hObject, eventdata, handles);


% --- Executes on button press in pushbutton15.
function pushbutton15_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
mesher();


% --- Executes on button press in pushbutton16.
function pushbutton16_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global sel_mode_d;
global X_point_sel;
global cut_sel;
global X_points;
global r2D;
global z2D;
global flux2D;
waitfor(msgbox('choose branch to move'))
axes(handles.axes1);
pt=ginput(1);
%find closest segment
dmin=1e10;
zone_sel=0;
side_sel=0;
for k=1:X_points.num
    for n=1:4
        vx=X_points.cut(k).arc(n).R(2:end)-X_points.cut(k).arc(n).R(1:end-1);
        vy=X_points.cut(k).arc(n).Z(2:end)-X_points.cut(k).arc(n).Z(1:end-1);
        vx1=pt(1)-X_points.cut(k).arc(n).R(1:end-1);
        vy1=pt(2)-X_points.cut(k).arc(n).Z(1:end-1);
        vx2=pt(1)-X_points.cut(k).arc(n).R(2:end);
        vy2=pt(2)-X_points.cut(k).arc(n).Z(2:end);
        d=zeros(size(vx));
        for m=1:length(vx)
            p=(vx(m)*vx1(m)+vy(m)*vy1(m))/(vx(m)^2+vy(m)^2);
            if(p>0&&p<1)
                d(m)=abs(vx(m).*vy1(m)-vy(m).*vx1(m))./sqrt(vx(m).^2+vy(m).^2);
            else
                d(m)=min(sqrt(vx1(m)^2+vy1(m)^2),sqrt(vx2(m)^2+vy2(m)^2));
            end
        end
        d=min(d);
        %     d=min(sqrt((zones.zone(k).north.R-pt(1)).^2+(zones.zone(k).north.Z-pt(2)).^2));
        if(d<dmin)
            dmin=d;
            X_point_sel=k;
            cut_sel=n;
            sel_mode_d=true;
        end
    end
end
update_figure1(hObject, eventdata, handles)
waitfor(msgbox('choose new branch end position'))
pt=ginput(1);
vec=[X_points.R(X_point_sel) X_points.Z(X_point_sel);pt(1) pt(2)];
arc=bezier(vec);
plot(arc(:,1),arc(:,2),'k-')
plot(vec(:,1),vec(:,2),'b-')
done=false;
while(~done)
    choice = questdlg('Add control point?', ...
        'X point limit mover', ...
        'Yes','No','No');
    switch choice
        case 'Yes'
            ok=false;
            while(~ok)
                nb = inputdlg('Position from X point');
                nb=str2num(nb{1});
                if((nb>=1)&&(nb<=length(vec')))
                    ok=true;
                else
                   waitfor(errordlg(['Position must be between 1 and ', num2str(length(vec'))]));
                end
            end
            pt=ginput(1);
            vec=[vec(1:nb,:);pt(1) pt(2);vec(nb+1:end,:)];
            arc=bezier(vec);
            update_figure1(hObject, eventdata, handles);
            plot(arc(:,1),arc(:,2),'k-')
            plot(vec(:,1),vec(:,2),'bo-')
            done=false;
        case 'No'
            done=true;
    end
end
c1.num=1;
c1.arc(1).x=arc(:,1);
c1.arc(1).y=arc(:,2);
c2=contour_better(r2D,z2D,flux2D,X_points.cut(X_point_sel).psilim(cut_sel));
X=intersect_contour(c1,c2);
cin.x=c1.arc(1).x';
cin.y=c1.arc(1).y';
p1.x=X.x(1);
p1.y=X.y(1);
cout=part_contour(cin,p1);
plot(cout.arc(1).x,cout.arc(1).y,'g-','Linewidth',3)
X_points.cut(X_point_sel).arc(cut_sel).R=cout.arc(1).x;
X_points.cut(X_point_sel).arc(cut_sel).Z=cout.arc(1).y;
X_points.cut(X_point_sel).arc(cut_sel).type=3;
sel_mode_d=false;
update_figure1(hObject, eventdata, handles);



function edit6_Callback(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global raise_value;
raise_value=str2num(get(hObject,'String'));
% Hints: get(hObject,'String') returns contents of edit6 as text
%        str2double(get(hObject,'String')) returns contents of edit6 as a double


% --- Executes during object creation, after setting all properties.
function edit6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
