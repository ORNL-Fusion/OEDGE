function varargout = custom_mesh_interface(varargin)
% CUSTOM_MESH_INTERFACE MATLAB code for custom_mesh_interface.fig
%      CUSTOM_MESH_INTERFACE, by itself, creates a new CUSTOM_MESH_INTERFACE or raises the existing
%      singleton*.
%
%      H = CUSTOM_MESH_INTERFACE returns the handle to a new CUSTOM_MESH_INTERFACE or the handle to
%      the existing singleton*.
%
%      CUSTOM_MESH_INTERFACE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in CUSTOM_MESH_INTERFACE.M with the given input arguments.
%
%      CUSTOM_MESH_INTERFACE('Property','Value',...) creates a new CUSTOM_MESH_INTERFACE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before custom_mesh_interface_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to custom_mesh_interface_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help custom_mesh_interface

% Last Modified by GUIDE v2.5 01-Dec-2015 10:25:37

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @custom_mesh_interface_OpeningFcn, ...
    'gui_OutputFcn',  @custom_mesh_interface_OutputFcn, ...
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


% --- Executes just before custom_mesh_interface is made visible.
function custom_mesh_interface_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to custom_mesh_interface (see VARARGIN)

% Choose default command line output for custom_mesh_interface
handles.output = hObject;
global npt_MI;
npt_MI=varargin{1};
global R_MI;
R_MI=varargin{2};
global Z_MI;
global dist_MI;
Z_MI=varargin{3};
global ctr_MI;
ctr_MI=[[1;npt_MI],[0;1]];
dist_MI=zeros(size(R_MI));
for k=2:length(R_MI)
    dist_MI(k)=dist_MI(k-1)+sqrt((R_MI(k)-R_MI(k-1))^2+(Z_MI(k)-Z_MI(k-1))^2);
end
dist_MI=dist_MI/dist_MI(end);
compute_meshpoint();
update_figure1(handles);
update_figure2(handles)
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes custom_mesh_interface wait for user response (see UIRESUME)
% uiwait(handles.figure1);

function compute_meshpoint()
global d_MI
global dist_MI
global R_MI
global Z_MI
global Rp_MI
global Zp_MI
global bez_MI
global ctr_MI
global npt_MI
bez_MI=bezier(ctr_MI);
d=zeros(1,npt_MI);
d(1)=0;
d(end)=1;
for k=2:npt_MI-1
    d(k)=interp1(bez_MI(:,1),bez_MI(:,2),k);
end
d_MI=d;
Rp_MI=zeros(1,length(d_MI));
Zp_MI=zeros(1,length(d_MI));
Rp_MI(1)=R_MI(1);
Rp_MI(end)=R_MI(end);
Zp_MI(1)=Z_MI(1);
Zp_MI(end)=Z_MI(end);
for k=2:length(d_MI)-1
    Rp_MI(k)=interp1(dist_MI,R_MI,d_MI(k));
    Zp_MI(k)=interp1(dist_MI,Z_MI,d_MI(k));
end


function update_figure1(handles)
global npt_MI
global ctr_MI
global bez_MI
global d_MI
axes(handles.axes1)
cla
hold on
plot(bez_MI(:,1),bez_MI(:,2),'k-')
plot(1:npt_MI,d_MI,'ro')
plot(ctr_MI(:,1),ctr_MI(:,2),'bs-')
xlim([1 npt_MI])
ylim([0 1])
% shading flat
set( gcf, 'toolbar', 'none' )
set( gcf, 'toolbar', 'figure' )

function update_figure2(handles)
global npt_MI
global ctr_MI
axes(handles.axes2)
Xlims=get(gca,'xlim');
Ylims=get(gca,'ylim');
cla reset
hold on
global Rwall
global Zwall
global r2D
global z2D
global flux2D
global X_points
global R_MI
global Z_MI
global Rp_MI
global Zp_MI
plot(Rwall,Zwall,'k-','Linewidth',2)
for k=1:X_points.num
    contour(r2D,z2D,flux2D,[X_points.psi(k) X_points.psi(k)],'r-')
end
plot(R_MI,Z_MI,'g-');
text(R_MI(1),Z_MI(1),'0','Fontsize',16);
plot(Rp_MI,Zp_MI,'k.');
text(R_MI(end),Z_MI(end),'1','Fontsize',16);
if(get(handles.checkbox1,'Value'))
    for k=1:npt_MI
        psi=interp2(r2D,z2D,flux2D,Rp_MI(k),Zp_MI(k));
        contour(r2D,z2D,flux2D,[psi psi],'b-')
    end
end
axis equal
if(~(sum(Xlims==[0 1])==2))
    xlim(Xlims)
    ylim(Ylims)
end
% shading flat
% set( gcf, 'toolbar', 'none' )
% set( gcf, 'toolbar', 'figure' )

% --- Outputs from this function are returned to the command line.
function varargout = custom_mesh_interface_OutputFcn(hObject, eventdata, handles)
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
global ctr_MI
axes(handles.axes1)
[x,y]=ginput(1)
newx=[ctr_MI(:,1);x];
newy=[ctr_MI(:,2);y];
[a,b]=sort(newx);
newx=newx(b);
newy=newy(b);
ctr_MI=[newx,newy];
compute_meshpoint();
update_figure1(handles);
update_figure2(handles);



% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
update_figure2(handles)
% Hint: get(hObject,'Value') returns toggle state of checkbox1


% --- Executes on button press in checkbox2.
function checkbox2_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox2


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
delete(handles.figure1);

% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global ctr_MI
global npt_MI
axes(handles.axes1)
[x,y]=ginput(1)
d=zeros(1,length(ctr_MI));
for k=1:length(ctr_MI)
    d(k)=sqrt((ctr_MI(k,1)-x)^2+(ctr_MI(k,2)-y)^2);
end
a=find(d==min(d));
if((a~=1)&&(a~=npt_MI))
    hold on
    plot(ctr_MI(a,1),ctr_MI(a,2),'ms','Linewidth',3)
    [x,y]=ginput(1)
    ctr_MI(a,1)=x;
    ctr_MI(a,2)=y;
    newx=[ctr_MI(:,1)];
    newy=[ctr_MI(:,2)];
    [a,b]=sort(newx);
    newx=newx(b);
    newy=newy(b);
    ctr_MI=[newx,newy];
end
compute_meshpoint();
update_figure1(handles);
update_figure2(handles);

% --- Executes during object deletion, before destroying properties.
function figure1_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hint: delete(hObject) closes the figure
delete(hObject);
