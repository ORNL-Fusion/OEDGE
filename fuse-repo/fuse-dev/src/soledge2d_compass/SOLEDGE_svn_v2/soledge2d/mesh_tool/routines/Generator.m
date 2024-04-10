function varargout = Generator(varargin)
% GENERATOR MATLAB code for Generator.fig
%      GENERATOR, by itself, creates a new GENERATOR or raises the existing
%      singleton*.
%
%      H = GENERATOR returns the handle to a new GENERATOR or the handle to
%      the existing singleton*.
%
%      GENERATOR('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GENERATOR.M with the given input arguments.
%
%      GENERATOR('Property','Value',...) creates a new GENERATOR or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before Generator_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to Generator_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help Generator

% Last Modified by GUIDE v2.5 19-Jan-2016 19:03:43

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @Generator_OpeningFcn, ...
    'gui_OutputFcn',  @Generator_OutputFcn, ...
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


% --- Executes just before Generator is made visible.
function Generator_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to Generator (see VARARGIN)

% Choose default command line output for Generator
handles.output = hObject;
% Update handles structure
guidata(hObject, handles);
generate_soledge_mesh;
generate_xz;
detect_core;
neigh_props;
update_figure1(handles);
update_figure2(handles);
update_tableMat(handles);
update_tablePump(handles);
update_tablePuff(handles)
% UIWAIT makes Generator wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = Generator_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

function update_figure1(handles)
global zones;
global Rwall;
global Zwall;
axes(handles.axes1)
% Xlims=get(gca,'xlim');
% Ylims=get(gca,'ylim');
cla reset
hold on
x=get(handles.uipanel1,'SelectedObject');
switch(get(x,'Tag'))
    case 'rb_grid'
        for k=1:zones.num
            plot(zones.zone(k).gridR,zones.zone(k).gridZ,'b-');
            plot(zones.zone(k).gridR',zones.zone(k).gridZ','b-');
            plot(zones.zone(k).gridRc,zones.zone(k).gridZc,'r.');
        end
    case 'rb_chi'
        for k=1:zones.num
            chi=[zones.zone(k).chi,zones.zone(k).chi(:,end)];
            chi=[chi;chi(end,:)];
            pcolor(zones.zone(k).gridR,zones.zone(k).gridZ,chi)
        end
    case 'rb_BR'
        for k=1:zones.num
            Br=[zones.zone(k).Br,zones.zone(k).Br(:,end)];
            Br=[Br;Br(end,:)];
            pcolor(zones.zone(k).gridR,zones.zone(k).gridZ,Br)
        end
        shading flat
    case 'rb_BZ'
        for k=1:zones.num
            Bz=[zones.zone(k).Bz,zones.zone(k).Bz(:,end)];
            Bz=[Bz;Bz(end,:)];
            pcolor(zones.zone(k).gridR,zones.zone(k).gridZ,Bz)
        end
        shading flat
    case 'rb_Bphi'
        for k=1:zones.num
            Bphi=[zones.zone(k).Bphi,zones.zone(k).Bphi(:,end)];
            Bphi=[Bphi;Bphi(end,:)];
            pcolor(zones.zone(k).gridR,zones.zone(k).gridZ,Bphi)
        end
        shading flat
    case 'rb_x'
        for k=1:zones.num
            xc=[zones.zone(k).x,zones.zone(k).x(:,end)];
            xc=[xc;xc(end,:)];
            pcolor(zones.zone(k).gridR,zones.zone(k).gridZ,xc)
        end
        shading flat
    case 'rb_z'
        for k=1:zones.num
            zc=[zones.zone(k).z,zones.zone(k).z(:,end)];
            zc=[zc;zc(end,:)];
            pcolor(zones.zone(k).gridR,zones.zone(k).gridZ,zc)
        end
        shading flat
end
plot(Rwall,Zwall,'k-','Linewidth',2);

axis equal
% if(~(sum(Xlims==[0 1])==2))
%     xlim(Xlims)
%     ylim(Ylims)
% end
% shading flat
switch(get(x,'Tag'))
    case 'rb_grid'
    case 'rb_chi'
    case 'rb_BR'
        colorbar
    case 'rb_BZ'
        colorbar
    case 'rb_Bphi'
        colorbar
    case 'rb_x'
        colorbar
        caxis([0 1])
    case 'rb_z'
        colorbar
        caxis([0 1])
end
set( gcf, 'toolbar', 'none' )
set( gcf, 'toolbar', 'figure' )



% --- Executes when selected object is changed in uipanel1.
function uipanel1_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel1
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
update_figure1(handles);


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global zones;
global megazone;
global use_penalization;
use_penalization=get(handles.checkbox3,'Value');
save_soledge_mesh_file;


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
generate_eirene_grid;
global eirene
eirene.nummaterial=1;
eirene.numpump=1;
eirene.material(1).name='W';
eirene.pump(1).name='Default';
update_figure2(handles);
update_tableMat(handles);
update_tablePump(handles);


function update_figure2(handles)
global zones;
global Rwall;
global Zwall;
global eirene;
global eireneOK
global X_points
global r2D
global z2D
global flux2D
axes(handles.axes2)
% Xlims=get(gca,'xlim');
% Ylims=get(gca,'ylim');
cla reset
hold on
Colorlistmat={'Black','Blue','Green','Cyan','Grey','Purple'};
Colors=[0 0 0;0 0 1;0 1 0;0 1 1;1 1 0 ; 0.3 0.3 0.3 ;  0.5 0 0.5];
Colorlistpump={'Red','Magenta','Yellow','Orange','Pink'};
Colorsp=[1 0 0;1 0 1;1 0.5 0 ;  1 0.75 0.8 ];
if(eireneOK)
    if(~get(handles.checkbox2,'Value'))
        plot([eirene.triR, eirene.triR(:,1)]'/100,[eirene.triZ, eirene.triZ(:,1)]'/100,'b-');
        plot(Rwall,Zwall,'k-','Linewidth',2);
    else
        list=eirene.wall.sequence;
        for n=1:length(list)
            ntri=eirene.triangle(list(n)).ntri;
            if(eirene.triangles_(ntri).BC1==1)
                plot(eirene.triR(ntri,[1,2])/100,eirene.triZ(ntri,[1,2])/100,'g-','LineWidth');
            else
                if(eirene.triangles_(ntri).BC1>1)
                    if(~eirene.wall.ispump(n))
                        plot(eirene.triR(ntri,[1,2])/100,eirene.triZ(ntri,[1,2])/100,'-','LineWidth',2,'Color',Colors(eirene.wall.typemat(n),:));
                        plot(mean(eirene.triR(ntri,[1,2])/100),mean(eirene.triZ(ntri,[1,2])/100),'.','Color',Colors(eirene.wall.typemat(n),:));
                    else
                        plot(eirene.triR(ntri,[1,2])/100,eirene.triZ(ntri,[1,2])/100,'-','LineWidth',2,'Color',Colorsp(eirene.wall.typemat(n),:));
                        plot(mean(eirene.triR(ntri,[1,2])/100),mean(eirene.triZ(ntri,[1,2])/100),'.','Color',Colorsp(eirene.wall.typemat(n),:));
                    end
                end
            end
            if(eirene.triangles_(ntri).BC2==1)
                plot(eirene.triR(ntri,[2,3])/100,eirene.triZ(ntri,[2,3])/100,'g-','LineWidth',1);
            else
                if(eirene.triangles_(ntri).BC2>1)
                    if(~eirene.wall.ispump(n))
                        plot(eirene.triR(ntri,[2,3])/100,eirene.triZ(ntri,[2,3])/100,'-','LineWidth',2,'Color',Colors(eirene.wall.typemat(n),:));
                        plot(mean(eirene.triR(ntri,[2,3])/100),mean(eirene.triZ(ntri,[2,3])/100),'.','Color',Colors(eirene.wall.typemat(n),:));
                    else
                        plot(eirene.triR(ntri,[2,3])/100,eirene.triZ(ntri,[2,3])/100,'-','LineWidth',2,'Color',Colorsp(eirene.wall.typemat(n),:));
                        plot(mean(eirene.triR(ntri,[2,3])/100),mean(eirene.triZ(ntri,[2,3])/100),'.','Color',Colorsp(eirene.wall.typemat(n),:));
                    end
                end
            end
            if(eirene.triangles_(ntri).BC3==1)
                plot(eirene.triR(ntri,[3,1])/100,eirene.triZ(ntri,[3,1])/100,'g-','LineWidth',1);
            else
                if(eirene.triangles_(ntri).BC3>1)
                    if(~eirene.wall.ispump(n))
                        plot(eirene.triR(ntri,[3,1])/100,eirene.triZ(ntri,[3,1])/100,'-','LineWidth',2,'Color',Colors(eirene.wall.typemat(n),:));
                        plot(mean(eirene.triR(ntri,[3,1])/100),mean(eirene.triZ(ntri,[3,1])/100),'.','Color',Colors(eirene.wall.typemat(n),:));
                    else
                        plot(eirene.triR(ntri,[3,1])/100,eirene.triZ(ntri,[3,1])/100,'-','LineWidth',2,'Color',Colorsp(eirene.wall.typemat(n),:));
                        plot(mean(eirene.triR(ntri,[3,1])/100),mean(eirene.triZ(ntri,[3,1])/100),'.','Color',Colorsp(eirene.wall.typemat(n),:));
                    end
                end
            end
        end
        
        for k=1:X_points.num
            contour(r2D,z2D,flux2D,[X_points.psi(k) X_points.psi(k)],'r--','LineWidth',2)
        end
    end
    
    for n=1:eirene.npuff
       plot(eirene.puff(n).R,eirene.puff(n).Z,'m.','MarkerSize',8,'Linewidth',3)
       text(eirene.puff(n).R,eirene.puff(n).Z,num2str(n),'Color','m','FontSize',18)
    end
end
xlim([min(Rwall) max(Rwall)])
ylim([min(Zwall) max(Zwall)])
axis equal


% --- Executes on button press in checkbox2.
function checkbox2_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
update_figure2(handles);
update_tableMat(handles);
update_tablePump(handles);
% Hint: get(hObject,'Value') returns toggle state of checkbox2


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global eirene
choice = questdlg('Roughly select points and press Enter when done. Pump or wall?','Choose if you want to set pump or wall',...
    'Pump','Wall','Wall');
% Handle response
global selwall
axes(handles.axes2);
switch choice
    case 'Pump'
        selwall=false;
        [x,y]=ginput();
        change_wall(x,y)
    case 'Wall'
        selwall=true;
        [x,y]=ginput();
        change_wall(x,y)
end
update_figure2(handles)


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global eirene
eirene.nummaterial=eirene.nummaterial+1;
eirene.material(eirene.nummaterial).name = inputdlg('Name for the material:');
update_tableMat(handles);


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global eirene
eirene.numpump=eirene.numpump+1;
eirene.pump(eirene.numpump).name = inputdlg('Name for the material:');
update_tablePump(handles);

function update_tableMat(handles)
global eirene
Colorlistmat={'Black','Blue','Green','Cyan','Grey','Purple'};
Colors=[0 0 0;0 0 1;0 1 0;0 1 1;1 1 0 ; 0.3 0.3 0.3 ;  0.5 0 0.5];
clear H;
for k=1:eirene.nummaterial
    if(iscellstr(eirene.material(k).name))
        H{k,1}=eirene.material(k).name{1};
    else
        H{k,1}=eirene.material(k).name;
    end
    H{k,2}=Colorlistmat{k};
end
set(handles.uitable1,'Data',H);

function update_tablePump(handles)
global eirene
Colorlistpump={'Red','Magenta','Yellow','Orange','Pink'};
Colorsp=[1 0 0;1 0 1;1 0.5 0 ;  1 0.75 0.8 ];
clear H;
for k=1:eirene.numpump
    if(iscellstr(eirene.pump(k).name))
        H{k,1}=eirene.pump(k).name{1};
    else
        H{k,1}=eirene.pump(k).name;
    end
    H{k,2}=Colorlistpump{k};
end
set(handles.uitable2,'Data',H);


% --- Executes when selected cell(s) is changed in uitable1.
function uitable1_CellSelectionCallback(hObject, eventdata, handles)
% hObject    handle to uitable1 (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) currently selecteds
% handles    structure with handles and user data (see GUIDATA)
ind=eventdata.Indices;
global eirene;
for k=1:eirene.nummaterial
    if(size(find(k==ind(:,1)))>0)
        eirene.material(k).sel=true;
    else
        eirene.material(k).sel=false;
    end
end


% --- Executes when selected cell(s) is changed in uitable2.
function uitable2_CellSelectionCallback(hObject, eventdata, handles)
% hObject    handle to uitable2 (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) currently selecteds
% handles    structure with handles and user data (see GUIDATA)
ind=eventdata.Indices;
global eirene;
for k=1:eirene.numpump
    if(size(find(k==ind(:,1)))>0)
        eirene.pump(k).sel=true;
    else
        eirene.pump(k).sel=false;
    end
end



function change_wall(x,y)
global selwall;
global eirene;
if(selwall)
    L=[];
    for k=1:eirene.nummaterial
        if(eirene.material(k).sel)
            L=[L,k];
        end
    end
    if(~isempty(L))
        mat=L(1);
        %find closest to the begining
        list=eirene.wall.sequence;
        ptR=zeros(size(list));
        ptZ=zeros(size(list));
        for n=1:length(list)
            ntri=eirene.triangle(list(n)).ntri;
            if(eirene.triangles_(ntri).BC1>1)
                ptR(n)=mean(eirene.triR(ntri,[1,2])/100);
                ptZ(n)=mean(eirene.triZ(ntri,[1,2])/100);
            end
            if(eirene.triangles_(ntri).BC2>1)
                ptR(n)=mean(eirene.triR(ntri,[2,3])/100);
                ptZ(n)=mean(eirene.triZ(ntri,[2,3])/100);
            end
            if(eirene.triangles_(ntri).BC3>1)
                ptR(n)=mean(eirene.triR(ntri,[1,3])/100);
                ptZ(n)=mean(eirene.triZ(ntri,[1,3])/100);
            end
        end
        d=sqrt((ptR-x(1)).^2+(ptZ-y(1)).^2);
        a=find(d==min(d));
        num=1:length(list);
        ptRn=circshift(ptR',-a)';
        ptZn=circshift(ptZ',-a)';
        numn=circshift(num',-a)';
        d=sqrt((ptRn-x(end)).^2+(ptZn-y(end)).^2);
        b=find(d==min(d));
        c=numn(b);
        %[a,b] or [b,a] ?
        d1=sqrt((x-ptRn(floor(b/2))).^2+(y-ptZn(floor(b/2))).^2);
        d1=min(d1);
        d2=sqrt((x-ptRn(floor((length(list)+b)/2))).^2+(y-ptZn(floor((length(list)+b)/2))).^2);
        d2=min(d2);
        if(d1<d2) %ab
            eirene.wall.typemat(numn(1:b))=mat;
            eirene.wall.ispump(numn(1:b))=false;
        else %ba
            eirene.wall.typemat(numn(b:end))=mat;
            eirene.wall.ispump(numn(b:end))=false;
        end
    end
else
    L=[];
    for k=1:eirene.numpump
        if(eirene.pump(k).sel)
            L=[L,k];
        end
    end
    if(~isempty(L))
        mat=L(1);
        %find closest to the begining
        list=eirene.wall.sequence;
        ptR=zeros(size(list));
        ptZ=zeros(size(list));
        for n=1:length(list)
            ntri=eirene.triangle(list(n)).ntri;
            if(eirene.triangles_(ntri).BC1>1)
                ptR(n)=mean(eirene.triR(ntri,[1,2])/100);
                ptZ(n)=mean(eirene.triZ(ntri,[1,2])/100);
            end
            if(eirene.triangles_(ntri).BC2>1)
                ptR(n)=mean(eirene.triR(ntri,[2,3])/100);
                ptZ(n)=mean(eirene.triZ(ntri,[2,3])/100);
            end
            if(eirene.triangles_(ntri).BC3>1)
                ptR(n)=mean(eirene.triR(ntri,[1,3])/100);
                ptZ(n)=mean(eirene.triZ(ntri,[1,3])/100);
            end
        end
        d=sqrt((ptR-x(1)).^2+(ptZ-y(1)).^2);
        a=find(d==min(d));
        num=1:length(list);
        ptRn=circshift(ptR',-a)';
        ptZn=circshift(ptZ',-a)';
        numn=circshift(num',-a)';
        d=sqrt((ptRn-x(end)).^2+(ptZn-y(end)).^2);
        b=find(d==min(d));
        c=numn(b);
        %[a,b] or [b,a] ?
        d1=sqrt((x-ptRn(floor(b/2))).^2+(y-ptZn(floor(b/2))).^2);
        d1=min(d1);
        d2=sqrt((x-ptRn(floor((length(list)+b)/2))).^2+(y-ptZn(floor((length(list)+b)/2))).^2);
        d2=min(d2);
        if(d1<d2) %ab
            eirene.wall.typemat(numn(1:b))=mat;
            eirene.wall.ispump(numn(1:b))=true;
        else %ba
            eirene.wall.typemat(numn(b:end))=mat;
            eirene.wall.ispump(numn(b:end))=true;
        end
    end
end


% --- Executes when entered data in editable cell(s) in uitable1.
function uitable1_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to uitable1 (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)
ind=eventdata.Indices;
global eirene;
for k=1:eirene.nummaterial
    if(size(find(k==ind(:,1)))>0)
        eirene.material(k).name=eventdata.NewData;
    end
end


% --- Executes when entered data in editable cell(s) in uitable2.
function uitable2_CellEditCallback(hObject, eventdata, handles)
% hObject    handle to uitable2 (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) edited
%	PreviousData: previous data for the cell(s) edited
%	EditData: string(s) entered by the user
%	NewData: EditData or its converted form set on the Data property. Empty if Data was not changed
%	Error: error string when failed to convert EditData to appropriate value for Data
% handles    structure with handles and user data (see GUIDATA)
ind=eventdata.Indices;
global eirene;
for k=1:eirene.numpump
    if(size(find(k==ind(:,1)))>0)
        eirene.pump(k).name=eventdata.NewData;
    end
end


% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global eirene
axes(handles.axes2);
[x,y]=ginput(1);
list=eirene.wall.sequence;
ptR=zeros(size(list));
ptZ=zeros(size(list));
for n=1:length(list)
    ntri=eirene.triangle(list(n)).ntri;
    if(eirene.triangles_(ntri).BC1>1)
        ptR(n)=mean(eirene.triR(ntri,[1,2])/100);
        ptZ(n)=mean(eirene.triZ(ntri,[1,2])/100);
    end
    if(eirene.triangles_(ntri).BC2>1)
        ptR(n)=mean(eirene.triR(ntri,[2,3])/100);
        ptZ(n)=mean(eirene.triZ(ntri,[2,3])/100);
    end
    if(eirene.triangles_(ntri).BC3>1)
        ptR(n)=mean(eirene.triR(ntri,[1,3])/100);
        ptZ(n)=mean(eirene.triZ(ntri,[1,3])/100);
    end
end
d=sqrt((ptR-x).^2+(ptZ-y).^2);
a=find(d==min(d));
eirene.npuff=eirene.npuff+1;
eirene.puff(eirene.npuff).ntri=eirene.triangle(list(a)).ntri;
ntri=eirene.puff(eirene.npuff).ntri;
if(eirene.triangles_(ntri).BC1>1)
    eirene.puff(eirene.npuff).side=1;
end
if(eirene.triangles_(ntri).BC2>1)
   eirene.puff(eirene.npuff).side=2;
end
if(eirene.triangles_(ntri).BC3>1)
    eirene.puff(eirene.npuff).side=3;
end
eirene.puff(eirene.npuff).R=ptR(a);
eirene.puff(eirene.npuff).Z=ptZ(a);
update_figure2(handles);
update_tablePuff(handles);


function update_tablePuff(handles)
global eirene;
if(eirene.npuff>0)
    D=zeros(eirene.npuff,2);
    for n=1:eirene.npuff
        D(n,1)=eirene.puff(n).R;
        D(n,2)=eirene.puff(n).Z;
    end
    set(handles.uitable3,'Data',D);
end


% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton8.
function pushbutton8_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global eirene
for k=1:eirene.numpump
    eirene.pump(k).number=1+k;
end
for k=1:eirene.nummaterial
    eirene.material(k).number=1+eirene.numpump+k;
end
list=eirene.wall.sequence;
for n=1:length(list)
    ntri=eirene.triangle(list(n)).ntri;
    if(eirene.triangles_(ntri).BC1>1)
        if(eirene.wall.ispump(n))
            pump=eirene.wall.typemat(n);
            eirene.triangles_(ntri).BC1=eirene.pump(pump).number;
        else
            mat=eirene.wall.typemat(n);
            eirene.triangles_(ntri).BC1=eirene.material(mat).number;
        end
    end
    if(eirene.triangles_(ntri).BC2>1)
        if(eirene.wall.ispump(n))
            pump=eirene.wall.typemat(n);
            eirene.triangles_(ntri).BC2=eirene.pump(pump).number;
        else
            mat=eirene.wall.typemat(n);
            eirene.triangles_(ntri).BC2=eirene.material(mat).number;
        end
    end
    if(eirene.triangles_(ntri).BC3>1)
        if(eirene.wall.ispump(n))
            pump=eirene.wall.typemat(n);
            eirene.triangles_(ntri).BC3=eirene.pump(pump).number;
        else
            mat=eirene.wall.typemat(n);
            eirene.triangles_(ntri).BC3=eirene.material(mat).number;
        end
    end
end
write_neighbors_file;
write_knot_files;
fid=fopen('puffs','w');
for n=1:eirene.npuff
    fprintf(fid,'%d\t',n);
    fprintf(fid,'%d\t',eirene.puff(n).ntri);
    fprintf(fid,'%d\t',eirene.puff(n).side);
    fprintf(fid,'%f\t',eirene.puff(n).R);
    fprintf(fid,'%f\n',eirene.puff(n).Z);
end
fclose(fid);
disp('saving done')


% --- Executes on button press in checkbox3.
function checkbox3_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox3
