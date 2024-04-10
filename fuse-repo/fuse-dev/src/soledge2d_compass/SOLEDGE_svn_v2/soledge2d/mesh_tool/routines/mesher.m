function varargout = mesher(varargin)
% MESHER M-file for mesher.fig
%      MESHER, by itself, creates a new MESHER or raises the existing
%      singleton*.
%
%      H = MESHER returns the handle to a new MESHER or the handle to
%      the existing singleton*.
%
%      MESHER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in MESHER.M with the given input arguments.
%
%      MESHER('Property','Value',...) creates a new MESHER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before mesher_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to mesher_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help mesher

% Last Modified by GUIDE v2.5 18-Mar-2016 10:03:21

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @mesher_OpeningFcn, ...
    'gui_OutputFcn',  @mesher_OutputFcn, ...
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


% --- Executes just before mesher is made visible.
function mesher_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to mesher (see VARARGIN)

% Choose default command line output for mesher
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
global sel_mode;
global modeAjust;
global align_psimode;
global alpha_ad;
alpha_ad=0.1;
global range_ad;
range_ad=0.2;
global length_ad;
length_ad=1000;
global orthoptim;
orthoptim=0.5;
global dmin_ad;
dmin_ad=2e-3;
global target_ad;
target_ad=0.75;
global drmin_ad;
drmin_ad=0.1;
align_psimode=false;
global newgraph2;
newgraph2=true;
modeAdjust=0;
sel_mode=false;
update_figure1(hObject, eventdata, handles);
% UIWAIT makes mesher wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = mesher_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function update_figure1(hObject, eventdata, handles)
global zones;
global Rwall;
global Zwall;
global equ_OK;
global wall_OK;
global zone_sel;
global side_sel;
global sel_mode;
global sub_sel;
global megazone;
global Pmegazone;
global r2D;
global z2D;
global flux2D;
global newgraph2;
axes(handles.axes1)
Xlims=get(gca,'xlim');
Ylims=get(gca,'ylim');
set( gcf, 'toolbar', 'none' )
set( gcf, 'toolbar', 'figure' )
cla reset
hold on
if(wall_OK)
    plot(Rwall,Zwall,'k-','LineWidth',2)
end
if((get(handles.checkbox1,'Value')))
    if((~get(handles.checkbox2,'Value'))&&(~get(handles.checkbox6,'Value')))
        for k=1:zones.num
            if(zones.zone(k).meshortho)
                plot(zones.zone(k).gridR,zones.zone(k).gridZ,'g-');
                plot(zones.zone(k).gridR',zones.zone(k).gridZ','g-');
            else
                plot(zones.zone(k).gridR,zones.zone(k).gridZ,'-','Color',[1.,0.8,0.2]);
                plot(zones.zone(k).gridR',zones.zone(k).gridZ','-','Color',[1.,0.8,0.2]);
            end
        end
        set(handles.text16,'String','Info:');
    end
    if((get(handles.checkbox2,'Value'))&&(~get(handles.checkbox6,'Value')))
        qualmin=1e10;
        for k=1:zones.num
            [m,p]=size(zones.zone(k).gridR);
            if(m~=0)
                vx=zones.zone(k).gridR(1,2)-zones.zone(k).gridR(1,1);
                vy=zones.zone(k).gridZ(1,2)-zones.zone(k).gridZ(1,1);
                vx1=zones.zone(k).gridR(2,1)-zones.zone(k).gridR(1,1);
                vy1=zones.zone(k).gridZ(2,1)-zones.zone(k).gridZ(1,1);
                qual=(vx1*vy-vy1*vx)/(sqrt(vx^2+vy^2)*sqrt(vx1^2+vy1^2));
                sens=sign(qual);
            else
                sens=1;
            end
            qual=zeros(m,p);
            for i=1:m-1
                for j=1:p-1
                    vx=zones.zone(k).gridR(i,j+1)-zones.zone(k).gridR(i,j);
                    vy=zones.zone(k).gridZ(i,j+1)-zones.zone(k).gridZ(i,j);
                    vx1=zones.zone(k).gridR(i+1,j)-zones.zone(k).gridR(i,j);
                    vy1=zones.zone(k).gridZ(i+1,j)-zones.zone(k).gridZ(i,j);
                    qual(i,j)=(vx1*vy-vy1*vx)/(sqrt(vx^2+vy^2)*sqrt(vx1^2+vy1^2))*sens;
                    if(qual(i,j)<qualmin)
                        qualmin=qual(i,j);
                    end
                    %                     patch([zones.zone(k).gridR(i,j), zones.zone(k).gridR(i+1,j), ...
                    %                         zones.zone(k).gridR(i+1,j+1), zones.zone(k).gridR(i,j+1)],...
                    %                         [zones.zone(k).gridZ(i,j), zones.zone(k).gridZ(i+1,j), ...
                    %                         zones.zone(k).gridZ(i+1,j+1), zones.zone(k).gridZ(i,j+1)], ...
                    %                         [qual, qual, qual, qual]*sens)
                end
            end
            pcolor(zones.zone(k).gridR,zones.zone(k).gridZ,qual)
        end
        colorbar
        caxis([0 1])
        shading flat
        load('redgreencmap.mat');
        colormap(cmap);
        set(handles.text16,'String',['Info: Minimun quality factor = ',num2str(qualmin)]);
    end
    if((get(handles.checkbox6,'Value')))
        qualmin=1e10;
        for k=1:zones.num
            [m,p]=size(zones.zone(k).gridR);
            qual=zeros(m,p);
            for i=1:m-1
                for j=1:p-1
                    vx=zones.zone(k).gridR(i,j+1)-zones.zone(k).gridR(i,j);
                    vy=zones.zone(k).gridZ(i,j+1)-zones.zone(k).gridZ(i,j);
                    vx1=zones.zone(k).gridR(i+1,j)-zones.zone(k).gridR(i,j);
                    vy1=zones.zone(k).gridZ(i+1,j)-zones.zone(k).gridZ(i,j);
                    qualr1=abs((vx1*vy-vy1*vx)/(sqrt(vx^2+vy^2)))*1000;
                    qualz1=abs((vx1*vy-vy1*vx)/(sqrt(vx1^2+vy1^2)))*1000;
                    qual(i,j)=min(qualr1,qualz1);
                    if(qual(i,j)<qualmin)
                        qualmin=qual(i,j);
                    end
                    %                     patch([zones.zone(k).gridR(i,j), zones.zone(k).gridR(i+1,j), ...
                    %                         zones.zone(k).gridR(i+1,j+1), zones.zone(k).gridR(i,j+1)],...
                    %                         [zones.zone(k).gridZ(i,j), zones.zone(k).gridZ(i+1,j), ...
                    %                         zones.zone(k).gridZ(i+1,j+1), zones.zone(k).gridZ(i,j+1)], ...
                    %                         [qual, qual, qual, qual])
                end
            end
            pcolor(zones.zone(k).gridR,zones.zone(k).gridZ,qual)
        end
        colorbar
        caxis([0 1])
        shading flat
        colormap('hot');
        set(handles.text16,'String',['Info: Minimun distance = ',num2str(qualmin),' mm.']);
    end
else
    set(handles.text16,'String','Info:');
    for k=1:zones.num
        R=(zones.zone(k).south.R(floor(end/2))+zones.zone(k).north.R(floor(end/2)))*0.5;
        Z=(zones.zone(k).south.Z(floor(end/2))+zones.zone(k).north.Z(floor(end/2)))*0.5;
        if(zones.zone(k).south.ismeshed)
            if(zones.zone(k).southaligned)
                plot(zones.zone(k).south.R,zones.zone(k).south.Z,'g-','LineWidth',3)
            else
                plot(zones.zone(k).south.R,zones.zone(k).south.Z,'g-','LineWidth',1)
            end
        else
            if(zones.zone(k).southaligned)
                plot(zones.zone(k).south.R,zones.zone(k).south.Z,'r-','LineWidth',3)
            else
                plot(zones.zone(k).south.R,zones.zone(k).south.Z,'r-','LineWidth',1)
            end
        end
        if(zones.zone(k).north.ismeshed)
            if(zones.zone(k).northaligned)
                plot(zones.zone(k).north.R,zones.zone(k).north.Z,'g-','LineWidth',3)
            else
                plot(zones.zone(k).north.R,zones.zone(k).north.Z,'g-','LineWidth',1)
            end
        else
            if(zones.zone(k).northaligned)
                plot(zones.zone(k).north.R,zones.zone(k).north.Z,'r-','LineWidth',3)
            else
                plot(zones.zone(k).north.R,zones.zone(k).north.Z,'r-','LineWidth',1)
            end
        end
        if(zones.zone(k).east.ismeshed)
            plot(zones.zone(k).east.R,zones.zone(k).east.Z,'g-','LineWidth',1)
        else
            plot(zones.zone(k).east.R,zones.zone(k).east.Z,'r-','LineWidth',1)
        end
        if(zones.zone(k).west.ismeshed)
            plot(zones.zone(k).west.R,zones.zone(k).west.Z,'g-','LineWidth',1)
        else
            plot(zones.zone(k).west.R,zones.zone(k).west.Z,'r-','LineWidth',1)
        end
        text(R,Z,num2str(k))
    end
    if(sel_mode)
        switch(side_sel)
            case 1
                if(zones.zone(zone_sel).northaligned==false)
                    plot(zones.zone(zone_sel).north.R,zones.zone(zone_sel).north.Z,'b--','LineWidth',2)
                else
                    psimin=Pmegazone.mz(zones.zone(zone_sel).pmz).align_psimin;
                    psimax=Pmegazone.mz(zones.zone(zone_sel).pmz).align_psimax;
                    contour(r2D,z2D,flux2D,[psimin psimin],'m--');
                    contour(r2D,z2D,flux2D,[psimax psimax],'m--');
                    if(~isempty(sub_sel))
                        plot(zones.zone(zone_sel).subNorth(sub_sel).R,zones.zone(zone_sel).subNorth(sub_sel).Z,'b--','LineWidth',2)
                    end
                end
            case 2
                if(zones.zone(zone_sel).southaligned==false)
                    plot(zones.zone(zone_sel).south.R,zones.zone(zone_sel).south.Z,'b--','LineWidth',2)
                else
                    psimin=Pmegazone.mz(zones.zone(zone_sel).pmz).align_psimin;
                    psimax=Pmegazone.mz(zones.zone(zone_sel).pmz).align_psimax;
                    contour(r2D,z2D,flux2D,[psimin psimin],'m--');
                    contour(r2D,z2D,flux2D,[psimax psimax],'m--');
                    if(~isempty(sub_sel))
                        plot(zones.zone(zone_sel).subSouth(sub_sel).R,zones.zone(zone_sel).subSouth(sub_sel).Z,'b--','LineWidth',2)
                    end
                end
            case 3
                plot(zones.zone(zone_sel).east.R,zones.zone(zone_sel).east.Z,'b--','LineWidth',2)
            case 4
                plot(zones.zone(zone_sel).west.R,zones.zone(zone_sel).west.Z,'b--','LineWidth',2)
            otherwise
        end
    end
    for k=1:megazone.num
        if(megazone.mz(k).ismeshed)
            plot(megazone.mz(k).refpoints.R,megazone.mz(k).refpoints.Z,'k.')
        end
    end
    for k=1:Pmegazone.num
        if(~isempty(Pmegazone.mz(k).refpoints))
            if(~isempty(Pmegazone.mz(k).refpoints.R))
                plot(Pmegazone.mz(k).refpoints.R,Pmegazone.mz(k).refpoints.Z,'k.')
            end
        end
    end
    global align_psimode;
    if(align_psimode)
        contour(r2D,z2D,flux2D,[str2num(get(handles.edit4,'String')) str2num(get(handles.edit4,'String'))],'m--');
        contour(r2D,z2D,flux2D,[str2num(get(handles.edit5,'String')) str2num(get(handles.edit5,'String'))],'m--');
    end
end
axis equal
if(~(newgraph2))
    xlim(Xlims)
    ylim(Ylims)
end
% shading flat
newgraph2=false;



% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global zones;
global side_sel;
global zone_sel;
global sub_sel;
global sel_mode;
axes(handles.axes1);
pt=ginput(1);
%find closest edge
dmin=1e10;
zone_sel=0;
side_sel=0;
for k=1:zones.num
    %north
    if(zones.zone(k).northaligned==false)
        vx=zones.zone(k).north.R(2:end)-zones.zone(k).north.R(1:end-1);
        vy=zones.zone(k).north.Z(2:end)-zones.zone(k).north.Z(1:end-1);
        vx1=pt(1)-zones.zone(k).north.R(1:end-1);
        vy1=pt(2)-zones.zone(k).north.Z(1:end-1);
        vx2=pt(1)-zones.zone(k).north.R(2:end);
        vy2=pt(2)-zones.zone(k).north.Z(2:end);
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
            zone_sel=k;
            side_sel=1;
            sel_mode=true;
        end
    else
        for n=1:2
            vx=zones.zone(k).subNorth(n).R(2:end)-zones.zone(k).subNorth(n).R(1:end-1);
            vy=zones.zone(k).subNorth(n).Z(2:end)-zones.zone(k).subNorth(n).Z(1:end-1);
            vx1=pt(1)-zones.zone(k).subNorth(n).R(1:end-1);
            vy1=pt(2)-zones.zone(k).subNorth(n).Z(1:end-1);
            vx2=pt(1)-zones.zone(k).subNorth(n).R(2:end);
            vy2=pt(2)-zones.zone(k).subNorth(n).Z(2:end);
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
                zone_sel=k;
                side_sel=1;
                sub_sel=n;
                sel_mode=true;
            end
        end
    end
    %south
    if(zones.zone(k).southaligned==false)
        vx=zones.zone(k).south.R(2:end)-zones.zone(k).south.R(1:end-1);
        vy=zones.zone(k).south.Z(2:end)-zones.zone(k).south.Z(1:end-1);
        vx1=pt(1)-zones.zone(k).south.R(1:end-1);
        vy1=pt(2)-zones.zone(k).south.Z(1:end-1);
        vx2=pt(1)-zones.zone(k).south.R(2:end);
        vy2=pt(2)-zones.zone(k).south.Z(2:end);
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
        %     d=min(sqrt((zones.zone(k).south.R-pt(1)).^2+(zones.zone(k).south.Z-pt(2)).^2));
        if(d<dmin)
            dmin=d;
            zone_sel=k;
            side_sel=2;
            sel_mode=true;
        end
    else
        for n=1:2
            vx=zones.zone(k).subSouth(n).R(2:end)-zones.zone(k).subSouth(n).R(1:end-1);
            vy=zones.zone(k).subSouth(n).Z(2:end)-zones.zone(k).subSouth(n).Z(1:end-1);
            vx1=pt(1)-zones.zone(k).subSouth(n).R(1:end-1);
            vy1=pt(2)-zones.zone(k).subSouth(n).Z(1:end-1);
            vx2=pt(1)-zones.zone(k).subSouth(n).R(2:end);
            vy2=pt(2)-zones.zone(k).subSouth(n).Z(2:end);
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
            %     d=min(sqrt((zones.zone(k).south.R-pt(1)).^2+(zones.zone(k).south.Z-pt(2)).^2));
            if(d<dmin)
                dmin=d;
                zone_sel=k;
                side_sel=2;
                sel_mode=true;
                sub_sel=n;
            end
        end
    end
    %east
    vx=zones.zone(k).east.R(2:end)-zones.zone(k).east.R(1:end-1);
    vy=zones.zone(k).east.Z(2:end)-zones.zone(k).east.Z(1:end-1);
    vx1=pt(1)-zones.zone(k).east.R(1:end-1);
    vy1=pt(2)-zones.zone(k).east.Z(1:end-1);
    vx2=pt(1)-zones.zone(k).east.R(2:end);
    vy2=pt(2)-zones.zone(k).east.Z(2:end);
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
    %     d=min(sqrt((zones.zone(k).east.R-pt(1)).^2+(zones.zone(k).east.Z-pt(2)).^2));
    if(d<dmin)
        dmin=d;
        zone_sel=k;
        side_sel=3;
        sel_mode=true;
    end
    vx=zones.zone(k).west.R(2:end)-zones.zone(k).west.R(1:end-1);
    vy=zones.zone(k).west.Z(2:end)-zones.zone(k).west.Z(1:end-1);
    vx1=pt(1)-zones.zone(k).west.R(1:end-1);
    vy1=pt(2)-zones.zone(k).west.Z(1:end-1);
    vx2=pt(1)-zones.zone(k).west.R(2:end);
    vy2=pt(2)-zones.zone(k).west.Z(2:end);
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
    %     d=min(sqrt((zones.zone(k).west.R-pt(1)).^2+(zones.zone(k).west.Z-pt(2)).^2));
    if(d<dmin)
        dmin=d;
        zone_sel=k;
        side_sel=4;
        sel_mode=true;
    end
end
if(side_sel==1)
    psival=zones.zone(zone_sel).pB.coord(1);
    set(handles.edit4,'String',num2str(psival));
    set(handles.edit5,'String',num2str(psival));
end
if(side_sel==2)
    psival=zones.zone(zone_sel).pA.coord(1);
    set(handles.edit4,'String',num2str(psival));
    set(handles.edit5,'String',num2str(psival));
end
if(side_sel>=3)
    set(handles.edit4,'String','');
    set(handles.edit5,'String','');
end
update_figure1(hObject, eventdata, handles)


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global zone_sel;
global side_sel;
global sel_mode;
sel_mode=false;
side_sel=0;
zone_sel=0;
update_figure1(hObject, eventdata, handles)



function edit1_Callback(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit1 as text
%        str2double(get(hObject,'String')) returns contents of edit1 as a double


% --- Executes during object creation, after setting all properties.
function edit1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit1 (see GCBO)
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
global zones;
global zone_sel;
global side_sel;
global megazone;
global modeAdjust;
mesh_pt_num=str2num(get(handles.edit1,'String'));
x=get(handles.uipanel2,'SelectedObject');
switch(get(x,'Tag'))
    case 'NoneT'
        RefineType=0;
    case 'LinearT'
        RefineType=1;
    case 'CustomT'
        RefineType=2;
end
x=get(handles.uipanel5,'SelectedObject');
switch(get(x,'Tag'))
    case 'LeftS'
        RefineSide=0;
    case 'RightS'
        RefineSide=1;
    case 'BothS'
        RefineSide=2;
end
paramL=str2num(get(handles.edit2,'String'));
paramR=str2num(get(handles.edit3,'String'));
mesh_segment(mesh_pt_num,RefineType,RefineSide,modeAdjust,paramL,paramR);
update_figure1(hObject, eventdata, handles)


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global megazone;
global Pmegazone;
global zones;
global align_psimode;
align_psimode=false;
for k=1:megazone.num
    megazone.mz(k).ismeshed=false;
end
for k=1:Pmegazone.num
    Pmegazone.mz(k).ismeshed=false;
    Pmegazone.mz(k).isaligned=false;
    Pmegazone.mz(k).refpoints.R=[];
    Pmegazone.mz(k).refpoints.Z=[];
    Pmegazone.mz(k).subrefpoints(1).R=[];
    Pmegazone.mz(k).subrefpoints(1).Z=[];
    Pmegazone.mz(k).subrefpoints(2).R=[];
    Pmegazone.mz(k).subrefpoints(2).Z=[];
    Pmegazone.mz(k).meshchanged=true;
end
for k=1:zones.num
    zones.zone(k).north.ismeshed=false;
    zones.zone(k).south.ismeshed=false;
    zones.zone(k).east.ismeshed=false;
    zones.zone(k).west.ismeshed=false;
    zones.zone(k).northaligned=false;
    zones.zone(k).southaligned=false;
    zones.zone(k).gridR=[];
    zones.zone(k).gridZ=[];
    zones.zone(k).orthomeshchanged=true;
end
update_figure1(hObject, eventdata, handles)


% --- Executes on button press in BothS.
function BothS_Callback(hObject, eventdata, handles)
% hObject    handle to BothS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of BothS


% --- Executes on button press in RightS.
function RightS_Callback(hObject, eventdata, handles)
% hObject    handle to RightS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of RightS


% --- Executes on button press in LeftS.
function LeftS_Callback(hObject, eventdata, handles)
% hObject    handle to LeftS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of LeftS


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.edit2,'String',num2str(get(hObject,'Value')));
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider2_Callback(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.edit3,'String',num2str(get(hObject,'Value')));
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


% --- Executes when selected object is changed in uipanel7.
function uipanel7_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel7
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)
global modeAdjust;
x=get(handles.uipanel7,'SelectedObject');
if(strcmp(get(x,'Tag'),'AbsoluteM'))
    set(handles.text3,'String','cm');
    set(handles.text4,'String','cm');
    modeAdjust=1;
else
    set(handles.text3,'String','');
    set(handles.text4,'String','');
    modeAdjust=0;
end


% --- Executes on button press in pushbutton5.
function pushbutton5_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global megazone;
global Pmegazone;
mesh_all_zones;
doit=true;
for k=1:megazone.num
    doit=doit&&megazone.mz(k).ismeshed;
end
for k=1:Pmegazone.num
    doit=doit&&Pmegazone.mz(k).ismeshed;
end
if(doit)
    conform_mesh;
end


% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
update_figure1(hObject, eventdata, handles);
% Hint: get(hObject,'Value') returns toggle state of checkbox1


% --- Executes on button press in checkbox2.
function checkbox2_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
update_figure1(hObject, eventdata, handles);
% Hint: get(hObject,'Value') returns toggle state of checkbox2


% --- Executes on button press in checkbox5.
function checkbox5_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if(get(hObject,'Value'))
    hsub=zones_prop();
    winontop(hsub);
end

% Hint: get(hObject,'Value') returns toggle state of checkbox5


% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
global zone_sel;
global zones;
global Pmegazone;
subdivide_selected_segment();
pmz=zones.zone(zone_sel).pmz;
Pmegazone.mz(pmz).isaligned=true;
Pmegazone.mz(pmz).align_psimin=str2num(get(handles.edit4,'String'));
Pmegazone.mz(pmz).align_psimax=str2num(get(handles.edit5,'String'));
Pmegazone.mz(pmz).ismeshed=false;
Pmegazone.mz(pmz).subrefpoints(1).R=[];
Pmegazone.mz(pmz).subrefpoints(1).Z=[];
Pmegazone.mz(pmz).subrefpoints(2).R=[];
Pmegazone.mz(pmz).subrefpoints(2).Z=[];
Pmegazone.mz(pmz).refpoints.R=[];
Pmegazone.mz(pmz).refpoints.Z=[];
for k=1:length(Pmegazone.mz(pmz).list)
    zones.zone(Pmegazone.mz(pmz).list(k)).north.ismeshed=false;
    zones.zone(Pmegazone.mz(pmz).list(k)).south.ismeshed=false;
end
global align_psimode;
align_psimode=false;
update_figure1(hObject, eventdata, handles);
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on slider movement.
function slider3_Callback(hObject, eventdata, handles)
% hObject    handle to slider3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global zones;
global zone_sel;
global side_sel;
global Pmegazone;
v=get(hObject,'Value');
if(side_sel==1)
    psival=zones.zone(zone_sel).pB.coord(1);
    pmz=zones.zone(zone_sel).pmz;
    nz=Pmegazone.mz(pmz).list(1);
    psimin=zones.zone(nz).pA.coord(1);
    psi=psival+v*(psimin-psival);
    set(handles.edit4,'String',num2str(psi));
end
if(side_sel==2)
    psival=zones.zone(zone_sel).pA.coord(1);
    pmz=zones.zone(zone_sel).pmz;
    nz=Pmegazone.mz(pmz).list(1);
    psimin=zones.zone(nz).pA.coord(1);
    psi=psival+v*(psimin-psival);
    set(handles.edit4,'String',num2str(psi));
end
global align_psimode;
align_psimode=true;
update_figure1(hObject, eventdata, handles);
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


% --- Executes on slider movement.
function slider4_Callback(hObject, eventdata, handles)
% hObject    handle to slider4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global zones;
global zone_sel;
global side_sel;
global Pmegazone;
v=get(hObject,'Value');
if(side_sel==1)
    psival=zones.zone(zone_sel).pB.coord(1);
    pmz=zones.zone(zone_sel).pmz;
    nz=Pmegazone.mz(pmz).list(end);
    psimax=zones.zone(nz).pB.coord(1);
    psi=psival+v*(psimax-psival);
    set(handles.edit5,'String',num2str(psi));
end
if(side_sel==2)
    psival=zones.zone(zone_sel).pA.coord(1);
    pmz=zones.zone(zone_sel).pmz;
    nz=Pmegazone.mz(pmz).list(end);
    psimax=zones.zone(nz).pB.coord(1);
    psi=psival+v*(psimax-psival);
    set(handles.edit5,'String',num2str(psi));
end
global align_psimode;
align_psimode=true;
update_figure1(hObject, eventdata, handles);
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function edit4_Callback(hObject, eventdata, handles)
% hObject    handle to edit4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
update_figure1(hObject, eventdata, handles);
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



function edit5_Callback(hObject, eventdata, handles)
% hObject    handle to edit5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
update_figure1(hObject, eventdata, handles);
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


% --- Executes on button press in pushbutton7.
function pushbutton7_Callback(hObject, eventdata, handles)
global zone_sel;
global zones;
global Pmegazone;
pmz=zones.zone(zone_sel).pmz;
Pmegazone.mz(pmz).isaligned=false;
for k=1:length(Pmegazone.mz(pmz).list)
    zones.zone(Pmegazone.mz(pmz).list(k)).northaligned=false;
    zones.zone(Pmegazone.mz(pmz).list(k)).southaligned=false;
end
% hObject    handle to pushbutton7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in pushbutton8.
function pushbutton8_Callback(hObject, eventdata, handles)
global zone_sel;
global zones;
global Pmegazone;
pmz=zones.zone(zone_sel).pmz;
Pmegazone.mz(pmz).align_psimin=str2num(get(handles.edit4,'String'));
Pmegazone.mz(pmz).align_psimax=str2num(get(handles.edit5,'String'));
global align_psimode;
align_psimode=false;
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on slider movement.
function slider5_Callback(hObject, eventdata, handles)
% hObject    handle to slider5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global orthoptim;
v=get(hObject,'Value');
orthoptim=(1-v);
global zones;
for k=1:zones.num
    zones.zone(k).orthomeshchanged=true;
end
% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider5_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function edit6_Callback(hObject, eventdata, handles)
% hObject    handle to edit6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global dmin_ad;
dmin_ad=str2num(get(hObject,'String'));
global zones;
for k=1:zones.num
    zones.zone(k).orthomeshchanged=true;
end
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



function edit8_Callback(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global range_ad;
range_ad=str2num(get(hObject,'String'));
global zones;
for k=1:zones.num
    zones.zone(k).orthomeshchanged=true;
end
% Hints: get(hObject,'String') returns contents of edit8 as text
%        str2double(get(hObject,'String')) returns contents of edit8 as a double


% --- Executes during object creation, after setting all properties.
function edit8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit9_Callback(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global length_ad;
length_ad=floor(str2num(get(hObject,'String')));
global zones;
for k=1:zones.num
    zones.zone(k).orthomeshchanged=true;
end
% Hints: get(hObject,'String') returns contents of edit9 as text
%        str2double(get(hObject,'String')) returns contents of edit9 as a double


% --- Executes during object creation, after setting all properties.
function edit9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit10_Callback(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global alpha_ad;
alpha_ad=str2num(get(hObject,'String'));
global zones;
for k=1:zones.num
    zones.zone(k).orthomeshchanged=true;
end
% Hints: get(hObject,'String') returns contents of edit10 as text
%        str2double(get(hObject,'String')) returns contents of edit10 as a double


% --- Executes during object creation, after setting all properties.
function edit10_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit11_Callback(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global target_ad;
target_ad=str2num(get(hObject,'String'));
global zones;
for k=1:zones.num
    zones.zone(k).orthomeshchanged=true;
end
% Hints: get(hObject,'String') returns contents of edit11 as text
%        str2double(get(hObject,'String')) returns contents of edit11 as a double


% --- Executes during object creation, after setting all properties.
function edit11_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton9.
function pushbutton9_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
tweak_mesh;
conform_mesh;
update_figure1(hObject, eventdata, handles);


function edit12_Callback(hObject, eventdata, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global drmin_ad
drmin_ad=str2num(get(hObject,'String'));
% Hints: get(hObject,'String') returns contents of edit12 as text
%        str2double(get(hObject,'String')) returns contents of edit12 as a double


% --- Executes during object creation, after setting all properties.
function edit12_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox6.
function checkbox6_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
update_figure1(hObject, eventdata, handles);
% Hint: get(hObject,'Value') returns toggle state of checkbox6


% --- Executes on button press in pushbutton10.
function pushbutton10_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Generator();


% --- Executes on button press in pushbutton11.
function pushbutton11_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global Pmegazone
global megazone
global zones
for k=1:Pmegazone.num
    Pmegazone.mz(k).meshchanged=true;
end
for k=1:zones.num
    zones.zone(k).north.ismeshed=false;
    zones.zone(k).south.ismeshed=false;
    zones.zone(k).east.ismeshed=false;
    zones.zone(k).west.ismeshed=false;
    zones.zone(k).northaligned=false;
    zones.zone(k).southaligned=false;
    zones.zone(k).gridR=[];
    zones.zone(k).gridZ=[];
    zones.zone(k).orthomeshchanged=true;
end


% --- Executes on button press in pushbutton12.
function pushbutton12_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton12 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
regularize_mesh_size;
conform_mesh;
update_figure1(hObject, eventdata, handles);
