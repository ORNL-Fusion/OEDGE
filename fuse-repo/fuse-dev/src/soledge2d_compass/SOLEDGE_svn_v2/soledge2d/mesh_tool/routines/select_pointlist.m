function varargout = select_pointlist(varargin)
% SELECT_POINTLIST M-file for select_pointlist.fig
%      SELECT_POINTLIST, by itself, creates a new SELECT_POINTLIST or raises the existing
%      singleton*.
%
%      H = SELECT_POINTLIST returns the handle to a new SELECT_POINTLIST or the handle to
%      the existing singleton*.
%
%      SELECT_POINTLIST('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SELECT_POINTLIST.M with the given input arguments.
%
%      SELECT_POINTLIST('Property','Value',...) creates a new SELECT_POINTLIST or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before select_pointlist_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to select_pointlist_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help select_pointlist

% Last Modified by GUIDE v2.5 02-Nov-2015 16:23:16

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @select_pointlist_OpeningFcn, ...
    'gui_OutputFcn',  @select_pointlist_OutputFcn, ...
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


% --- Executes just before select_pointlist is made visible.
function select_pointlist_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to select_pointlist (see VARARGIN)

% Choose default command line output for select_pointlist
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

global list_res;
list_res=[];
global parent_fig;
parent_fig=varargin{1};
% UIWAIT makes select_pointlist wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = select_pointlist_OutputFcn(hObject, eventdata, handles)
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


% --- Executes on button press in pushbutton2.
function pushbutton2_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global frontiers;
global list_res;
global parent_fig;
global X_points;
global r2D;
global z2D;
global flux2D;
%look for intersections with X_point radial limits (should be two)
P1=[];
for k=1:X_points.num
    for n=1:4
        c1.num=1;
        c1.arc(1).x=X_points.cut(k).arc(n).psiR;
        c1.arc(1).y=X_points.cut(k).arc(n).psiZ;
        c2.num=1;
        c2.arc(1).x=list_res(:,1);
        c2.arc(1).y=list_res(:,2);
        inter=intersect_contour(c1,c2);
        if(inter.num==1)
            %check if psi surface already tested
            in=false;
            [m,p]=size(P1);
            for k1=1:m
                if(P1(k1,3)==X_points.cut(k).psilim(n))
                    in=true;
                end
            end
            if(~in) %if no doublon
                P1=[P1;inter.x(1),inter.y(1),X_points.cut(k).psilim(n)];
            end
        end
    end
end
frontiers.num=frontiers.num+1;
[m,p]=size(list_res);
psipt=zeros(1,m);
for k=1:m
   psipt(k)=interp2(r2D,z2D,flux2D,list_res(k,1),list_res(k,2)); 
end
[a,b]=sort(psipt);
frontiers.lim(frontiers.num).R=list_res(b,1);
frontiers.lim(frontiers.num).Z=list_res(b,2);
frontiers.lim(frontiers.num).P1=P1;
frontiers.lim(frontiers.num).sel=false;
frontiers.lim(frontiers.num).psimin=min(P1(:,3));
frontiers.lim(frontiers.num).psimax=max(P1(:,3));
uiresume(parent_fig);
delete(handles.figure1);


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global selectedlines;
global list_res;
[m,p]=size(selectedlines);
[m2,p2]=size(list_res);
listnew=[];
for i=1:m2
    in=true;
    for i2=1:m
        if(selectedlines(i2,1)==i)
            in=false;
        end
    end
    if(in)
        listnew=[listnew;list_res(i,:)];
    end
end
list_res=listnew;
update_table(hObject, eventdata, handles);

function update_table(hObject, eventdata, handles)
global list_res;
global parent_fig;
set(handles.uitable1,'Data',list_res);
figure(parent_fig);
plot(list_res(:,1),list_res(:,2),'b.-')


% --- Executes on button press in pushbutton4.
function pushbutton4_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global parent_fig;
global list_res;
figure(parent_fig);
pt=ginput(1);
list_res=[list_res;pt(1),pt(2)];
update_table(hObject, eventdata, handles);


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


% --- Executes when selected cell(s) is changed in uitable1.
function uitable1_CellSelectionCallback(hObject, eventdata, handles)
% hObject    handle to uitable1 (see GCBO)
% eventdata  structure with the following fields (see UITABLE)
%	Indices: row and column indices of the cell(s) currently selecteds
% handles    structure with handles and user data (see GUIDATA)
global selectedlines;
selectedlines=eventdata.Indices;