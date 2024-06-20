function varargout = zones_prop(varargin)
% ZONES_PROP MATLAB code for zones_prop.fig
%      ZONES_PROP, by itself, creates a new ZONES_PROP or raises the existing
%      singleton*.
%
%      H = ZONES_PROP returns the handle to a new ZONES_PROP or the handle to
%      the existing singleton*.
%
%      ZONES_PROP('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ZONES_PROP.M with the given input arguments.
%
%      ZONES_PROP('Property','Value',...) creates a new ZONES_PROP or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before zones_prop_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to zones_prop_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help zones_prop

% Last Modified by GUIDE v2.5 11-Nov-2015 19:06:46

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @zones_prop_OpeningFcn, ...
                   'gui_OutputFcn',  @zones_prop_OutputFcn, ...
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


% --- Executes just before zones_prop is made visible.
function zones_prop_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to zones_prop (see VARARGIN)

% Choose default command line output for zones_prop
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
update_tableZ(hObject, eventdata, handles);
% UIWAIT makes zones_prop wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = zones_prop_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

function update_tableZ(hObject, eventdata, handles)
global zones;
d=zeros(zones.num,4);
for i=1:zones.num
   d(i,1)=zones.zone(i).Neighbour.north;
   d(i,2)=zones.zone(i).Neighbour.south;
   d(i,3)=zones.zone(i).Neighbour.east;
   d(i,4)=zones.zone(i).Neighbour.west;
end
set(handles.uitable1,'Data',d);
