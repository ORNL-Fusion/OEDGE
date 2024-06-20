function WasOnTop = WinOnTop( FigureHandle, IsOnTop )
%WINONTOP allows to trigger figure's "Always On Top" state
% INPUT ARGUMENTS:
% * FigureHandle - Matlab's figure handle, scalar
% * IsOnTop      - logical scalar or empty array
%
% USAGE:
% * WinOnTop( hfigure, bool );
% * WinOnTop( hfigure );            - equal to WinOnTop( hfigure,true);
% * WinOnTop();                     - equal to WinOnTop( gcf, true);
% * WasOnTop = WinOnTop(...);       - returns boolean value "if figure WAS on top"
% * IsOnTop  = WinOnTop(hfigure,[]) - gets "if figure is on top" property
% 
% LIMITATIONS:
% * java enabled
% * figure must be visible
% * figure's "WindowStyle" should be "normal"
%
%
% Written by Igor
% i3v@mail.ru
%
% 16 June 2013 - Initial version
% 27 June 2013 - removed custom "ishandle_scalar" function call
%

%% Parse Inputs

if ~exist('FigureHandle','var');FigureHandle = gcf; end
assert(...
          isscalar( FigureHandle ) && ishandle( FigureHandle ) &&  strcmp(get(FigureHandle,'Type'),'figure'),...
          'WinOnTop:Bad_FigureHandle_input',...
          '%s','Provided FigureHandle input is not a figure handle'...
       );

assert(...
            strcmp('on',get(FigureHandle,'Visible')),...
            'WinOnTop:FigInisible',...
            '%s','Figure Must be Visible'...
       );

assert(...
            strcmp('normal',get(FigureHandle,'WindowStyle')),...
            'WinOnTop:FigWrongWindowStyle',...
            '%s','WindowStyle Must be Normal'...
       );
   
if ~exist('IsOnTop','var');IsOnTop=true;end
assert(...
          islogical( IsOnTop ) &&  isscalar( IsOnTop) || isempty( IsOnTop ), ...
          'WinOnTop:Bad_IsOnTop_input',...
          '%s','Provided IsOnTop input is neither boolean, nor empty'...
      );
%% Pre-checks

error(javachk('swing',mfilename)) % Swing components must be available.
  
  
%% Action

% Flush the Event Queue of Graphic Objects and Update the Figure Window.
drawnow expose

warnStruct=warning('off','MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');
jFrame = get(handle(FigureHandle),'JavaFrame');
warning(warnStruct.state,'MATLAB:HandleGraphics:ObsoletedProperty:JavaFrame');

drawnow

v=ver('Matlab');


if str2double(v.Version) >= 8.4
    jFrame_fHGxClient = jFrame.fHG2Client;
else
    jFrame_fHGxClient = jFrame.fHG1Client;
end


WasOnTop = jFrame_fHGxClient.getWindow.isAlwaysOnTop;

if ~isempty(IsOnTop)
    jFrame_fHGxClient.getWindow.setAlwaysOnTop(IsOnTop);
end

end

