function varargout = retMenu(varargin)
% retMenu - gui for retinotopic mapping program parameters.
%
% SOD 10/2005: created it, consolidating several existing gui's.

% RETMENU M-file for retMenu.fig
%      RETMENU, by itself, creates a new RETMENU or raises the existing
%      singleton*.
%
%      H = RETMENU returns the handle to a new RETMENU or the handle to
%      the existing singleton*.
%
%      RETMENU('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in RETMENU.M with the given input arguments.
%
%      RETMENU('Property','Value',...) creates a new RETMENU or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before retMenu_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to retMenu_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help retMenu

% Last Modified by GUIDE v2.5 27-Oct-2005 22:38:25

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @retMenu_OpeningFcn, ...
                   'gui_OutputFcn',  @retMenu_OutputFcn, ...
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

%----------------------------------------------------------
% --- Executes just before retMenu is made visible.
%----------------------------------------------------------
function retMenu_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to retMenu (see VARARGIN)

% Choose default menu params for retMenu
set(handles.experiment,     'String',setRetinotopyParams);
fixString  = {'disk', 'small cross +','double disk'};
set(handles.fixation(1),       'String',fixString);
set(handles.savestimparams, 'Value',1); 
set(handles.repetitions,    'String','1');    %#
set(handles.runPriority,    'String','2');    %# 7
set(handles.skipCycleFrames,'String','0');    %#
set(handles.prescanDuration,'String','12');   %seconds
set(handles.period,         'String','24');   %seconds
set(handles.numCycles,      'String','6');    %#
set(handles.interleaves,    'String','N/A');    %#
set(handles.tr,             'String','1.5');  %seconds
% all .mat files where retMenu lives are considered potentials
tmpdir = which(mfilename);
tmpdir = fileparts(tmpdir);
tmp    = dir([tmpdir filesep '*.mat']);
if ~isempty(tmp),
    tmp2{1}   = 'None'; 
    for n=1:length(tmp),
        tmp2{n+1} = tmp(n).name;
    end;
    set(handles.loadMatrix,     'String',tmp2);  
else,
    set(handles.loadMatrix,     'String','None');
end;
set(handles.saveMatrix,     'String','');  
% all directories in matlabroot displays are considered 
% calibration directories
try,
    mydir{1} = 'None';
    tmpdir = dir(getDisplayPath);
    count = 2;
    for n=1:length(tmpdir), % stupid loop
        if tmpdir(n).isdir & ~strcmp(tmpdir(n).name(1),'.'),
            mydir{count} = tmpdir(n).name;
            count = count+1;
        end;
    end;
    set(handles.calibration,'String',mydir);
catch,
    set(handles.calibration,'String','None');
end;
set(handles.stimSize,      'String','9');%'max');

% default command line output for retMenu = menuparams
tmp      = get(handles.experiment,'String');
data.experiment      = tmp(get(handles.experiment,'Value'));      
tmp      = get(handles.fixation(1),'String');
data.fixation      = tmp(get(handles.fixation(1),'Value'));      
data.savestimparams  = get(handles.savestimparams,'Value');;
data.repetitions     = str2double(get(handles.repetitions,'String'));
data.runPriority     = str2double(get(handles.runPriority,'String'));
data.skipCycleFrames = str2double(get(handles.skipCycleFrames,'String'));
data.prescanDuration = str2double(get(handles.prescanDuration,'String'));
data.period          = str2double(get(handles.period,'String'));
data.numCycles       = str2double(get(handles.numCycles,'String'));
data.interleaves     = str2double(get(handles.interleaves,'String'));
data.tr              = str2double(get(handles.tr,'String'));
data.loadMatrix      = 'None';
data.saveMatrix      = [];
data.calibration     = 'Skyra'; %'None';
data.stimSize        = str2double(get(handles.stimSize,'String'));
if isnan(data.stimSize),
    data.stimSize    = 'max';
end;

% store in handles structure
handles.data = data;

% Update handles structure
guidata(hObject, handles);


% UIWAIT makes retMenu wait for user response (see UIRESUME)
uiwait(handles.figure1);
return;

%----------------------------------------------------------
% --- Outputs from this function are returned to the command line.
%----------------------------------------------------------
function varargout = retMenu_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);

% we want strings from these not cells
if iscell(handles.data.experiment),
    handles.data.experiment = handles.data.experiment{1};
end;
if iscell(handles.data.fixation),
    handles.data.fixation = handles.data.fixation{1};
end;
if iscell(handles.data.loadMatrix),
    handles.data.loadMatrix = handles.data.loadMatrix{1};
end;
if iscell(handles.data.calibration),
    handles.data.calibration = handles.data.calibration{1};
end;

% directory this file lives in
tmpdir = which(mfilename);
tmpdir = fileparts(tmpdir);
tmpdir = [tmpdir filesep];

% if we save the image matrix we put it in the directory this file lives in
if ~isempty(handles.data.saveMatrix),
    handles.data.saveMatrix = [tmpdir handles.data.saveMatrix];
end;

% if no image matrix are found to be loaded -> empty
if strcmp(lower(handles.data.loadMatrix),'none')
    handles.data.loadMatrix = [];
else,
    handles.data.loadMatrix = [tmpdir handles.data.loadMatrix];
end;

if strcmp(lower(handles.data.calibration),'none'),
    handles.data.calibration = [];
end;
    
% output
varargout{1} = handles.data;

% quit
delete(handles.figure1);
return;

%----------------------------------------------------------
% --- These are only executed when changed by user
%----------------------------------------------------------
function experiment_Callback(hObject, eventdata, handles)
% returns experiment contents as cell array
contents = get(hObject,'String');
% returns selected item from experiment
handles.data.experiment = contents{get(hObject,'Value')};
% Set default options for some experiments
% This is convenient (for me) but may be tricky if you are not aware other
% settings are changing as well! So print a warning message at least.
switch handles.data.experiment,
     case {'Translating Bars 1', 'Translating Bars 2', 'Translating Bars 3'},
        handles.data.savestimparams = 1;
        set(handles.savestimparams, 'Value',1); 
        handles.data.tr = 1.5;
        set(handles.tr, 'String','1.5');  
        handles.data.numCycles = 1;
        set(handles.numCycles, 'String','1');    %#
        handles.data.period    = 180;
        set(handles.period,    'String','180');   %seconds
        handles.data.prescanDuration = 12;
        set(handles.prescanDuration,    'String','12');   %seconds
        
        disp(sprintf('[%s]:WARNING: when setting experiment to %s, num Cycles (%d) and period (%d) and others have been changed as well!',...
            mfilename,handles.data.experiment,handles.data.numCycles,handles.data.period));
        
    case {'Translating Bars 8 Pass (Dumoulin)'},
        % period should be:
        % 20 bar positions * 8 directions * TR + 4 blanks period * ceil(10.5/TR) * TR
        % 20*8*TR + 4*ceil(10.5/TR)*TR
        handles.data.savestimparams = 1;
        set(handles.savestimparams, 'Value',1); 
        handles.data.tr = 1.8; %1.5
        set(handles.tr, 'String','1.8');
        handles.data.numCycles = 1;
        set(handles.numCycles, 'String','1');    %#
        handles.data.period    = 331.2; %282
        set(handles.period,    'String','331.2');   %seconds
        handles.data.prescanDuration = 12.6; %12;
        set(handles.prescanDuration,    'String','12.6');   %seconds
        
        disp(sprintf('[%s]:WARNING: when setting experiment to %s, num Cycles (%d) and period (%d) and others have been changed as well!',...
            mfilename,handles.data.experiment,handles.data.numCycles,handles.data.period));
    otherwise,
        % do nothing
end;
    
guidata(hObject,handles);
return;

function experiment_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
return;


%----------------------------------------------------------
function tr_Callback(hObject, eventdata, handles)
handles.data.tr = str2double(get(hObject,'String'));
guidata(hObject,handles);
return;

function tr_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
return;

%----------------------------------------------------------
function interleaves_Callback(hObject, eventdata, handles)
handles.data.interleaves = str2double(get(hObject,'String'));
guidata(hObject,handles);
return;

function interleaves_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
return;

%----------------------------------------------------------
function numCycles_Callback(hObject, eventdata, handles)
handles.data.numCycles = str2double(get(hObject,'String'));
guidata(hObject,handles);
return;

function numCycles_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
return;

%----------------------------------------------------------
function period_Callback(hObject, eventdata, handles)
handles.data.period=str2double(get(hObject,'String'));
guidata(hObject,handles);
return;

function period_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
return;

%----------------------------------------------------------
function prescanDuration_Callback(hObject, eventdata, handles)
handles.data.prescanDuration = str2double(get(hObject,'String'));
guidata(hObject,handles);
return;

function prescanDuration_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
return;

%----------------------------------------------------------
function skipCycleFrames_Callback(hObject, eventdata, handles)
handles.data.skipCycleFrames = str2double(get(hObject,'String'));
guidata(hObject,handles);
return;

function skipCycleFrames_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
return;

%----------------------------------------------------------
function repetitions_Callback(hObject, eventdata, handles)
handles.data.repetitions = str2double(get(hObject,'String'));
guidata(hObject,handles);
return;

function repetitions_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
return;

%----------------------------------------------------------
function runPriority_Callback(hObject, eventdata, handles)
handles.data.runPriority=str2double(get(hObject,'String'));
guidata(hObject,handles);
return;

function runPriority_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
return;

%----------------------------------------------------------
function fixation_Callback(hObject, eventdata, handles)
contents = get(hObject,'String');
handles.data.fixation=contents{get(hObject,'Value')};
guidata(hObject,handles);
return;

function fixation_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
return;

%----------------------------------------------------------
% --- Executes on button press in go.
function go_Callback(hObject, eventdata, handles)
% now continue and finish
uiresume(handles.figure1);
return;

%----------------------------------------------------------
function Untitled_1_Callback(hObject, eventdata, handles)
return;


%----------------------------------------------------------
function saveMatrix_Callback(hObject, eventdata, handles)
handles.data.saveMatrix=get(hObject,'String');
guidata(hObject,handles);
return;

function saveMatrix_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
return;

%----------------------------------------------------------
% not sure what this one does - i should remove it at some point
function popupmenu4_Callback(hObject, eventdata, handles)
% Hints: contents = get(hObject,'String') returns popupmenu4 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu4
return;

function popupmenu4_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
return

%----------------------------------------------------------
function loadMatrix_Callback(hObject, eventdata, handles)
contents = get(hObject,'String');
if iscell(contents),
    handles.data.loadMatrix=contents{get(hObject,'Value')};
else,
    handles.data.loadMatrix=contents;
end;
guidata(hObject,handles);
return

function loadMatrix_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
return;


%----------------------------------------------------------
function calibration_Callback(hObject, eventdata, handles)
contents = get(hObject,'String');
if iscell(contents),
    handles.data.calibration=contents{get(hObject,'Value')};
else,
    handles.data.calibration=contents;
end;
guidata(hObject,handles);
return

function calibration_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
return;


%----------------------------------------------------------
function stimSize_Callback(hObject, eventdata, handles)
handles.data.stimSize=str2double(get(hObject,'String'));
if isnan(handles.data.stimSize),
    handles.data.stimSize = 'max';
end;
guidata(hObject,handles);
return;


function stimSize_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
return;


%----------------------------------------------------------
function savestimparams_Callback(hObject, eventdata, handles)
handles.data.savestimparams = get(hObject,'Value');
guidata(hObject,handles);
return;


