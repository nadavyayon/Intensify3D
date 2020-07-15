function varargout = User_GUI_Intensify3D(varargin)

% USER_GUI_INTENSIFY3D MATLAB code for User_GUI_Intensify3D.fig
%      USER_GUI_INTENSIFY3D, by itself, creates a new USER_GUI_INTENSIFY3D or raises the existing
%      singleton*.
%
%      H = USER_GUI_INTENSIFY3D returns the handle to a new USER_GUI_INTENSIFY3D or the handle to
%      the existing singleton*.
%
%      USER_GUI_INTENSIFY3D('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in USER_GUI_INTENSIFY3D.M with the given input arguments.
%
%      USER_GUI_INTENSIFY3D('Property','Value',...) creates a new USER_GUI_INTENSIFY3D or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before User_GUI_Intensify3D_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to User_GUI_Intensify3D_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help User_GUI_Intensify3D

% Last Modified by GUIDE v2.5 14-Jul-2020 13:55:02

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @User_GUI_Intensify3D_OpeningFcn, ...
                   'gui_OutputFcn',  @User_GUI_Intensify3D_OutputFcn, ...
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


% --- Executes just before User_GUI_Intensify3D is made visible.
function User_GUI_Intensify3D_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to User_GUI_Intensify3D (see VARARGIN)

% Choose default command line output for User_GUI_Intensify3D
handles.output = hObject;
addpath(pwd);
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes User_GUI_Intensify3D wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = User_GUI_Intensify3D_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function StackFolder_Callback(hObject, eventdata, handles)
% hObject    handle to StackFolder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of StackFolder as text
%        str2double(get(hObject,'String')) returns contents of StackFolder as a double


% --- Executes during object creation, after setting all properties.
function StackFolder_CreateFcn(hObject, eventdata, handles)
% hObject    handle to StackFolder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ParentFolderBrowse.
function ParentFolderBrowse_Callback(hObject, eventdata, handles)
% hObject    handle to ParentFolderBrowse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Folder = uigetdir(cd,'MATLAB Root Directory');
set(handles.StackFolder,'string',Folder);
import java.lang.*;
r=Runtime.getRuntime;
ncpu=r.availableProcessors;
Directories = dir(Folder);
DirectoryIndex = find([Directories.isdir])'; 
for i = 1:length(DirectoryIndex)-2
    StringArray{i,1} = [num2str(i),') ',Directories(DirectoryIndex(i+2)).name];
end
set(handles.listbox1, 'string', string(StringArray));
set(handles.text47, 'String', num2str(ncpu/2));

% --- Executes on button press in FileBrowse.
function FileBrowse_Callback(hObject, eventdata, handles)
% hObject    handle to FileBrowse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
addpath(pwd);
[MultiTiffFile, MultiTiffFileFolder]  = uigetfile({'*.tif;'},'Select Multi Tiff file',' ');
cd(MultiTiffFileFolder);
mkdir(MultiTiffFile(1:end-4));
Ti = imfinfo(MultiTiffFile);
cd(MultiTiffFile(1:end-4));
h = waitbar(0,'Creating image sequence, Please wait...');
for i = 1:length(Ti)
    T = imread([MultiTiffFileFolder MultiTiffFile],i);
    imwrite(T,[MultiTiffFile(1:end-4),'_', sprintf('%3.4d',i),'.tif']);
    waitbar(i/length(Ti));
end
% close(h);

set(handles.StackFolder,'string',[MultiTiffFileFolder MultiTiffFile(1:end-4)]);
import java.lang.*;
r=Runtime.getRuntime;
ncpu=r.availableProcessors;
set(handles.text15, 'String', num2str(ncpu/2));

function StartIm_Callback(hObject, eventdata, handles)
% hObject    handle to StartIm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of StartIm as text
%        str2double(get(hObject,'String')) returns contents of StartIm as a double


% --- Executes during object creation, after setting all properties.
function StartIm_CreateFcn(hObject, eventdata, handles)
% hObject    handle to StartIm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function EndIm_Callback(hObject, eventdata, handles)
% hObject    handle to EndIm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of EndIm as text
%        str2double(get(hObject,'String')) returns contents of EndIm as a double


% --- Executes during object creation, after setting all properties.
function EndIm_CreateFcn(hObject, eventdata, handles)
% hObject    handle to EndIm (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function MBI_Callback(hObject, eventdata, handles)
% hObject    handle to MBI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MBI as text
%        str2double(get(hObject,'String')) returns contents of MBI as a double


% --- Executes during object creation, after setting all properties.
function MBI_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MBI (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function MSFS_Callback(hObject, eventdata, handles)
% hObject    handle to MSFS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of MSFS as text
%        str2double(get(hObject,'String')) returns contents of MSFS as a double


% --- Executes during object creation, after setting all properties.
function MSFS_CreateFcn(hObject, eventdata, handles)
% hObject    handle to MSFS (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in Execute.
function Execute_Callback(hObject, eventdata, handles)

try close f 
catch
end

TissueDetection = find([get(handles.TissueDetection_1,'Value'), get(handles.TissueDetection_2,'Value'), get(handles.TissueDetection_3,'Value'), get(handles.TissueDetection_4,'Value')]);
NormType = find([get(handles.ZType_1,'Value'), get(handles.ZType_2,'Value'), get(handles.ZType_3,'Value'), get(handles.ZType_4,'Value')]);
XYNorm = find([get(handles.XYType_2,'Value'), get(handles.XYType_1,'Value')]);
AcrossNormType = find([get(handles.Across_1,'Value'), get(handles.Across_2,'Value'), get(handles.Across_3,'Value'), get(handles.Across_4,'Value')]);

UserSelectedThrehold = str2num(get(handles.edit26,'String'));
TissueSmoothing = str2num(get(handles.edit21,'String'));
ExistingSupportFilesFolder = get(handles.checkbox1,'Value');
Folder = get(handles.StackFolder,'string');
Sensitivity = str2num(get(handles.edit20,'string'));
MBI = str2num(get(handles.edit24,'string'));

ImageToMeasure = str2num(get(handles.ImageToShow,'string'));
FolderToMeasure = str2num(get(handles.edit35,'string'));
clc
addpath(pwd);
MSFS = str2num(get(handles.edit25,'string'));
Threads = str2num(get(handles.edit23,'string'));
idx = mod(MSFS,2)<1;
MSFSo = floor(MSFS);
MSFSo(idx) = MSFSo(idx)+1;
try
    Intensify3D_Plus(Folder,Sensitivity,MSFSo,MBI,TissueDetection-1,NormType,Threads,ImageToMeasure,FolderToMeasure,handles,UserSelectedThrehold,TissueSmoothing,ExistingSupportFilesFolder,XYNorm-1,AcrossNormType)
catch me
    set(handles.Execute,'FontSize',20);
    set(handles.Execute,'ForegroundColor','black');
    set(handles.Execute,'String','Start');
    msgbox(['Something went wrong, Check that all parameters are filled in correctly. Error in function - ',me.stack(1).name,'.    Line - ',num2str(me.stack(1).line),'.   ',me.message])
end

% --- Executes on button press in togglebutton1.
function togglebutton1_Callback(hObject, eventdata, handles)
% hObject    handle to togglebutton1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of toTogglegglebutton1



function NSTD_Callback(hObject, eventdata, handles)
% hObject    handle to NSTD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of NSTD as text
%        str2double(get(hObject,'String')) returns contents of NSTD as a double


% --- Executes during object creation, after setting all properties.
function NSTD_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NSTD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function ZType_1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ZType_4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function ZType_2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ZType_2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function ZType_3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ZType_3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function ZType_4_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ZType_4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function ParentFolderBrowse_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ParentFolderBrowse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in ImageShow.
function ImageShow_Callback(hObject, eventdata, handles)
% hObject    handle to ImageShow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Folder = get(handles.StackFolder,'string');
Directories = dir(Folder);
DirectoryIndex = find([Directories.isdir])'; 
DirectoryIndex = DirectoryIndex(3:end);
StringDir = [Directories(DirectoryIndex(str2num(get(handles.edit35,'string')))).name]; 
cd([Folder,'\',StringDir]);
FolderInfo = dir('*.tif');
FirstImageFileIndex = find([FolderInfo.bytes]>10*1000,1,'first');
Image2Show = str2num(get(handles.ImageToShow,'string'));
set(handles.ImageName, 'String', FolderInfo(FirstImageFileIndex+Image2Show-1).name);
LastImageFileIndex = find([FolderInfo.bytes]>10*1000,1,'last');
axes(handles.axes1);
I = imread([Folder,'\',StringDir,'' filesep '',FolderInfo(FirstImageFileIndex+Image2Show-1).name]);
RI = imref2d(size(I));
handles.I=I;
guidata(hObject,handles)
MaxBsckgroundGuess = quantile(I(:),0.995);
ImagingMediaint = quantile(I(:),0.5);

set(handles.edit26,'string',num2str(ImagingMediaint));
set(handles.edit24,'string',num2str(MaxBsckgroundGuess))
set(handles.edit25,'string',num2str(round(max(size(I))/5)));
imshow(I,RI);

set(handles.MaxInt,'string',num2str(max(I(:))));
set(handles.MinInt,'string',num2str(min(I(:))));
Factor = single(MaxBsckgroundGuess)/single(max(I(:)));
handles.Factor=Factor;
guidata(hObject,handles)
set(handles.sliderThreshold,'Value',Factor);
imcontrast
ThresholdRed = round(256*Factor);
axes(handles.axes1);
NewMap=colormap(gray(256));
handles.CM=colormap(gray(256));
guidata(hObject,handles)
NewMap(ThresholdRed+1:end,:) = repmat([1 0 0],256-ThresholdRed,1);
colormap(handles.axes1,NewMap);


function Threads_Callback(hObject, eventdata, handles)
% hObject    handle to Threads (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Threads as text
%        str2double(get(hObject,'String')) returns contents of Threads as a double


% --- Executes during object creation, after setting all properties.
function Threads_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Threads (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on key press with focus on ZType_2 and none of its controls.
function ZType_2_KeyPressFcn(hObject, eventdata, handles)
% hObject    handle to ZType_2 (see GCBO)
% eventdata  structure with the following fields (see UICONTROL)
%	Key: name of the key that was pressed, in lower case
%	Character: character interpretation of the key(s) that was pressed
%	Modifier: name(s) of the modifier key(s) (i.e., control, shift) pressed
% handles    structure with handles and user data (see GUIDATA)

function ImageToShow_Callback(hObject, eventdata, handles)
% hObject    handle to ImageToShow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ImageToShow as text
%        str2double(get(hObject,'String')) returns contents of ImageToShow as a double

% --- Executes during object creation, after setting all properties.
function ImageToShow_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ImageToShow (see GCBO)
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


% --- Executes during object creation, after setting all properties.
function uipanel1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uipanel1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function uipanel3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uipanel3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object creation, after setting all properties.
function uipanel6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to uipanel6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes during object deletion, before destroying properties.
function uipanel1_DeleteFcn(hObject, eventdata, handles)
% hObject    handle to uipanel1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on slider movement.
function sliderThreshold_Callback(hObject, eventdata, handles)
% hObject    handle to sliderThreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
usermax = (get(handles.sliderThreshold, 'value'));
% Max = str2num(get(handles.MaxInt,'string'));
MinMax = get(handles.axes1,'CLim');
set(handles.MaxInt,'string',MinMax(2));
set(handles.MinInt,'string',MinMax(1));

% Min = str2num(get(handles.MinInt,'string'));
usermaxtrue = round(MinMax(1) + usermax*(MinMax(2) - MinMax(1)));
set(handles.edit24,'string',num2str(usermaxtrue))
Factor = usermax;
if Factor == 1; Factor = 255/256; end
if Factor == 0; Factor = 1/256; end

ThresholdRed = round(256*Factor);
axes(handles.axes1);
NewMap=handles.CM;
colormap(NewMap);
NewMap=colormap;
NewMap(ThresholdRed+1:end,:) = repmat([1 0 0],256-ThresholdRed,1);
colormap(handles.axes1,NewMap);



% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function sliderThreshold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliderThreshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --------------------------------------------------------------------
function About1_ClickedCallback(hObject, eventdata, handles)

% hObject    handle to About1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
website_link = 'https://github.com/nadavyayon/Intensify3D'
rs = questdlg({'Intensify3D - Fluorescent Image Normalization Tool,  2017 v1.0';
['More information can be found at ',website_link];'Nadav Yayon, Soreq Lab, Hebrew University of Jerusalem'},'About', 'Visit Website','Ok')

if rs == 'Visit Website'
web(website_link,'-browser')
end


% --- Executes on slider movement.
function slider5_Callback(hObject, eventdata, handles)
% hObject    handle to slider5 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

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


% --- Executes on button press in pushbutton8.
function pushbutton8_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit15_Callback(hObject, eventdata, handles)
% hObject    handle to edit15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit15 as text
%        str2double(get(hObject,'String')) returns contents of edit15 as a double


% --- Executes during object creation, after setting all properties.
function edit15_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit15 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function slider6_Callback(hObject, eventdata, handles)
% hObject    handle to slider6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider6_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in pushbutton9.
function pushbutton9_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit16_Callback(hObject, eventdata, handles)
% hObject    handle to edit16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit16 as text
%        str2double(get(hObject,'String')) returns contents of edit16 as a double


% --- Executes during object creation, after setting all properties.
function edit16_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit16 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit18_Callback(hObject, eventdata, handles)
% hObject    handle to edit18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit18 as text
%        str2double(get(hObject,'String')) returns contents of edit18 as a double


% --- Executes during object creation, after setting all properties.
function edit18_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit18 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit17_Callback(hObject, eventdata, handles)
% hObject    handle to edit17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit17 as text
%        str2double(get(hObject,'String')) returns contents of edit17 as a double


% --- Executes during object creation, after setting all properties.
function edit17_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit17 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton10.
function pushbutton10_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton10 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
FolderSupport = uigetdir(cd,'MATLAB Root Directory');
set(handles.edit19,'string',FolderSupport);




function edit19_Callback(hObject, eventdata, handles)
% hObject    handle to edit19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit19 as text
%        str2double(get(hObject,'String')) returns contents of edit19 as a double


% --- Executes during object creation, after setting all properties.
function edit19_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit19 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit20_Callback(hObject, eventdata, handles)
% hObject    handle to edit20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit20 as text
%        str2double(get(hObject,'String')) returns contents of edit20 as a double


% --- Executes during object creation, after setting all properties.
function edit20_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit20 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit21_Callback(hObject, eventdata, handles)
% hObject    handle to edit21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit21 as text
%        str2double(get(hObject,'String')) returns contents of edit21 as a double


% --- Executes during object creation, after setting all properties.
function edit21_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit21 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton11.
function pushbutton11_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton11 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function edit22_Callback(hObject, eventdata, handles)
% hObject    handle to edit22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit22 as text
%        str2double(get(hObject,'String')) returns contents of edit22 as a double


% --- Executes during object creation, after setting all properties.
function edit22_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit22 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit23_Callback(hObject, eventdata, handles)
% hObject    handle to edit23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit23 as text
%        str2double(get(hObject,'String')) returns contents of edit23 as a double


% --- Executes during object creation, after setting all properties.
function edit23_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit23 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit24_Callback(hObject, eventdata, handles)
% hObject    handle to edit24 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit24 as text
%        str2double(get(hObject,'String')) returns contents of edit24 as a double


% --- Executes during object creation, after setting all properties.
function edit24_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit24 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit25_Callback(hObject, eventdata, handles)
% hObject    handle to edit25 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit25 as text
%        str2double(get(hObject,'String')) returns contents of edit25 as a double


% --- Executes during object creation, after setting all properties.
function edit25_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit25 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit26_Callback(hObject, eventdata, handles)
% hObject    handle to edit26 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit26 as text
%        str2double(get(hObject,'String')) returns contents of edit26 as a double


% --- Executes during object creation, after setting all properties.
function edit26_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit26 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit35_Callback(hObject, eventdata, handles)
% hObject    handle to edit35 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit35 as text
%        str2double(get(hObject,'String')) returns contents of edit35 as a double


% --- Executes during object creation, after setting all properties.
function edit35_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit35 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------
function uibuttongroup1_ButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to uibuttongroup1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(uibuttongroup1_ButtonDownFcn,'SelectedObject',ZType_2);  % Set the object.



% --- Executes on selection change in listbox1.
function listbox1_Callback(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox1


% --- Executes during object creation, after setting all properties.
function listbox1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox1.
function checkbox1_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox1
