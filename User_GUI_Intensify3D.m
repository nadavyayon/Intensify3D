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

% Last Modified by GUIDE v2.5 13-Feb-2018 12:38:31

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


% --- Executes on button press in StackFolderBrowse.
function StackFolderBrowse_Callback(hObject, eventdata, handles)
% hObject    handle to StackFolderBrowse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Folder = uigetdir(cd,'MATLAB Root Directory');
set(handles.StackFolder,'string',Folder);
import java.lang.*;
r=Runtime.getRuntime;
ncpu=r.availableProcessors;
set(handles.text15, 'String', num2str(ncpu/2));

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
close(h);

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
TissueDetection = find([get(handles.TissueDetection_1,'Value'), get(handles.TissueDetection_2,'Value'), get(handles.TissueDetection_3,'Value')]);
NormType = find([get(handles.ZType_1,'Value'), get(handles.ZType_2,'Value'), get(handles.ZType_3,'Value'), get(handles.ZType_4,'Value')]);

clc
addpath(pwd);
MSFS = str2num(get(handles.MSFS,'string'));
Threads = str2num(get(handles.Threads,'string'));
if mod(MSFS,2)==0; set(handles.MSFS,'string',num2str(MSFS+1)); end;
try
    Intensify3D(get(handles.StackFolder,'string'),...
    str2num(get(handles.StartIm,'string')),str2num(get(handles.EndIm,'string')),...
    str2num(get(handles.NSTD,'string')),str2num(get(handles.MBI,'string')),str2num(get(handles.MSFS,'string')),...
    TissueDetection-1,NormType,Threads,str2num(get(handles.ImageToShow,'string')),handles);
catch me
    set(handles.Execute,'FontSize',20);
    set(handles.Execute,'ForegroundColor','black');
    set(handles.Execute,'String','Start');
    msgbox('Something went wrong, Check that all parameters are filled in correctly.')
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
% hObject    handle to ZType_1 (see GCBO)
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
function StackFolderBrowse_CreateFcn(hObject, eventdata, handles)
% hObject    handle to StackFolderBrowse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in ImageShow.
function ImageShow_Callback(hObject, eventdata, handles)
% hObject    handle to ImageShow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Folder = get(handles.StackFolder,'string');
cd(Folder);
FolderInfo = dir('*.tif');
FirstImageFileIndex = find([FolderInfo.bytes]>10*1000,1,'first');
set(handles.StartIm,'String','1');
Image2Show = str2num(get(handles.ImageToShow,'string'));
set(handles.ImageName, 'String', FolderInfo(FirstImageFileIndex+Image2Show-1).name);
LastImageFileIndex = find([FolderInfo.bytes]>10*1000,1,'last');
set(handles.EndIm,'String',num2str(LastImageFileIndex-FirstImageFileIndex+1));
axes(handles.axes1);
I = imread([Folder,'' filesep '',FolderInfo(FirstImageFileIndex+Image2Show-1).name]);
RI = imref2d(size(I));
handles.I=I;
guidata(hObject,handles)
MaxBsckgroundGuess = quantile(I(:),0.995);
set(handles.MBI,'string',num2str(MaxBsckgroundGuess))
set(handles.MSFS,'string',num2str(round(max(size(I))/5)));
imshow(I,RI);

set(handles.Threads,'string',num2str(str2num(get(handles.text15,'string'))-1));
set(handles.MaxInt,'string',num2str(max(I(:))));
set(handles.MinInt,'string',num2str(min(I(:))));
Factor = single(MaxBsckgroundGuess)/single(max(I(:)));
handles.Factor=Factor;
guidata(hObject,handles)
set(handles.sliderThreshold,'Value',Factor);
imcontrast
ThresholdRed = round(256*Factor);
axes(handles.axes1);
NewMap=colormap;
handles.CM=colormap;
guidata(hObject,handles)
NewMap(ThresholdRed+1:end,:) = repmat([1 0 0],256-ThresholdRed,1);
colormap(NewMap);


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
set(handles.MBI,'string',num2str(usermaxtrue))
Factor = usermax;
if Factor == 1; Factor = 255/256; end
if Factor == 0; Factor = 1/256; end

ThresholdRed = round(256*Factor);
axes(handles.axes1);
NewMap=handles.CM;
colormap(NewMap);
NewMap=colormap;
NewMap(ThresholdRed+1:end,:) = repmat([1 0 0],256-ThresholdRed,1);
colormap(NewMap);



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

% 
% % --- Executes on slider movement.
% function MinBrighness_Callback(hObject, eventdata, handles)
% % hObject    handle to MinBrighness (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% % Min = str2num(get(handles.MinInt,'string'));
% MinFactor = get(hObject,'Value');
% MinBrighness = round(256*MinFactor);
% NewMap=handles.CM;
% NewMap(1:MinBrighness,:) = repmat([0 0 0],MinBrighness,1);
% colormap(NewMap)
% 
% % Hints: get(hObject,'Value') returns position of slider
% %        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
% 
% 
% % --- Executes during object creation, after setting all properties.
% function MinBrighness_CreateFcn(hObject, eventdata, handles)
% % hObject    handle to MinBrighness (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    empty - handles not created until after all CreateFcns called
% 
% % Hint: slider controls usually have a light gray background.
% if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
%     set(hObject,'BackgroundColor',[.9 .9 .9]);
% end
% 
% 
% % --- Executes on slider movement.
% function MaxBrighness_Callback(hObject, eventdata, handles)
% % hObject    handle to MaxBrighness (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    structure with handles and user data (see GUIDATA)
% MaxFactor = get(hObject,'Value');
% MaxBrighness = round(256*MaxFactor);
% MaxNewMap = gray(MaxBrighness);
% MaxNewMap = [MaxNewMap; ones(256-MaxBrighness,3);]
% colormap(MaxNewMap)
% 
% % Hints: get(hObject,'Value') returns position of slider
% %        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
% 
% 
% % --- Executes during object creation, after setting all properties.
% function MaxBrighness_CreateFcn(hObject, eventdata, handles)
% % hObject    handle to MaxBrighness (see GCBO)
% % eventdata  reserved - to be defined in a future version of MATLAB
% % handles    empty - handles not created until after all CreateFcns called
% 
% % Hint: slider controls usually have a light gray background.
% if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
%     set(hObject,'BackgroundColor',[.9 .9 .9]);
% end


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
