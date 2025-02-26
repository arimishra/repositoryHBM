function varargout = vuOnePaneBlockROIDrawer(varargin)
% VUONEPANEBLOCKROIDRAWER M-file for vuOnePaneBlockROIDrawer.fig
%      VUONEPANEBLOCKROIDRAWER, by itself, creates a new VUONEPANEBLOCKROIDRAWER or raises the existing
%      singleton*.
%
%      H = VUONEPANEBLOCKROIDRAWER returns the handle to a new VUONEPANEBLOCKROIDRAWER or the handle to
%      the existing singleton*.
%
%      VUONEPANEBLOCKROIDRAWER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in VUONEPANEBLOCKROIDRAWER.M with the given input arguments.
%
%      VUONEPANEBLOCKROIDRAWER('Property','Value',...) creates a new VUONEPANEBLOCKROIDRAWER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before vuOnePaneROIDrawer_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to vuOnePaneBlockROIDrawer_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help vuOnePaneBlockROIDrawer

% Last Modified by GUIDE v2.5 26-Aug-2010 12:26:34

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @vuOnePaneBlockROIDrawer_OpeningFcn, ...
                   'gui_OutputFcn',  @vuOnePaneBlockROIDrawer_OutputFcn, ...
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


% --- Executes just before vuOnePaneBlockROIDrawer is made visible.
function vuOnePaneBlockROIDrawer_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to vuOnePaneBlockROIDrawer (see VARARGIN)

% Choose default command line output for vuOnePaneBlockROIDrawer
handles.output = hObject;

% Check/Make input a Meta Image
if length(varargin) < 1
	error('MATLAB:vuOnePaneROIDrawer:NotEnoughInputs', 'Not enough input arguments.');
end

image = varargin{1};
if (~isStructure(image))
    image = vuGenerateMetaImage(image);
end

setappdata(gcf,'measure',0);

% Default
pane = 3;

% Which viewer
if length(varargin) >= 2
    try
        pane = varargin{2};
        if pane < 1 || pane > 3
            pane = 3;
        end
    catch
        pane = 3;
    end
else
    pane = 3;
end

% Permute our image based on which pane is viewed
if pane == 1
    image.Data = permute(image.Data,[1 3 2]);
    image.Dims = [image.Dims(2) image.Dims(3) image.Dims(1)];
    image.Spc = [image.Spc(2) image.Spc(3) image.Spc(1)];
    image.Origin = [image.Origin(2) image.Origin(3) image.Origin(1)];
elseif pane == 2
    image.Data = permute(image.Data,[2 3 1]);
    image.Dims = [image.Dims(1) image.Dims(3) image.Dims(2)];
    image.Spc = [image.Spc(1) image.Spc(3) image.Spc(2)];
    image.Origin = [image.Origin(1) image.Origin(3) image.Origin(2)];
end

setappdata(gcf,'image',image);
setappdata(gcf,'view',1);

% Calculate middle slices
midSlice = floor(image.Dims/2);
slice = midSlice(3);

% Which slice
if length(varargin) == 3
    try
        slice = varargin{3};
        if slice < 1 || slice > image.Dims(3)
            slice = midSlice(3);
        end
    catch
        slice = midSlice(3);
    end
else
    slice = midSlice(3);
end

% Setup the sliders
set(handles.sliceSlider,'Min',0.999)
set(handles.sliceSlider,'Max',image.Dims(3))
set(handles.sliceSlider,'SliderStep',[1./(image.Dims(3)-0.999) 1./(image.Dims(3)-0.999)]);
set(handles.sliceSlider,'Value',slice);

setappdata(gcf,'ROI',[]);
setappdata(gcf,'numberOfROIs',0);
setappdata(gcf,'ROIDepth',3);
setappdata(gcf,'showROI',0);

% Check Depth
if (image.Dims(3) < 3)
    set(handles.maxRadio,'Value',1)
    set(handles.threeRadio,'Enable','off')
    set(handles.fiveRadio,'Enable','off')
elseif (image.Dims(3) < 5)
    set(handles.maxRadio,'Value',1)
    set(handles.threeRadio,'Enable','on')
    set(handles.fiveRadio,'Enable','off')
else
    set(handles.threeRadio,'Enable','on')
    set(handles.fiveRadio,'Enable','on')
end

% Show the image
axes(handles.mainFig)
imshow(squeeze(image.Data(:,:,slice,1)),[])

% Windowing Sliders
minValue = min(image.Data(:));
maxValue = max(image.Data(:));
if ((minValue-maxValue)==0)
    minValue = minValue-0.001;
end
setappdata(gcf,'curMin',minValue);
setappdata(gcf,'curMax',maxValue);

axes(handles.mainFig)
caxis([minValue maxValue]);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes vuOnePaneBlockROIDrawer wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = vuOnePaneBlockROIDrawer_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on slider movement.
function sliceSlider_Callback(hObject, eventdata, handles)
% hObject    handle to sliceSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

curVal = round(get(handles.sliceSlider,'Value'));

image = getappdata(gcf,'image');

set(handles.sliceSlider,'Value',curVal)
axes(handles.mainFig)
imshow(squeeze(image.Data(:,:,curVal,1)),[])
curMin = getappdata(gcf,'curMin');
curMax = getappdata(gcf,'curMax');
caxis([curMin curMax]);

mapString = get(handles.mapMenu,'String');
mapValue = get(handles.mapMenu,'Value');
eval(['colormap ' cell2mat(mapString(mapValue))]);

showROI = getappdata(gcf,'showROI');
if (showROI)
    showROIFunction(handles)
end

% --- Executes during object creation, after setting all properties.
function sliceSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliceSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% Function to Check the Structure of the image
function isStruct = isStructure(X)

if(isstruct(X))
  %Check for meta image structure
  if(isfield(X,'Data') &&  isfield(X,'Dims') && ...
          isfield(X,'Spc') && isfield(X,'Origin'))
      isStruct = true;
      if (length(X.Dims) < 3 || length(X.Dims) > 4)
          error('MATLAB:vuThreePaneViewer:UnknownDims', 'vuThreePaneViewer can only handle images of 3 or 4 dimensions.');
      end
  else
      error('MATLAB:vuThreePaneViewer:InvalidStruct', 'The input image structure is not valid.');
  end
else
    if (ndims(X)~=3 && ndims(X)~=4)
          error('MATLAB:vuThreePaneViewer:UnknownDims', 'vuThreePaneViewer can only handle images of 3 or 4 dimensions.');
    end
    isStruct = false;
end

return


% --- Executes on button press in noneRadio.
function noneRadio_Callback(hObject, eventdata, handles)
% hObject    handle to noneRadio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of noneRadio

axes(handles.mainFig)
datacursormode off
zoom off

% --- Executes on button press in dataRadio.
function dataRadio_Callback(hObject, eventdata, handles)
% hObject    handle to dataRadio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of dataRadio

axes(handles.mainFig)
zoom off
datacursormode on


% --- Executes on button press in zoomRadio.
function zoomRadio_Callback(hObject, eventdata, handles)
% hObject    handle to zoomRadio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of zoomRadio

axes(handles.mainFig)
zoom on
datacursormode off

% --- Executes on selection change in mapMenu.
function mapMenu_Callback(hObject, eventdata, handles)
% hObject    handle to mapMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns mapMenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from mapMenu
mapString = get(handles.mapMenu,'String');
mapValue = get(handles.mapMenu,'Value');
axes(handles.mainFig)
eval(['colormap ' cell2mat(mapString(mapValue))]);

% --- Executes during object creation, after setting all properties.
function mapMenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to mapMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in measureRadio.
function measureRadio_Callback(hObject, eventdata, handles)
% hObject    handle to measureRadio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of measureRadio

axes(handles.mainFig)
[x,y] = ginput(2);
h = imdistline(gca,x,y);


% --- Executes on button press in windowLevelButton.
function windowLevelButton_Callback(hObject, eventdata, handles)
% hObject    handle to windowLevelButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

image = getappdata(gcf,'image');
curVal = round(get(handles.sliceSlider,'Value'));

% Windowing Sliders
slice = image.Data(:,:,curVal);
minValue = min(slice(:));
maxValue = max(slice(:));
if ((minValue-maxValue)==0)
    minValue = minValue-0.001;
end

% Compare to current Value
minValue = max(minValue,getappdata(gcf,'curMin'));
maxValue = min(maxValue,getappdata(gcf,'curMax'));

axes(handles.mainFig)
caxis([minValue maxValue]);
waitfor(imcontrast(gca))
newValues = get(gca,'CLim');
setappdata(gcf,'curMin',newValues(1));
setappdata(gcf,'curMax',newValues(2));


% --- Executes during object creation, after setting all properties.
function windowLevelButton_CreateFcn(hObject, eventdata, handles)
% hObject    handle to windowLevelButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in xyRadio.
function xyRadio_Callback(hObject, eventdata, handles)
% hObject    handle to xyRadio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of xyRadio
image = getappdata(gcf,'image');
ROI = getappdata(gcf,'ROI');
if (getappdata(gcf,'view')==2)
    image.Data = permute(image.Data,[3 1 2]);
    if (~isempty(ROI))
        ROI = permute(ROI,[3 1 2]);
    end
elseif (getappdata(gcf,'view')==3)
    image.Data = permute(image.Data,[1 3 2]);
    if (~isempty(ROI))
        ROI = permute(ROI,[1 3 2]);
    end
end
setappdata(gcf,'image',image);
setappdata(gcf,'ROI',ROI);
setappdata(gcf,'view',1);
midSlice = floor(image.Dims/2);
slice = midSlice(3);

% Setup the sliders
set(handles.sliceSlider,'Min',0.999)
set(handles.sliceSlider,'Max',image.Dims(3))
set(handles.sliceSlider,'SliderStep',[1./(image.Dims(3)-0.999) 1./(image.Dims(3)-0.999)]);
set(handles.sliceSlider,'Value',slice);

sliceSlider_Callback(hObject, eventdata, handles)

% Check Depth
if (image.Dims(3) < 3)
    set(handles.maxRadio,'Value',1)
    set(handles.threeRadio,'Enable','off')
    set(handles.fiveRadio,'Enable','off')
elseif (image.Dims(3) < 5)
    set(handles.maxRadio,'Value',1)
    set(handles.threeRadio,'Enable','on')
    set(handles.fiveRadio,'Enable','off')
else
    set(handles.threeRadio,'Enable','on')
    set(handles.fiveRadio,'Enable','on')
end

% Check ROI Depth
if (get(handles.maxRadio,'Value'))
    maxRadio_Callback(hObject, eventdata, handles)
elseif (get(handles.otherRadio,'Value'))
    otherRadio_Callback(hObject, eventdata, handles)    
end

% --- Executes on button press in zxRadio.
function zxRadio_Callback(hObject, eventdata, handles)
% hObject    handle to zxRadio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of zxRadio
image = getappdata(gcf,'image');
ROI = getappdata(gcf,'ROI');
if (getappdata(gcf,'view')==1)
    image.Data = permute(image.Data,[2 3 1]);
    if (~isempty(ROI))
        ROI = permute(ROI,[2 3 1]);
    end
elseif (getappdata(gcf,'view')==3)
    image.Data = permute(image.Data,[3 2 1]);
    if (~isempty(ROI))
        ROI = permute(ROI,[3 2 1]);
    end
end
setappdata(gcf,'image',image);
setappdata(gcf,'ROI',ROI);
setappdata(gcf,'view',2);
midSlice = floor(image.Dims/2);
slice = midSlice(2);

% Setup the sliders
set(handles.sliceSlider,'Min',0.999)
set(handles.sliceSlider,'Max',image.Dims(2))
set(handles.sliceSlider,'SliderStep',[1./(image.Dims(2)-0.999) 1./(image.Dims(2)-0.999)]);
set(handles.sliceSlider,'Value',slice);

sliceSlider_Callback(hObject, eventdata, handles)

% Check Depth
if (image.Dims(2) < 3)
    set(handles.maxRadio,'Value',1)
    set(handles.threeRadio,'Enable','off')
    set(handles.fiveRadio,'Enable','off')
elseif (image.Dims(2) < 5)
    set(handles.maxRadio,'Value',1)
    set(handles.threeRadio,'Enable','on')
    set(handles.fiveRadio,'Enable','off')
else
    set(handles.threeRadio,'Enable','on')
    set(handles.fiveRadio,'Enable','on')
end

% Check ROI Depth
if (get(handles.maxRadio,'Value'))
    maxRadio_Callback(hObject, eventdata, handles)
elseif (get(handles.otherRadio,'Value'))
    otherRadio_Callback(hObject, eventdata, handles)    
end

% --- Executes on button press in zyRadio.
function zyRadio_Callback(hObject, eventdata, handles)
% hObject    handle to zyRadio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of zyRadio
image = getappdata(gcf,'image');
ROI = getappdata(gcf,'ROI');
if (getappdata(gcf,'view')==1)
    image.Data = permute(image.Data,[1 3 2]);
    if (~isempty(ROI))
        ROI = permute(ROI,[1 3 2]);
    end
elseif (getappdata(gcf,'view')==2)
    image.Data = permute(image.Data,[3 2 1]);
    if (~isempty(ROI))
        ROI = permute(ROI,[3 2 1]);
    end
end
setappdata(gcf,'image',image);
setappdata(gcf,'ROI',ROI);
setappdata(gcf,'view',3);

midSlice = floor(image.Dims/2);
slice = midSlice(1);



% Setup the sliders
set(handles.sliceSlider,'Min',0.999)
set(handles.sliceSlider,'Max',image.Dims(1))
set(handles.sliceSlider,'SliderStep',[1./(image.Dims(1)-0.999) 1./(image.Dims(1)-0.999)]);
set(handles.sliceSlider,'Value',slice);

sliceSlider_Callback(hObject, eventdata, handles)

% Check Depth
if (image.Dims(1) < 3)
    set(handles.maxRadio,'Value',1)
    set(handles.threeRadio,'Enable','off')
    set(handles.fiveRadio,'Enable','off')
elseif (image.Dims(1) < 5)
    set(handles.maxRadio,'Value',1)
    set(handles.threeRadio,'Enable','on')
    set(handles.fiveRadio,'Enable','off')
else
    set(handles.threeRadio,'Enable','on')
    set(handles.fiveRadio,'Enable','on')
end

% Check ROI Depth
if (get(handles.maxRadio,'Value'))
    maxRadio_Callback(hObject, eventdata, handles)
elseif (get(handles.otherRadio,'Value'))
    otherRadio_Callback(hObject, eventdata, handles)    
end


function otherEdit_Callback(hObject, eventdata, handles)
% hObject    handle to otherEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of otherEdit as text
%        str2double(get(hObject,'String')) returns contents of otherEdit as a double
value = get(handles.otherEdit,'String');
oldValue = getappdata(gcf,'ROIDepth');
if (isempty(str2num(value)))
    set(handles.otherEdit,'String',oldValue);
elseif (str2num(value)<1)
    set(handles.otherEdit,'String',oldValue);
elseif (str2num(value)>get(handles.sliceSlider,'Max'))
    set(handles.otherEdit,'String',oldValue);
else
    set(handles.otherEdit,'String',round(str2num(value)));
    setappdata(gcf,'ROIDepth',str2num(value))
end

% --- Executes during object creation, after setting all properties.
function otherEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to otherEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in threeRadio.
function threeRadio_Callback(hObject, eventdata, handles)
% hObject    handle to threeRadio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of threeRadio
setappdata(gcf,'ROIDepth',3);
set(handles.otherEdit,'Visible','off');

% --- Executes on button press in fiveRadio.
function fiveRadio_Callback(hObject, eventdata, handles)
% hObject    handle to fiveRadio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of fiveRadio
setappdata(gcf,'ROIDepth',5);
set(handles.otherEdit,'Visible','off');

% --- Executes on button press in maxRadio.
function maxRadio_Callback(hObject, eventdata, handles)
% hObject    handle to maxRadio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of maxRadio
setappdata(gcf,'ROIDepth',get(handles.sliceSlider,'Max'));
set(handles.otherEdit,'Visible','off');

% --- Executes on button press in otherRadio.
function otherRadio_Callback(hObject, eventdata, handles)
% hObject    handle to otherRadio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of otherRadio
set(handles.otherEdit,'Visible','on');
set(handles.otherEdit,'String','3');
setappdata(gcf,'ROIDepth',3);


% --- Executes on button press in createROIButton.
function createROIButton_Callback(hObject, eventdata, handles)
% hObject    handle to createROIButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get our points
k = waitforbuttonpress;
point1 = get(gca,'CurrentPoint');    
finalRect = rbbox;                  
point2 = get(gca,'CurrentPoint');   
point1 = point1(1,1:2);             
point2 = point2(1,1:2);
p1 = round(min(point1,point2)); 
p2 = round(max(point1,point2));

% Create ROI
image = getappdata(gcf,'image');
ROI = zeros(size(image.Data));
ROIDepth = getappdata(gcf,'ROIDepth');

% If Max use all slices
if (get(handles.maxRadio,'Value'))
    ROI(p1(2):p2(2),p1(1):p2(1),:) = 1;
else
    % Current Slice
    curVal = round(get(handles.sliceSlider,'Value'));
    
    % Find our plus/minus slices
    if (mod(ROIDepth,2))
        firstValue = (ROIDepth-1)/2;
        lastValue = firstValue;
    else
        firstValue = ROIDepth/2;
        lastValue = firstValue-1;
    end
    
    ROI(p1(2):p2(2),p1(1):p2(1),max(curVal-firstValue,1):min(curVal+lastValue,end)) = 1;
end

setappdata(gcf,'ROI',ROI);
setappdata(gcf,'numberOfROIs',getappdata(gcf,'numberOfROIs')+1);
setappdata(gcf,'showROI',1);
set(handles.showROIButton,'Visible','on')
set(handles.showROIButton,'String','Hide ROI')
showROIFunction(handles);
assignin('base',['Created_ROI_' num2str(getappdata(gcf,'numberOfROIs'))],ROI);

% --- Executes on button press in showROIButton.
function showROIButton_Callback(hObject, eventdata, handles)
% hObject    handle to showROIButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
showROI = getappdata(gcf,'showROI');

% Invert show ROI flag
showROI = ~showROI;
setappdata(gcf,'showROI',showROI)

% Change Button Text
if (showROI)
    set(handles.showROIButton,'String','Hide ROI')

    % Call show ROI function
    showROIFunction(handles);
    
else
    set(handles.showROIButton,'String','Show ROI')
    
    % Redraw slice
    sliceSlider_Callback(hObject, eventdata, handles)
end

function showROIFunction(handles)

% Get ROI
ROI = getappdata(gcf,'ROI');
ROI(ROI==0) = 0.4;

axes(handles.mainFig);
curVal = round(get(handles.sliceSlider,'Value'));
alpha(squeeze(ROI(:,:,curVal)))
