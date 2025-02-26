function varargout = vuOnePaneViewer(varargin)
% VUONEPANEVIEWER M-file for vuOnePaneViewer.fig
%      VUONEPANEVIEWER, by itself, creates a new VUONEPANEVIEWER or raises the existing
%      singleton*.
%
%      H = VUONEPANEVIEWER returns the handle to a new VUONEPANEVIEWER or the handle to
%      the existing singleton*.
%
%      VUONEPANEVIEWER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in VUONEPANEVIEWER.M with the given input arguments.
%
%      VUONEPANEVIEWER('Property','Value',...) creates a new VUONEPANEVIEWER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before vuOnePaneViewer_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to vuOnePaneViewer_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help vuOnePaneViewer

% Last Modified by GUIDE v2.5 28-Oct-2010 15:14:01

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @vuOnePaneViewer_OpeningFcn, ...
                   'gui_OutputFcn',  @vuOnePaneViewer_OutputFcn, ...
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


% --- Executes just before vuOnePaneViewer is made visible.
function vuOnePaneViewer_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to vuOnePaneViewer (see VARARGIN)

% Choose default command line output for vuOnePaneViewer
handles.output = hObject;

% Check/Make input a Meta Image
if length(varargin) < 1
	error('MATLAB:vuOnePaneViewer:NotEnoughInputs', 'Not enough input arguments.');
end

image = varargin{1};
if (~isStructure(image))
    image = vuGenerateMetaImage(image);
end

% Calculate middle slices
midSlice = floor(image.Dims/2);

setappdata(gcf,'measure',0);

% Get and parse options
p = inputParser;
p.addParamValue('pane',3,@(x) isa(x,'double'));
p.addParamValue('slice',midSlice(3),@(x) isa(x,'double'));
p.addParamValue('dynamic',1,@(x) isa(x,'double'));
p.FunctionName='vuOnePaneViewer';
p.parse(varargin{2:end});

pane = p.Results.pane;
slice = p.Results.slice;
dynamic = p.Results.dynamic;

% Check Pane
if pane < 1 || pane > 3
    pane = 3;
end

% Check Slice
if slice < 1 || slice > image.Dims(3)
    slice = midSlice(3);
end

% Check Dynamic
if length(image.Dims) >= 4 && (dynamic < 1 || dynamic > image.Dims(4))
    dynamic = 1;
end

% Permute our image based on which pane is viewed
if pane == 1
    image.Data = permute(image.Data,[1 3 2 4]);
    image.Dims = [image.Dims(2) image.Dims(3) image.Dims(1)];
    image.Spc = [image.Spc(2) image.Spc(3) image.Spc(1)];
    image.Origin = [image.Origin(2) image.Origin(3) image.Origin(1)];
elseif pane == 2
    image.Data = permute(image.Data,[2 3 1 4]);
    image.Dims = [image.Dims(1) image.Dims(3) image.Dims(2)];
    image.Spc = [image.Spc(1) image.Spc(3) image.Spc(2)];
    image.Origin = [image.Origin(1) image.Origin(3) image.Origin(2)];
end

setappdata(gcf,'image',image);
setappdata(gcf,'dynamic',dynamic);

% Setup the sliders
set(handles.sliceSlider,'Min',0.999)
set(handles.sliceSlider,'Max',image.Dims(3))
set(handles.sliceSlider,'SliderStep',[1./(image.Dims(3)-0.999) 1./(image.Dims(3)-0.999)]);
set(handles.sliceSlider,'Value',slice);

% Show the image
axes(handles.mainFig)
hImg = imshow(squeeze(image.Data(:,:,slice,dynamic)),[]);
setappdata(gcf,'hImg',hImg);

% Windowing Sliders
minValue = min(image.Data(:));
maxValue = max(image.Data(:));
if ((minValue-maxValue)==0)
    minValue = minValue-0.001;
end
setappdata(gcf,'curMin',minValue);
setappdata(gcf,'curMax',maxValue);
set(handles.lowerSlide,'Min', minValue);
set(handles.lowerSlide,'Max', maxValue-0.0001);
set(handles.upperSlide,'Min', minValue+0.0001);
set(handles.upperSlide,'Max', maxValue);

% Set slider to 2%-98%
sortedImage = sort(image.Data(:));
lowerValue = sortedImage(round(0.02*length(sortedImage)));
upperValue = sortedImage(round(0.98*length(sortedImage)));
if(lowerValue == upperValue)
    sortedImage(sortedImage==upperValue) = [];
    if (~isempty(sortedImage))
        upperValue = sortedImage(1);
    else
        lowerValue = minValue;
        upperValue = maxValue;
    end
end
setappdata(gcf,'curMin',lowerValue);
setappdata(gcf,'curMax',upperValue);
set(handles.lowerSlide,'Value', lowerValue);
set(handles.upperSlide,'Value', upperValue);
set(handles.lowerText,'String',num2str(lowerValue));
set(handles.upperText,'String',num2str(upperValue));

axes(handles.mainFig)
caxis([lowerValue upperValue]);

% Initialize axes limits
axes(handles.mainFig);
set(handles.mainFig,'xlim',[0.5 image.Dims(1)+0.5],'ylim',[0.5 image.Dims(2)+0.5]);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes vuOnePaneViewer wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = vuOnePaneViewer_OutputFcn(hObject, eventdata, handles) 
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
dynamic = getappdata(gcf,'dynamic');
xlimits = xlim(handles.mainFig);
ylimits = ylim(handles.mainFig);

set(handles.sliceSlider,'Value',curVal)
axes(handles.mainFig)
hImg = imshow(squeeze(image.Data(:,:,curVal,dynamic)),[]);
set(handles.mainFig,'xlim',xlimits,'ylim',ylimits);
setappdata(gcf,'hImg',hImg);
curMin = getappdata(gcf,'curMin');
curMax = getappdata(gcf,'curMax');
caxis([curMin curMax]);

mapString = get(handles.mapMenu,'String');
mapValue = get(handles.mapMenu,'Value');
eval(['colormap ' cell2mat(mapString(mapValue))]);

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
pan off
datacursormode off
zoom off

% --- Executes on button press in dataRadio.
function dataRadio_Callback(hObject, eventdata, handles)
% hObject    handle to dataRadio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of dataRadio

axes(handles.mainFig)
pan off
zoom off
datacursormode on


% --- Executes on button press in zoomRadio.
function zoomRadio_Callback(hObject, eventdata, handles)
% hObject    handle to zoomRadio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of zoomRadio

axes(handles.mainFig)
pan off
zoom on
datacursormode off


% --- Executes on slider movement.
function upperSlide_Callback(hObject, eventdata, handles)
% hObject    handle to upperSlide (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

curMin = getappdata(gcf,'curMin');
curMax = getappdata(gcf,'curMax');
usrMax = get(handles.upperSlide,'Value');
usrMax = round(usrMax);
if (usrMax <= curMin)
    set(handles.upperSlide,'Value',curMax);
else
    setappdata(gcf,'curMax',usrMax);
    set(handles.upperText,'String',num2str(usrMax));
end
curMin = getappdata(gcf,'curMin');
curMax = getappdata(gcf,'curMax');
axes(handles.mainFig)
caxis([curMin curMax]);

% --- Executes during object creation, after setting all properties.
function upperSlide_CreateFcn(hObject, eventdata, handles)
% hObject    handle to upperSlide (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function lowerSlide_Callback(hObject, eventdata, handles)
% hObject    handle to lowerSlide (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
curMin = getappdata(gcf,'curMin');
curMax = getappdata(gcf,'curMax');
usrMin = get(handles.lowerSlide,'Value');
usrMin = round(usrMin);
if (usrMin >= curMax)
    set(handles.lowerSlide,'Value',curMin);
else
    setappdata(gcf,'curMin',usrMin);
    set(handles.lowerText,'String',num2str(usrMin));
end
curMin = getappdata(gcf,'curMin');
curMax = getappdata(gcf,'curMax');
axes(handles.mainFig)
caxis([curMin curMax]);

% --- Executes during object creation, after setting all properties.
function lowerSlide_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lowerSlide (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function upperText_Callback(hObject, eventdata, handles)
% hObject    handle to upperText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of upperText as text
%        str2double(get(hObject,'String')) returns contents of upperText as a double
curMin = getappdata(gcf,'curMin');
curMax = getappdata(gcf,'curMax');
usrMax = get(handles.upperText,'String');
usrMax = str2num(usrMax);
if (~isempty(usrMax))
    if (usrMax < curMin)
        usrMax = curMin+1;
        set(handles.upperText,'String',num2str(curMin+1));
    end
    setappdata(gcf,'curMax',usrMax);
    set(handles.upperSlide,'Value',usrMax);
else
    set(handles.upperText,'String',num2str(curMax));
end
curMin = getappdata(gcf,'curMin');
curMax = getappdata(gcf,'curMax');
axes(handles.mainFig)
caxis([curMin curMax]);

% --- Executes during object creation, after setting all properties.
function upperText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to upperText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function lowerText_Callback(hObject, eventdata, handles)
% hObject    handle to lowerText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lowerText as text
%        str2double(get(hObject,'String')) returns contents of lowerText as a double
curMin = getappdata(gcf,'curMin');
curMax = getappdata(gcf,'curMax');
usrMin = get(handles.lowerText,'String');
usrMin = str2num(usrMin);
if (~isempty(usrMin))
    if (usrMin > curMax)
        usrMin = curMax-1;
        set(handles.lowerText,'String',num2str(curMax-1));
    end
    setappdata(gcf,'curMin',usrMin);
    set(handles.lowerSlide,'Value',usrMin);
else
    set(handles.lowerText,'String',num2str(curMin));
end
curMin = getappdata(gcf,'curMin');
curMax = getappdata(gcf,'curMax');
axes(handles.mainFig)
caxis([curMin curMax]);

% --- Executes during object creation, after setting all properties.
function lowerText_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lowerText (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


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
curVal = round(get(handles.sliceSlider,'Value'));
image = getappdata(gcf,'image');
dynamic = getappdata(gcf,'dynamic');
axes(handles.mainFig)
hImg = getappdata(gcf,'hImg');
pan off
zoom off
[x,y] = ginput(2);

hline = imdistline(gca,x,y);
api = iptgetapi(hline);
api.setLabelTextFormatter('%02.1f'); 

% Uncomment below for spacing application -- comment above
% % Convert XData and YData to meters using conversion factor.
% XDataInMeters = get(hImg,'XData')*image.Spc(1); 
% YDataInMeters = get(hImg,'YData')*image.Spc(2);
%      
% % Set XData and YData of image to reflect desired units.    
% set(hImg,'XData',XDataInMeters,'YData',YDataInMeters);    
% set(gca,'XLim',XDataInMeters,'YLim',YDataInMeters);
% 
% % Specify initial position of distance tool on Harvard Bridge.
% hline = imdistline(gca,x*image.Spc(1),y*image.Spc(2));
% api = iptgetapi(hline);
% api.setLabelTextFormatter('%02.0fmm'); 

% --- Executes on button press in panRadio.
function panRadio_Callback(hObject, eventdata, handles)
% hObject    handle to panRadio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of panRadio
axes(handles.mainFig)
zoom off
datacursormode off
pan on
