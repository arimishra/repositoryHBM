function varargout = vuOnePaneViewerPointSelector(varargin)
% VUONEPANEVIEWERPOINTSELECTOR M-file for vuOnePaneViewerPointSelector.fig
%      VUONEPANEVIEWERPOINTSELECTOR, by itself, creates a new VUONEPANEVIEWERPOINTSELECTOR or raises the existing
%      singleton*.
%
%      H = VUONEPANEVIEWERPOINTSELECTOR returns the handle to a new VUONEPANEVIEWERPOINTSELECTOR or the handle to
%      the existing singleton*.
%
%      VUONEPANEVIEWERPOINTSELECTOR('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in VUONEPANEVIEWERPOINTSELECTOR.M with the given input arguments.
%
%      VUONEPANEVIEWERPOINTSELECTOR('Property','Value',...) creates a new VUONEPANEVIEWERPOINTSELECTOR or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before vuOnePaneViewerPointSelector_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to vuOnePaneViewerPointSelector_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help vuOnePaneViewerPointSelector

% Last Modified by GUIDE v2.5 29-Jun-2009 13:15:31

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @vuOnePaneViewerPointSelector_OpeningFcn, ...
                   'gui_OutputFcn',  @vuOnePaneViewerPointSelector_OutputFcn, ...
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


% --- Executes just before vuOnePaneViewerPointSelector is made visible.
function vuOnePaneViewerPointSelector_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to vuOnePaneViewerPointSelector (see VARARGIN)

% Choose default command line output for vuOnePaneViewerPointSelector
handles.output = hObject;

% Check/Make input a Meta Image
if length(varargin) < 1
	error('MATLAB:vuOnePaneViewer:NotEnoughInputs', 'Not enough input arguments.');
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

% Vessel and Tumor Point array
setappdata(gcf,'vesselPoints',[]);
setappdata(gcf,'tumorPoints',[]);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes vuOnePaneViewerPointSelector wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = vuOnePaneViewerPointSelector_OutputFcn(hObject, eventdata, handles) 
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

% Display points
tumorPoints = getappdata(gcf,'tumorPoints');
vesselPoints = getappdata(gcf,'vesselPoints');
if (~isempty(tumorPoints))
    for i = find(tumorPoints(:,3)==curVal)'
        h = impoint(gca,tumorPoints(i,1),tumorPoints(i,2));
        api = iptgetapi(h);
        api.setColor('red');
    end
end
if (~isempty(vesselPoints))
    for i = find(vesselPoints(:,3)==curVal)'
        h = impoint(gca,vesselPoints(i,1),vesselPoints(i,2));
        api = iptgetapi(h);
        api.setColor('green');
    end
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

% --- Executes on mouse press over figure background, over a disabled or
% --- inactive control, or over an axes background.
function figure1_WindowButtonDownFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

cp = get(gca,'CurrentPoint');
curVal = round(get(handles.sliceSlider,'Value'));

% Make sure within current axis
if (cp(1,1)>=1) && (cp(1,2)>=1) && (cp(1,1)<=max(get(gca,'XLim'))) && (cp(1,2)<=max(get(gca,'YLim')))

    % Adding a tumor point
    if (get(handles.tumorPointRadio,'Value'))
        tumorPoints = getappdata(gcf,'tumorPoints');
        tumorPoints = [tumorPoints ; cp(1,1) cp(1,2) curVal];
        setappdata(gcf,'tumorPoints',tumorPoints);
        h = impoint(gca,cp(1,1),cp(1,2));
        api = iptgetapi(h);
        api.setColor('red');
    % Adding a vessel point
    elseif (get(handles.vesselPointRadio,'Value'))
        vesselPoints = getappdata(gcf,'vesselPoints');
        vesselPoints = [vesselPoints ; cp(1,1) cp(1,2) curVal];
        setappdata(gcf,'vesselPoints',vesselPoints);
        h = impoint(gca,cp(1,1),cp(1,2));
        api = iptgetapi(h);
        api.setColor('green');
    end

end


% --- Executes on button press in clearAllButton.
function clearAllButton_Callback(hObject, eventdata, handles)
% hObject    handle to clearAllButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Clear our points
setappdata(gcf,'tumorPoints',[]);
setappdata(gcf,'vesselPoints',[]);

% Redraw
sliceSlider_Callback(hObject, eventdata, handles)

% --- Executes on button press in clearLastButton.
function clearLastButton_Callback(hObject, eventdata, handles)
% hObject    handle to clearLastButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Check which radio is selected
if (get(handles.tumorPointRadio,'Value'))
    
    % Get the tumor points
    tumorPoints = getappdata(gcf,'tumorPoints');
    % Take off last point
    tumorPoints = tumorPoints(1:end-1,:);    
    % Save
    setappdata(gcf,'tumorPoints',tumorPoints);
    % Redraw
    sliceSlider_Callback(hObject, eventdata, handles)
    
elseif (get(handles.vesselPointRadio,'Value'))
    
    % Get the tumor points
    vesselPoints = getappdata(gcf,'vesselPoints');
    % Take off last point
    vesselPoints = vesselPoints(1:end-1,:);    
    % Save
    setappdata(gcf,'vesselPoints',vesselPoints);
    % Redraw
    sliceSlider_Callback(hObject, eventdata, handles)
    
end


% --- Executes on button press in oneQuarterButton.
function oneQuarterButton_Callback(hObject, eventdata, handles)
% hObject    handle to oneQuarterButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get the tumor points
tumorPoints = getappdata(gcf,'tumorPoints');

% If more than 2 points go to 1/4 mark
if (size(tumorPoints,1)>1)
    
    % Get minimum and maximum slice
    minSlice = min(tumorPoints(:,3));
    maxSlice = max(tumorPoints(:,3));
    
    % Calculate 1/4 mark
    curVal = round(minSlice + (maxSlice - minSlice)/4);
    
    % Set slider and redraw
    set(handles.sliceSlider,'Value',curVal);
    sliceSlider_Callback(hObject, eventdata, handles)
    
end

% --- Executes on button press in twoQuarterButton.
function twoQuarterButton_Callback(hObject, eventdata, handles)
% hObject    handle to twoQuarterButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get the tumor points
tumorPoints = getappdata(gcf,'tumorPoints');

% If more than 2 points go to 1/4 mark
if (size(tumorPoints,1)>1)
    
    % Get minimum and maximum slice
    minSlice = min(tumorPoints(:,3));
    maxSlice = max(tumorPoints(:,3));
    
    % Calculate 1/4 mark
    curVal = round(minSlice + 2*(maxSlice - minSlice)/4);
    
    % Set slider and redraw
    set(handles.sliceSlider,'Value',curVal);
    sliceSlider_Callback(hObject, eventdata, handles)
    
end

% --- Executes on button press in threeQuarterButton.
function threeQuarterButton_Callback(hObject, eventdata, handles)
% hObject    handle to threeQuarterButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get the tumor points
tumorPoints = getappdata(gcf,'tumorPoints');

% If more than 2 points go to 1/4 mark
if (size(tumorPoints,1)>1)
    
    % Get minimum and maximum slice
    minSlice = min(tumorPoints(:,3));
    maxSlice = max(tumorPoints(:,3));
    
    % Calculate 1/4 mark
    curVal = round(minSlice + 3*(maxSlice - minSlice)/4);
    
    % Set slider and redraw
    set(handles.sliceSlider,'Value',curVal);
    sliceSlider_Callback(hObject, eventdata, handles)
    
end


% --- Executes on button press in savePointsButton.
function savePointsButton_Callback(hObject, eventdata, handles)
% hObject    handle to savePointsButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

csvwrite('tumorPoints.dat',getappdata(gcf,'tumorPoints'));
csvwrite('vesselPoints.dat',getappdata(gcf,'vesselPoints'));
