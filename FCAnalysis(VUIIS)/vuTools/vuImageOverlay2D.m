function varargout = vuImageOverlay2D(varargin)
% VUIMAGEOVERLAY2D M-file for vuImageOverlay2D.fig
%      VUIMAGEOVERLAY2D, by itself, creates a new VUIMAGEOVERLAY2D or raises the existing
%      singleton*.
%
%      H = VUIMAGEOVERLAY2D returns the handle to a new VUIMAGEOVERLAY2D or the handle to
%      the existing singleton*.
%
%      VUIMAGEOVERLAY2D('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in VUIMAGEOVERLAY2D.M with the given input arguments.
%
%      VUIMAGEOVERLAY2D('Property','Value',...) creates a new VUIMAGEOVERLAY2D or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before vuImageOverlay2D_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to vuImageOverlay2D_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help vuImageOverlay2D

% Last Modified by GUIDE v2.5 12-Nov-2010 13:06:53

% Begin initialization code - DO NOT EDIT
gui_Singleton = 0;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @vuImageOverlay2D_OpeningFcn, ...
                   'gui_OutputFcn',  @vuImageOverlay2D_OutputFcn, ...
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


% --- Executes just before vuImageOverlay2D is made visible.
function vuImageOverlay2D_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to vuImageOverlay2D (see VARARGIN)

% Choose default command line output for vuImageOverlay2D
handles.output = hObject;

% Read input
im1 = varargin{1};
im2 = varargin{2};

% Scale images to 1-256
im1 = im1-min(im1(:));
im1 = im1.*(255/max(im1(:)))+1;
im2 = im2-min(im2(:));
im2 = im2.*(255/max(im2(:)))+1;

% Save images
setappdata(gcf,'im1',im1);
setappdata(gcf,'im2',im2);

% Create RGB Images
createRGBImage(1, hObject, eventdata, handles)
createRGBImage(0, hObject, eventdata, handles)

% View data
viewData(hObject, eventdata, handles)

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes vuImageOverlay2D wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = vuImageOverlay2D_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function topLowerEdit_Callback(hObject, eventdata, handles)
% hObject    handle to topLowerEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of topLowerEdit as text
%        str2double(get(hObject,'String')) returns contents of topLowerEdit as a double

% Get new value
winLowerValue = str2num(get(handles.topLowerEdit,'String'));
winUpperValue = str2num(get(handles.topUpperEdit,'String'));

% Is it valid
if ((~isempty(winLowerValue))&&(winLowerValue>=1)&&(winLowerValue<=winUpperValue))
    
    % Recalculate and view images
    createRGBImage(1, hObject, eventdata, handles)
    viewData(hObject, eventdata, handles)
else
    % Reset
    set(handles.topLowerEdit,'String','1')
    topUpperEdit_Callback(hObject, eventdata, handles)
end


% --- Executes during object creation, after setting all properties.
function topLowerEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to topLowerEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function topUpperEdit_Callback(hObject, eventdata, handles)
% hObject    handle to topUpperEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of topUpperEdit as text
%        str2double(get(hObject,'String')) returns contents of topUpperEdit as a double

% Get new value
winLowerValue = str2num(get(handles.topLowerEdit,'String'));
winUpperValue = str2num(get(handles.topUpperEdit,'String'));

% Is it valid
if ((~isempty(winUpperValue))&&(winUpperValue>=winLowerValue)&&(winUpperValue<=256))
    
    % Recalculate and view images
    createRGBImage(1, hObject, eventdata, handles)
    viewData(hObject, eventdata, handles)
else
    % Reset
    set(handles.topUpperEdit,'String','256')
    topUpperEdit_Callback(hObject, eventdata, handles)
end

% --- Executes during object creation, after setting all properties.
function topUpperEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to topUpperEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function bottomLowerEdit_Callback(hObject, eventdata, handles)
% hObject    handle to bottomLowerEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bottomLowerEdit as text
%        str2double(get(hObject,'String')) returns contents of bottomLowerEdit as a double

% Get new value
winLowerValue = str2num(get(handles.bottomLowerEdit,'String'));
winUpperValue = str2num(get(handles.bottomUpperEdit,'String'));

% Is it valid
if ((~isempty(winLowerValue))&&(winLowerValue>=1)&&(winLowerValue<=winUpperValue))
    
    % Recalculate and view images
    createRGBImage(0, hObject, eventdata, handles)
    viewData(hObject, eventdata, handles)
else
    % Reset
    set(handles.bottomLowerEdit,'String','1')
    bottomLowerEdit_Callback(hObject, eventdata, handles)
end

% --- Executes during object creation, after setting all properties.
function bottomLowerEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bottomLowerEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function bottomUpperEdit_Callback(hObject, eventdata, handles)
% hObject    handle to bottomUpperEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of bottomUpperEdit as text
%        str2double(get(hObject,'String')) returns contents of bottomUpperEdit as a double

% Get new value
winLowerValue = str2num(get(handles.bottomLowerEdit,'String'));
winUpperValue = str2num(get(handles.bottomUpperEdit,'String'));

% Is it valid
if ((~isempty(winUpperValue))&&(winUpperValue>=winLowerValue)&&(winUpperValue<=256))
    
    % Recalculate and view images
    createRGBImage(0, hObject, eventdata, handles)
    viewData(hObject, eventdata, handles)
else
    % Reset
    set(handles.bottomUpperEdit,'String','256')
    bottomUpperEdit_Callback(hObject, eventdata, handles)
end

% --- Executes during object creation, after setting all properties.
function bottomUpperEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bottomUpperEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in topMapMenu.
function topMapMenu_Callback(hObject, eventdata, handles)
% hObject    handle to topMapMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns topMapMenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from topMapMenu

createRGBImage(1, hObject, eventdata, handles)
viewData(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function topMapMenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to topMapMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in bottomMapMenu.
function bottomMapMenu_Callback(hObject, eventdata, handles)
% hObject    handle to bottomMapMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns bottomMapMenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from bottomMapMenu

createRGBImage(0, hObject, eventdata, handles)
viewData(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function bottomMapMenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bottomMapMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in overlayCheckbox.
function overlayCheckbox_Callback(hObject, eventdata, handles)
% hObject    handle to overlayCheckbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of overlayCheckbox

viewData(hObject, eventdata, handles)

% --- Executes on slider movement.
function alphaSlider_Callback(hObject, eventdata, handles)
% hObject    handle to alphaSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

viewData(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function alphaSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to alphaSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function createRGBImage(topFlag, hObject, eventdata, handles)

if (topFlag)
    
    % Get image
    im = getappdata(gcf,'im1');
    
    % Get window limits
    minWin = round(str2double(get(handles.topLowerEdit,'String')));
    maxWin = round(str2double(get(handles.topUpperEdit,'String')));
    
    % Get map info
    mapString = get(handles.topMapMenu,'String');
    mapValue = get(handles.topMapMenu,'Value');
    map = eval([cell2mat(mapString(mapValue)) '(' int2str(maxWin-minWin+1) ');']);
    
    % Make map 256 length
    fullMap = zeros(256,3);
    fullMap(minWin:maxWin,:) = map;
    fullMap(1:minWin,:) = repmat(fullMap(minWin,:),[minWin 1]);
    fullMap(maxWin:end,:) = repmat(fullMap(maxWin,:),[256-maxWin+1 1]);
    
    % Create RGB
    imRGB = ind2rgb(round(im),fullMap);
    
    % Save RGB
    setappdata(gcf,'im1RGB',imRGB);
    
else
    
    % Get image
    im = getappdata(gcf,'im2');
    
    % Get window limits
    minWin = round(str2double(get(handles.bottomLowerEdit,'String')));
    maxWin = round(str2double(get(handles.bottomUpperEdit,'String')));
    
    % Get map info
    mapString = get(handles.bottomMapMenu,'String');
    mapValue = get(handles.bottomMapMenu,'Value');
    map = eval([cell2mat(mapString(mapValue)) '(' int2str(maxWin-minWin+1) ');']);
    
    % Make map 256 length
    fullMap = zeros(256,3);
    fullMap(minWin:maxWin,:) = map;
    fullMap(1:minWin,:) = repmat(fullMap(minWin,:),[minWin 1]);
    fullMap(maxWin:end,:) = repmat(fullMap(maxWin,:),[256-maxWin+1 1]);
    
    % Create RGB
    imRGB = ind2rgb(round(im),fullMap);
    
    % Save RGB
    setappdata(gcf,'im2RGB',imRGB);
end

function viewData(hObject, eventdata, handles)

% Get RGB images
im1RGB = getappdata(gcf,'im1RGB');
im2RGB = getappdata(gcf,'im2RGB');

im2 = getappdata(gcf,'im2');

% Show images
axes(handles.topFig)
imshow(im1RGB);
axes(handles.bottomFig)
imshow(im2RGB);
axes(handles.overlayFig)
imshow(im1RGB);
if (get(handles.overlayCheckbox,'Value'))
    maxAlpha = get(handles.alphaSlider,'Value');
    maxWin = round(str2double(get(handles.bottomUpperEdit,'String')));
    alphaData = double((im2-1)./maxWin)*maxAlpha;
    alphaData(alphaData>1) = 1;
    hold on
    image(im2RGB,'AlphaData',alphaData)
    hold off
end
