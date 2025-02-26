function varargout = vuOnePaneROIDrawer(varargin)
% VUONEPANEROIDRAWER M-file for vuOnePaneROIDrawer.fig
%      VUONEPANEROIDRAWER, by itself, creates a new VUONEPANEROIDRAWER or raises the existing
%      singleton*.
%
%      H = VUONEPANEROIDRAWER returns the handle to a new VUONEPANEROIDRAWER or the handle to
%      the existing singleton*.
%
%      VUONEPANEROIDRAWER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in VUONEPANEROIDRAWER.M with the given input arguments.
%
%      VUONEPANEROIDRAWER('Property','Value',...) creates a new VUONEPANEROIDRAWER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before vuOnePaneROIDrawer_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to vuOnePaneROIDrawer_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help vuOnePaneROIDrawer

% Last Modified by GUIDE v2.5 19-Mar-2012 13:31:24

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @vuOnePaneROIDrawer_OpeningFcn, ...
                   'gui_OutputFcn',  @vuOnePaneROIDrawer_OutputFcn, ...
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


% --- Executes just before vuOnePaneROIDrawer is made visible.
function vuOnePaneROIDrawer_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to vuOnePaneROIDrawer (see VARARGIN)

% Choose default command line output for vuOnePaneROIDrawer
handles.output = hObject;

% Check/Make input a Meta Image
if length(varargin) < 1
	error('MATLAB:vuOnePaneROIDrawer:NotEnoughInputs', 'Not enough input arguments.');
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
p.addParamValue('rgb',false,@(x) isa(x,'double'));
p.addParamValue('createRGB',false,@(x) isa(x,'double'));
p.addParamValue('pane',3,@(x) isa(x,'double'));
p.addParamValue('slice',midSlice(3),@(x) isa(x,'double'));
p.addParamValue('dynamic',1,@(x) isa(x,'double'));
p.addParamValue('roiNumber',0,@(x) isa(x,'double'));
p.FunctionName='vuOnePaneViewer';
p.parse(varargin{2:end});

pane = p.Results.pane;
slice = p.Results.slice;
dynamic = p.Results.dynamic;
rgbImage = p.Results.rgb;
createRGB = p.Results.createRGB;
roiNumber = p.Results.roiNumber;

% Check Pane
if pane < 1 || pane > 3
    pane = 3;
end

% Check Dynamic
if length(image.Dims) >= 4 && (dynamic < 1 || dynamic > image.Dims(4))
    dynamic = 1;
end

% Permute our image based on which pane is viewed
if pane == 1
    image.Data = permute(image.Data,[1 3 2 4]);
    image.Dims(1:3) = [image.Dims(2) image.Dims(3) image.Dims(1)];
    image.Spc(1:3) = [image.Spc(2) image.Spc(3) image.Spc(1)];
    image.Origin(1:3) = [image.Origin(2) image.Origin(3) image.Origin(1)];
elseif pane == 2
    image.Data = permute(image.Data,[2 3 1 4]);
    image.Dims(1:3) = [image.Dims(1) image.Dims(3) image.Dims(2)];
    image.Spc(1:3) = [image.Spc(1) image.Spc(3) image.Spc(2)];
    image.Origin(1:3) = [image.Origin(1) image.Origin(3) image.Origin(2)];
end


% 3 color image?
if rgbImage
    
    % Already 3 color?
    if (~createRGB)
                
        % Save RGB image
        setappdata(gcf,'image',image.Data);
        setappdata(gcf,'rgbFlag',1);
        
        roi = image;
        roi.Dims(3) = [];
        roi.Data = zeros(roi.Dims);
    else
        % Create 3 color RGB image
        rgbMap = [0 0 0;0 0 1;0 1 0;1 0 0];
        imageRGB = zeros([image.Dims(1:2) 3 image.Dims(3)]);
        for i = 1:image.Dims(3)
            imageRGB(:,:,:,i) = ind2rgb(image.Data(:,:,i), rgbMap);
        end

        % Save RGB image
        setappdata(gcf,'image',imageRGB);
        setappdata(gcf,'rgbFlag',1);
        
        roi = image;
        roi.Data(:) = 0;
        
        image.Dims(4) = image.Dims(3);
        image.Dims(3) = 3;
    end
    
    % Check Slice
    if slice < 1 || slice > image.Dims(4)
        slice = midSlice(3);
    end
    
    % Setup the sliders
    set(handles.sliceSlider,'Min',0.999)
    set(handles.sliceSlider,'Max',image.Dims(4))
    set(handles.sliceSlider,'SliderStep',[1./(image.Dims(4)-0.999) 1./(image.Dims(4)-0.999)]);
    set(handles.sliceSlider,'Value',slice);
    
    set(handles.lowerSlide,'Visible','off');
    set(handles.upperSlide,'Visible','off');
    
else
    
    % Save image
    setappdata(gcf,'image',image);
    setappdata(gcf,'rgbFlag',0);

    % Set slider to 2%-98%
    sortedImage = sort(image.Data(:));
    lowerValue = sortedImage(round(0.02*length(sortedImage)));
    upperValue = sortedImage(round(0.98*length(sortedImage)));
    setappdata(gcf,'curMin',lowerValue);
    setappdata(gcf,'curMax',upperValue);
    if (sortedImage(1) == sortedImage(end))
        sortedImage(end) = sortedImage(1) + 0.0002;
    end
    set(handles.lowerSlide,'Min', sortedImage(1));
    set(handles.lowerSlide,'Max', sortedImage(end)-0.0001);
    set(handles.upperSlide,'Min', sortedImage(1)+0.0001);
    set(handles.upperSlide,'Max', sortedImage(end));
    set(handles.lowerSlide,'Value', lowerValue);
    set(handles.upperSlide,'Value', upperValue);
    
    roi = image;
    roi.Data(:) = 0;
    
    % Check Slice
    if slice < 1 || slice > image.Dims(3)
        slice = midSlice(3);
    end    
    
    % Setup the sliders
    set(handles.sliceSlider,'Min',0.999)
    set(handles.sliceSlider,'Max',image.Dims(3))
    set(handles.sliceSlider,'SliderStep',[1./(image.Dims(3)-0.999) 1./(image.Dims(3)-0.999)]);
    set(handles.sliceSlider,'Value',slice);
    
end

% Save options
setappdata(gcf,'dynamic',dynamic);
setappdata(gcf,'roiNumber',roiNumber)

setappdata(gcf,'ROI',roi);
setappdata(gcf,'roiPoints',[]);

% Initialize axes limits
axes(handles.mainFig);
set(handles.mainFig,'xlim',[0.5 image.Dims(1)+0.5],'ylim',[0.5 image.Dims(2)+0.5]);

viewData(hObject, eventdata, handles)

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes vuOnePaneROIDrawer wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = vuOnePaneROIDrawer_OutputFcn(hObject, eventdata, handles) 
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
roiPoints = getappdata(gcf,'roiPoints');

% Check if roi present
if (length(roiPoints) >= curVal && ~isempty(roiPoints{curVal}))

    set(handles.roiButton,'String','Clear ROI');
    set(handles.copyPasteROIButton,'String','Copy ROI');
    
else
    
    set(handles.roiButton,'String','New ROI Slice');
    set(handles.copyPasteROIButton,'String','Paste ROI');

    
end

viewData(hObject, eventdata, handles);


% --- Executes during object creation, after setting all properties.
function sliceSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliceSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in roiButton.
function roiButton_Callback(hObject, eventdata, handles)
% hObject    handle to roiButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

curVal = round(get(handles.sliceSlider,'Value'));
roiPoints = getappdata(gcf,'roiPoints');

set(handles.noneButton,'Value',1);
zoom off
pan off

% Check if roi present
if (length(roiPoints) >= curVal && ~isempty(roiPoints{curVal}))

    set(handles.roiButton,'String','New ROI Slice');
    set(handles.copyPasteROIButton,'String','Paste ROI');
    
    % Save ROI
    roi = getappdata(gcf,'ROI');
    roi.Data(:,:,curVal) = zeros(size(roi.Data(:,:,curVal)));
    setappdata(gcf,'ROI',roi);

    % Save ROI points
    roiPoints = getappdata(gcf,'roiPoints');
    roiPoints{curVal} = [];
    setappdata(gcf,'roiPoints',roiPoints);
    
    viewData(hObject, eventdata, handles);

else

    % Disable buttons
    set(handles.sliceSlider,'Enable','off');
    set(handles.roiButton,'Enable','off');
    set(handles.finishButton,'Enable','off');
    set(handles.uipanel2,'visible','off');
    set(handles.copyPasteROIButton,'Enable','off');
    set(handles.editROIButton,'Enable','off');
    
    % Draw ROI
    axes(handles.mainFig);
    [BW, xi, yi] = roipoly;
    
    % Enable buttons
    set(handles.sliceSlider,'Enable','on');
    set(handles.roiButton,'Enable','on');
    set(handles.finishButton,'Enable','on');
    set(handles.uipanel2,'visible','on');
    set(handles.copyPasteROIButton,'Enable','on');
    set(handles.editROIButton,'Enable','on');
    set(handles.roiButton,'String','Clear ROI'); 
    set(handles.copyPasteROIButton,'String','Copy ROI');

    % Save ROI
    roi = getappdata(gcf,'ROI');
    roi.Data(:,:,curVal) = single(BW);
    setappdata(gcf,'ROI',roi);

    % Save ROI points
    roiPoints = getappdata(gcf,'roiPoints');
    roiPoints{curVal} = [xi';yi'];
    setappdata(gcf,'roiPoints',roiPoints);

    line(roiPoints{curVal}(1,:),roiPoints{curVal}(2,:),'Marker','s','MarkerSize',4,'MarkerEdgeColor','m','LineWidth',1.5,'Color','white')
    
end

% --- Executes on button press in finishButton.
function finishButton_Callback(hObject, eventdata, handles)
% hObject    handle to finishButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ROI = getappdata(gcf,'ROI');
roiNumber = getappdata(gcf,'roiNumber');
assignin('base',sprintf('vuOnePaneROIDrawer_ROI_%i',roiNumber),ROI);
close gcf

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

% View data function
function viewData(hObject, eventdata, handles)

image = getappdata(gcf,'image');
rgbFlag = getappdata(gcf,'rgbFlag');
curVal = round(get(handles.sliceSlider,'Value'));
dynamic = getappdata(gcf,'dynamic');
roiPoints = getappdata(gcf,'roiPoints');
xlimits = xlim(handles.mainFig);
ylimits = ylim(handles.mainFig);

% 3 color image
if (rgbFlag)
    
    % Show image
    axes(handles.mainFig)
    imshow(squeeze(image(:,:,:,curVal,dynamic)))
    set(handles.mainFig,'xlim',xlimits,'ylim',ylimits);

% Grayscale image
else
    
    % Window levels
    curMin = getappdata(gcf,'curMin');
    curMax = getappdata(gcf,'curMax');
    

    % Show images
    axes(handles.mainFig)
    imshow(image.Data(:,:,curVal,dynamic),[curMin curMax])
    set(handles.mainFig,'xlim',xlimits,'ylim',ylimits);

end

% ROI 
if (length(roiPoints) >= curVal && ~isempty(roiPoints{curVal}))

    line(roiPoints{curVal}(1,:),roiPoints{curVal}(2,:),'Marker','s','MarkerSize',4,'MarkerEdgeColor','m','LineWidth',1.5,'Color','white')
    
end


% --- Executes when selected object is changed in uipanel2.
function uipanel2_SelectionChangeFcn(hObject, eventdata, handles)
% hObject    handle to the selected object in uipanel2 
% eventdata  structure with the following fields (see UIBUTTONGROUP)
%	EventName: string 'SelectionChanged' (read only)
%	OldValue: handle of the previously selected object or empty if none was selected
%	NewValue: handle of the currently selected object
% handles    structure with handles and user data (see GUIDATA)

noneValue = get(handles.noneButton,'Value');
zoomValue = get(handles.zoomButton,'Value');
panValue = get(handles.panButton,'Value');

if (noneValue)
    zoom off
    pan off
elseif (zoomValue)
    pan off
    zoom on
elseif (panValue)
    zoom off
    pan on
end


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
usrMax = round(usrMax*1000)/1000;
if (usrMax <= curMin)
    usrMax = curMin+(get(handles.upperSlide,'Max')-get(handles.lowerSlide,'Min'))*0.001;
    if (usrMax > get(handles.upperSlide,'Max'));
        usrMax = get(handles.upperSlide,'Max');
    end
    set(handles.upperSlide,'Value',usrMax);
end
setappdata(gcf,'curMax',usrMax);
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
usrMin = round(usrMin*1000)/1000;
if (usrMin >= curMax)
    usrMin = curMax-(get(handles.upperSlide,'Max')-get(handles.lowerSlide,'Min'))*0.001;
    if (usrMin < get(handles.lowerSlide,'Min'))
        usrMin = get(handles.lowerSlide,'Min');
    end
    set(handles.lowerSlide,'Value',usrMin);
end
setappdata(gcf,'curMin',usrMin);
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


% --- Executes on button press in editROIButton.
function editROIButton_Callback(hObject, eventdata, handles)
% hObject    handle to editROIButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

curVal = round(get(handles.sliceSlider,'Value'));
roiPoints = getappdata(gcf,'roiPoints');

set(handles.noneButton,'Value',1);
zoom off
pan off

% Check if roi present
if (length(roiPoints) >= curVal && ~isempty(roiPoints{curVal}))

    % Disable buttons
    set(handles.sliceSlider,'Enable','off');
    set(handles.roiButton,'Enable','off');
    set(handles.finishButton,'Enable','off');
    set(handles.uipanel2,'visible','off');
    set(handles.copyPasteROIButton,'Enable','off');
    set(handles.editROIButton,'Enable','off');
    
    % Create editable poly
    points = roiPoints{curVal}';
    h = impoly(gca, points(1:end-1,:));
    newROIPoints = wait(h);
    newROIPoints(end+1,:) = newROIPoints(1,:);
    
    % Save ROI
    roi = getappdata(gcf,'ROI');
    roi.Data(:,:,curVal) = single(roipoly(getimage(gca),newROIPoints(:,1)',newROIPoints(:,2)'));
    setappdata(gcf,'ROI',roi);
    
    % Save new ROI points
    roiPoints = getappdata(gcf,'roiPoints');
    roiPoints{curVal} = newROIPoints';
    setappdata(gcf,'roiPoints',roiPoints);
    
    % Enable buttons
    set(handles.sliceSlider,'Enable','on');
    set(handles.roiButton,'Enable','on');
    set(handles.finishButton,'Enable','on');
    set(handles.uipanel2,'visible','on');
    set(handles.copyPasteROIButton,'Enable','on');
    set(handles.editROIButton,'Enable','on');
    set(handles.roiButton,'String','Clear ROI'); 
    set(handles.copyPasteROIButton,'String','Copy ROI');
    
    viewData(hObject, eventdata, handles);
end


% --- Executes on button press in copyPasteROIButton.
function copyPasteROIButton_Callback(hObject, eventdata, handles)
% hObject    handle to copyPasteROIButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


curVal = round(get(handles.sliceSlider,'Value'));
roiPoints = getappdata(gcf,'roiPoints');

set(handles.noneButton,'Value',1);
zoom off
pan off

% Check if roi present (COPY)
if (length(roiPoints) >= curVal && ~isempty(roiPoints{curVal}))
    
    % Save ROI on "clipboard"
    setappdata(gcf,'roiClipboard',roiPoints{curVal});
    
% (PASTE)
else 
    
    % Disable buttons
    set(handles.sliceSlider,'Enable','off');
    set(handles.roiButton,'Enable','off');
    set(handles.finishButton,'Enable','off');
    set(handles.uipanel2,'visible','off');
    set(handles.copyPasteROIButton,'Enable','off');
    set(handles.editROIButton,'Enable','off');
    
    % Get the clipboard
    roiClipboard = getappdata(gcf,'roiClipboard');
    
    % Create editable poly
    h = impoly(gca, roiClipboard(:,1:end-1)');
    newROIPoints = wait(h);
    newROIPoints(end+1,:) = newROIPoints(1,:);
    
    % Save ROI
    roi = getappdata(gcf,'ROI');
    roi.Data(:,:,curVal) = single(roipoly(getimage(gca),newROIPoints(:,1)',newROIPoints(:,2)'));
    setappdata(gcf,'ROI',roi);
    
    % Save new ROI points
    roiPoints = getappdata(gcf,'roiPoints');
    roiPoints{curVal} = newROIPoints';
    setappdata(gcf,'roiPoints',roiPoints);
    
    % Enable buttons
    set(handles.sliceSlider,'Enable','on');
    set(handles.roiButton,'Enable','on');
    set(handles.finishButton,'Enable','on');
    set(handles.uipanel2,'visible','on');
    set(handles.copyPasteROIButton,'Enable','on');
    set(handles.editROIButton,'Enable','on');
    set(handles.roiButton,'String','Clear ROI'); 
    set(handles.copyPasteROIButton,'String','Copy ROI'); 
    
    viewData(hObject, eventdata, handles);
    
end
