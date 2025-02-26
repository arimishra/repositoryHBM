function varargout = vuOrthogonalOverlay(varargin)
% VUORTHOGONALOVERLAY allows for two 3D images and their corresponding
% orthogonal views to be seen simultaneously with an option to overlay
% each.
%      Usage: vuorthogonaloverlay(IM1,IM2)
%               IM1 = Bottom Layer Image Array 
%               IM2 = Top Layer Image Array (overlaid image)
%               Note:  This program requires a large amount of RAM, and
%               is designed to accept Meta Image structures, if input is
%               not a structure, vuGenerateMetaImage will convert it,
%               however, the dimensions may be incorrect.
%               To prevent running out of RAM, close all other programs,
%               and clear variables for successive runs.
%               Recommended for use in Core Computer Lab, but will run on
%               any computer with at least 1 GB of RAM.
%
%Author: S.Narayan - 10/5/2007.
%Version 1.1(10/10/2007): Increased to six panels with independent colormap and window
%intensity sliders,slice synchronization, overlay on/off,more robust code
%Version 1.2(10/15/2007): Streamlined code for faster performance, re-nested slice
%callbacks to eliminate redundant code.
%Version 1.3(10/20/2007: Added uicontrol object Update Image to prevent excessive
%updateImage calls. Introduced Recursion Limit set to 200, to prevent,
%running out of RAM resulting in freezing or Matlab crash, also clears
%defined variables after initialization and updateImage function, appdata
%remains defined.
%
%
%      H = VUORTHOGONALOVERLAY returns the handle to a new VUORTHOGONALOVERLAY or the handle to
%      the existing singleton*.
%
%      VUORTHOGONALOVERLAY('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in VUORTHOGONALOVERLAY.M with the given input arguments.
%
%      VUORTHOGONALOVERLAY('Property','Value',...) creates a new VUORTHOGONALOVERLAY or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before vuOrthogonalOverlay_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to vuOrthogonalOverlay_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help vuOrthogonalOverlay

% Last Modified by GUIDE v2.5 19-Oct-2007 16:08:48

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @vuOrthogonalOverlay_OpeningFcn, ...
                   'gui_OutputFcn',  @vuOrthogonalOverlay_OutputFcn, ...
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


% --- Executes just before vuOrthogonalOverlay is made visible.
function vuOrthogonalOverlay_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to vuOrthogonalOverlay (see VARARGIN)

% Choose default command line output for vuOrthogonalOverlay
handles.output = hObject;

% Get default command line output from handles structure
varargout{1} = handles.output;

% Assign Data to cell array of structures
imageData=cell(1,2);
% Check/Make input a Meta Image
if (~isstruct(varargin{1}))  | (~isstruct(varargin{2}))
    imageData{1,1}= vuGenerateMetaImage(varargin{1});
    imageData{1,2}= vuGenerateMetaImage(varargin{2});
else
    imageData{1,1}= varargin{1};
    imageData{1,2}= varargin{2};
end
setappdata(handles.axes1,'image',imageData{1,1})
setappdata(handles.axes7,'image',imageData{1,2});

midSlice = floor(imageData{1,1}.Dims/2);
midSlice2 = floor(imageData{1,2}.Dims/2);

%Set limits for slice sliders
set(handles.sliceX,'Min',0.999)
set(handles.sliceX,'Max',imageData{1,1}.Dims(2))
set(handles.sliceX,'SliderStep',[1./(imageData{1,1}.Dims(2)-0.999) 1./(imageData{1,1}.Dims(2)-0.999)]);
set(handles.sliceX,'Value',midSlice(2));
set(handles.sliceX2,'Min',0.999)
set(handles.sliceX2,'Max',imageData{1,2}.Dims(2))
set(handles.sliceX2,'SliderStep',[1./(imageData{1,2}.Dims(2)-0.999) 1./(imageData{1,2}.Dims(2)-0.999)]);
set(handles.sliceX2,'Value',midSlice2(2));

set(handles.sliceY,'Min',0.999)
set(handles.sliceY,'Max',imageData{1,1}.Dims(1))
set(handles.sliceY,'SliderStep',[1./(imageData{1,1}.Dims(1)-0.999) 1./(imageData{1,1}.Dims(1)-0.999)]);
set(handles.sliceY,'Value',midSlice(1));
set(handles.sliceY2,'Min',0.999)
set(handles.sliceY2,'Max',imageData{1,2}.Dims(1))
set(handles.sliceY2,'SliderStep',[1./(imageData{1,2}.Dims(1)-0.999) 1./(imageData{1,2}.Dims(1)-0.999)]);
set(handles.sliceY2,'Value',midSlice(1));

set(handles.sliceZ,'Min',0.999)
set(handles.sliceZ,'Max',imageData{1,1}.Dims(3))
set(handles.sliceZ,'SliderStep',[1./(imageData{1,1}.Dims(3)-0.999) 1./(imageData{1,1}.Dims(3)-0.999)]);
set(handles.sliceZ,'Value',midSlice(3));
set(handles.sliceZ2,'Min',0.999)
set(handles.sliceZ2,'Max',imageData{1,2}.Dims(3))
set(handles.sliceZ2,'SliderStep',[1./(imageData{1,2}.Dims(3)-0.999) 1./(imageData{1,2}.Dims(3)-0.999)]);
set(handles.sliceZ2,'Value',midSlice(3));

%Set limits for window sliders
maxValue = max(max(max(imageData{1,1}.Data(:))));
minValue = min(min(min(imageData{1,1}.Data(:))));

maxValue2 = max(max(max(imageData{1,2}.Data(:))));
minValue2 = min(min(min(imageData{1,2}.Data(:))));

if ((minValue-maxValue)==0)
    minValue = maxValue-0.001;
elseif ((minValue2-maxValue2)==0)
    minValue2 = maxValue2-0.001;
end

setappdata(handles.minValue,'defMin',minValue);
setappdata(handles.maxValue,'defMax',maxValue);
setappdata(handles.minValue,'curMin',minValue);
setappdata(handles.maxValue,'curMax',maxValue);
setappdata(handles.minValue2,'defMin',minValue2);
setappdata(handles.maxValue2,'defMax',maxValue2);
setappdata(handles.minValue2,'curMin',minValue2);
setappdata(handles.maxValue2,'curMax',maxValue2);

set(handles.minValue,'Min', minValue);
set(handles.minValue,'Max', maxValue-0.0001);
set(handles.maxValue,'Min', minValue+0.0001);
set(handles.maxValue,'Max', maxValue);
set(handles.minValue2,'Min', minValue2);
set(handles.minValue2,'Max', maxValue2-0.0001);
set(handles.maxValue2,'Min', minValue2+0.0001);
set(handles.maxValue2,'Max', maxValue2);

set(handles.minValue,'Value', minValue);
set(handles.maxValue,'Value', maxValue);
set(handles.minValue2,'Value', minValue2);
set(handles.maxValue2,'Value', maxValue2);

set(handles.minEdit,'String',num2str(minValue));
set(handles.maxEdit,'String',num2str(maxValue));
set(handles.minEdit2,'String',num2str(minValue2));
set(handles.maxEdit2,'String',num2str(maxValue2));

set(handles.maxValue,'Enable','off')
set(handles.maxEdit,'Enable','off')
set(handles.minValue,'Enable','off')
set(handles.minEdit,'Enable','off')
set(handles.maxValue2,'Enable','off')
set(handles.maxEdit2,'Enable','off')
set(handles.minValue2,'Enable','off')
set(handles.minEdit2,'Enable','off')


%Set initial values/operating status of opacity sliders
set(handles.opaque1,'Value',0.5)
set(handles.opaque2,'Value',0.5)
set(handles.opaque3,'Value',0.5)
set(handles.opaque1,'Enable','off')
set(handles.opaque2,'Enable','off')
set(handles.opaque3,'Enable','off')

% Calculate and apply current colormap & scaling
ind_im.Index=uint16((imageData{1,1}.Data-min(imageData{1,1}.Data(:)))./(max(imageData{1,1}.Data(:))-min(imageData{1,1}.Data(:))).*65535+1);
ind_im2.Index=uint16((imageData{1,2}.Data-min(imageData{1,2}.Data(:)))./(max(imageData{1,2}.Data(:))-min(imageData{1,2}.Data(:))).*65535+1);
ind_im.Map=gray(65536);
ind_im2.Map=gray(65536);
ind_im.MapName='gray';
ind_im2.MapName='gray';
imageData{2,1}=ind_im;
imageData{2,2}=ind_im2;
setappdata(handles.axes1,'index',imageData{2,1})
setappdata(handles.axes7,'index',imageData{2,2})


% Display Images
axes(handles.axes1)
imagesc(ind2rgb(squeeze(ind_im.Index(midSlice(2),:,:)),ind_im.Map));
axes(handles.axes3)
imagesc(ind2rgb(squeeze(ind_im.Index(:,midSlice(1),:)),ind_im.Map));
axes(handles.axes5)
imagesc(ind2rgb(squeeze(ind_im.Index(:,:,midSlice(3))),ind_im.Map));

axes(handles.axes7)
imagesc(ind2rgb(squeeze(ind_im2.Index(midSlice2(2),:,:)),ind_im2.Map));
axes(handles.axes8)
imagesc(ind2rgb(squeeze(ind_im2.Index(:,midSlice2(1),:)),ind_im2.Map));
axes(handles.axes9)
imagesc(ind2rgb(squeeze(ind_im2.Index(:,:,midSlice2(3))),ind_im2.Map));


% Image Number Labels
set(handles.imageNum,'String',num2str(midSlice(2)));
set(handles.imageNum2,'String',num2str(midSlice(1)));
set(handles.imageNum3,'String',num2str(midSlice(3)));
set(handles.imageNum7,'String',num2str(midSlice2(2)));
set(handles.imageNum8,'String',num2str(midSlice2(1)));
set(handles.imageNum9,'String',num2str(midSlice2(3)));

% --- Outputs from this function are returned to the command line.
function varargout = vuOrthogonalOverlay_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes vuOrthogonalOverlay wait for user response (see UIRESUME)
% uiwait(handles.figure1);

clear imageData
clear midSlice
clear midSlice2
clear minValue
clear maxValue
clear minValue2
clear maxValue
clear ind_im
clear ind_im2

function updateImage(handles)
          
mapString = get(handles.cMapmenu,'String');
mapValue = get(handles.cMapmenu,'Value');
map=cell2mat(mapString(mapValue));
% map = cell2mat(map)
eval(['m1=' map '(65536);']);
setappdata(handles.cMapmenu,'map',m1)
mapString = get(handles.cMapmenu2,'String');
mapValue = get(handles.cMapmenu2,'Value');
map=cell2mat(mapString(mapValue));
% map = cell2mat(map)
eval(['m2=' map '(65536);']);
setappdata(handles.cMapmenu2,'map',m2)

imageData{1,1}=getappdata(handles.axes1,'image');
imageData{1,2}=getappdata(handles.axes7,'image');
maxValue=getappdata(handles.maxValue,'curMax');
minValue=getappdata(handles.minValue,'curMin');
maxValue2=getappdata(handles.maxValue2,'curMax');
minValue2=getappdata(handles.minValue2,'curMin');

ind_im.Index=uint16((imageData{1,1}.Data-minValue)./(maxValue-minValue).*65535+1);
ind_im2.Index=uint16((imageData{1,2}.Data-minValue2)./(maxValue2-minValue2).*65535+1);
ind_im.Map=getappdata(handles.cMapmenu,'map');
ind_im2.Map=getappdata(handles.cMapmenu2,'map');
% if isempty(ind_im.Map)
%     ind_im.Map=gray(65536)

imageData{2,1}=ind_im;
imageData{2,2}=ind_im2;
setappdata(handles.axes1,'index',imageData{2,1})
setappdata(handles.axes7,'index',imageData{2,2})

overlay=get(handles.overlayOn,'Value');
sync=get(handles.synchronize,'Value');
islice{1,1}=round(get(handles.sliceX,'Value'));
islice{1,2}=round(get(handles.sliceY,'Value'));
islice{1,3}=round(get(handles.sliceZ,'Value'));
islice{2,1}=round(get(handles.sliceX2,'Value'));
islice{2,2}=round(get(handles.sliceY2,'Value'));
islice{2,3}=round(get(handles.sliceZ2,'Value'));
alpha_val=cell(1,3);
alpha_val{1,1}=get(handles.opaque1,'Value');
alpha_val{1,2}=get(handles.opaque2,'Value');
alpha_val{1,3}=get(handles.opaque3,'Value');

if sync==1;
        axes(handles.axes7)
        imagesc(ind2rgb(squeeze(imageData{2,2}.Index(islice{1,1},:,:)),imageData{2,2}.Map));
        axes(handles.axes8)
        imagesc(ind2rgb(squeeze(imageData{2,2}.Index(:,islice{1,2},:)),imageData{2,2}.Map));
        axes(handles.axes9)
        imagesc(ind2rgb(squeeze(imageData{2,2}.Index(:,:,islice{1,3})),imageData{2,2}.Map));        
        axes(handles.axes1)
        imagesc(ind2rgb(squeeze(imageData{2,1}.Index(islice{1,1},:,:)),imageData{2,1}.Map));
    if overlay==1
        hold on;
        alpha(imagesc(ind2rgb(squeeze(imageData{2,2}.Index(islice{1,1},:,:)),imageData{2,2}.Map)),alpha_val{1,1})
    end
        axes(handles.axes3)
        imagesc(ind2rgb(squeeze(imageData{2,1}.Index(:,islice{1,2},:)),imageData{2,1}.Map));
    if overlay
        hold on;
        alpha(imagesc(ind2rgb(squeeze(imageData{2,2}.Index(:,islice{1,2},:)),imageData{2,2}.Map)),alpha_val{1,2})
    end        
        axes(handles.axes5)
        imagesc(ind2rgb(squeeze(imageData{2,1}.Index(:,:,islice{1,3})),imageData{2,1}.Map));
    if overlay
        hold on;
        alpha(imagesc(ind2rgb(squeeze(imageData{2,2}.Index(:,:,islice{1,3})),imageData{2,2}.Map)),alpha_val{1,3})
    end
    set(handles.imageNum7,'String',islice{1,1})
    set(handles.sliceX2,'Value',islice{1,1})
    set(handles.imageNum8,'String',islice{1,2})
    set(handles.sliceY2,'Value',islice{1,2})
    set(handles.imageNum9,'String',islice{1,3})
    set(handles.sliceZ2,'Value',islice{1,3})
elseif sync==0
        axes(handles.axes1)
        imagesc(ind2rgb(squeeze(imageData{2,1}.Index(islice{1,1},:,:)),imageData{2,1}.Map));
    if overlay==1
        hold on;
        alpha(imagesc(ind2rgb(squeeze(imageData{2,2}.Index(islice{1,1},:,:)),imageData{2,2}.Map)),alpha_val{1,1})
    end
        axes(handles.axes3)
        imagesc(ind2rgb(squeeze(imageData{2,1}.Index(:,islice{1,2},:)),imageData{2,1}.Map));
    if overlay
        hold on;
        alpha(imagesc(ind2rgb(squeeze(imageData{2,2}.Index(:,islice{1,2},:)),imageData{2,2}.Map)),alpha_val{1,2})
    end        
        axes(handles.axes5)
        imagesc(ind2rgb(squeeze(imageData{2,1}.Index(:,:,islice{1,3})),imageData{2,1}.Map));
    if overlay
        hold on;
        alpha(imagesc(ind2rgb(squeeze(imageData{2,2}.Index(:,:,islice{1,3})),imageData{2,2}.Map)),alpha_val{1,3})
    end
        axes(handles.axes7)
        imagesc(ind2rgb(squeeze(imageData{2,2}.Index(islice{2,1},:,:)),imageData{2,2}.Map));
        axes(handles.axes8)
        imagesc(ind2rgb(squeeze(imageData{2,2}.Index(:,islice{2,2},:)),imageData{2,2}.Map));
        axes(handles.axes9)
        imagesc(ind2rgb(squeeze(imageData{2,2}.Index(:,:,islice{2,3})),imageData{2,2}.Map));
end
clear imageData
clear islice
clear alpha_val
clear ind_im
clear ind_im2

% --- Executes on button press in synchronize.
function synchronize_Callback(hObject, eventdata, handles)
% hObject    handle to synchronize (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of synchronize
sync=get(hObject,'Value');
imageData{1,1}=getappdata(handles.axes1,'index');
imageData{1,2}=getappdata(handles.axes7,'index');
islice{1,1}=round(get(handles.sliceX,'Value'));
islice{1,2}=round(get(handles.sliceY,'Value'));
islice{1,3}=round(get(handles.sliceZ,'Value'));
if sync==1;
    %sync image 2 sliders with image 1 sliders
    cur_val=get(handles.imageNum,'String');
    set(handles.sliceX2,'Value',str2num(cur_val))
    set(handles.imageNum7,'String',cur_val)
    cur_val2=get(handles.imageNum2,'String');
    set(handles.sliceY2,'Value',str2num(cur_val2))
    set(handles.imageNum8,'String',cur_val2)
    cur_val3=get(handles.imageNum3,'String');
    set(handles.sliceZ2,'Value',str2num(cur_val3))
    set(handles.imageNum9,'String',cur_val3)
    %update images
    axes(handles.axes7)
    imagesc(ind2rgb(squeeze(imageData{1,2}.Index(islice{1,1},:,:)),imageData{1,2}.Map));
    axes(handles.axes8)
    imagesc(ind2rgb(squeeze(imageData{1,2}.Index(:,islice{1,2},:)),imageData{1,2}.Map));
    axes(handles.axes9)
    imagesc(ind2rgb(squeeze(imageData{1,2}.Index(:,:,islice{1,3})),imageData{1,2}.Map));        
    %lock sliders
    set(handles.sliceX2,'Enable','off')
    set(handles.sliceY2,'Enable','off')
    set(handles.sliceZ2,'Enable','off')
    set(handles.imageNum7,'Enable','off')
    set(handles.imageNum8,'Enable','off')
    set(handles.imageNum9,'Enable','off')
else
    set(handles.sliceX2,'Enable','on')
    set(handles.imageNum7,'Enable','on')
    set(handles.sliceY2,'Enable','on')
    set(handles.imageNum8,'Enable','on')
    set(handles.sliceZ2,'Enable','on')
    set(handles.imageNum9,'Enable','on')
end

clear imageData
clear islice
% --- Executes on button press in autoWindow.
function autoWindow_Callback(hObject, eventdata, handles)
% hObject    handle to autoWindow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of autoWindow
auto_on = get(hObject,'Value');

if auto_on==0.0
    set(handles.maxValue, 'Enable', 'on');
    set(handles.minValue, 'Enable', 'on');
    set(handles.maxEdit, 'Enable', 'on');
    set(handles.minEdit, 'Enable', 'on');
    maxValue=get(handles.maxValue,'Value');
    minValue=get(handles.minValue,'Value');
    setappdata(handles.maxValue,'curMax',maxValue)
    setappdata(handles.minValue,'curMin',minValue)
else    
    minValue = getappdata(handles.maxValue,'defMin');
    maxValue = getappdata(handles.minValue,'defMax');
    setappdata(handles.maxValue,'curMax',maxValue)
    setappdata(handles.minValue,'curMin',minValue)
    set(handles.minValue,'Min', minValue);
    set(handles.maxValue,'Max', maxValue);
    set(handles.minValue,'Min', minValue);
    set(handles.maxValue,'Max', maxValue);
    set(handles.minValue,'Value', minValue);
    set(handles.maxValue,'Value', maxValue);
    set(handles.minEdit,'String',num2str(minValue));
    set(handles.maxEdit,'String',num2str(maxValue));
    set(handles.maxValue, 'Enable', 'off');
    set(handles.minValue, 'Enable', 'off');
    set(handles.maxEdit, 'Enable', 'off');
    set(handles.minEdit, 'Enable', 'off');
end

% --- Executes on button press in overlayOn.
function overlayOn_Callback(hObject, eventdata, handles)
% hObject    handle to overlayOn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of overlayOn
imageData{1,1}=getappdata(handles.axes1,'index');
imageData{1,2}=getappdata(handles.axes7,'index');
islice{1,1}=round(get(handles.sliceX,'Value'));
islice{1,2}=round(get(handles.sliceY,'Value'));
islice{1,3}=round(get(handles.sliceZ,'Value'));
alpha_val=get(handles.opaque1,'Value');
overlay=get(hObject,'Value');
if overlay == 1
    set(handles.opaque1,'Enable','on')
    set(handles.opaque2,'Enable','on')
    set(handles.opaque3,'Enable','on')
    hold(handles.axes1);
    alpha(imagesc(ind2rgb(squeeze(imageData{1,2}.Index(islice{1,1},:,:)),imageData{1,2}.Map),'Parent',handles.axes1),alpha_val);
    hold(handles.axes3);
    alpha(imagesc(ind2rgb(squeeze(imageData{1,2}.Index(:,islice{1,2},:)),imageData{1,2}.Map),'Parent',handles.axes3),alpha_val)    
    hold(handles.axes5);
    alpha(imagesc(ind2rgb(squeeze(imageData{1,2}.Index(:,:,islice{1,3})),imageData{1,2}.Map),'Parent',handles.axes5),alpha_val)
else
    set(handles.opaque1,'Enable','off')
    set(handles.opaque2,'Enable','off')
    set(handles.opaque3,'Enable','off')
    axes(handles.axes1)
    imagesc(ind2rgb(squeeze(imageData{1,1}.Index(islice{1,1},:,:)),imageData{1,1}.Map));
    hold off
    axes(handles.axes3)
    imagesc(ind2rgb(squeeze(imageData{1,1}.Index(:,islice{1,2},:)),imageData{1,1}.Map));
    hold off
    axes(handles.axes5)
    imagesc(ind2rgb(squeeze(imageData{1,1}.Index(:,:,islice{1,3})),imageData{1,1}.Map));    
    hold off
end

clear imageData
clear islice
% --- Executes on selection change in cMapmenu.
function cMapmenu_Callback(hObject, eventdata, handles)
% hObject    handle to cMapmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns cMapmenu contents as cell array
%        contents{get(hObject,'Value')} returns selected item from cMapmenu
updateImage(handles)    
    
% --- Executes during object creation, after setting all properties.
function cMapmenu_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cMapmenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on slider movement.
function sliceX_Callback(hObject, eventdata, handles)
% hObject    handle to sliceX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

overlay=get(handles.overlayOn,'Value');
sync=get(handles.synchronize,'Value');
islice=round(get(handles.sliceX,'Value'));
alpha_val=get(handles.opaque1,'Value');
imageData=cell(1,2);
imageData{1,1}=getappdata(handles.axes1,'index');
imageData{1,2}=getappdata(handles.axes7,'index');

if sync==1;
    axes(handles.axes1)
    imagesc(ind2rgb(squeeze(imageData{1,1}.Index(islice,:,:)),imageData{1,1}.Map));
    if overlay==1
        hold on;
        alpha(imagesc(ind2rgb(squeeze(imageData{1,2}.Index(islice,:,:)),imageData{1,2}.Map),'Parent',handles.axes1),alpha_val);
    end
    axes(handles.axes7)
    imagesc(ind2rgb(squeeze(imageData{1,2}.Index(islice,:,:)),imageData{1,2}.Map));
    set(handles.sliceX2,'Value',islice)
    islice=num2str(islice);
    set(handles.imageNum7,'String',islice)
    set(handles.imageNum,'String',islice)
elseif sync==0
    axes(handles.axes1)
    imagesc(ind2rgb(squeeze(imageData{1,1}.Index(islice,:,:)),imageData{1,1}.Map));
    if overlay==1
        hold on;
        alpha(imagesc(ind2rgb(squeeze(imageData{1,2}.Index(islice,:,:)),imageData{1,2}.Map),'Parent',handles.axes1),alpha_val);
    end
    islice=num2str(islice);
    set(handles.imageNum,'String',islice)
end    


% --- Executes during object creation, after setting all properties.
function sliceX_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliceX (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function opaque1_Callback(hObject, eventdata, handles)
% hObject    handle to opaque1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
overlay=get(handles.overlayOn,'Value');
sync=get(handles.synchronize,'Value');
islice=round(get(handles.sliceX,'Value'));
alpha_val=get(handles.opaque1,'Value');
imageData=cell(1,2);
imageData{1,1}=getappdata(handles.axes1,'index');
imageData{1,2}=getappdata(handles.axes7,'index');

axes(handles.axes1);
imagesc(ind2rgb(squeeze(imageData{1,1}.Index(islice,:,:)),imageData{1,1}.Map));
hold on;
alpha(imagesc(ind2rgb(squeeze(imageData{1,2}.Index(islice,:,:)),imageData{1,2}.Map),'Parent',handles.axes1),alpha_val);


% --- Executes during object creation, after setting all properties.
function opaque1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to opaque1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function sliceY_Callback(hObject, eventdata, handles)
% hObject    handle to sliceY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
overlay=get(handles.overlayOn,'Value');
sync=get(handles.synchronize,'Value');
islice=round(get(handles.sliceY,'Value'));
alpha_val=get(handles.opaque2,'Value');
imageData=cell(1,2);
imageData{1,1}=getappdata(handles.axes1,'index');
imageData{1,2}=getappdata(handles.axes7,'index');

if sync==1;
    axes(handles.axes3)
    imagesc(ind2rgb(squeeze(imageData{1,1}.Index(:,islice,:)),imageData{1,1}.Map));
    if overlay==1
        hold on
        alpha(imagesc(ind2rgb(squeeze(imageData{1,2}.Index(:,islice,:)),imageData{1,2}.Map),'Parent',handles.axes3),alpha_val)
    end
    axes(handles.axes8)
    imagesc(ind2rgb(squeeze(imageData{1,2}.Index(:,islice,:)),imageData{1,2}.Map));
    set(handles.sliceY2,'Value',islice)
    islice=num2str(islice);
    set(handles.imageNum8,'String',islice)
    set(handles.imageNum2,'String',islice)
elseif sync==0;
    axes(handles.axes3)
    imagesc(ind2rgb(squeeze(imageData{1,1}.Index(:,islice,:)),imageData{1,1}.Map));
    if overlay==1
        hold on
        alpha(imagesc(ind2rgb(squeeze(imageData{1,2}.Index(:,islice,:)),imageData{1,2}.Map),'Parent',handles.axes3),alpha_val)
    end
    islice=num2str(islice);
    set(handles.imageNum2,'String',islice)
end    

% --- Executes during object creation, after setting all properties.
function sliceY_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliceY (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function opaque2_Callback(hObject, eventdata, handles)
% hObject    handle to opaque2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
overlay=get(handles.overlayOn,'Value');
sync=get(handles.synchronize,'Value');
islice=round(get(handles.sliceY,'Value'));
alpha_val=get(handles.opaque2,'Value');
imageData=cell(1,2);
imageData{1,1}=getappdata(handles.axes1,'index');
imageData{1,2}=getappdata(handles.axes7,'index');

axes(handles.axes3)
imagesc(ind2rgb(squeeze(imageData{1,1}.Index(:,islice,:)),imageData{1,1}.Map));
hold on;
alpha(imagesc(ind2rgb(squeeze(imageData{1,2}.Index(:,islice,:)),imageData{1,2}.Map),'Parent',handles.axes3),alpha_val);

% --- Executes during object creation, after setting all properties.
function opaque2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to opaque2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function sliceZ_Callback(hObject, eventdata, handles)
% hObject    handle to sliceZ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
overlay=get(handles.overlayOn,'Value');
sync=get(handles.synchronize,'Value');
islice=round(get(handles.sliceZ,'Value'));
alpha_val=get(handles.opaque3,'Value');
imageData=cell(1,2);
imageData{1,1}=getappdata(handles.axes1,'index');
imageData{1,2}=getappdata(handles.axes7,'index');


if sync==1;
    axes(handles.axes5)
    imagesc(ind2rgb(squeeze(imageData{1,1}.Index(:,:,islice)),imageData{1,1}.Map));
    if overlay==1;
        hold on;
        alpha(imagesc(ind2rgb(squeeze(imageData{1,2}.Index(:,:,islice)),imageData{1,2}.Map),'Parent',handles.axes5),alpha_val)
    end
    axes(handles.axes9)
    imagesc(ind2rgb(squeeze(imageData{1,2}.Index(:,:,islice)),imageData{1,2}.Map));
    set(handles.sliceZ2,'Value',islice)
    islice=num2str(islice);
    set(handles.imageNum9,'String',islice)
    set(handles.imageNum3,'String',islice)    
elseif sync==0;
    axes(handles.axes5)
    imagesc(ind2rgb(squeeze(imageData{1,1}.Index(:,:,islice)),imageData{1,1}.Map));
    if overlay==1;
        hold on;
        alpha(imagesc(ind2rgb(squeeze(imageData{1,2}.Index(:,:,islice)),imageData{1,2}.Map),'Parent',handles.axes5),alpha_val)
    end
    islice=num2str(islice);
    set(handles.imageNum3,'String',islice)
end    


% --- Executes during object creation, after setting all properties.
function sliceZ_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliceZ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function opaque3_Callback(hObject, eventdata, handles)
% hObject    handle to opaque3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
overlay=get(handles.overlayOn,'Value');
sync=get(handles.synchronize,'Value');
islice=round(get(handles.sliceZ,'Value'));
alpha_val=get(handles.opaque3,'Value');
imageData=cell(1,2);
imageData{1,1}=getappdata(handles.axes1,'index');
imageData{1,2}=getappdata(handles.axes7,'index');

axes(handles.axes5)
imagesc(ind2rgb(squeeze(imageData{1,1}.Index(:,:,islice)),imageData{1,1}.Map));
hold on;
alpha(imagesc(ind2rgb(squeeze(imageData{1,2}.Index(:,:,islice)),imageData{1,2}.Map),'Parent',handles.axes5),alpha_val);


% --- Executes during object creation, after setting all properties.
function opaque3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to opaque3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function sliceX2_Callback(hObject, eventdata, handles)
% hObject    handle to sliceX2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

islice=round(get(hObject,'Value'));
imageData=getappdata(handles.axes7,'index');
axes(handles.axes7)
imagesc(ind2rgb(squeeze(imageData.Index(islice,:,:)),imageData.Map));
islice=num2str(islice);
set(handles.imageNum7,'String',islice)

% --- Executes during object creation, after setting all properties.
function sliceX2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliceX2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function sliceY2_Callback(hObject, eventdata, handles)
% hObject    handle to sliceY2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

islice=round(get(hObject,'Value'));
imageData=getappdata(handles.axes7,'index');
axes(handles.axes8)
imagesc(ind2rgb(squeeze(imageData.Index(:,islice,:)),imageData.Map));
islice=num2str(islice);
set(handles.imageNum8,'String',islice)
      
% --- Executes during object creation, after setting all properties.
function sliceY2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliceY2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function sliceZ2_Callback(hObject, eventdata, handles)
% hObject    handle to sliceZ2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

islice=round(get(hObject,'Value'));
imageData=getappdata(handles.axes7,'index');
axes(handles.axes9)
imagesc(ind2rgb(squeeze(imageData.Index(:,:,islice)),imageData.Map));
islice=num2str(islice);
set(handles.imageNum9,'String',islice)

% --- Executes during object creation, after setting all properties.
function sliceZ2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to sliceZ2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function maxValue_Callback(hObject, eventdata, handles)
% hObject    handle to maxValue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

maxValue = round(get(hObject,'Value'));
setappdata(handles.maxValue,'curMax',maxValue)
set(handles.maxEdit,'String',num2str(maxValue))

% --- Executes during object creation, after setting all properties.
function maxValue_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maxValue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

% --- Executes on slider movement.
function minValue_Callback(hObject, eventdata, handles)
% hObject    handle to minValue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
minValue = round(get(hObject,'Value'));
setappdata(handles.minValue,'curMin',minValue)
set(handles.minEdit,'String',num2str(minValue))

% --- Executes during object creation, after setting all properties.
function minValue_CreateFcn(hObject, eventdata, handles)
% hObject    handle to minValue (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function maxEdit_Callback(hObject, eventdata, handles)
% hObject    handle to maxEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of maxEdit as text
%        str2double(get(hObject,'String')) returns contents of maxEdit as a double
maxValue = get(hObject,'String');
maxValue = str2num(maxValue);
defMax = getappdata(handles.maxValue,'defMax');
curMin = getappdata(handles.minValue,'defMin');
if maxValue>=defMax
    maxValue = num2str(defMax);
    set(handles.maxEdit,'String',maxValue)
    set(handles.maxValue,'Value',defMax)
    setappdata(handles.maxValue,'curMax',defMax)
elseif maxValue<=curMin
    maxValue =num2str(curMin);
    set(handles.maxEdit,'String',minValue)
    set(handles.maxValue,'Value',curMin)
    setappdata(handles.maxValue,'curMax',curMin)
else
    set(handles.maxValue,'Value',maxValue)
    setappdata(handles.maxValue,'curMax',maxValue)
    maxValue=num2str(maxValue);
    set(handles.maxEdit,'String',maxValue)
end

% --- Executes during object creation, after setting all properties.
function maxEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maxEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function minEdit_Callback(hObject, eventdata, handles)
% hObject    handle to minEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of minEdit as text
%        str2double(get(hObject,'String')) returns contents of minEdit as a double
minValue = get(hObject,'String');
minValue = str2num(minValue);
curMax = round(getappdata(handles.maxValue,'curMax'));
defMin = round(getappdata(handles.minValue,'defMin'));
if minValue>=curMax
    minValue = num2str(curMax);
    set(handles.minEdit,'String',minValue)
    set(handles.minValue,'Value',curMax)
    setappdata(handles.minValue,'curMin',curMax)
elseif minValue<=defMin
    minValue =num2str(defMin);
    set(handles.minEdit,'String',minValue)
    set(handles.minValue,'Value',defMin)
    setappdata(handles.minValue,'curMin',defMin)
else
    set(handles.minValue,'Value',minValue)
    setappdata(handles.minValue,'curMin',minValue)
    minValue=num2str(minValue);
    set(handles.minEdit,'String',minValue)
end

% --- Executes during object creation, after setting all properties.
function minEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to minEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function imageNum_Callback(hObject, eventdata, handles)
% hObject    handle to imageNum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of imageNum as text
%        str2double(get(hObject,'String')) returns contents of imageNum as a double
islice=get(hObject,'String');
islice=round(str2num(islice));
overlay=get(handles.overlayOn,'Value');
sync=get(handles.synchronize,'Value');
alpha_val=get(handles.opaque1,'Value');
maxValue = get(handles.sliceX,'Max');
minValue = get(handles.sliceX,'Min');
if islice>=maxValue
    set(handles.sliceX,'Value',maxValue)
    islice=maxValue;
    maxValue = num2str(maxValue);
    set(handles.imageNum,'String',maxValue)
elseif islice<=minValue
    set(handles.sliceX,'Value',minValue)
    islice=minValue;
    minValue =num2str(minValue);
    set(handles.imageNum,'String',minValue)
else
    set(handles.sliceX,'Value',islice)
    islice=num2str(islice);
    set(handles.imageNum,'String',islice)
    islice=str2num(islice);
end

imageData=cell(1,2);
imageData{1,1}=getappdata(handles.axes1,'index');
imageData{1,2}=getappdata(handles.axes7,'index');

if sync==1;
    axes(handles.axes7)
    imagesc(ind2rgb(squeeze(imageData{1,2}.Index(islice,:,:)),imageData{1,2}.Map));
    axes(handles.axes1)
    imagesc(ind2rgb(squeeze(imageData{1,1}.Index(islice,:,:)),imageData{1,1}.Map));
    if overlay==1
        hold on;
        alpha(imagesc(ind2rgb(squeeze(imageData{1,2}.Index(islice,:,:)),imageData{1,2}.Map),'Parent',handles.axes1),alpha_val);
    end
    set(handles.sliceX2,'Value',islice)
    islice=num2str(islice);
    set(handles.imageNum7,'String',islice)
elseif sync==0
    axes(handles.axes1)
    imagesc(ind2rgb(squeeze(imageData{1,1}.Index(islice,:,:)),imageData{1,1}.Map));
    if overlay==1
        hold on;
        alpha(imagesc(ind2rgb(squeeze(imageData{1,2}.Index(islice,:,:)),imageData{1,2}.Map),'Parent',handles.axes1),alpha_val);
    end
end    
    
% --- Executes during object creation, after setting all properties.
function imageNum_CreateFcn(hObject, eventdata, handles)
% hObject    handle to imageNum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function imageNum2_Callback(hObject, eventdata, handles)
% hObject    handle to imageNum2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of imageNum2 as text
%        str2double(get(hObject,'String')) returns contents of imageNum2 as a double
islice=get(hObject,'String');
islice=round(str2num(islice));
overlay=get(handles.overlayOn,'Value');
sync=get(handles.synchronize,'Value');
alpha_val=get(handles.opaque2,'Value');
maxValue = get(handles.sliceY,'Max');
minValue = get(handles.sliceY,'Min');
if islice>=maxValue
    set(handles.sliceY,'Value',maxValue)
    islice=maxValue;
    maxValue = num2str(maxValue);
    set(handles.imageNum2,'String',maxValue)
elseif islice<=minValue
    set(handles.sliceY,'Value',minValue)
    islice=minValue;
    minValue =num2str(minValue);
    set(handles.imageNum2,'String',minValue)
else
    set(handles.sliceY,'Value',islice)
    islice=num2str(islice);
    set(handles.imageNum2,'String',islice)
    islice=str2num(islice);
end

imageData=cell(1,2);
imageData{1,1}=getappdata(handles.axes1,'index');
imageData{1,2}=getappdata(handles.axes7,'index');

if sync==1;
    axes(handles.axes8)
    imagesc(ind2rgb(squeeze(imageData{1,2}.Index(:,islice,:)),imageData{1,2}.Map));
    axes(handles.axes3)
    imagesc(ind2rgb(squeeze(imageData{1,1}.Index(:,islice,:)),imageData{1,1}.Map));
    if overlay==1
        hold on
        alpha(imagesc(ind2rgb(squeeze(imageData{1,2}.Index(:,islice,:)),imageData{1,2}.Map),'Parent',handles.axes3),alpha_val)
    end
    set(handles.sliceY2,'Value',islice)
    islice=num2str(islice);
    set(handles.imageNum8,'String',islice)
elseif sync==0;
    axes(handles.axes3)
    imagesc(ind2rgb(squeeze(imageData{1,1}.Index(:,islice,:)),imageData{1,1}.Map));
    if overlay==1
        hold on
        alpha(imagesc(ind2rgb(squeeze(imageData{1,2}.Index(:,islice,:)),imageData{1,2}.Map),'Parent',handles.axes3),alpha_val)
    end
end    

% --- Executes during object creation, after setting all properties.
function imageNum2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to imageNum2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function imageNum3_Callback(hObject, eventdata, handles)
% hObject    handle to imageNum3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of imageNum3 as text
%        str2double(get(hObject,'String')) returns contents of imageNum3 as a double
islice=get(hObject,'String');
islice=round(str2num(islice));
overlay=get(handles.overlayOn,'Value');
sync=get(handles.synchronize,'Value');
alpha_val=get(handles.opaque3,'Value');
maxValue = get(handles.sliceZ,'Max');
minValue = get(handles.sliceZ,'Min');
if islice>=maxValue;
    set(handles.sliceZ,'Value',maxValue)
    islice=maxValue;
    maxValue = num2str(maxValue);
    set(handles.imageNum3,'String',maxValue)
elseif islice<=minValue;
    set(handles.sliceZ,'Value',minValue)
    islice=minValue;
    minValue =num2str(minValue);
    set(handles.imageNum3,'String',minValue)
else
    set(handles.sliceZ,'Value',islice)
    islice=num2str(islice);
    set(handles.imageNum3,'String',islice)
    islice=str2num(islice);
end

imageData=cell(1,2);
imageData{1,1}=getappdata(handles.axes1,'index');
imageData{1,2}=getappdata(handles.axes7,'index');

if sync==1;
    axes(handles.axes9)
    imagesc(ind2rgb(squeeze(imageData{1,2}.Index(:,:,islice)),imageData{1,2}.Map));
    axes(handles.axes5)
    imagesc(ind2rgb(squeeze(imageData{1,1}.Index(:,:,islice)),imageData{1,1}.Map));
    if overlay==1;
        hold on;
        alpha(imagesc(ind2rgb(squeeze(imageData{1,2}.Index(:,:,islice)),imageData{1,2}.Map),'Parent',handles.axes5),alpha_val)
    end
    set(handles.sliceZ2,'Value',islice)
    islice=num2str(islice);
    set(handles.imageNum9,'String',islice)
elseif sync==0;
    axes(handles.axes5)
    imagesc(ind2rgb(squeeze(imageData{1,1}.Index(:,:,islice)),imageData{1,1}.Map));
    if overlay==1;
        hold on;
        alpha(imagesc(ind2rgb(squeeze(imageData{1,2}.Index(:,:,islice)),imageData{1,2}.Map),'Parent',handles.axes5),alpha_val)
    end
end    


% --- Executes during object creation, after setting all properties.
function imageNum3_CreateFcn(hObject, eventdata, handles)
% hObject    handle to imageNum3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function imageNum7_Callback(hObject, eventdata, handles)
% hObject    handle to imageNum7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of imageNum7 as text
%        str2double(get(hObject,'String')) returns contents of imageNum7 as a double
islice=get(hObject,'String');
islice=round(str2num(islice));
maxValue = get(handles.sliceX2,'Max');
minValue = get(handles.sliceX2,'Min');
if islice>=maxValue
    set(handles.sliceX2,'Value',maxValue)
    islice=maxValue;
    maxValue = num2str(maxValue);
    set(handles.imageNum7,'String',maxValue)
elseif islice<=minValue
    set(handles.sliceX2,'Value',minValue)
    islice=minValue;
    minValue =num2str(minValue);
    set(handles.imageNum7,'String',minValue)
else
    set(handles.sliceX2,'Value',islice)
    islice=num2str(islice);
    set(handles.imageNum7,'String',islice)
    islice=str2num(islice);
end

imageData=getappdata(handles.axes7,'index');
axes(handles.axes7)
imagesc(ind2rgb(squeeze(imageData.Index(islice,:,:)),imageData.Map));


% --- Executes during object creation, after setting all properties.
function imageNum7_CreateFcn(hObject, eventdata, handles)
% hObject    handle to imageNum7 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function imageNum8_Callback(hObject, eventdata, handles)
% hObject    handle to imageNum8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of imageNum8 as text
%        str2double(get(hObject,'String')) returns contents of imageNum8 as a double
islice=get(hObject,'String');
islice=round(str2num(islice));
maxValue = get(handles.sliceY2,'Max');
minValue = get(handles.sliceY2,'Min');
if islice>=maxValue
    set(handles.sliceY2,'Value',maxValue)
    islice=maxValue;
    maxValue = num2str(maxValue);
    set(handles.imageNum8,'String',maxValue)
elseif islice<=minValue
    set(handles.sliceY2,'Value',minValue)
    islice=minValue;
    minValue =num2str(minValue);
    set(handles.imageNum8,'String',minValue)
else
    set(handles.sliceY2,'Value',islice)
    islice=num2str(islice);
    set(handles.imageNum8,'String',islice)
    islice=str2num(islice);
end

imageData=getappdata(handles.axes7,'index');
axes(handles.axes8)
imagesc(ind2rgb(squeeze(imageData.Index(:,islice,:)),imageData.Map));

% --- Executes during object creation, after setting all properties.
function imageNum8_CreateFcn(hObject, eventdata, handles)
% hObject    handle to imageNum8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function imageNum9_Callback(hObject, eventdata, handles)
% hObject    handle to imageNum9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of imageNum9 as text
%        str2double(get(hObject,'String')) returns contents of imageNum9 as a double
islice=get(hObject,'String');
islice=round(str2num(islice));
maxValue = get(handles.sliceZ2,'Max');
minValue = get(handles.sliceZ2,'Min');
if islice>=maxValue
    set(handles.sliceZ2,'Value',maxValue)
    islice=maxValue;
    maxValue = num2str(maxValue);
    set(handles.imageNum9,'String',maxValue)
elseif islice<=minValue
    set(handles.sliceZ2,'Value',minValue)
    islice=minValue;
    minValue =num2str(minValue);
    set(handles.imageNum9,'String',minValue)
else
    set(handles.sliceZ2,'Value',islice)
    islice=num2str(islice);
    set(handles.imageNum9,'String',islice)
    islice=str2num(islice);
end

imageData=getappdata(handles.axes7,'index');
axes(handles.axes9)
imagesc(ind2rgb(squeeze(imageData.Index(:,:,islice)),imageData.Map));

% --- Executes during object creation, after setting all properties.
function imageNum9_CreateFcn(hObject, eventdata, handles)
% hObject    handle to imageNum9 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




% --- Executes on selection change in cMapmenu2.
function cMapmenu2_Callback(hObject, eventdata, handles)
% hObject    handle to cMapmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns cMapmenu2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from cMapmenu2
updateImage(handles)



% --- Executes during object creation, after setting all properties.
function cMapmenu2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to cMapmenu2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function maxValue2_Callback(hObject, eventdata, handles)
% hObject    handle to maxValue2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

maxValue = round(get(hObject,'Value'));
setappdata(handles.maxValue2,'curMax',maxValue)
set(handles.maxEdit2,'String',num2str(maxValue))


% --- Executes during object creation, after setting all properties.
function maxValue2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maxValue2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function minValue2_Callback(hObject, eventdata, handles)
% hObject    handle to minValue2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
minValue = round(get(hObject,'Value'));
setappdata(handles.minValue2,'curMin',minValue)
set(handles.minEdit2,'String',num2str(minValue))


% --- Executes during object creation, after setting all properties.
function minValue2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to minValue2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function maxEdit2_Callback(hObject, eventdata, handles)
% hObject    handle to maxEdit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of maxEdit2 as text
%        str2double(get(hObject,'String')) returns contents of maxEdit2 as a double
maxValue = get(hObject,'String');
maxValue = str2num(maxValue);
defMax = getappdata(handles.maxValue2,'defMax');
curMin = getappdata(handles.minValue2,'defMin');
if maxValue>=defMax
    maxValue = num2str(defMax);
    set(handles.maxEdit2,'String',maxValue)
    set(handles.maxValue2,'Value',defMax)
    setappdata(handles.maxValue2,'curMax',defMax)
elseif maxValue<=curMin
    maxValue =num2str(curMin);
    set(handles.maxEdit2,'String',minValue)
    set(handles.maxValue2,'Value',curMin)
    setappdata(handles.maxValue2,'curMax',curMin)
else
    set(handles.maxValue2,'Value',maxValue)
    setappdata(handles.maxValue2,'curMax',maxValue)
    maxValue=num2str(maxValue);
    set(handles.maxEdit2,'String',maxValue)
end

% --- Executes during object creation, after setting all properties.
function maxEdit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to maxEdit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function minEdit2_Callback(hObject, eventdata, handles)
% hObject    handle to minEdit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of minEdit2 as text
%        str2double(get(hObject,'String')) returns contents of minEdit2 as a double
minValue = get(hObject,'String');
minValue = str2num(minValue);
curMax = round(getappdata(handles.maxValue,'curMax'));
defMin = round(getappdata(handles.minValue2,'defMin'));
if minValue>=curMax
    minValue = num2str(curMax-0.0001);
    set(handles.minEdit2,'String',minValue)
    set(handles.minValue2,'Value',curMax)
    setappdata(handles.minValue2,'curMin',curMax)
elseif minValue<=defMin
    minValue =num2str(defMin);
    set(handles.minEdit2,'String',minValue)
    set(handles.minValue2,'Value',defMin)
    setappdata(handles.minValue2,'curMin',defMin)
else
    set(handles.minValue2,'Value',minValue)
    setappdata(handles.minValue2,'curMin',minValue)
    minValue=num2str(minValue);
    set(handles.minEdit2,'String',minValue)
end

% --- Executes during object creation, after setting all properties.
function minEdit2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to minEdit2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in autoWindow2.
function autoWindow2_Callback(hObject, eventdata, handles)
% hObject    handle to autoWindow2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of autoWindow2
auto_on = get(hObject,'Value');

if auto_on==0.0
    set(handles.maxValue2, 'Enable', 'on');
    set(handles.minValue2, 'Enable', 'on');
    set(handles.maxEdit2, 'Enable', 'on');
    set(handles.minEdit2, 'Enable', 'on');
    maxValue=round(get(handles.maxValue2,'Value'));
    minValue=round(get(handles.minValue2,'Value'));
    setappdata(handles.maxValue2,'curMax',maxValue)
    setappdata(handles.minValue2,'curMin',minValue)
else    
    minValue = round(getappdata(handles.maxValue2,'defMin'));
    maxValue = round(getappdata(handles.minValue2,'defMax'));
    setappdata(handles.maxValue2,'curMax',maxValue)
    setappdata(handles.minValue2,'curMin',minValue)
    set(handles.minValue2,'Min', minValue);
    set(handles.maxValue2,'Max', maxValue);
    set(handles.minValue2,'Min', minValue);
    set(handles.maxValue2,'Max', maxValue);
    set(handles.minValue2,'Value', minValue);
    set(handles.maxValue2,'Value', maxValue);
    set(handles.minEdit2,'String',num2str(minValue));
    set(handles.maxEdit2,'String',num2str(maxValue));
    set(handles.maxValue2, 'Enable', 'off');
    set(handles.minValue2, 'Enable', 'off');
    set(handles.maxEdit2, 'Enable', 'off');
    set(handles.minEdit2, 'Enable', 'off');
end


% --- Executes on button press in updateAll.
function updateAll_Callback(hObject, eventdata, handles)
% hObject    handle to updateAll (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
updateImage(handles)

