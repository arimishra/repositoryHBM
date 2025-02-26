function varargout = vuSingleSpectraViewer(varargin)
% VUSINGLESPECTRAVIEWER M-file for vuSingleSpectraViewer.fig
%      VUSINGLESPECTRAVIEWER, by itself, creates a new VUSINGLESPECTRAVIEWER or raises the existing
%      singleton*.
%
%      H = VUSINGLESPECTRAVIEWER returns the handle to a new VUSINGLESPECTRAVIEWER or the handle to
%      the existing singleton*.
%
%      VUSINGLESPECTRAVIEWER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in VUSINGLESPECTRAVIEWER.M with the given input arguments.
%
%      VUSINGLESPECTRAVIEWER('Property','Value',...) creates a new VUSINGLESPECTRAVIEWER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before vuSingleSpectraViewer_OpeningFunction gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to vuSingleSpectraViewer_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help vuSingleSpectraViewer

% Last Modified by GUIDE v2.5 02-Dec-2009 15:23:52

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @vuSingleSpectraViewer_OpeningFcn, ...
                   'gui_OutputFcn',  @vuSingleSpectraViewer_OutputFcn, ...
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


% --- Executes just before vuSingleSpectraViewer is made visible.
function vuSingleSpectraViewer_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to vuSingleSpectraViewer (see VARARGIN)

% Choose default command line output for vuSingleSpectraViewer
handles.output = hObject;

% Check inputs
if length(varargin) < 1
	error('MATLAB:vuSpectraViewer:NotEnoughInputs', 'Not enough input arguments.');
end
if length(varargin) > 1
	error('MATLAB:vuSpectraViewer:TooManyInputs', 'Too many input arguments.');
end

% Save our inputs
spectra = varargin{1};
setappdata(gcf,'spectra',spectra);

% Process our data
processSpectra(handles);
processedSpec = getappdata(gcf,'processedSpec');

% Initial values
zeroFreqStart = 4.65;
setappdata(gcf,'zeroFreqValue',zeroFreqStart);
setappdata(gcf,'lineBroadeningValue',0);

% Calculate range
if isfield(processedSpec.Parms,'sample_frequency')
    
    % Bandwidth
    bandwidth = processedSpec.Parms.sample_frequency;
    
    if isfield(processedSpec.Parms,'synthesizer_frequency')
        
        % Scale
        bandwidth = bandwidth/(processedSpec.Parms.synthesizer_frequency*1e-6);
        
    end
    
else
    
    % Default
    bandwidth = 1;
    
end
xRange = linspace(-processedSpec.Dims(end)/2,processedSpec.Dims(end)/2 - 1,processedSpec.Dims(end));
xRange = xRange*(bandwidth/processedSpec.Dims(end)) + zeroFreqStart;
setappdata(gcf,'xRange',xRange)

% Setup Zoom Sliders
zoomLeft = 1;
zoomRight = processedSpec.Dims(end);
set(handles.zoomLeftSlider,'Min',-zoomRight)
set(handles.zoomLeftSlider,'Max',-zoomLeft-1)
set(handles.zoomLeftSlider,'SliderStep',[1/(zoomRight-zoomLeft-1) 50/(zoomRight-zoomLeft-1)]);
set(handles.zoomLeftSlider,'Value',-zoomRight);
set(handles.zoomLeftEdit,'String',round(xRange(zoomRight)*100)/100);
set(handles.zoomRightSlider,'Min',-zoomRight+1)
set(handles.zoomRightSlider,'Max',-zoomLeft)
set(handles.zoomRightSlider,'SliderStep',[1/(zoomRight-zoomLeft-1) 50/(zoomRight-zoomLeft-1)]);
set(handles.zoomRightSlider,'Value',-zoomLeft);
set(handles.zoomRightEdit,'String',round(xRange(1)*100)/100);

% Setup Phase correction
set(handles.phaseCorrectSlider,'Min',-180)
set(handles.phaseCorrectSlider,'Max',180)
set(handles.phaseCorrectSlider,'SliderStep',[1/360 10/360]);
set(handles.phaseCorrectSlider,'Value',0);
set(handles.phaseEdit,'String',0);

% Setup Zero Frequency
set(handles.zeroFreqSlider,'Min',-10)
set(handles.zeroFreqSlider,'Max',10)
set(handles.zeroFreqSlider,'SliderStep',[0.01 0.1]);
set(handles.zeroFreqSlider,'Value',zeroFreqStart);
set(handles.zeroFreqEdit,'String',zeroFreqStart);

% Plot
plotSpectrum(handles);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes vuSingleSpectraViewer wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = vuSingleSpectraViewer_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function phaseCorrectEdit_Callback(hObject, eventdata, handles)
% hObject    handle to phaseCorrectEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of phaseCorrectEdit as text
%        str2double(get(hObject,'String')) returns contents of phaseCorrectEdit as a double


% --- Executes during object creation, after setting all properties.
function phaseCorrectEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to phaseCorrectEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function zoomLeftSlider_Callback(hObject, eventdata, handles)
% hObject    handle to zoomLeftSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

% Slider movement - get new value
zoomLeftValue = -1*get(handles.zoomLeftSlider,'Value');
zoomLeftValue = round(zoomLeftValue);
zoomRightValue = -1*get(handles.zoomRightSlider,'Value');
xRange = getappdata(gcf,'xRange');

% Check if it's more than zoom right
if (zoomLeftValue <= zoomRightValue)
    zoomLeftValue = zoomRightValue - 1;
end

% Find ppm value
ppmValue = round(xRange(zoomLeftValue)*100)/100;

% Store value and set edit box
set(handles.zoomLeftSlider,'Value',-1*zoomLeftValue);
set(handles.zoomLeftEdit,'String',ppmValue);

% Plot current spectrum
plotSpectrum(handles);

% --- Executes during object creation, after setting all properties.
function zoomLeftSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to zoomLeftSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function zoomRightSlider_Callback(hObject, eventdata, handles)
% hObject    handle to zoomRightSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

% Slider movement - get new value
zoomRightValue = -1*get(handles.zoomRightSlider,'Value');
zoomRightValue = round(zoomRightValue);
zoomLeftValue = -1*get(handles.zoomLeftSlider,'Value');
xRange = getappdata(gcf,'xRange');

% Check if it's more than zoom right
if (zoomLeftValue <= zoomRightValue)
    zoomRightValue = zoomLeftValue + 1;
end

% Find ppm value
ppmValue = round(xRange(zoomRightValue)*100)/100;

% Store value and set edit box
set(handles.zoomRightSlider,'Value',-1*zoomRightValue);
set(handles.zoomRightEdit,'String',ppmValue);

% Plot current spectrum
plotSpectrum(handles);

% --- Executes during object creation, after setting all properties.
function zoomRightSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to zoomRightSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function zoomLeftEdit_Callback(hObject, eventdata, handles)
% hObject    handle to zoomLeftEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of zoomLeftEdit as text
%        str2double(get(hObject,'String')) returns contents of zoomLeftEdit as a double

% Edit box changed
zoomLeftText = get(handles.zoomLeftEdit,'String');
zoomLeftSliderValue = -1*get(handles.zoomLeftSlider,'Value');
zoomLeftMin = -1*get(handles.zoomLeftSlider,'Min');
zoomRightValue = -1*get(handles.zoomRightSlider,'Value');
xRange = getappdata(gcf,'xRange');

% Try to convert to number
[zoomLeftText isNum] = str2num(zoomLeftText);

% Check if it's a valid number
if isNum
    
    % Get the index value
    [tmp,zoomLeftValue] = min(abs(xRange-zoomLeftText));
    
    % Round 
    zoomLeftValue = round(zoomLeftValue);
    
    % Compare range
    if ~(zoomLeftValue > zoomRightValue && zoomLeftValue < zoomLeftMin)
       
        % Put back to original value
        zoomLeftValue = zoomLeftSliderValue;
        
    end
        
else

    % Put back to original value
    zoomLeftValue = zoomLeftSliderValue;
    
end

% Find ppm value
ppmValue = round(xRange(zoomLeftValue)*100)/100;

% Store value and set edit box
set(handles.zoomLeftSlider,'Value',-1*zoomLeftValue);
set(handles.zoomLeftEdit,'String',ppmValue);

% Plot current spectrum
plotSpectrum(handles);

% --- Executes during object creation, after setting all properties.
function zoomLeftEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to zoomLeftEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function zoomRightEdit_Callback(hObject, eventdata, handles)
% hObject    handle to zoomRightEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of zoomRightEdit as text
%        str2double(get(hObject,'String')) returns contents of zoomRightEdit as a double

% Edit box changed
zoomRightText = get(handles.zoomRightEdit,'String');
zoomRightSliderValue = -1*get(handles.zoomRightSlider,'Value');
zoomLeftValue = -1*get(handles.zoomLeftSlider,'Value');
xRange = getappdata(gcf,'xRange');

% Try to convert to number
[zoomRightText isNum] = str2num(zoomRightText);

% Check if it's a valid number
if isNum
    
    % Get the index value
    [tmp,zoomRightValue] = min(abs(xRange-zoomRightText));
    
    % Round 
    zoomRightValue = round(zoomRightValue);
    
    % Compare range
    if ~(zoomLeftValue > zoomRightValue && zoomRightValue > 0 )
       
        % Put back to original value
        zoomRightValue = zoomRightSliderValue;
        
    end
        
else

    % Put back to original value
    zoomRightValue = zoomRightSliderValue;
    
end

% Find ppm value
ppmValue = round(xRange(zoomRightValue)*100)/100;

% Store value and set edit box
set(handles.zoomRightSlider,'Value',-1*zoomRightValue);
set(handles.zoomRightEdit,'String',ppmValue);
  
% Plot current spectrum
plotSpectrum(handles);

% --- Executes during object creation, after setting all properties.
function zoomRightEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to zoomRightEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on slider movement.
function phaseCorrectSlider_Callback(hObject, eventdata, handles)
% hObject    handle to phaseCorrectSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

% Get value
phaseValue = get(handles.phaseCorrectSlider,'Value');

% Round
phaseValue = fix(phaseValue*100)/100;
set(handles.phaseCorrectSlider,'Value',phaseValue);

% Display
set(handles.phaseEdit,'String',phaseValue);

% Replot
plotSpectrum(handles);

% --- Executes during object creation, after setting all properties.
function phaseCorrectSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to phaseCorrectSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end

function plotSpectrum(handles)
% Function to plot the current selected spectrum

% Get the data
processedSpec = getappdata(gcf,'processedSpec');
zoomLeftValue = -1*get(handles.zoomLeftSlider,'Value');
zoomRightValue = -1*get(handles.zoomRightSlider,'Value');
xRange = getappdata(gcf,'xRange');

% Phase adjustment
% Get data
phaseValue = get(handles.phaseCorrectSlider,'Value');

% Multiply by phase value
processedSpec.Data = processedSpec.Data.*exp(-i.*(phaseValue*pi/180));

% Determine type of plot
realData = get(handles.realButton,'Value');
absData = get(handles.absButton,'Value');
powData = get(handles.powerButton,'Value');
if (realData)
    plotData = squeeze(real(processedSpec.Data));
elseif (absData)
    plotData = squeeze(abs(real(processedSpec.Data)));
else
    plotData = squeeze(abs(processedSpec.Data));
end

% Plot full range
axes(handles.spectraFig)
plot(xRange,plotData);
set(gca, 'XDir', 'reverse')
axis tight

% Get the min and max of the y-axis values
yMinMax = get(gca,'YLim');
ylim(yMinMax);

% Check zoom in range
xMinMax = get(gca,'XLim');
if (get(handles.zoomRangeCheckBox,'Value'))
    
    % Get x limits
    xMinMaxIdx(2) = -get(handles.zoomLeftSlider,'Value');
    xMinMaxIdx(1) = -get(handles.zoomRightSlider,'Value');
    xlim([xRange(xMinMaxIdx(1)) xRange(xMinMaxIdx(2))])
    
end

% Draw zoom lines
line([xRange(zoomLeftValue) xRange(zoomLeftValue)],yMinMax,'Color','r');
line([xRange(zoomRightValue) xRange(zoomRightValue)],yMinMax,'Color','r');

% Check for zero line
if get(handles.zeroLineCheckBox,'Value')
    
    % Draw zero line
    h = line(xMinMax,[0 0],'LineStyle','--','Color','g');
    
    % Save handle
    setappdata(gcf,'zeroLineHandle',h)
    
end

function processSpectra(handles)
% Function to plot the current selected spectrum

% Get raw data
spectra = getappdata(gcf,'spectra');

% Copy parmeters
processedSpec = spectra;

% Line broading
if (get(handles.lineBroadeningCheckBox,'Value'))
    
    % Get value
    lineBroadValue = getappdata(gcf,'lineBroadeningValue');
    
    % Determine type
    if get(handles.lbExpRadio,'Value')

        % 1D Lorentzian apodization with Lorentzian linebroadening
        processedSpec.Data = processedSpec.Data.*exp(-lineBroadValue*(0:(spectra.Dims(2)-1))*2*pi/(2*spectra.Parms.sample_frequency))';
        
    else
           
        % 1D Gaussian
        processedSpec.Data = processedSpec.Data.*exp(-(lineBroadValue*(0:(spectra.Dims(2)-1))*2*pi/(2*spectra.Parms.sample_frequency)).^2/(2*log(2)))';
        
    end

end

% Fourier transform
processedSpec.Data = fftshift(fft(processedSpec.Data));


% Save processsed data
setappdata(gcf,'processedSpec',processedSpec);


% --- Executes on button press in realButton.
function realButton_Callback(hObject, eventdata, handles)
% hObject    handle to realButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of realButton

% Plot
plotSpectrum(handles);


% --- Executes on button press in absButton.
function absButton_Callback(hObject, eventdata, handles)
% hObject    handle to absButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of absButton

% Plot
plotSpectrum(handles);


% --- Executes on button press in powerButton.
function powerButton_Callback(hObject, eventdata, handles)
% hObject    handle to powerButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of powerButton

% Plot
plotSpectrum(handles);


% --- Executes on slider movement.
function zeroFreqSlider_Callback(hObject, eventdata, handles)
% hObject    handle to zeroFreqSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider

% Get value
zeroFreqValue = get(handles.zeroFreqSlider,'Value');

% Round
zeroFreqValue = fix(zeroFreqValue*1000)/1000;
set(handles.zeroFreqSlider,'Value',zeroFreqValue);

% Display
set(handles.zeroFreqEdit,'String',zeroFreqValue);

% Edit range
xRange = getappdata(gcf,'xRange');
zeroFreqOld = getappdata(gcf,'zeroFreqValue');
xRange = xRange - zeroFreqOld + zeroFreqValue;
setappdata(gcf,'xRange',xRange);
setappdata(gcf,'zeroFreqValue',zeroFreqValue);

% Replot
plotSpectrum(handles);

% --- Executes during object creation, after setting all properties.
function zeroFreqSlider_CreateFcn(hObject, eventdata, handles)
% hObject    handle to zeroFreqSlider (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end



function zeroFreqEdit_Callback(hObject, eventdata, handles)
% hObject    handle to zeroFreqEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of zeroFreqEdit as text
%        str2double(get(hObject,'String')) returns contents of zeroFreqEdit as a double

% Edit box changed
zeroFreqText = get(handles.zeroFreqEdit,'String');
zeroFreqSliderValue = get(handles.zeroFreqSlider,'Value');
zeroFreqSliderMin = get(handles.zeroFreqSlider,'Min');
zeroFreqSliderMax = get(handles.zeroFreqSlider,'Max');

% Try to convert to number
[zeroFreqText isNum] = str2num(zeroFreqText);

% Check if it's a valid number
if isNum
    
    % Round
    zeroFreqText = fix(zeroFreqText*1000)/1000;
    
    % Compare range
    if (zeroFreqText < zeroFreqSliderMin) || (zeroFreqText > zeroFreqSliderMax)
       
        % Put back to original value
        zeroFreqText = zeroFreqSliderValue;
        
    end
        
else

    % Put back to original value
    zeroFreqText = zeroFreqSliderValue;
    
end

% Store value and set edit box
set(handles.zeroFreqSlider,'Value',zeroFreqText);
set(handles.zeroFreqEdit,'String',zeroFreqText);
  
% Edit range
xRange = getappdata(gcf,'xRange');
zeroFreqOld = getappdata(gcf,'zeroFreqValue');
xRange = xRange - zeroFreqOld + zeroFreqText;
setappdata(gcf,'xRange',xRange);
setappdata(gcf,'zeroFreqValue',zeroFreqText);

% Plot current spectrum
plotSpectrum(handles);

% --- Executes during object creation, after setting all properties.
function zeroFreqEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to zeroFreqEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in lineBroadeningCheckBox.
function lineBroadeningCheckBox_Callback(hObject, eventdata, handles)
% hObject    handle to lineBroadeningCheckBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of lineBroadeningCheckBox

processSpectra(handles);

% Plot current spectrum
plotSpectrum(handles);

function lineBroadeningEdit_Callback(hObject, eventdata, handles)
% hObject    handle to lineBroadeningEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lineBroadeningEdit as text
%        str2double(get(hObject,'String')) returns contents of lineBroadeningEdit as a double

% Get data
lineBroadText = get(handles.lineBroadeningEdit,'String');

% Try to convert to number
[lineBroadText isNum] = str2num(lineBroadText);

if (isNum)
    
    % Save Value
    setappdata(gcf,'lineBroadeningValue',lineBroadText);
    processSpectra(handles);
    
    % Plot current spectrum
    plotSpectrum(handles);
    
else
    
    % Reset
    lineBroadValue = getappdata(gcf,'lineBroadeningValue');
    set(handles.lineBroadeningEdit,'String',lineBroadValue);

end


% --- Executes during object creation, after setting all properties.
function lineBroadeningEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lineBroadeningEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in analyzeButton.
function analyzeButton_Callback(hObject, eventdata, handles)
% hObject    handle to analyzeButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get the data
processedSpec = getappdata(gcf,'processedSpec');
zoomLeftValue = -1*get(handles.zoomLeftSlider,'Value');
zoomRightValue = -1*get(handles.zoomRightSlider,'Value');
xRange = getappdata(gcf,'xRange');

% Phase adjustment
% Get data
phaseValue = get(handles.phaseCorrectSlider,'Value');

% Multiply by phase value
processedSpec.Data = processedSpec.Data.*exp(-i.*phaseValue);

% Determine type of plot
realData = get(handles.realButton,'Value');
absData = get(handles.absButton,'Value');
powData = get(handles.powerButton,'Value');
if (realData)
    plotData = squeeze(real(processedSpec.Data));
elseif (absData)
    plotData = squeeze(abs(real(processedSpec.Data)));
else
    plotData = squeeze(abs(processedSpec.Data));
end

% Find maximum in range
[maxValue,idx] = max(plotData(zoomRightValue:zoomLeftValue));
idx = idx - 1 + zoomRightValue;

% Half
halfMaxValue = maxValue/2;

% Find Left/Right index where goes below half
leftIdx = idx;
while (plotData(leftIdx) > halfMaxValue)
    leftIdx = leftIdx + 1;
end
rightIdx = idx;
while (plotData(rightIdx) > halfMaxValue)
    rightIdx = rightIdx - 1;
end


% Linearly interpret
leftIdxFrac = (plotData(leftIdx-1) - halfMaxValue) / (plotData(leftIdx-1) - plotData(leftIdx));
rightIdxFrac = (halfMaxValue - plotData(rightIdx)) / (plotData(rightIdx+1) - plotData(rightIdx));

% Get values in ppm
leftPPMValue = xRange(leftIdx-1) + leftIdxFrac*(abs(xRange(leftIdx-1)-xRange(leftIdx)));
rightPPMValue = xRange(rightIdx) + rightIdxFrac*(abs(xRange(rightIdx+1)-xRange(rightIdx)));

% FWHM
FWHM = leftPPMValue - rightPPMValue;
FWHM_Hz = FWHM*processedSpec.Parms.synthesizer_frequency*1e-6;

% Area
area = trapz(xRange(zoomRightValue:zoomLeftValue),abs(plotData(zoomRightValue:zoomLeftValue)));

% Set labels
set(handles.peakLabel,'String',['Peak Value : ' sprintf('%0.2f',maxValue)])
set(handles.fwhmLabel,'String',['FWHM : ' sprintf('%0.2f',FWHM_Hz) 'Hz'])
set(handles.areaLabel,'String',['Area : ' sprintf('%0.2f',area)])

function phaseEdit_Callback(hObject, eventdata, handles)
% hObject    handle to phaseEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of phaseEdit as text
%        str2double(get(hObject,'String')) returns contents of phaseEdit as a double
% Edit box changed
phaseText = get(handles.phaseEdit,'String');
phaseSliderValue = get(handles.phaseCorrectSlider,'Value');
phaseSliderMin = get(handles.phaseCorrectSlider,'Min');
phaseSliderMax = get(handles.phaseCorrectSlider,'Max');

% Try to convert to number
[phaseText isNum] = str2num(phaseText);

% Check if it's a valid number
if isNum
    
    % Round
    phaseText = fix(phaseText*1000)/1000;
    
    % Compare range
    if (phaseText < phaseSliderMin) || (phaseText > phaseSliderMax)
       
        % Put back to original value
        phaseText = phaseSliderValue;
        
    end
        
else

    % Put back to original value
    phaseText = phaseSliderValue;
    
end

% Store value and set edit box
set(handles.phaseCorrectSlider,'Value',phaseText);
set(handles.phaseEdit,'String',phaseText);

% Plot current spectrum
plotSpectrum(handles);

% --- Executes during object creation, after setting all properties.
function phaseEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to phaseEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in zoomRangeCheckBox.
function zoomRangeCheckBox_Callback(hObject, eventdata, handles)
% hObject    handle to zoomRangeCheckBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of zoomRangeCheckBox

% Plot current spectrum
plotSpectrum(handles);


% --- Executes on button press in noneRadio.
function noneRadio_Callback(hObject, eventdata, handles)
% hObject    handle to noneRadio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of noneRadio

axes(handles.spectraFig)
datacursormode off
zoom off

% Plot current spectrum
plotSpectrum(handles);


% --- Executes on button press in zoomRadio.
function zoomRadio_Callback(hObject, eventdata, handles)
% hObject    handle to zoomRadio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of zoomRadio

axes(handles.spectraFig)
datacursormode off
h = zoom;
set(h,'Enable','on');


% --- Executes on button press in dataCursorRadio.
function dataCursorRadio_Callback(hObject, eventdata, handles)
% hObject    handle to dataCursorRadio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of dataCursorRadio

axes(handles.spectraFig)
zoom off
datacursormode on


% --- Executes on button press in lbExpRadio.
function lbExpRadio_Callback(hObject, eventdata, handles)
% hObject    handle to lbExpRadio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of lbExpRadio

processSpectra(handles);

% Plot current spectrum
plotSpectrum(handles);

% --- Executes on button press in lbGaussRadio.
function lbGaussRadio_Callback(hObject, eventdata, handles)
% hObject    handle to lbGaussRadio (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of lbGaussRadio

processSpectra(handles);

% Plot current spectrum
plotSpectrum(handles);


% --- Executes on button press in zeroLineCheckBox.
function zeroLineCheckBox_Callback(hObject, eventdata, handles)
% hObject    handle to zeroLineCheckBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of zeroLineCheckBox

axes(handles.spectraFig)
xMinMax = get(gca,'XLim');

% Check for zero line
if get(handles.zeroLineCheckBox,'Value')
    
    % Draw zero line
    h = line(xMinMax,[0 0],'LineStyle','--','Color','g');
    
    % Save handle
    setappdata(gcf,'zeroLineHandle',h)
    
else
    
    % Delete zero line
    delete(getappdata(gcf,'zeroLineHandle'))
    
end

