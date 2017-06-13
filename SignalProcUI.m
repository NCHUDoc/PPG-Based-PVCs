function varargout = SignalProcUI(varargin)
% SIGNALPROCUI MATLAB code for SignalProcUI.fig
%      SIGNALPROCUI, by itself, creates a new SIGNALPROCUI or raises the existing
%      singleton*.
%
%      H = SIGNALPROCUI returns the handle to a new SIGNALPROCUI or the handle to
%      the existing singleton*.
%
%      SIGNALPROCUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in SIGNALPROCUI.M with the given input arguments.
%
%      SIGNALPROCUI('Property','Value',...) creates a new SIGNALPROCUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before SignalProcUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to SignalProcUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help SignalProcUI

% Last Modified by GUIDE v2.5 09-Jun-2017 19:15:51

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @SignalProcUI_OpeningFcn, ...
                   'gui_OutputFcn',  @SignalProcUI_OutputFcn, ...
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


% --- Executes just before SignalProcUI is made visible.
function SignalProcUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to SignalProcUI (see VARARGIN)

% Choose default command line output for SignalProcUI
handles.output = hObject;

handles.line_signal=line('parent',handles.axes_sig,...
    'xdata',[],...
    'ydata',[],...
    'linewidth',2,'color',[1 0 0],...
    'linestyle','-');
% handles.line_signal2=line('parent',handles.axes_sig,...%     'xdata',[],...%     'ydata',[],...
%     'linewidth',2,'color',[1 0 0],...
%     'linestyle','-');

set(handles.slider_sig,'value',1);
set(handles.pop_filter,'value',3);
handles.File={}; 
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes SignalProcUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = SignalProcUI_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on selection change in listbox_subject.
function listbox_subject_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_subject (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_subject contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_subject


% --- Executes during object creation, after setting all properties.
function listbox_subject_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_subject (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listbox_channel.
function listbox_channel_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_channel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_channel contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_channel
% ch=get(handles.listbox_channel,'value');
% 
% set(handles.line_signal,'xdata',[1:size(handles.Data.signal,1)],...
%                         'ydata',handles.Data.signal(:,ch))

pop_dispmode_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function listbox_channel_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_channel (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on button press in load.
function load_Callback(hObject, eventdata, handles)
% hObject    handle to load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[handles.Data.filename, handles.Data.filepath]=uigetfile('*.dat','Please select a data file');
% add by mingfan 20170614
if handles.Data.filename == 0
    return
else
 handles.Data.signal = load(handles.Data.filename, handles.Data.filepath,'-ascii');
end
% [filename filepath]=uigetfile('*.dat');
% fullpath = [filepath filename];
% signal = load(fullpath);
% handles.Data.signal = load([filepath filename]); % 40000*40


subjnum = length(handles.File)+1;
for i=1:subjnum
    liststr_subj{i}=['Subject #' num2str(i)];
end
set(handles.listbox_subject,'string',liststr_subj,'value',subjnum)

chnum = size(handles.Data.signal,2);
                   
for i=1:chnum
    liststr_ch{i}=['Channel #' num2str(i)];
end
set(handles.listbox_channel,'string',liststr_ch,'value',1);

set(handles.line_signal,'xdata',[1:size(handles.Data.signal,1)],...
                        'ydata',handles.Data.signal(:,1))


guidata(hObject, handles);

% --- Executes on button press in print.
function print_Callback(hObject, eventdata, handles)
% hObject    handle to print (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[filename filepath]=uiputfile('*.jpeg');
fullpath=[filepath filename];

set(gcf,'PaperUnits','normalized',...
              'PaperType','A4',...
              'PaperPositionMode','auto');

print(gcf,'-djpeg',fullpath);
% --- Executes on button press in save.
function save_Callback(hObject, eventdata, handles)
% hObject    handle to save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
Data=handles.Data;

[filename filepath]=uiputfile('*.mat');
fullpath=[filepath filename];

save(fullpath,'Data')
% --- Executes on button press in reset.
function reset_Callback(hObject, eventdata, handles)
% hObject    handle to reset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close gcf
SignalProcUI

% --- Executes on button press in close.
function close_Callback(hObject, eventdata, handles)
% hObject    handle to close (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
close gcf



% --- Executes on slider movement.
function slider_sig_Callback(hObject, eventdata, handles)
% hObject    handle to slider_sig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
pop_dispmode_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function slider_sig_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_sig (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end




% pop_dispmode_Callback(hObject,eventdata,handles)







% --- Executes on button press in filter.
function filter_Callback(hObject, eventdata, handles)
% hObject    handle to filter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
filtertype = get(handles.pop_filter,'value');
cutoff=str2num(get(handles.edit_filterpara,'string'));
samplerate=str2num(get(handles.edit_samplerate,'string'));

if filtertype == 1
    type = 'low';
elseif filtertype == 2
    type = 'high';
elseif filtertype == 3
    type = 'bandpass';
elseif filtertype ==4
    type = 'stop';
end

handles.Data.signal = filter_2sIIR(handles.Data.signal', [cutoff(1) cutoff(end)], samplerate, 6, type);

guidata(hObject,handles)

pop_dispmode_Callback(hObject,eventdata,handles)


% --- Executes on selection change in pop_filter.
function pop_filter_Callback(hObject, eventdata, handles)
% hObject    handle to pop_filter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pop_filter contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pop_filter


% --- Executes during object creation, after setting all properties.
function pop_filter_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pop_filter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function edit_onset_Callback(hObject, eventdata, handles)
% hObject    handle to edit_onset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_onset as text
%        str2double(get(hObject,'String')) returns contents of edit_onset as a double


% --- Executes during object creation, after setting all properties.
function edit_onset_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_onset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_period_Callback(hObject, eventdata, handles)
% hObject    handle to edit_period (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_period as text
%        str2double(get(hObject,'String')) returns contents of edit_period as a double


% --- Executes during object creation, after setting all properties.
function edit_period_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_period (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
% --- Executes on selection change in pop_dispmode.
function pop_dispmode_Callback(hObject, eventdata, handles)
% hObject    handle to pop_dispmode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pop_dispmode contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pop_dispmode
ch=get(handles.listbox_channel,'value');
dispmode = get(handles.pop_dispmode,'value');

if dispmode == 1  % whole series
    set(handles.axes_sig,'xlim',[1 size(handles.Data.signal,1)])
	set(handles.line_signal,'xdata',[1:size(handles.Data.signal,1)],...
       'ydata',handles.Data.signal(:,ch))
    set(handles.slider_sig,'enable','off');
elseif dispmode ==2 % segment
	displength = round(str2num(get(handles.edit_displength,'string')));
	datalength=size(handles.Data.signal,1);
    
    
	set(handles.slider_sig,'enable','on',...
						'min',1,'max',ceil(datalength/displength),...
                        'sliderstep',[1/ceil(datalength/displength) 3/ceil(datalength/displength)]);
	segment = round(get(handles.slider_sig,'value'));
	set(handles.axes_sig,'xlim',[displength*(segment-1)+1  displength*segment])
	set(handles.line_signal,'xdata',[1:size(handles.Data.signal,1)],...
                       'ydata',handles.Data.signal(:,ch))
elseif dispmode ==3 % event avg
    period=round(str2num(get(handles.edit_period,'string')));
    
    set(handles.axes_sig,'xlim',[period(1) period(end)])
    set(handles.line_signal,'xdata',[period(1):period(end)],...
        'ydata',handles.Data.evtaverage(:,ch))
    
    set(handles.slider_sig,'enable','off') 
end


check_onset = get(handles.checkbox_onset,'value');
if isfield(handles,'line_onset')
    set(handles.line_onset(:),'xdata',[],'ydata',[])
end

if check_onset ==1 & dispmode~=3
    onset = round(str2num(get(handles.edit_onset,'string')));
    
    for i=1:length(onset)
        handles.line_onset(i)=line('parent',handles.axes_sig,...
            'xdata',[onset(i) onset(i)],...
            'ydata',[min(handles.Data.signal(:,ch)) max(handles.Data.signal(:,ch))],...
            'color','b','linewidth',1.5,...
            'linestyle','--');
    end
end  

guidata(hObject,handles)
% --- Executes during object creation, after setting all properties.
function pop_dispmode_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pop_dispmode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function edit_displength_Callback(hObject, eventdata, handles)
% hObject    handle to edit_displength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_displength as text
%        str2double(get(hObject,'String')) returns contents of edit_displength as a double
set(handles.slider_sig,'value',1)  % prevent out-of-range value
pop_dispmode_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function edit_displength_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_displength (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit_filterpara_Callback(hObject, eventdata, handles)
% hObject    handle to edit_filterpara (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_filterpara as text
%        str2double(get(hObject,'String')) returns contents of edit_filterpara as a double


% --- Executes during object creation, after setting all properties.
function edit_filterpara_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_filterpara (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end





% --- Executes on button press in evtaverage.
function evtaverage_Callback(hObject, eventdata, handles)
% hObject    handle to evtaverage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
onset=round(str2num(get(handles.edit_onset,'string')));
period=round(str2num(get(handles.edit_period,'string')));

epoch=zeros(length([period(1):period(end)]),size(handles.Data.signal,2),length(onset));
for i=1:length(onset)
    for j=1:size(handles.Data.signal,2)
        epoch(:,j,i)=handles.Data.signal(onset(i)+period(1):onset(i)+period(end),j);
    end
end
handles.Data.evtaverage=mean(epoch,3);

guidata(hObject,handles)
set(handles.pop_dispmode,'value',3)
pop_dispmode_Callback(hObject,eventdata,handles);




function edit_samplerate_Callback(hObject, eventdata, handles)
% hObject    handle to edit_samplerate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_samplerate as text
%        str2double(get(hObject,'String')) returns contents of edit_samplerate as a double


% --- Executes during object creation, after setting all properties.
function edit_samplerate_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_samplerate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in checkbox_onset.
function checkbox_onset_Callback(hObject, eventdata, handles)
% hObject    handle to checkbox_onset (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of checkbox_onset
pop_dispmode_Callback(hObject, eventdata, handles)
