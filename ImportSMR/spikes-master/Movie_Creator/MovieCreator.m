function fig = MovieCreator()
% This is the machine-generated representation of a Handle Graphics object
% and its children.  Note that handle values may change when these objects
% are re-created. This may cause problems with any callbacks written to
% depend on the value of the handle at the time the object was saved.
% This problem is solved by saving the output as a FIG-file.
%
% To reopen this object, just type the name of the M-file at the MATLAB
% prompt. The M-file and its associated MAT-file must be on your path.
% 
% NOTE: certain newer features in MATLAB may not have been saved in this
% M-file due to limitations of this format, which has been superseded by
% FIG-files.  Figures which have been annotated using the plot editor tools
% are incompatible with the M-file/MAT-file format, and should be saved as
% FIG-files.

load MovieCreator

h0 = figure('Color',[0.8 0.8 0.8], ...
	'Colormap',mat0, ...
	'FileName','D:\matlab\RowlandRoutines\MovieCreator\MovieCreator.m', ...
	'MenuBar','none', ...
	'Name','Movie Creator', ...
	'NumberTitle','off', ...
	'PaperPosition',[18 180 576 432], ...
	'PaperUnits','points', ...
	'Position',[715 268 164 447], ...
	'Resize','off', ...
	'Tag','MovieToolFig', ...
	'ToolBar','none');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[0 0 0], ...
	'FontName','Courier', ...
	'FontSize',14, ...
	'FontWeight','bold', ...
	'ForegroundColor',[1 1 1], ...
	'ListboxTop',0, ...
	'Position',[4.5 297.75 111 29.25], ...
	'String','Temporal Movie Toolbox', ...
	'Style','text', ...
	'Tag','StaticText1');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[0.752941176470588 0.752941176470588 0.752941176470588], ...
	'ListboxTop',0, ...
	'Position',[5.25 214.5 110.25 80.25], ...
	'Style','frame', ...
	'Tag','Frame1');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[0.752941176470588 0.752941176470588 0.752941176470588], ...
	'FontSize',10, ...
	'FontWeight','light', ...
	'ListboxTop',0, ...
	'Position',[8.25 273 51.75 16.5], ...
	'String','Timeslice', ...
	'Style','text', ...
	'Tag','StaticText2');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[0.752941176470588 0.752941176470588 0.752941176470588], ...
	'FontSize',10, ...
	'FontWeight','light', ...
	'ListboxTop',0, ...
	'Position',[12.75 247.5 51.75 16.5], ...
	'String','Overlap', ...
	'Style','text', ...
	'Tag','StaticText2');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[1 1 1], ...
	'Position',[63 275.25 45 15], ...
	'String',' ', ...
	'Style','listbox', ...
	'Tag','TimeSliceMenu', ...
	'Value',1);
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[1 1 1], ...
	'Position',[63 248.25 45 15], ...
	'Style','listbox', ...
	'Tag','OverlapMenu', ...
	'Value',1);
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[0.752941176470588 0.752941176470588 0.752941176470588], ...
	'ListboxTop',0, ...
	'Position',[57.75 263.25 55.5 12], ...
	'String','(Milliseconds)', ...
	'Style','text', ...
	'Tag','StaticText3');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[0.752941176470588 0.752941176470588 0.752941176470588], ...
	'FontWeight','bold', ...
	'ListboxTop',0, ...
	'Position',[13.5 220.5 95.25 20.25], ...
	'String','Create Movie', ...
	'Tag','CreateMovieButton');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[0.752941176470588 0.752941176470588 0.752941176470588], ...
	'ListboxTop',0, ...
	'Position',[6 134.25 108.75 77.25], ...
	'Style','frame', ...
	'Tag','Frame2');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[0.752941176470588 0.752941176470588 0.752941176470588], ...
	'FontWeight','bold', ...
	'ListboxTop',0, ...
	'Position',[12.75 140.25 95.25 20.25], ...
	'String','Play Movie', ...
	'Tag','PlayMovieButton');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[0.752941176470588 0.752941176470588 0.752941176470588], ...
	'FontSize',10, ...
	'FontWeight','light', ...
	'ListboxTop',0, ...
	'Position',[10.5 189 51.75 16.5], ...
	'String','Speed', ...
	'Style','text', ...
	'Tag','StaticText2');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[0.752941176470588 0.752941176470588 0.752941176470588], ...
	'FontSize',10, ...
	'FontWeight','light', ...
	'ListboxTop',0, ...
	'Position',[11.25 169.5 51.75 16.5], ...
	'String','Loops', ...
	'Style','text', ...
	'Tag','StaticText2');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[1 1 1], ...
	'Callback',mat1, ...
	'ListboxTop',0, ...
	'Position',[60.75 192.75 45 15], ...
	'Style','edit', ...
	'Tag','EditSpeed');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[1 1 1], ...
	'Callback',mat2, ...
	'ListboxTop',0, ...
	'Position',[60.75 171.75 45 15], ...
	'Style','edit', ...
	'Tag','EditLoops');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[0.752941176470588 0.752941176470588 0.752941176470588], ...
	'ListboxTop',0, ...
	'Position',[6.75 44.25 108.75 87], ...
	'Style','frame', ...
	'Tag','Frame3');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[0.752941176470588 0.752941176470588 0.752941176470588], ...
	'Callback',mat3, ...
	'ListboxTop',0, ...
	'Position',[15 108.75 93 16.5], ...
	'String','Microsoft Video (.avi)', ...
	'Style','radiobutton', ...
	'Tag','ChooseAvi', ...
	'Value',1);
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[0.752941176470588 0.752941176470588 0.752941176470588], ...
	'Callback',mat4, ...
	'ListboxTop',0, ...
	'Position',[14.25 93 93 16.5], ...
	'String','MPEG Movie (.mpg)', ...
	'Style','radiobutton', ...
	'Tag','ChooseMpeg');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[0.752941176470588 0.752941176470588 0.752941176470588], ...
	'Callback',mat5, ...
	'ListboxTop',0, ...
	'Position',[14.25 77.25 93 16.5], ...
	'String','Matlab MAT file (.mat)', ...
	'Style','radiobutton', ...
	'Tag','ChooseMat');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[0.752941176470588 0.752941176470588 0.752941176470588], ...
	'FontWeight','demi', ...
	'ListboxTop',0, ...
	'Position',[13.5 51 95.25 20.25], ...
	'String','Export Movie', ...
	'Tag','ExportMovieButton');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[0.752941176470588 0.752941176470588 0.752941176470588], ...
	'ListboxTop',0, ...
	'Position',[7.5 7.5 108.75 32.25], ...
	'Style','frame', ...
	'Tag','Frame4');
h1 = uicontrol('Parent',h0, ...
	'Units','points', ...
	'BackgroundColor',[0.752941176470588 0.752941176470588 0.752941176470588], ...
	'FontWeight','bold', ...
	'ListboxTop',0, ...
	'Position',[15 13.5 95.25 20.25], ...
	'String','Load .MAT File', ...
	'Tag','LoadMatButton');
if nargout > 0, fig = h0; end
