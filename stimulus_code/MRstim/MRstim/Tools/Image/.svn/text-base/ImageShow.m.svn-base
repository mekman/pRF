function h = ImageShow(Object, FigTitle, bVal, wVal, Map, Mark, markPixels)
%  ImageShow.m	Function to show an image in its true pixel size with axis off.
%				imshow is used in this function.
%
%  By YZW, 04/02/97
%
%  Usage:  h = ImageShow(Object, FigTitle, bVal, wVal, Map, Mark, markPixels);
%			Object:		input image to be shown
%			FigTitle:	string of figure name
%			bVal:		pixel value displays as black
%			wVal:		pixel value displays as white
%						Pixel values between bVal and wVal display as
%						intermediate shades of gray
%			Map:		colormap for image display
%			Mark:		string of marks
%			markPixels:	number of pixels corresponding to the mark
%			h:			returned handle to the figure
%
%  Default values:
%			Object:		required
%			FigTitle:	= ''
%			bVal		= min(min(Object))
%			wVal		= max(max(Object))
%			Map			= gray(256)
%			Mark		= ''
%			markPixels	= []
%

% default values
if nargin==0						% no input parameter
	help ImageShow;
	return;
elseif nargin==1
	FigTitle=''; bVal=[]; wVal=[]; Map=[]; Mark=''; markPixels=[];
elseif nargin == 2
	bVal=[]; wVal=[]; Map=[]; Mark=''; markPixels=[];
elseif nargin==3
	wVal=[]; Map=[]; Mark=''; markPixels=[];
elseif nargin==4
	Map=[]; Mark=''; markPixels=[];
elseif nargin==5
	Mark=''; markPixels=[];
elseif nargin==6
	disp('Mark must be paired with markPixels.');
	return;
end;

if (~ischar(FigTitle))
	disp('FigTitle must be a string');
	FigTitle = '';
end

if (isempty(bVal))	bVal = min(min(Object));	end
if (isempty(wVal))	wVal = max(max(Object));	end
if (isempty(Map))	Map  = gray(256);			end

if (bVal==wVal)		bVal=0;		wVal=1;			end

if (~ischar(Mark))
	disp('Mark must be a string');
	Mark = '';
end

h = figure('visible', 'off');

imshow(Object, [bVal wVal]);
colormap(Map);
axis('off');
truesize;

if (~isempty(FigTitle))			% put on figure title
	set(h, 'Name', FigTitle);
	set(h, 'NumberTitle', 'off');
end

if (~isempty(markPixels))			% put on scale mark
	[m, n] = size(Object);
	h_line = line([n-markPixels-10 n-10], [m-10 m-10]);
	set(h_line, 'Color', 'yellow');
	h_line = line([n-markPixels-10 n-markPixels-10], [m-10 m-10-5]);
	set(h_line, 'Color', 'yellow');
	h_line = line([n-10 n-10], [m-10 m-10-5]);
	set(h_line, 'Color', 'yellow');
	h_text = text(n-markPixels-15,m-20,Mark);
	set(h_text, 'Color', 'yellow');
end

set(h, 'visible', 'on');			% show figure

% end of function IMAGESHOW
