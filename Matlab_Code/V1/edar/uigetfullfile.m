function f = uigetfullfile( FilterSpec, DialogTitle, DefaultFile )
%UIGETFULLFILE   Standard open file dialog box returning full file names
%   F = UIGETFULLFILE displays a dialog box for the user to fill in, and
%   returns a full file specification F.
%
%   F = UIGETFULLFILE(FILTERSPEC)              See help UIGETFILE
%   F = UIGETFULLFILE(FILTERSPEC, TITLE)       See help UIGETFILE
%   F = UIGETFULLFILE(FILTERSPEC, TITLE, FILE) See help UIGETFILE
%
%   Multiple file selection is enabled.
%
%   See also UIPUTFULLFILE, UIGETFILE, FULLFILE

% Copyright 2013-2020 Sigma Connectivity AB
% Author: Claes Hovmalm
% $Id: uigetfullfile.m 8 2020-03-31 12:08:26Z hovm1011 $
% $HeadURL: file:///C:/Users/hovm1011/OneDrive%20-%20Sigma%20AB/Tools/MATLAB_repository/trunk/chutilities/uigetfullfile.m $

narginchk(0,3)

switch nargin
    case 0
        [FileName,PathName] = uigetfile('*.*','MultiSelect','on') ;
    case 1
        [FileName,PathName] = uigetfile( FilterSpec, 'MultiSelect','on' ) ;
    case 2
        [FileName,PathName] = uigetfile( FilterSpec, DialogTitle, 'MultiSelect','on' ) ;
    case 3
        [FileName,PathName] = uigetfile( FilterSpec, DialogTitle, DefaultFile, 'MultiSelect','on' ) ;
    otherwise
        error([NameOfCurrentFunction ': Incorrect number of input arguments.'])
end

% if ~isempty(FileName)
%     f = fullfile(PathName,FileName) ;
% else
%     f = 0 ;
% end

% if FileName
%     f = fullfile(PathName,FileName) ;
% else
%     f = 0 ;
% end

if ischar(FileName)||iscellstr(FileName)
    f = fullfile(PathName,FileName) ;
else
    f = 0 ;
end
   
end

