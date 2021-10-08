function f = uiputfullfile( FilterSpec, DialogTitle, DefaultFile )
%UIPUTFULLFILE   Standard save file dialog box returning full file name
%   F = UIGETFULLFILE displays a dialog box for the user to fill in, and
%   returns a full file specification F.
%
%   F = UIGETFULLFILE(FILTERSPEC)              See help UIPUTFILE
%   F = UIGETFULLFILE(FILTERSPEC, TITLE)       See help UIPUTFILE
%   F = UIGETFULLFILE(FILTERSPEC, TITLE, FILE) See help UIPUTFILE
%
%   See also UIGETFULLFILE, UIPUTFILE, FULLFILE

% Copyright 2013-2020 Sigma Connectivity AB
% Author: Claes Hovmalm
% $Id: uiputfullfile.m 8 2020-03-31 12:08:26Z hovm1011 $
% $HeadURL: file:///C:/Users/hovm1011/OneDrive%20-%20Sigma%20AB/Tools/MATLAB_repository/trunk/chutilities/uiputfullfile.m $


narginchk(0,3)

switch nargin
    case 0
        [FileName,PathName] = uiputfile ;
    case 1
        [FileName,PathName] = uiputfile( FilterSpec ) ;
    case 2
        [FileName,PathName] = uiputfile( FilterSpec, DialogTitle ) ;
    case 3
        [FileName,PathName] = uiputfile( FilterSpec, DialogTitle, DefaultFile ) ;
    otherwise
        error([NameOfCurrentFunction ': Incorrect number of input arguments.'])
end

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
