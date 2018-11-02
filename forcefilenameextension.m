function [ newextensionFilenameStr, baseFilenameStr ] = forcefilenameextension( sourceFilenameStr, extensionStr, varargin )
%
%   [ newextensionFilenameStr, baseFilenameStr ] = forceFilenameExtension( sourceFilenameStr, extensionStr [, filenameSuffixStr ] );
%
%   This function performs two functions. 
%
%   First, it takes an input filename (sourceFilenameStr), which may
%   include path information and a file extension, and creates a new
%   filename string (newextensionFilenameStr) that is the same as the input
%   filename but with a new extension (extensionStr).
%
%
%   Second, it returned the base (baseFilenameStr) of the input filename
%   string (sourceFilenameStr) where the base if the filename stripped of
%   any path and extension information. This string is useful for titling
%   in plots and naming results files.
%
%   If an optional, third argument filenameSuffixStr is provided, the
%   filename of both the outputs will have a suffix of filenameSuffixStr
%
%   Example
%
%       [ a, b ] = forceFilenameExtension ('c:/foo/foobar.txt', 'pdf', '_suffix' );
%
%       will result in a='c:/foo/foobar_suffix.pdf' and b='foobar_suffix'
%
%   The code is a simple application of fileparts and fullfile
%
%   $author$, University of Utah
%   $date$
%

%% Default outputs
newextensionFilenameStr = '';numel( newextensionFilenameStr );
baseFilenameStr = '';numel( baseFilenameStr );

%% Check input number and types
pInput = inputParser();
addRequired( pInput, 'sourceFilenameStr', ...
    @(x)( ischar(x) ) );
addRequired( pInput, 'extensionStr', ...
    @(x)( ischar(x) ) );
addOptional( pInput, 'filenameSuffixStr', '', ...
    @(x)( ischar(x) ) );
try
    pInput.parse( sourceFilenameStr, extensionStr, varargin{:} );
catch mExp
    error( 'forcefilenameextension:invalidInputParameter', ...
        'Error: %s', mExp.message );
end%% Extensive error checking
filenameSuffixStr = pInput.Results.filenameSuffixStr;
clear pInput varargin

%% Make sure first character of extensionStr is .
if( extensionStr(1) ~= '.' )
    extensionStr = [ '.' extensionStr ];
end

[ filePathStr, baseFilenameStr  ] = ...
    fileparts( sourceFilenameStr );
baseFilenameStr = [ baseFilenameStr, filenameSuffixStr ];
newextensionFilenameStr = ...
    fullfile( filePathStr, [ baseFilenameStr, extensionStr ] );

end



