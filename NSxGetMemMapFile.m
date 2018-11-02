function [ NSxFileArray, NSxbasicHeader, NSxchannelHeader, NSxTimeStamps ] = NSxGetMemMapFile( NSxFilename, varargin )
% NSxGetMemMapFile - Create memory map object for reading NSx files
%
%   function [ NSxFileArray, NSxbasicHeader, NSxchannelHeader, NSxTimeStamps ] = NSxGetMemMapFile( NSxFilename, varargin )
%
%   This function creates a memory map object to the data contained in an
%   NSx file, which is a very memory and time efficient method of reading
%   and writing such files.
%
%   Only NEV specification 2.1 and later are supported.
%
%   Inputs
%
%       NSxFilename         -   Name of NSx file to create a memory map to.
%
%       DataPacketIndex     -   Optional parameter pair of
%                               'DataPacketIndex' and a value of which data
%                               packet to create the map to. Starting with
%                               NEV Specification 2.2, NSx files can have
%                               multiple data packets. Nevertheless,
%                               generally only one data packet exists. If
%                               not provided, the default practice is to
%                               use the first data packet.
%
%       SampleOffset        -   Optional parameter pair of 'SampleOffset'
%                               and a scalar offset of the number of time
%                               samples to skip at the start of the file.
%                               If the offset exceeds the size of the
%                               buffer, an error condition will occur. If
%                               note provided, the default practice is to
%                               use an offset of zero.
%
%       WriteableFlag       -   Optional parameter pair of 'WriteableFlag'
%                               and a Boolean to indicate if the memory map
%                               should have the ability to write to the
%                               file (WriteableFlag=true) or not write to
%                               the file (WriteableFlag=false). If not
%                               provided, the default practice is to assume
%                               it is not writable. Note: one can modify
%                               the content of an existing memory map file
%                               but one cannot extend (i.e., append) to a
%                               memory map file. Further, one cannot
%                               trunctate the file.
%
%   Outputs
%
%       NSxFileArray        -   Memory map object to data of NSx file. The
%                               data are accessable via
%                               NSxFileArray.Data.NSxData, but use the name
%                               of NSxFileArray provided in the calling
%                               program.
%
%       NSxbasicHeader      -   The NSx basic header in a structure format
%                               from the NSx file. To get this header, the
%                               subroutine NSxGetHeaders is called and one
%                               should refer to that function for details
%                               on the meaning of the structure entries, or
%                               refer to the NEV Spec.
%
%       NSxchannelHeader    -   The NSx extended headers in an array of
%                               structures format from the NSx file. To get
%                               this header, the subroutine NSxGetHeaders
%                               is called and one should refer to that
%                               function for details on the meaning of the
%                               structure entries, or refer to the NEV
%                               Spec.
%
%       NSxTimeStamps       -   A unsigned 32-bit interger vector of the
%                               time stamps associated with the data in
%                               file. The least significant bit has a value
%                               of 1/30000th of a second.
%
%   $Author: WarrenDJ $, University of Utah
%   $DateTime: 2016/12/02 12:05:25 $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Magic Numbers
DataPacketIndexDefault = 1;
WriteableFlagDefault = false;
SampleOffsetDefault = 0;
SizeOfInt16 = 2;

%% Force outputs to exist
NSxFileArray = [];numel(NSxFileArray);
NSxbasicHeader = [];numel(NSxbasicHeader);
NSxchannelHeader = [];numel(NSxchannelHeader);
NSxTimeStamps = [];numel(NSxTimeStamps);

%% Parse inputs
pInput = inputParser();
pInput.addRequired( 'NSxFilename', ...
    @(x) ( ischar(x) && isvector(x) && exist( x, 'file' ) ) );
pInput.addParamValue( 'DataPacketIndex', DataPacketIndexDefault, ...
    @(x) ( isscalar(x) && isnumeric(x) && ( x > 0 ) ) ); %#ok<NVREPL>
pInput.addParamValue( 'SampleOffset', SampleOffsetDefault, ...
    @(x) ( isscalar(x) && isnumeric(x) && ( x >= 0 ) ) ); %#ok<NVREPL>
pInput.addParamValue( 'WriteableFlag', WriteableFlagDefault, ...
    @(x) ( isscalar(x) && islogical(x) ) ); %#ok<NVREPL>
try
    pInput.parse( NSxFilename, varargin{:} );
catch mExp
    error( 'NSxGetMemMapFile:invalidInputParameter', ...
        'Error: %s', mExp.message );
end
DataPacketIndex = pInput.Results.DataPacketIndex;
SampleOffset = pInput.Results.SampleOffset;
WriteableFlag = pInput.Results.WriteableFlag;
clear pInput

%% Get file headers
fidRead = fopen( NSxFilename, 'rb' );
[ NSxbasicHeader, NSxchannelHeader ] = NSxGetHeaders( fidRead );
fidRead = fclose( fidRead );clear fidRead %#ok<NASGU>
if( isempty( NSxbasicHeader ) || isempty( NSxchannelHeader ) )
    error( 'NSxGetMemMapFile:invalidNSxFile', ...
        'Error: The file %s does not appear to be a NSx file', ...
        NSxFilename );
end

%% Verify data packet exists
if( DataPacketIndex > NSxbasicHeader.NumDataPackets )
    error( 'NSxGetMemMapFile:invalidNSxDataPacket', ...
        'Error: Request to map data packet %d when file only has %d packets', ...
        DataPacketIndex, NSxbasicHeader.NumDataPackets );
end

%% Calculate offset and size of mapped matrix
NSxOffset = (NSxbasicHeader.NumChannels*SizeOfInt16) * SampleOffset;
NSxNumberSamples = NSxbasicHeader.NumTimeSamples( DataPacketIndex ) - ...
    SampleOffset;
if( NSxNumberSamples <= 0 )
    error( 'NSxGetMemMapFile:invalidOffset', ...
        'Error: Request to map data packet sample offset of %d when file only has %d samples', ...
        SampleOffset, NSxbasicHeader.NumTimeSamples( DataPacketIndex ) );
end

%% Create memory map object
NSxFileArray = memmapfile( NSxFilename, ...
    'Offset', (NSxbasicHeader.dataMatrixPtr(DataPacketIndex) + NSxOffset), ...
    'Format', { 'int16', ...
    [ NSxbasicHeader.NumChannels, ...
    NSxNumberSamples ], ...
    'NSxData' }, ...
    'Writable', WriteableFlag, ...
    'Repeat', 1 );

if( nargout >= 4 )
    NSxTimeStamps = ...
        uint32( (NSxbasicHeader.TimeStamp( DataPacketIndex )):...
        (NSxbasicHeader.Period):...
        (NSxbasicHeader.TimeStampEnd( DataPacketIndex )) );
end

return

