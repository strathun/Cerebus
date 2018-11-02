function [ NEVbasicHeader, NEVwaveformHeader, NEVexternalIOHeader, NEVarrayHeader, NEVcommentsText, NEVVideoSync, NEVTrackObj, NEVOtherHeaders ] = NEVGetHeaders( fidRead, varargin )
% NEVGetHeaders - Reads headers of a NEV file into basic, waveform, external I/O and other structures.
%
%	[ NEVbasicHeader, NEVwaveformHeader, NEVexternalIOHeader, NEVarrayHeader, NEVcommentsText, NEVVideoSync, NEVTrackObj, NEVOtherHeaders ] = NEVGetHeaders( fidRead, varargin );
%
%  Reads the header data from a nev file - see the nev spec for interpretation
%
%   Inputs
%
%       fidRead             -   File ID of previously opened NEV file or a
%                               string containing the file to open
%
%                               Although not necessary in all cases, it is
%                               recommended that the file be opened with
%                               'rb' permissions, for read only and binary
%                               format(binary format is default, anyway).
%                           	Although not necessary in all cases, it is
%                               recommended that the file be opened with
%                               'ieee-le' machine format, for little-end
%                               byte ordering that is typically the case
%                               for Windows machines. Although not
%                               necessary in all cases, it is recommended
%                               that the file be opened with 'windows-1252'
%                               encoding so that the extended characters in
%                               text strings are properly interpreted. This
%                               later item is very important for Unix
%                               machines.
%
%       fileSource          -   Parameter pair to indicated neural
%                               recording source of file. Use the parameter
%                               'fileSource' and the parameter value of
%                               'Cerebus' or 'Trellis'. If none provided,
%                               then neither is assumed. Note: This will
%                               override any file source parameter that is
%                               extracted from the application name in the
%                               header.
%
%   Outputs (See most recent NEV specificiation for details
%
%       NEVbasicHeader      -	Header data implemented as structure scalar
%                               detailing information common to entire NEV
%                               file.
%
%       NEVwaveformHeader   -	Header data implemented as structure vector
%                               detailing information unique to particular
%                               electrode channels, including input channel
%                               information, label information, and
%                               filtering information.
%
%       NEVexternalIOHeader -	Header data implemented as structure scalar
%                               detailing information unique to external
%                               experimental data, including input channel
%                               information and label information.
%
%       NEVarrayHeader      -	Header data implemented as structure scalar
%                               detailing information unique to particular
%                               array used, including array name and map
%                               file used in creation of data.
%
%       NEVcommentsText     -	Header data implemented as structure vector
%                               with the text from any comments in the
%                               headers.
%
%       NEVVideoSync        -	Header data implemented as structure scalar
%                               detailing information related to video
%                               synchronization events in the headers.
%
%       NEVTrackObj         -	Header data implemented as structure scalar
%                               detailing information related to motion
%                               tracking events in the headers.
%
%       NEVOtherHeaders     -   Headers that don't fall into any other
%                               defined structure.
%
%   $author$, University of Utah
%   $date$
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Magic Numbers
maximumChannelNumberCerebus = 255;  % 06-Apr-2015: largest neural channel number in BR specification
maximumChannelNumberTrellis = 10240;  % 06-Apr-2015: largest neural channel number in Ripple specification
basicHeaderPacketSize = 8 + 2 + 2 + 4 + 4 + 4 + 4 + 8*2 + 32 + 256 + 4;
extentedHeaderPacketSize = 32;
sizeOfInt8 = 1;numel( sizeOfInt8 ); %  bytes
sizeOfInt16 = 2;numel( sizeOfInt16 ); %  bytes
sizeOfInt32 = 4;numel( sizeOfInt32 ); % % bytes
sizeOfInt64 = 8;numel( sizeOfInt64 ); % % bytes
sizeOfSingle = 4;numel( sizeOfSingle ); % bytes
sizeOfDouble = 8;numel( sizeOfDouble ); % bytes

%% Assure some output
NEVbasicHeader = [];
NEVwaveformHeader= [];
NEVexternalIOHeader= [];
NEVarrayHeader = [];
NEVcommentsText = [];
NEVVideoSync = [];
NEVTrackObj = [];
NEVOtherHeaders = [];

%% extensive checking of inputs
pInput = inputParser();
pInput.addRequired( 'fidRead', ...
    @(x)( ( isscalar( x ) && ( x ~= -1 ) ) || ( ischar( x ) && exist( x, 'file' ) ) ) );
pInput.addParamValue( 'fileSource', [], ...
    @(x)( ischar(x) ) ); %#ok<NVREPL>
try
    pInput.parse( fidRead, varargin{:} );
catch mExp
    error( 'NEVGetHeaders:invalidInputParameter', ...
        'Error: %s', mExp.message );
end%% Extensive error checking
fileSource = pInput.Results.fileSource;
clear pInput

%%	Open file if string
if( ischar( fidRead ) && exist( fidRead, 'file' ) )
    fidReadIn = fidRead;
    fidRead = fopen( fidReadIn, 'rb', 'ieee-le', 'windows-1252' );
    if( fidRead < 0 )
        warning( 'NSxGetHeaders:FileNameError', ...
            'Unable to open file\n' );
        return
    end
else
    fidReadIn = [];
end

%%	Get details of file
[ filename, permissions, machineformat, encoding ] = fopen(fidRead); 
if( isempty( filename ) )
    warning( 'NEVGetHeaders:FileNameError', ...
        'Unable to get filename of open file\n' );
    if(~isempty(fidReadIn));fclose(fidRead);end
    return;
end;
if( ~isempty( strfind( lower( permissions ), 't' ) ) )
    warning( 'NEVGetHeaders:FilePermissionsError', ...
        'File %s opened in text mode, may result in problems interpreting strings\n', filename );
end;
if( ~isempty( strfind( lower( permissions ), 'w' ) ) )
    warning( 'NEVGetHeaders:FilePermissionsError', ...
        'File %s opened for writing, may result in problems\n', filename );
end;
if( ~isempty( strfind( lower( permissions ), 'a' ) ) )
    warning( 'NEVGetHeaders:FilePermissionsError', ...
        'File %s opened for appending, may result in problems\n', filename );
end;
if( ~strcmpi( machineformat, 'ieee-le' ) )
    warning( 'NEVGetHeaders:FileFormatError', ...
        'File %s not opened in little-endian mode, may result in problems interpreting numbers\n', filename );
end;
if( ~strcmpi( encoding, 'windows-1252' ) )
    warning( 'NEVGetHeaders:FileEncodingError', ...
        'File %s not opened with correct text encoding, may result in problems interpreting strings\n', filename );
end;
clear permissions machineformat encoding

%%	Get details of file
fileDir = dir( filename ); 
if( isempty( fileDir ) )
    warning( 'NEVGetHeaders:FileNameError', ...
        'Unable to get details of open file\n' );
    if(~isempty(fidReadIn));fclose(fidRead);end
    return;
end;

%%	Position file to beginning
if( fseek( fidRead, 0, 'bof' ) == -1 )
    warning( 'NEVGetHeaders:FilePositioningError', ...
        'Invalid file positioning with message %s\n', ferror(fidRead,'clear') );
    if(~isempty(fidReadIn));fclose(fidRead);end;
    return;
end;

%% Precreate basic header
NEVbasicHeader = struct( ...
    'id', blanks(8), ...
    'filespecMajor', 0, ...
    'filespecMinor', 0, ...
    'fileSource', '', ...
    'fileformat', 0, ...
    'dataptr', 0, ...
    'datasize', 0, ...
    'TimeRes', 0, ...
    'SampleRes', 0, ...
    'FileTime', zeros( 8, 1), ...
    'AppName', blanks(32), ...
    'Comment', blanks(256), ...
    'StartTimeNIPTimeStamp', uint32( 0 ), ...
    'NumHeaders', 0, ...
    'maximumChannelNumber', 0, ...
    'filespecDouble', 0, ...
    'fileformat16Bit', 0, ...
    'NumPackets', 0, ...
    'SerialDateNumber', 0, ...
    'NumBytes', 0, ...
    'Filename', blanks(1) ...
    );

%%	Read Basic Headers
ncountTest = 8;[ temp, ncount ] = fread( fidRead, [1,ncountTest], 'char*1' );
if( ncount ~= ncountTest );warning( 'NEVGetHeaders:readCountError', 'Unable to read correct number of elements' );if(~isempty(fidReadIn));fclose(fidRead);end;return;
else
    NEVbasicHeader.id = char( temp );
end

ncountTest = 1;[ NEVbasicHeader.filespecMajor, ncount ] = fread( fidRead, [1,ncountTest], 'uint8' );
if( ncount ~= ncountTest );warning( 'NEVGetHeaders:readCountError', 'Unable to read correct number of elements' );if(~isempty(fidReadIn));fclose(fidRead);end;return;end

ncountTest = 1;[ NEVbasicHeader.filespecMinor, ncount ] = fread( fidRead, [1,ncountTest], 'uint8' );
if( ncount ~= ncountTest );warning( 'NEVGetHeaders:readCountError', 'Unable to read correct number of elements' );if(~isempty(fidReadIn));fclose(fidRead);end;return;end

ncountTest = 1;[ NEVbasicHeader.fileformat, ncount ] = fread( fidRead, [1,ncountTest], 'uint16' );
if( ncount ~= ncountTest );warning( 'NEVGetHeaders:readCountError', 'Unable to read correct number of elements' );if(~isempty(fidReadIn));fclose(fidRead);end;return;end

ncountTest = 1;[ NEVbasicHeader.dataptr, ncount ] = fread( fidRead, [1,ncountTest], 'uint32' );
if( ncount ~= ncountTest );warning( 'NEVGetHeaders:readCountError', 'Unable to read correct number of elements' );if(~isempty(fidReadIn));fclose(fidRead);end;return;end

ncountTest = 1;[ NEVbasicHeader.datasize, ncount ] = fread( fidRead, [1,ncountTest], 'uint32' );
if( ncount ~= ncountTest );warning( 'NEVGetHeaders:readCountError', 'Unable to read correct number of elements' );if(~isempty(fidReadIn));fclose(fidRead);end;return;end

ncountTest = 1;[ NEVbasicHeader.TimeRes, ncount ] = fread( fidRead, [1,ncountTest], 'uint32' );
if( ncount ~= ncountTest );warning( 'NEVGetHeaders:readCountError', 'Unable to read correct number of elements' );if(~isempty(fidReadIn));fclose(fidRead);end;return;end

ncountTest = 1;[ NEVbasicHeader.SampleRes, ncount ] = fread( fidRead, [1,ncountTest], 'uint32' );
if( ncount ~= ncountTest );warning( 'NEVGetHeaders:readCountError', 'Unable to read correct number of elements' );if(~isempty(fidReadIn));fclose(fidRead);end;return;end

ncountTest = 8;[ NEVbasicHeader.FileTime, ncount ] = fread( fidRead, [1,ncountTest], 'uint16' );
if( ncount ~= ncountTest );warning( 'NEVGetHeaders:readCountError', 'Unable to read correct number of elements' );if(~isempty(fidReadIn));fclose(fidRead);end;return;end

ncountTest = 32;[ temp, ncount ] = fread( fidRead, [1,ncountTest], 'char*1' );
if( ncount ~= ncountTest );warning( 'NEVGetHeaders:readCountError', 'Unable to read correct number of elements' );if(~isempty(fidReadIn));fclose(fidRead);end;return;
else
    NEVbasicHeader.AppName = char( temp );
end

ncountTest = 256;[ temp, ncount ] = fread( fidRead, [1,ncountTest], 'char*1' );
if( ncount ~= ncountTest );warning( 'NEVGetHeaders:readCountError', 'Unable to read correct number of elements' );if(~isempty(fidReadIn));fclose(fidRead);end;return;
else
    if( ...
            ( temp( 252 ) == 0 ) && ...
            any( temp( 253:256 ) ~= 0 ) ...
            )
        % Ripple Trellis Version 1.4.1 and later
        NEVbasicHeader.Comment = char( [ temp(1:252) zeros(1,4,'uint8') ] );
        NEVbasicHeader.StartTimeNIPTimeStamp = typecast( uint8( temp( 253:256 ) ), 'uint32' );
    else
        NEVbasicHeader.Comment = char( temp(1:end) );
        NEVbasicHeader.StartTimeNIPTimeStamp = uint32( 0 );
    end
end

ncountTest = 1;[ NEVbasicHeader.NumHeaders, ncount ] = fread( fidRead, [1,ncountTest], 'uint32' );
if( ncount ~= ncountTest );warning( 'NEVGetHeaders:readCountError', 'Unable to read correct number of elements' );if(~isempty(fidReadIn));fclose(fidRead);end;return;end

NEVbasicHeader.filespecDouble = fix( NEVbasicHeader.filespecMajor ) + ...
    fix( NEVbasicHeader.filespecMinor ) / 10;
if( NEVbasicHeader.filespecDouble < 1 ) % handle really old NSAS files
    NEVbasicHeader.filespecDouble = 10 * ...
        NEVbasicHeader.filespecDouble;
end
NEVbasicHeader.fileformat16Bit = mod( NEVbasicHeader.fileformat, 2 );
NEVbasicHeader.NumPackets = ( fileDir.bytes - NEVbasicHeader.dataptr ) ...
    / ( NEVbasicHeader.datasize );
if( mod( NEVbasicHeader.NumPackets, 1 ) ~= 0 )
    warning( 'NEVGetHeaders:PacketCountError', ...
        'Non-integer number of data packet at %f\n', NEVbasicHeader.NumPackets );
     NEVbasicHeader.NumPackets = floor( NEVbasicHeader.NumPackets );
end;
dateVector = NEVbasicHeader.FileTime( [1 2 4 5 6 7] );
if( dateVector(1) < 50 )
    dateVector(1) = dateVector(1) + 2000;
else
    if( dateVector(1) < 100 )
        dateVector(1) = dateVector(1) + 1900;
    end
end
NEVbasicHeader.SerialDateNumber = datenum( dateVector );
NEVbasicHeader.Filename = filename;
NEVbasicHeader.NumBytes = fileDir.bytes;
clear filename fileDir dateVector

%% Figure out if is Cerebus or Trellis file
if( strcmpi( fileSource, 'Trellis' ) )
    maximumChannelNumber = maximumChannelNumberTrellis;
    NEVbasicHeader.fileSource = 'Trellis';
elseif( strcmpi( fileSource, 'Cerebus' ) )
    maximumChannelNumber = maximumChannelNumberCerebus;
    NEVbasicHeader.fileSource = 'Cerebus';
elseif( ~isempty( strfind( lower( NEVbasicHeader.AppName ), lower( 'trellis' ) ) ) )
    maximumChannelNumber = maximumChannelNumberTrellis;
    NEVbasicHeader.fileSource = 'Trellis';
elseif( ~isempty( strfind( lower( NEVbasicHeader.AppName ), lower( 'file dialog' ) ) ) )
    maximumChannelNumber = maximumChannelNumberCerebus;
    NEVbasicHeader.fileSource = 'Cerebus';
else
    maximumChannelNumber = maximumChannelNumberTrellis;
    NEVbasicHeader.fileSource = 'Trellis';
    warning( 'NEVGetHeaders:UnclearNEVfileSource', ...
        'It is unclear if this is a Cerebus or Trellis NEV file, assuming %s\n', ...
        NEVbasicHeader.fileSource );
end
NEVbasicHeader.maximumChannelNumber = maximumChannelNumber;

%% Precreate waveform header
NEVwaveformHeader = struct( ...
    'id', num2cell( zeros( maximumChannelNumber, 1) ), ...
    'existNEUEVWAV', num2cell( zeros( maximumChannelNumber, 1) ), ...
    'existNEUEVLBL', num2cell( zeros( maximumChannelNumber, 1) ), ...
    'existNEUEVFLT', num2cell( zeros( maximumChannelNumber, 1) ), ...
    'pinch', num2cell( zeros( maximumChannelNumber, 1) ), ...
    'pinnum', num2cell( zeros( maximumChannelNumber, 1) ), ...
    'sf', num2cell( zeros( maximumChannelNumber, 1) ), ...
    'energythreshold', num2cell( zeros( maximumChannelNumber, 1) ), ...
    'highthreshold', num2cell( zeros( maximumChannelNumber, 1) ), ...
    'lowthreshold', num2cell( zeros( maximumChannelNumber, 1) ), ...
    'numberunits', num2cell( zeros( maximumChannelNumber, 1) ), ...
    'numberbytes', num2cell( zeros( maximumChannelNumber, 1) ), ...
    'dummy', num2cell( zeros( maximumChannelNumber, 10 ), 2 ), ...
    'label', num2cell( zeros( maximumChannelNumber, 17 ), 2 ), ...
    'dummy1', num2cell( zeros( maximumChannelNumber, 6 ), 2 ), ...
    'highfreqcorner', num2cell( zeros( maximumChannelNumber, 1) ), ...
    'highfreqorder', num2cell( zeros( maximumChannelNumber, 1) ), ...
    'highfreqtype', num2cell( zeros( maximumChannelNumber, 1) ), ...
    'lowfreqcorner', num2cell( zeros( maximumChannelNumber, 1) ), ...
    'lowfreqorder', num2cell( zeros( maximumChannelNumber, 1) ), ...
    'lowfreqtype', num2cell( zeros( maximumChannelNumber, 1) ), ...
    'StimAmpsf', num2cell( zeros( maximumChannelNumber, 1) ), ...
    'dummy2', num2cell( zeros( maximumChannelNumber, 2 ), 2 ) ...
    );

%% Precreate external header
NEVexternalIOHeader.existNSASEXEV = 0;
NEVexternalIOHeader.existSTIMINFO = 0;
NEVexternalIOHeader.existDIGLABELSerial = 0;
NEVexternalIOHeader.existDIGLABELParallel = 0;
NEVexternalIOHeader.PeriodicRes = 0;
NEVexternalIOHeader.DigitalConfig = 0;
for m=1:5;
    NEVexternalIOHeader.AnalogConfig(m).Config = 0;
    NEVexternalIOHeader.AnalogConfig(m).Level = 0;
end;clear m
NEVexternalIOHeader.dummy = zeros(6,1);
NEVexternalIOHeader.SerialLabel = blanks( 16 );
NEVexternalIOHeader.SerialDummy = zeros(7,1);
NEVexternalIOHeader.ParallelLabel = blanks( 16 );
NEVexternalIOHeader.ParallelDummy = zeros(7,1);

%% Precreate array header
NEVarrayHeader.existARRAYNME = 0;
NEVarrayHeader.existMAPFILE = 0;
NEVarrayHeader.arrayName = blanks(24);
NEVarrayHeader.mapfile = blanks(24);

%% Precreate comments header
NEVcommentsText = struct( ...
    'existCCOMMENT', num2cell( zeros( 1, 1 ) ), ...
    'existECOMMENT', num2cell( zeros( 1, 1 ) ), ...
    'Comment', num2cell( blanks(24), 2 ) ...
    );

%% Precreate video sync header
NEVVideoSync.existVIDEOSYN = 0;
NEVVideoSync.videoSourceID = -1;
NEVVideoSync.videoSourceName = blanks(16);
NEVVideoSync.videoFrameRate = -1;
NEVVideoSync.videoDummy = zeros(2,1);

%% Precreate motion tracking header
NEVTrackObj.existTRACKOBJ = 0;
NEVTrackObj.TrackingType =0;
NEVTrackObj.TrackingID = 0;
NEVTrackObj.TrackingNumberPoints = 0;
NEVTrackObj.TrackingVideoSource = blanks(16);
NEVTrackObj.trackingDummy = zeros(2,1);

%% Precreate other header
NEVOtherHeaders = struct( ...
    'existOther', num2cell( zeros( 1, 1 ) ), ...
    'PackedIDString', num2cell( blanks(8), 2 ), ...
    'Payload', num2cell( blanks(24), 2 ) ...
    );

%%	Read Extended Headers
maximumChannelNumberObserved = 0;
for n=1:NEVbasicHeader.NumHeaders;
    
    % Check that file position is correct
    if( ftell( fidRead ) ~= basicHeaderPacketSize + (n-1)*extentedHeaderPacketSize )
        warning( 'NEVGetHeaders:badFilePosition', ...
            'Invalid file positioning\n' );
        if(~isempty(fidReadIn));fclose(fidRead);end;
        return
    end
    
    % Read extended header type string
    ncountTest = 8;[ temp, ncount ] = fread( fidRead, [1,ncountTest], 'char*1' );
    if( ncount ~= ncountTest );warning( 'NEVGetHeaders:readCountError', 'Unable to read correct number of elements' );if(~isempty(fidReadIn));fclose(fidRead);end;return;
    else
        HeaderName = upper( char( temp ) );
    end
    
    % handle based on type string
    switch HeaderName
        case { 'NEUEVWAV', 'NEUEVLBL', 'NEUEVFLT' } % Neural event waveform data
            
            ncountTest = 1;[ id, ncount ] = fread( fidRead, [1,ncountTest], 'uint16' );
            if( ncount ~= ncountTest );warning( 'NEVGetHeaders:readCountError', 'Unable to read correct number of elements' );if(~isempty(fidReadIn));fclose(fidRead);end;return;end
            if( ( id < 1 ) || ( id > NEVbasicHeader.maximumChannelNumber ) )
                warning( 'NEVGetHeaders:badChannelID', ...
                    'Invalid packet ID value (%f) found\n', id );
            end
            maximumChannelNumberObserved = ...
                max( maximumChannelNumberObserved, id );
            NEVwaveformHeader(id).id = id;
            
            switch HeaderName
                case 'NEUEVWAV' % Basic neural event waveform data
                    
                    NEVwaveformHeader(id).existNEUEVWAV = 1;
                    
                    ncountTest = 1;[ NEVwaveformHeader(id).pinch, ncount ] = fread( fidRead, [1,ncountTest], 'uint8' );
                    if( ncount ~= ncountTest );warning( 'NEVGetHeaders:readCountError', 'Unable to read correct number of elements' );if(~isempty(fidReadIn));fclose(fidRead);end;return;end
                    
                    ncountTest = 1;[ NEVwaveformHeader(id).pinnum, ncount ] = fread( fidRead, [1,ncountTest], 'uint8' );
                    if( ncount ~= ncountTest );warning( 'NEVGetHeaders:readCountError', 'Unable to read correct number of elements' );if(~isempty(fidReadIn));fclose(fidRead);end;return;end
                    
                    ncountTest = 1;[ NEVwaveformHeader(id).sf, ncount ] = fread( fidRead, [1,ncountTest], 'uint16' );
                    if( ncount ~= ncountTest );warning( 'NEVGetHeaders:readCountError', 'Unable to read correct number of elements' );if(~isempty(fidReadIn));fclose(fidRead);end;return;end
                    
                    ncountTest = 1;[ NEVwaveformHeader(id).energythreshold, ncount ] = fread( fidRead, [1,ncountTest], 'uint16' );
                    if( ncount ~= ncountTest );warning( 'NEVGetHeaders:readCountError', 'Unable to read correct number of elements' );if(~isempty(fidReadIn));fclose(fidRead);end;return;end
                    
                    ncountTest = 1;[ NEVwaveformHeader(id).highthreshold, ncount ] = fread( fidRead, [1,ncountTest], 'int16' );
                    if( ncount ~= ncountTest );warning( 'NEVGetHeaders:readCountError', 'Unable to read correct number of elements' );if(~isempty(fidReadIn));fclose(fidRead);end;return;end
                    
                    ncountTest = 1;[ NEVwaveformHeader(id).lowthreshold, ncount ] = fread( fidRead, [1,ncountTest], 'int16' );
                    if( ncount ~= ncountTest );warning( 'NEVGetHeaders:readCountError', 'Unable to read correct number of elements' );if(~isempty(fidReadIn));fclose(fidRead);end;return;end
                    
                    ncountTest = 1;[ NEVwaveformHeader(id).numberunits, ncount ] = fread( fidRead, [1,ncountTest], 'uint8' );
                    if( ncount ~= ncountTest );warning( 'NEVGetHeaders:readCountError', 'Unable to read correct number of elements' );if(~isempty(fidReadIn));fclose(fidRead);end;return;end
                    
                    ncountTest = 1;[ NEVwaveformHeader(id).numberbytes, ncount ] = fread( fidRead, [1,ncountTest], 'uint8' );
                    if( ncount ~= ncountTest );warning( 'NEVGetHeaders:readCountError', 'Unable to read correct number of elements' );if(~isempty(fidReadIn));fclose(fidRead);end;return;end
                    if( NEVbasicHeader.fileformat16Bit == 1 )
                        NEVwaveformHeader(id).numberbytes = 2;
                    else
                        if( NEVwaveformHeader(id).numberbytes == 0 )
                            NEVwaveformHeader(id).numberbytes = 1;
                        end
                    end
                    
                    if( strcmpi( NEVbasicHeader.fileSource, 'Trellis' ) )
                        ncountTest = 1;[ NEVwaveformHeader(id).StimAmpsf, ncount ] = fread( fidRead, [1,ncountTest], 'float32' );
                        if( ncount ~= ncountTest );warning( 'NEVGetHeaders:readCountError', 'Unable to read correct number of elements' );if(~isempty(fidReadIn));fclose(fidRead);end;return;end
                        ncountTest = 6;[ NEVwaveformHeader(id).dummy, ncount ] = fread( fidRead, [1,ncountTest], 'uint8' );
                        if( ncount ~= ncountTest );warning( 'NEVGetHeaders:readCountError', 'Unable to read correct number of elements' );if(~isempty(fidReadIn));fclose(fidRead);end;return;end
                    else
                        NEVwaveformHeader(id).StimAmpsf = single( 0.0 );
                        ncountTest = 10;[ NEVwaveformHeader(id).dummy, ncount ] = fread( fidRead, [1,ncountTest], 'uint8' );
                        if( ncount ~= ncountTest );warning( 'NEVGetHeaders:readCountError', 'Unable to read correct number of elements' );if(~isempty(fidReadIn));fclose(fidRead);end;return;end
                    end
                    
                case 'NEUEVLBL' % Neural event channel label
                    
                    NEVwaveformHeader(id).existNEUEVLBL = 1;
                    
                    ncountTest = 16;[ temp, ncount ] = fread( fidRead, [1,ncountTest], 'char*1' );
                    if( ncount ~= ncountTest );warning( 'NEVGetHeaders:readCountError', 'Unable to read correct number of elements' );if(~isempty(fidReadIn));fclose(fidRead);end;return;
                    else
                        NEVwaveformHeader(id).label = char( temp );
                    end
                    
                    ncountTest = 6;[ NEVwaveformHeader(id).dummy1, ncount ] = fread( fidRead, [1,ncountTest], 'uint8' );
                    if( ncount ~= ncountTest );warning( 'NEVGetHeaders:readCountError', 'Unable to read correct number of elements' );if(~isempty(fidReadIn));fclose(fidRead);end;return;end

                case 'NEUEVFLT' % Neural event channel filtering

                    NEVwaveformHeader(id).existNEUEVFLT = 1;

                    ncountTest = 1;[ NEVwaveformHeader(id).highfreqcorner, ncount ] = fread( fidRead, [1,ncountTest], 'uint32' );
                    if( ncount ~= ncountTest );warning( 'NEVGetHeaders:readCountError', 'Unable to read correct number of elements' );if(~isempty(fidReadIn));fclose(fidRead);end;return;end
                    
                    ncountTest = 1;[ NEVwaveformHeader(id).highfreqorder, ncount ] = fread( fidRead, [1,ncountTest], 'uint32' );
                    if( ncount ~= ncountTest );warning( 'NEVGetHeaders:readCountError', 'Unable to read correct number of elements' );if(~isempty(fidReadIn));fclose(fidRead);end;return;end
                    
                    ncountTest = 1;[ NEVwaveformHeader(id).highfreqtype, ncount ] = fread( fidRead, [1,ncountTest], 'uint16' );
                    if( ncount ~= ncountTest );warning( 'NEVGetHeaders:readCountError', 'Unable to read correct number of elements' );if(~isempty(fidReadIn));fclose(fidRead);end;return;end
                    
                    ncountTest = 1;[ NEVwaveformHeader(id).lowfreqcorner, ncount ] = fread( fidRead, [1,ncountTest], 'uint32' );
                    if( ncount ~= ncountTest );warning( 'NEVGetHeaders:readCountError', 'Unable to read correct number of elements' );if(~isempty(fidReadIn));fclose(fidRead);end;return;end
                    
                    ncountTest = 1;[ NEVwaveformHeader(id).lowfreqorder, ncount ] = fread( fidRead, [1,ncountTest], 'uint32' );
                    if( ncount ~= ncountTest );warning( 'NEVGetHeaders:readCountError', 'Unable to read correct number of elements' );if(~isempty(fidReadIn));fclose(fidRead);end;return;end
                    
                    ncountTest = 1;[ NEVwaveformHeader(id).lowfreqtype, ncount ] = fread( fidRead, [1,ncountTest], 'uint16' );
                    if( ncount ~= ncountTest );warning( 'NEVGetHeaders:readCountError', 'Unable to read correct number of elements' );if(~isempty(fidReadIn));fclose(fidRead);end;return;end
                    
                    ncountTest = 2;[ NEVwaveformHeader(id).dummy2, ncount ] = fread( fidRead, [1,ncountTest], 'uint8' );
                    if( ncount ~= ncountTest );warning( 'NEVGetHeaders:readCountError', 'Unable to read correct number of elements' );if(~isempty(fidReadIn));fclose(fidRead);end;return;end

            end;clear id
            
        case { 'STIMINFO', 'NSASEXEV' } % External I/O info
            
            switch HeaderName % External I/O info, note goes back to NSAS
                case 'STIMINFO'
                    NEVexternalIOHeader.existSTIMINFO = 1;
                case 'NSASEXEV' % External I/O info
                    NEVexternalIOHeader.existNSASEXEV = 1;
            end

            ncountTest = 1;[ NEVexternalIOHeader.PeriodicRes, ncount ] = fread( fidRead, [1,ncountTest], 'uint16' );
            if( ncount ~= ncountTest );warning( 'NEVGetHeaders:readCountError', 'Unable to read correct number of elements' );if(~isempty(fidReadIn));fclose(fidRead);end;return;end
            
            ncountTest = 1;[ NEVexternalIOHeader.DigitalConfig, ncount ] = fread( fidRead, [1,ncountTest], 'uint8' );
            if( ncount ~= ncountTest );warning( 'NEVGetHeaders:readCountError', 'Unable to read correct number of elements' );if(~isempty(fidReadIn));fclose(fidRead);end;return;end
            
            for m=1:5;
                
                ncountTest = 1;[ NEVexternalIOHeader.AnalogConfig(m).Config, ncount ] = fread( fidRead, [1,ncountTest], 'uint8' );
                if( ncount ~= ncountTest );warning( 'NEVGetHeaders:readCountError', 'Unable to read correct number of elements' );if(~isempty(fidReadIn));fclose(fidRead);end;return;end
                
                ncountTest = 1;[ NEVexternalIOHeader.AnalogConfig(m).Level, ncount ] = fread( fidRead, [1,ncountTest], 'int16' );
                if( ncount ~= ncountTest );warning( 'NEVGetHeaders:readCountError', 'Unable to read correct number of elements' );if(~isempty(fidReadIn));fclose(fidRead);end;return;end
                
            end;clear m
            
            ncountTest = 6;[ NEVexternalIOHeader.dummy, ncount ] = fread( fidRead, [1,ncountTest], 'uint8' );
            if( ncount ~= ncountTest );warning( 'NEVGetHeaders:readCountError', 'Unable to read correct number of elements' );if(~isempty(fidReadIn));fclose(fidRead);end;return;end
            
        case 'DIGLABEL' % External I/O label
            
            ncountTest = 16;[ temp, ncount ] = fread( fidRead, [1,ncountTest], 'char*1' );
            if( ncount ~= ncountTest );warning( 'NEVGetHeaders:readCountError', 'Unable to read correct number of elements' );if(~isempty(fidReadIn));fclose(fidRead);end;return;
            else
                tempLabel = char( temp );
            end
            
            ncountTest = 1;[ tempMode, ncount ] = fread( fidRead, [1,ncountTest], 'uint8' );
            if( ncount ~= ncountTest );warning( 'NEVGetHeaders:readCountError', 'Unable to read correct number of elements' );if(~isempty(fidReadIn));fclose(fidRead);end;return;end

            ncountTest = 7;[ tempDummy, ncount ] = fread( fidRead, [1,ncountTest], 'uint8' );
            if( ncount ~= ncountTest );warning( 'NEVGetHeaders:readCountError', 'Unable to read correct number of elements' );if(~isempty(fidReadIn));fclose(fidRead);end;return;end

            switch tempMode
                case 0
                    NEVexternalIOHeader.existDIGLABELSerial = 1;
                    NEVexternalIOHeader.SerialLabel = tempLabel;
                    NEVexternalIOHeader.SerialDummy = tempDummy;
                case 1
                    NEVexternalIOHeader.existDIGLABELParallel = 1;
                    NEVexternalIOHeader.ParallelLabel = tempLabel;
                    NEVexternalIOHeader.ParallelDummy = tempDummy;
                otherwise
            end;clear tempLabel tempMode tempDummy % Digital mode
        
        case 'ARRAYNME' % name of array

            NEVarrayHeader.existARRAYNME = 1;
            
            ncountTest = 24;[ temp, ncount ] = fread( fidRead, [1,ncountTest], 'char*1' );
            if( ncount ~= ncountTest );warning( 'NEVGetHeaders:readCountError', 'Unable to read correct number of elements' );if(~isempty(fidReadIn));fclose(fidRead);end;return;
            else
                NEVarrayHeader.arrayName = char( temp );
            end
            
        case 'MAPFILE' % name of array

            NEVarrayHeader.existMAPFILE = 1;
            ncountTest = 24;[ temp, ncount ] = fread( fidRead, [1,ncountTest], 'char*1' );
            if( ncount ~= ncountTest );warning( 'NEVGetHeaders:readCountError', 'Unable to read correct number of elements' );if(~isempty(fidReadIn));fclose(fidRead);end;return;
            else
                NEVarrayHeader.mapfile = char( temp );
            end
            
        case { 'ECOMMENT', 'CCOMMENT' } % name of array
            
            if( ...
                    ( NEVcommentsText(end).existECOMMENT ~=0 ) || ...
                    ( NEVcommentsText(end).existCCOMMENT ~=0 ) ...
                    )
                NEVcommentsText(end+1,1).existECOMMENT = 0; %#ok<AGROW>
                NEVcommentsText(end,1).existCCOMMENT = 0;
            end
            
            switch HeaderName
                case 'ECOMMENT'
                    NEVcommentsText(end,1).existECOMMENT = 1; 
                case 'CCOMMENT'
                    NEVcommentsText(end,1).existCCOMMENT = 1; 
            end;
            
            ncountTest = 24;[ temp, ncount ] = fread( fidRead, [1,ncountTest], 'char*1' );
            if( ncount ~= ncountTest );warning( 'NEVGetHeaders:readCountError', 'Unable to read correct number of elements' );if(~isempty(fidReadIn));fclose(fidRead);end;return;
            else
                NEVcommentsText(end,1).Comment = char( temp ); 
            end
            
        case 'VIDEOSYN' % Video sync
            
            ncountTest = 1;[ tempSourceID, ncount ] = fread( fidRead, [1,ncountTest], 'uint16' );
            if( ncount ~= ncountTest );warning( 'NEVGetHeaders:readCountError', 'Unable to read correct number of elements' );if(~isempty(fidReadIn));fclose(fidRead);end;return;end
            
            ncountTest = 16;[ temp, ncount ] = fread( fidRead, [1,ncountTest], 'char*1' );
            if( ncount ~= ncountTest );warning( 'NEVGetHeaders:readCountError', 'Unable to read correct number of elements' );if(~isempty(fidReadIn));fclose(fidRead);end;return;
            else
                tempLabel = char( temp );
            end
            
            ncountTest = 1;[ tempFrameRate, ncount ] = fread( fidRead, [1,ncountTest], 'float32' );
            if( ncount ~= ncountTest );warning( 'NEVGetHeaders:readCountError', 'Unable to read correct number of elements' );if(~isempty(fidReadIn));fclose(fidRead);end;return;end
            
            ncountTest = 2;[ tempDummy, ncount ] = fread( fidRead, [1,ncountTest], 'uint8' );
            if( ncount ~= ncountTest );warning( 'NEVGetHeaders:readCountError', 'Unable to read correct number of elements' );if(~isempty(fidReadIn));fclose(fidRead);end;return;end

            if( tempSourceID >= 0 )
                NEVVideoSync.existVIDEOSYN = 1;
                NEVVideoSync.videoSourceID = tempSourceID;
                NEVVideoSync.videoSourceName = tempLabel;
                NEVVideoSync.videoFrameRate = tempFrameRate;
                NEVVideoSync.videoDummy = tempDummy;
            end;clear tempSourceID tempLabel tempFrameRate tempDummy
        
        case 'TRACKOBJ' % Tracking
            
            ncountTest = 1;[ tempTrackingType, ncount ] = fread( fidRead, [1,ncountTest], 'uint16' );
            if( ncount ~= ncountTest );warning( 'NEVGetHeaders:readCountError', 'Unable to read correct number of elements' );if(~isempty(fidReadIn));fclose(fidRead);end;return;end
            
            ncountTest = 1;[ tempTrackingID, ncount ] = fread( fidRead, [1,ncountTest], 'uint16' );
            if( ncount ~= ncountTest );warning( 'NEVGetHeaders:readCountError', 'Unable to read correct number of elements' );if(~isempty(fidReadIn));fclose(fidRead);end;return;end
            
            ncountTest = 1;[ tempTrackingNumberPoints, ncount ] = fread( fidRead, [1,ncountTest], 'uint16' );
            if( ncount ~= ncountTest );warning( 'NEVGetHeaders:readCountError', 'Unable to read correct number of elements' );if(~isempty(fidReadIn));fclose(fidRead);end;return;end
            
            ncountTest = 16;[ temp, ncount ] = fread( fidRead, [1,ncountTest], 'char*1' );
            if( ncount ~= ncountTest );warning( 'NEVGetHeaders:readCountError', 'Unable to read correct number of elements' );if(~isempty(fidReadIn));fclose(fidRead);end;return;
            else
                tempLabel = char( temp );
            end
            
            ncountTest = 2;[ tempDummy, ncount ] = fread( fidRead, [1,ncountTest], 'uint8' );
            if( ncount ~= ncountTest );warning( 'NEVGetHeaders:readCountError', 'Unable to read correct number of elements' );if(~isempty(fidReadIn));fclose(fidRead);end;return;end

            if( tempTrackingType >= 0 )
                NEVTrackObj.TrackingType = tempTrackingType;
                NEVTrackObj.TrackingID = tempTrackingID;
                NEVTrackObj.TrackingNumberPoints = tempTrackingNumberPoints;
                NEVTrackObj.TrackingVideoSource = tempLabel;
                NEVTrackObj.trackingDummy = tempDummy;
            end;clear tempTrackingType tempTrackingID tempTrackingNumberPoints tempLabel tempDummy
        
        otherwise
            
            fprintf( 'Unhandled header packet type with code %s\n', HeaderName );
            NEVOtherHeaders.PackedIDString = HeaderName;
            ncountTest = 24;[ temp, ncount ] = fread( fidRead, [1,ncountTest], 'char*1' ); 
            if( ncount ~= ncountTest );warning( 'NEVGetHeaders:readCountError', 'Unable to read correct number of elements' );if(~isempty(fidReadIn));fclose(fidRead);end;return;
            else
                NEVOtherHeaders.Payload = temp;
            end;clear temp
            
    end;clear HeaderName% switch on packet id
end;clear n

%% Resize NEVwaveformHeader based on data available
if( maximumChannelNumberObserved == 0 )
    NEVwaveformHeader = [];
else
    if( maximumChannelNumberObserved < NEVbasicHeader.maximumChannelNumber )
        NEVwaveformHeader = NEVwaveformHeader(1:maximumChannelNumberObserved);
    end
end
NEVbasicHeader.maximumChannelNumber = maximumChannelNumberObserved;
NEVwaveformHeader = NEVwaveformHeader(...
    ( [ NEVwaveformHeader(:).existNEUEVWAV ] ~= 0 ) | ...
    ( [ NEVwaveformHeader(:).existNEUEVLBL ] ~= 0 ) | ...
    ( [ NEVwaveformHeader(:).existNEUEVFLT ] ~= 0 ) ...
    );
if( ...
        ( ~NEVexternalIOHeader.existNSASEXEV ) && ...
        ( ~NEVexternalIOHeader.existSTIMINFO ) && ...
        ( ~NEVexternalIOHeader.existDIGLABELSerial ) && ...
        ( ~NEVexternalIOHeader.existDIGLABELParallel ) ...
        )
    NEVexternalIOHeader = [];
end
if( ...
        ( ~NEVarrayHeader.existARRAYNME ) && ...
        ( ~NEVarrayHeader.existMAPFILE ) ...
        )
    NEVarrayHeader = [];
end
NEVcommentsText = NEVcommentsText(...
    ( [ NEVcommentsText(:).existCCOMMENT ] | ...
    [ NEVcommentsText(:).existECOMMENT ] ...
    ));
if( isempty( NEVcommentsText ) );NEVcommentsText = [];end
if( ...
        ( ~NEVVideoSync.existVIDEOSYN ) ...
        )
    NEVVideoSync = [];
end
if( ...
        ( ~NEVTrackObj.existTRACKOBJ ) ...
        )
    NEVTrackObj = [];
end
NEVOtherHeaders = NEVOtherHeaders(...
    ( [ NEVOtherHeaders(:).existOther ] ~= 0 ) ...
    );
if( isempty( NEVOtherHeaders ) );NEVOtherHeaders = [];end
if(~isempty(fidReadIn));fclose(fidRead);end;

return;




