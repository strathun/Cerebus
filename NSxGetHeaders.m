function [ NSxbasicHeader, NSxchannelHeader] = NSxGetHeaders( fidRead )
% NSxGetHeaders - Reads headers of a NSx file into basic and channel structures.
%
%	function [ NSxbasicHeader, NSxchannelHeader] = NSxGetHeaders( fidRead );
%
%  Reads the header data from a NSx file - see the nev spec for interpretation
%
%   Inputs
%
%       fidRead             -   File ID of previously opened NSx file or a
%                               string containing the file to open
%
%                               Although not necessary in all cases, it is
%                               recommended that the file be opened with
%                               'rb' permissions, for read only and binary
%                               format(binary format is default, anyway).
%                           	Althoughnot necessary in all cases, it is
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
%   Outputs (See most recent NEV specificiation for details
%
%       NSxbasicHeader      -	Header data implemented as structure scalar
%                               detailing information common to entire NSx
%                               file.
%
%       NSxchannelHeader   -	Header data implemented as structure vector
%                               detailing information unique to particular
%                               electrode channels, including input channel
%                               information, label information, and
%                               filtering information.
%
%   $Author: WarrenDJ $, University of Utah
%   $DateTime: 2017/03/20 14:46:35 $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Magic Numbers
FileIDSpecBySpec = [ ...
    2.1, ...
    2.2, ...
    2.3 ...
    ];
FileIDCodes  = { ...
    'NEURALSG', ...
    'NEURALCD', ...
    'NEURALCD' ...
    };
basicHeaderPacketSizeBySpec = [ ...
    8 + 16 + 4 + 4, ...
    8 + 2 + 4 + 16 + 256 + 4 + 4 + 16 + 4, ...
    8 + 2 + 4 + 16 + 256 + 4 + 4 + 16 + 4 ...
    ];
extentedHeaderPacketSizeBySpec = [ ...
    4, ...
    ( 2 + 2 + 16 + 1 + 1 + 2 + 2 + 2 + 2 + 16 + 4 + 4 + 2 + 4 + 4 + 2), ...
    ( 2 + 2 + 16 + 1 + 1 + 2 + 2 + 2 + 2 + 16 + 4 + 4 + 2 + 4 + 4 + 2) ...
    ];
dataPacketHeaderSizeBySpec = [ ...
    0, ...
    1 + 4 + 4, ...
    1 + 4 + 4 ...
    ];
defaultSampleRes = 30000;
defaulttempType = 'CC';
defaulttempHeader = 1;
maxDigitalValueDefault = 32764;
minDigitalValueDefault = -32764;
maxAnalogValueElectrodeDefault = 8191;
minAnalogValueElectrodeDefault = -8191;
unitsOFDataElectrodeDefault = [ 'uV' char( zeros(1,14) ) ];
maxAnalogValueExternalDefault = 5000;
minAnalogValueExternalDefault = -5000;
unitsOFDataExternalDefault = [ 'mV' char( zeros(1,14) ) ];
highfreqcornerElectrodeDefault = 300; % Analog filtering by Cerebus amp
highfreqorderElectrodeDefault = 1;
highfreqtypeElectrodeDefault = 1;
lowfreqcornerElectrodeDefault = 7500000;
lowfreqorderElectrodeDefault = 3;
lowfreqtypeElectrodeDefault = 1;
sizeOfInt8 = 1;numel( sizeOfInt8 ); %  bytes
sizeOfInt16 = 2;numel( sizeOfInt16 ); %  bytes
sizeOfInt32 = 4;numel( sizeOfInt32 ); % % bytes
sizeOfInt64 = 8;numel( sizeOfInt64 ); % % bytes
sizeOfSingle = 4;numel( sizeOfSingle ); % bytes
sizeOfDouble = 8;numel( sizeOfDouble ); % bytes
sampleSize = sizeOfInt16;

%% Assure some output
NSxbasicHeader = [];
NSxchannelHeader = [];

%% extensive checking of inputs
pInput = inputParser();
pInput.addRequired( 'fidRead', ...
    @(x)( ( isscalar( x ) && ( x ~= -1 ) ) || ( ischar( x ) && exist( x, 'file' ) ) ) );
try
    pInput.parse( fidRead );
catch mExp
    error( 'NSxGetHeaders:invalidInputParameter', ...
        'Error: %s', mExp.message );
end%% Extensive error checking
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
    warning( 'NSxGetHeaders:FileNameError', ...
        'Unable to get filename of open file\n' );
    if(~isempty(fidReadIn));fclose(fidRead);end
    return;
end;
if( ~isempty( strfind( lower( permissions ), 't' ) ) )
    warning( 'NSxGetHeaders:FilePermissionsError', ...
        'File %s opened in text mode, may result in problems interpreting strings\n', filename );
end;
if( ~isempty( strfind( lower( permissions ), 'w' ) ) )
    warning( 'NSxGetHeaders:FilePermissionsError', ...
        'File %s opened for writing, may result in problems\n', filename );
end;
if( ~isempty( strfind( lower( permissions ), 'a' ) ) )
    warning( 'NSxGetHeaders:FilePermissionsError', ...
        'File %s opened for appending, may result in problems\n', filename );
end;
if( ~strcmpi( machineformat, 'ieee-le' ) )
    warning( 'NSxGetHeaders:FileFormatError', ...
        'File %s not opened in little-endian mode, may result in problems interpreting numbers\n', filename );
end;
if( ~strcmpi( encoding, 'windows-1252' ) )
    warning( 'NSxGetHeaders:FileEncodingError', ...
        'File %s not opened with correct text encoding, may result in problems interpreting strings\n', filename );
end;
clear permissions machineformat encoding

%%	Get details of file
fileDir = dir( filename );
if( isempty( fileDir ) )
    warning( 'NSxGetHeaders:FileNameError', ...
        'Unable to get details of open file\n' );
    if(~isempty(fidReadIn));fclose(fidRead);end
    return;
end;

%%	Position file to beginning
if( fseek( fidRead, 0, 'bof' ) == -1 )
    warning( 'NSxGetHeaders:FilePositioningError', ...
        'Invalid file positioning with message %s\n', ferror(fidRead,'clear') );
    if(~isempty(fidReadIn));fclose(fidRead);end
    return;
end;

%% Precreate basic header
NSxbasicHeader = struct( ...
    'id', blanks(8), ...
    'filespecMajor', 0, ...
    'filespecMinor', 0, ...
    'dataptr', 0, ...
    'label', blanks(16), ...
    'Comment', blanks(256), ...
    'StartTimeNIPTimeStamp', uint32( 0 ), ...
    'Period', 0, ...
    'SampleRes', 0, ...
    'FileTime', zeros( 8, 1), ...
    'NumChannels', 0, ...
    'filespecDouble', 0, ...
    'SerialDateNumber', 0, ...
    'Filename', blanks(1), ...
    'NumBytes', 0, ...
    'NumDataPackets', 0, ...
    'dataMatrixPtr', 0, ...
    'TimeStamp', 0, ...
    'NumTimeSamples', 0, ...
    'TimeStampEnd', 0 ...
    );

%%	Read Basic Headers
ncountTest = 8;[ temp, ncount ] = fread( fidRead, [1,ncountTest], 'char*1' );
if( ncount ~= ncountTest );warning( 'NSxGetHeaders:readCountError', 'Unable to read correct number of elements' );return;
else
    NSxbasicHeader.id = char( temp );
end

zFile = find( strcmpi( NSxbasicHeader.id, FileIDCodes ), 1 );
if( isempty( zFile ) )
    warning( 'NSxGetHeaders:badFileTypeID', ...
        'Unexpected file type id code of %s\n', NSxbasicHeader.id );
    if( ~isempty( fidReadIn ) );fclose(fidRead);end
    return;
else
    
    FileIDSpec = FileIDSpecBySpec( zFile );
    switch FileIDSpec
        case 2.1
            
            % fill in neceesary data not contained in file header
            NSxbasicHeader.filespecMajor = 2;
            NSxbasicHeader.filespecMinor = 1;
            NSxbasicHeader.filespecDouble = FileIDSpec;
            NSxbasicHeader.filespecMajor = ...
                floor( NSxbasicHeader.filespecDouble );
            NSxbasicHeader.filespecMinor = ...
                floor( 10*mod( NSxbasicHeader.filespecDouble, 1 ) );
            
            warning( 'NSxGetHeaders:earlierFileTypeID', ...
                'NEV Specification %3.1f has limited support\n', NSxbasicHeader.filespecDouble );
            
            basicHeaderPacketSize = basicHeaderPacketSizeBySpec( zFile );
            extentedHeaderPacketSize = extentedHeaderPacketSizeBySpec( zFile );
            dataPacketHeaderSize = dataPacketHeaderSizeBySpec( zFile );
            clear zFile;
            
            ncountTest = 16;[ temp, ncount ] = fread( fidRead, [1,ncountTest], 'char*1' );
            if( ncount ~= ncountTest );warning( 'NSxGetHeaders:readCountError', 'Unable to read correct number of elements' );if(~isempty(fidReadIn));fclose(fidRead);end;return;
            else
                NSxbasicHeader.label = char( temp );
            end
            
            ncountTest = 1;[ NSxbasicHeader.Period, ncount ] = fread( fidRead, [1,ncountTest], 'uint32' );
            if( ncount ~= ncountTest );warning( 'NSxGetHeaders:readCountError', 'Unable to read correct number of elements' );if(~isempty(fidReadIn));fclose(fidRead);end;return;end
            
            ncountTest = 1;[ NSxbasicHeader.NumChannels, ncount ] = fread( fidRead, [1,ncountTest], 'uint32' );
            if( ncount ~= ncountTest );warning( 'NSxGetHeaders:readCountError', 'Unable to read correct number of elements' );if(~isempty(fidReadIn));fclose(fidRead);end;return;end
            
            NSxbasicHeader.dataptr(1) = basicHeaderPacketSize + ...
                NSxbasicHeader.NumChannels * extentedHeaderPacketSize;
            NSxbasicHeader.dataMatrixPtr(1) = NSxbasicHeader.dataptr(1) + ...
                dataPacketHeaderSize;
            NSxbasicHeader.SampleRes = defaultSampleRes;
            dateVector = datevec( fileDir(1).datenum );
            NSxbasicHeader.FileTime = [ dateVector(1:2) weekday( fileDir(1).datenum ) dateVector(3:end) ];
            clear dateVector
            
        case { 2.2, 2.3 }
            
            ncountTest = 1;[ NSxbasicHeader.filespecMajor, ncount ] = fread( fidRead, [1,ncountTest], 'uint8' );
            if( ncount ~= ncountTest );warning( 'NSxGetHeaders:readCountError', 'Unable to read correct number of elements' );if(~isempty(fidReadIn));fclose(fidRead);end;return;end
            
            ncountTest = 1;[ NSxbasicHeader.filespecMinor, ncount ] = fread( fidRead, [1,ncountTest], 'uint8' );
            if( ncount ~= ncountTest );warning( 'NSxGetHeaders:readCountError', 'Unable to read correct number of elements' );if(~isempty(fidReadIn));fclose(fidRead);end;return;end
            
            NSxbasicHeader.filespecDouble = NSxbasicHeader.filespecMajor + ...
                ( NSxbasicHeader.filespecMinor / 10 );
            switch ( NSxbasicHeader.filespecDouble)
                case { 2.2, 2.3 }
                    FileIDSpec = NSxbasicHeader.filespecDouble;
            end
            
            zFile = find( FileIDSpecBySpec == FileIDSpec );
            if( isempty( zFile ) )
                warning( 'NSxGetHeaders:badFileTypeID', ...
                    'Unexpected file type id code of %s\n', NSxbasicHeader.id );
                return;
            end
            
            basicHeaderPacketSize = basicHeaderPacketSizeBySpec( zFile );
            extentedHeaderPacketSize = extentedHeaderPacketSizeBySpec( zFile );
            dataPacketHeaderSize = dataPacketHeaderSizeBySpec( zFile );
            clear zFile;
            
            ncountTest = 1;[ NSxbasicHeader.dataptr(1), ncount ] = fread( fidRead, [1,ncountTest], 'uint32' );
            if( ncount ~= ncountTest );warning( 'NSxGetHeaders:readCountError', 'Unable to read correct number of elements' );if(~isempty(fidReadIn));fclose(fidRead);end;return;end
            NSxbasicHeader.dataMatrixPtr(1) = NSxbasicHeader.dataptr(1) + dataPacketHeaderSize;

            
            ncountTest = 16;[ temp, ncount ] = fread( fidRead, [1,ncountTest], 'char*1' );
            if( ncount ~= ncountTest );warning( 'NSxGetHeaders:readCountError', 'Unable to read correct number of elements' );if(~isempty(fidReadIn));fclose(fidRead);end;return;
            else
                NSxbasicHeader.label = char( temp );
            end
            
            ncountTest = 256;[ temp, ncount ] = fread( fidRead, [1,ncountTest], 'char*1' );
            if( ncount ~= ncountTest );warning( 'NSxGetHeaders:readCountError', 'Unable to read correct number of elements' );if(~isempty(fidReadIn));fclose(fidRead);end;return;
            else
                if( ...
                        ( temp( 252 ) == 0 ) && ...
                        any( temp( 253:256 ) ~= 0 ) ...
                        )
                    % Ripple Trellis Version 1.4.1 and later
                    NSxbasicHeader.Comment = char( [ temp(1:252) zeros(1,4,'uint8') ] );
                    NSxbasicHeader.StartTimeNIPTimeStamp = typecast( uint8( temp( 253:256 ) ), 'uint32' );
                else
                    NSxbasicHeader.Comment = char( temp(1:end) );
                    NSxbasicHeader.StartTimeNIPTimeStamp = uint32( 0 );
                end
            end
            
            ncountTest = 1;[ NSxbasicHeader.Period, ncount ] = fread( fidRead, [1,ncountTest], 'uint32' );
            if( ncount ~= ncountTest );warning( 'NSxGetHeaders:readCountError', 'Unable to read correct number of elements' );if(~isempty(fidReadIn));fclose(fidRead);end;return;end
            
            ncountTest = 1;[ NSxbasicHeader.SampleRes, ncount ] = fread( fidRead, [1,ncountTest], 'uint32' );
            if( ncount ~= ncountTest );warning( 'NSxGetHeaders:readCountError', 'Unable to read correct number of elements' );if(~isempty(fidReadIn));fclose(fidRead);end;return;end
            
            ncountTest = 8;[ NSxbasicHeader.FileTime, ncount ] = fread( fidRead, [1,ncountTest], 'uint16' );
            if( ncount ~= ncountTest );warning( 'NSxGetHeaders:readCountError', 'Unable to read correct number of elements' );if(~isempty(fidReadIn));fclose(fidRead);end;return;end
            
            ncountTest = 1;[ NSxbasicHeader.NumChannels, ncount ] = fread( fidRead, [1,ncountTest], 'uint32' );
            if( ncount ~= ncountTest );warning( 'NSxGetHeaders:readCountError', 'Unable to read correct number of elements' );if(~isempty(fidReadIn));fclose(fidRead);end;return;end
            
            if( ( ...
                    NSxbasicHeader.dataptr(1) ~= ...
                    ( basicHeaderPacketSize + ...
                    extentedHeaderPacketSize * NSxbasicHeader.NumChannels ) ) )
                warning( 'NSxGetHeaders:badPacketSize', ...
                    'First pointer to first packet is corrupt but now fixed, %d not equal expected %d\n', ...
                    NSxbasicHeader.dataptr(1), ...
                    ( basicHeaderPacketSize + ...
                    extentedHeaderPacketSize * NSxbasicHeader.NumChannels ) );
                NSxbasicHeader.dataptr(1) = ...
                ( basicHeaderPacketSize + ...
                    extentedHeaderPacketSize * NSxbasicHeader.NumChannels );
            end
            
        otherwise
            warning( 'NSxGetHeaders:badFileTypeID', ...
                'Unexpected file specification for %f\n', FileIDSpec );
            return;
    end%switch FileIDSpec
end

% Verify and save some header data
if( NSxbasicHeader.filespecDouble ~= FileIDSpec )
    warning( 'NSxGetHeaders:badFileTypeID', ...
        'Unexpected file specification for %f\n', FileIDSpec );
    if(~isempty(fidReadIn));fclose(fidRead);end;
    return;
end;clear FileIDSpec
dateVector = NSxbasicHeader.FileTime( [1 2 4 5 6 7] );
if( dateVector(1) < 50 )
    dateVector(1) = dateVector(1) + 2000;
else
    if( dateVector(1) < 100 )
        dateVector(1) = dateVector(1) + 1900;
    end
end
NSxbasicHeader.SerialDateNumber = datenum( dateVector );
NSxbasicHeader.Filename = filename;
NSxbasicHeader.NumBytes = fileDir.bytes;
clear filename fileDir dateVector

%% Precreate channel header
NSxchannelHeader = struct( ...
    'id', num2cell( zeros( NSxbasicHeader.NumChannels, 1) ), ...
    'existExtendedHeader', num2cell( zeros( NSxbasicHeader.NumChannels, 1) ), ...
    'Type', num2cell( zeros( NSxbasicHeader.NumChannels, 3 ), 2 ), ...
    'label', num2cell( zeros( NSxbasicHeader.NumChannels, 17 ), 2 ), ...
    'pinch', num2cell( zeros( NSxbasicHeader.NumChannels, 1) ), ...
    'pinnum', num2cell( zeros( NSxbasicHeader.NumChannels, 1) ), ...
    'minDigitalValue', num2cell( zeros( NSxbasicHeader.NumChannels, 1) ), ...
    'maxDigitalValue', num2cell( zeros( NSxbasicHeader.NumChannels, 1) ), ...
    'minAnalogValue', num2cell( zeros( NSxbasicHeader.NumChannels, 1) ), ...
    'maxAnalogValue', num2cell( zeros( NSxbasicHeader.NumChannels, 1) ), ...
    'unitsOfData', num2cell( zeros( NSxbasicHeader.NumChannels, 1) ), ...
    'highfreqcorner', num2cell( zeros( NSxbasicHeader.NumChannels, 1) ), ...
    'highfreqorder', num2cell( zeros( NSxbasicHeader.NumChannels, 1) ), ...
    'highfreqtype', num2cell( zeros( NSxbasicHeader.NumChannels, 1) ), ...
    'lowfreqcorner', num2cell( zeros( NSxbasicHeader.NumChannels, 1) ), ...
    'lowfreqorder', num2cell( zeros( NSxbasicHeader.NumChannels, 1) ), ...
    'lowfreqtype', num2cell( zeros( NSxbasicHeader.NumChannels, 1) ), ...
    'resolution', num2cell( zeros( NSxbasicHeader.NumChannels, 1) ) ...
    );

%%	Read Extended Headers
for n=1:NSxbasicHeader.NumChannels;
    
    % Check that file position is correct
    if( ftell( fidRead ) ~= basicHeaderPacketSize + (n-1)*extentedHeaderPacketSize )
        warning( 'NSxGetHeaders:badFilePosition', ...
            'Invalid file positioning\n' );
        return
    end
    
    switch NSxbasicHeader.filespecDouble
        case 2.1
            
            % Only Cerebus, some limits apply
            maximumChannelNumber = 144;
            nChannelsOfElectrodes = 128;
            numberExternalStop = 144;

            ncountTest = 1;[ id, ncount ] = fread( fidRead, [1,ncountTest], 'uint32' );
            if( ncount ~= ncountTest );warning( 'NSxGetHeaders:readCountError', 'Unable to read correct number of elements' );if(~isempty(fidReadIn));fclose(fidRead);end;return;end
            if( ( id < 1 ) || ( id > maximumChannelNumber ) )
                warning( 'NSxGetHeaders:badChannelID', ...
                    'Invalid packet ID found (%d)\n', id );
                return
            end
            NSxchannelHeader(n).id = id;clear id;
            
            NSxchannelHeader(n).existExtendedHeader = 1;
            NSxchannelHeader(n).Type = defaulttempType;
            NSxchannelHeader(n).minDigitalValue =  minDigitalValueDefault;
            NSxchannelHeader(n).maxDigitalValue =  maxDigitalValueDefault;
            if( NSxchannelHeader(n).id <= nChannelsOfElectrodes )
                NSxchannelHeader(n).minAnalogValue = minAnalogValueElectrodeDefault;
                NSxchannelHeader(n).maxAnalogValue = maxAnalogValueElectrodeDefault;
                NSxchannelHeader(n).unitsOfData = unitsOFDataElectrodeDefault;
            else
                if( ...
                        ( NSxchannelHeader(n).id >= (nChannelsOfElectrodes+1) ) && ...
                        ( NSxchannelHeader(n).id <= numberExternalStop ) ...
                        )
                    NSxchannelHeader(n).minAnalogValue = minAnalogValueExternalDefault;
                    NSxchannelHeader(n).maxAnalogValue = maxAnalogValueExternalDefault;
                    NSxchannelHeader(n).unitsOfData = unitsOFDataExternalDefault;
                else
                    warning( 'NSxGetHeaders:badChannelID', ...
                        'Invalid packet ID found (%d)\n', NSxchannelHeader(n).id );
                    NSxchannelHeader(n).minAnalogValue = intmin('int16');
                    NSxchannelHeader(n).maxAnalogValue = intmax('int16');
                    NSxchannelHeader(n).unitsOfData = [ 'uV' char( zeros(1,14) ) ];
                end
            end
            if( NSxchannelHeader(n).id <= nChannelsOfElectrodes )
                NSxchannelHeader(n).highfreqcorner = highfreqcornerElectrodeDefault;
                NSxchannelHeader(n).highfreqorder = highfreqorderElectrodeDefault;
                NSxchannelHeader(n).highfreqtype = highfreqtypeElectrodeDefault;
                NSxchannelHeader(n).lowfreqcorner = lowfreqcornerElectrodeDefault;
                NSxchannelHeader(n).lowfreqorder = lowfreqorderElectrodeDefault;
                NSxchannelHeader(n).lowfreqtype = lowfreqtypeElectrodeDefault;
            else
                if( ...
                        ( NSxchannelHeader(n).id >= (nChannelsOfElectrodes+1) ) ...
                        && ( NSxchannelHeader(n).id <= numberExternalStop ) ...
                        )
                    % No specific information
                else
                    warning( 'NSxGetHeaders:badChannelID', ...
                        'Invalid packet ID found (%d)\n', id );
                    NSxchannelHeader(n).highfreqcorner = highfreqcornerElectrodeDefault;
                    NSxchannelHeader(n).highfreqorder = highfreqorderElectrodeDefault;
                    NSxchannelHeader(n).highfreqtype = highfreqtypeElectrodeDefault;
                    NSxchannelHeader(n).lowfreqcorner = lowfreqcornerElectrodeDefault;
                    NSxchannelHeader(n).lowfreqorder = lowfreqorderElectrodeDefault;
                    NSxchannelHeader(n).lowfreqtype = lowfreqtypeElectrodeDefault;
                end
            end
            clear maximumChannelNumber nChannelsOfElectrodes numberExternalStop
            
        case { 2.2, 2.3 }
            
            ncountTest = 2;[ tempType, ncount ] = fread( fidRead, [1,ncountTest], 'char*1' );
            if( ncount ~= ncountTest );warning( 'NSxGetHeaders:readCountError', 'Unable to read correct number of elements' );if(~isempty(fidReadIn));fclose(fidRead);end;return;
            else
            end
            tempType = char( tempType );
            if( strcmp( tempType, defaulttempType ) == 0 )
                warning( 'NSxGetHeaders:badextenedHeaderType', ...
                    'Unexpected type code of %s\n', tempType );
                return;
            end
            
            ncountTest = 1;[ id, ncount ] = fread( fidRead, [1,ncountTest], 'uint16' );
            if( ncount ~= ncountTest );warning( 'NSxGetHeaders:readCountError', 'Unable to read correct number of elements' );if(~isempty(fidReadIn));fclose(fidRead);end;return;end
            NSxchannelHeader(n).Type = tempType;clear tempType
            NSxchannelHeader(n).id = id;clear id
            
            NSxchannelHeader(n).existExtendedHeader = 1;
            
            ncountTest = 16;[ temp, ncount ] = fread( fidRead, [1,ncountTest], 'char*1' );
            if( ncount ~= ncountTest );warning( 'NSxGetHeaders:readCountError', 'Unable to read correct number of elements' );if(~isempty(fidReadIn));fclose(fidRead);end;return;
            else
                NSxchannelHeader(n).label = char( temp );
            end
            
            ncountTest = 1;[ NSxchannelHeader(n).pinch, ncount ] = fread( fidRead, [1,ncountTest], 'uint8' );
            if( ncount ~= ncountTest );warning( 'NSxGetHeaders:readCountError', 'Unable to read correct number of elements' );if(~isempty(fidReadIn));fclose(fidRead);end;return;end
            
            ncountTest = 1;[ NSxchannelHeader(n).pinnum, ncount ] = fread( fidRead, [1,ncountTest], 'uint8' );
            if( ncount ~= ncountTest );warning( 'NSxGetHeaders:readCountError', 'Unable to read correct number of elements' );if(~isempty(fidReadIn));fclose(fidRead);end;return;end
            
            ncountTest = 1;[ NSxchannelHeader(n).minDigitalValue, ncount ] = fread( fidRead, [1,ncountTest], 'int16' );
            if( ncount ~= ncountTest );warning( 'NSxGetHeaders:readCountError', 'Unable to read correct number of elements' );if(~isempty(fidReadIn));fclose(fidRead);end;return;end
            
            ncountTest = 1;[ NSxchannelHeader(n).maxDigitalValue, ncount ] = fread( fidRead, [1,ncountTest], 'int16' );
            if( ncount ~= ncountTest );warning( 'NSxGetHeaders:readCountError', 'Unable to read correct number of elements' );if(~isempty(fidReadIn));fclose(fidRead);end;return;end
            
            ncountTest = 1;[ NSxchannelHeader(n).minAnalogValue, ncount ] = fread( fidRead, [1,ncountTest], 'int16' );
            if( ncount ~= ncountTest );warning( 'NSxGetHeaders:readCountError', 'Unable to read correct number of elements' );if(~isempty(fidReadIn));fclose(fidRead);end;return;end
            
            ncountTest = 1;[ NSxchannelHeader(n).maxAnalogValue, ncount ] = fread( fidRead, [1,ncountTest], 'int16' );
            if( ncount ~= ncountTest );warning( 'NSxGetHeaders:readCountError', 'Unable to read correct number of elements' );if(~isempty(fidReadIn));fclose(fidRead);end;return;end
            
            ncountTest = 16;[ temp, ncount ] = fread( fidRead, [1,ncountTest], 'char*1' );
            if( ncount ~= ncountTest );warning( 'NSxGetHeaders:readCountError', 'Unable to read correct number of elements' );if(~isempty(fidReadIn));fclose(fidRead);end;return;
            else
                NSxchannelHeader(n).unitsOfData = char( temp );
            end
            
            ncountTest = 1;[ NSxchannelHeader(n).highfreqcorner, ncount ] = fread( fidRead, [1,ncountTest], 'uint32' );
            if( ncount ~= ncountTest );warning( 'NSxGetHeaders:readCountError', 'Unable to read correct number of elements' );if(~isempty(fidReadIn));fclose(fidRead);end;return;end
            
            ncountTest = 1;[ NSxchannelHeader(n).highfreqorder, ncount ] = fread( fidRead, [1,ncountTest], 'uint32' );
            if( ncount ~= ncountTest );warning( 'NSxGetHeaders:readCountError', 'Unable to read correct number of elements' );if(~isempty(fidReadIn));fclose(fidRead);end;return;end
            
            ncountTest = 1;[ NSxchannelHeader(n).highfreqtype, ncount ] = fread( fidRead, [1,ncountTest], 'uint16' );
            if( ncount ~= ncountTest );warning( 'NSxGetHeaders:readCountError', 'Unable to read correct number of elements' );if(~isempty(fidReadIn));fclose(fidRead);end;return;end
            
            ncountTest = 1;[ NSxchannelHeader(n).lowfreqcorner, ncount ] = fread( fidRead, [1,ncountTest], 'uint32' );
            if( ncount ~= ncountTest );warning( 'NSxGetHeaders:readCountError', 'Unable to read correct number of elements' );if(~isempty(fidReadIn));fclose(fidRead);end;return;end
            
            ncountTest = 1;[ NSxchannelHeader(n).lowfreqorder, ncount ] = fread( fidRead, [1,ncountTest], 'uint32' );
            if( ncount ~= ncountTest );warning( 'NSxGetHeaders:readCountError', 'Unable to read correct number of elements' );if(~isempty(fidReadIn));fclose(fidRead);end;return;end
            
            ncountTest = 1;[ NSxchannelHeader(n).lowfreqtype, ncount ] = fread( fidRead, [1,ncountTest], 'uint16' );
            if( ncount ~= ncountTest );warning( 'NSxGetHeaders:readCountError', 'Unable to read correct number of elements' );if(~isempty(fidReadIn));fclose(fidRead);end;return;end
            
        otherwise
            warning( 'NSxGetHeaders:badFileTypeID', ...
                'Unexpected file specification for %f\n', NSxbasicHeader.filespecDouble );
            return;
    end%switch NSxbasicHeader.filespecDouble
    NSxchannelHeader(n).resolution = ...
        ( double( NSxchannelHeader(n).maxAnalogValue ) - double( NSxchannelHeader(n).minAnalogValue ) ) / ...
        ( double( NSxchannelHeader(n).maxDigitalValue ) - double( NSxchannelHeader(n).minDigitalValue ) );
end;clear n
clear maximumChannelNumber nChannelsOfElectrodes numberExternalStop


%% Resize NSxchannelHeader based on data available
NSxchannelHeader = ...
    NSxchannelHeader(( [  NSxchannelHeader(:).existExtendedHeader ] == 1 ));

%%	Read Data Packet Headers
switch NSxbasicHeader.filespecDouble
    case 2.1
        NSxbasicHeader.NumDataPackets = 1;
        NSxbasicHeader.TimeStamp(NSxbasicHeader.NumDataPackets) = 0;
        NSxbasicHeader.NumTimeSamples(NSxbasicHeader.NumDataPackets) = ...
            ( NSxbasicHeader.NumBytes - NSxbasicHeader.dataMatrixPtr(NSxbasicHeader.NumDataPackets) ) ...
            / ( sampleSize * NSxbasicHeader.NumChannels );
        NSxbasicHeader.TimeStampEnd(NSxbasicHeader.NumDataPackets) = ...
            NSxbasicHeader.TimeStamp(NSxbasicHeader.NumDataPackets) + NSxbasicHeader.Period * ...
            max( 0, (  NSxbasicHeader.NumTimeSamples(NSxbasicHeader.NumDataPackets) - 1 ) );
        if( mod( NSxbasicHeader.NumTimeSamples, 1 ) ~= 0 )
            warning( 'NSxGetHeaders:PacketSizeError', ...
                'Non-integer number of data packet at %f\n', NSxbasicHeader.NumTimeSamples );
            NSxbasicHeader.NumTimeSamples(NSxbasicHeader.NumDataPackets) = ...
                floor( NSxbasicHeader.NumTimeSamples(NSxbasicHeader.NumDataPackets) );
            return;
        end;
    case { 2.2, 2.3 }
        NSxbasicHeader.NumDataPackets = 0;
        while( ~feof( fidRead ) )
            
            % At least space for data header packet?
            if( ( ftell( fidRead )  + dataPacketHeaderSize ) <= NSxbasicHeader.NumBytes )
                
                ncountTest = 1;[ tempHeader, ncount ] = fread( fidRead, [1,ncountTest], 'uint8' );
                if( ncount ~= ncountTest );warning( 'NSxGetHeaders:readCountError', 'Unable to read correct number of elements' );if(~isempty(fidReadIn));fclose(fidRead);end;return;end
                if( tempHeader ~= defaulttempHeader )
                    warning( 'NSxGetHeaders:badPacketHeader', ...
                        'Invalid packet header found\n' );
                    break
                end;clear tempHeader
                
                NSxbasicHeader.NumDataPackets = ...
                    NSxbasicHeader.NumDataPackets + 1;
                
                NSxbasicHeader.dataptr(NSxbasicHeader.NumDataPackets) = ftell( fidRead ) - ncountTest;
                
                ncountTest = 1;[ NSxbasicHeader.TimeStamp(NSxbasicHeader.NumDataPackets), ncount ] = fread( fidRead, [1,ncountTest], 'uint32' );
                if( ncount ~= ncountTest );warning( 'NSxGetHeaders:readCountError', 'Unable to read correct number of elements' );if(~isempty(fidReadIn));fclose(fidRead);end;return;end
                
                ncountTest = 1;[ NSxbasicHeader.NumTimeSamples(NSxbasicHeader.NumDataPackets), ncount ] = fread( fidRead, [1,ncountTest], 'uint32' );
                if( ncount ~= ncountTest );warning( 'NSxGetHeaders:readCountError', 'Unable to read correct number of elements' );if(~isempty(fidReadIn));fclose(fidRead);end;return;end
                
                NSxbasicHeader.dataMatrixPtr(NSxbasicHeader.NumDataPackets) = ftell( fidRead );
                
                numberBytesInPacket = ( sampleSize * ...
                    NSxbasicHeader.NumChannels * ...
                    NSxbasicHeader.NumTimeSamples(NSxbasicHeader.NumDataPackets) );
                if( ( NSxbasicHeader.dataMatrixPtr(NSxbasicHeader.NumDataPackets) + numberBytesInPacket )...
                        > NSxbasicHeader.NumBytes )
                    warning( 'NSxGetHeaders:PacketSizeError', ...
                        'Incomplete data packet, truncating number of time samples to what is available\n' );
                    NSxbasicHeader.NumTimeSamples(NSxbasicHeader.NumDataPackets) = floor( ...
                        ( NSxbasicHeader.NumBytes - NSxbasicHeader.dataMatrixPtr(NSxbasicHeader.NumDataPackets) ) ...
                        / ( sampleSize * NSxbasicHeader.NumChannels ) );
                    numberBytesInPacket = ( sampleSize * ...
                        NSxbasicHeader.NumChannels * ...
                        NSxbasicHeader.NumTimeSamples(NSxbasicHeader.NumDataPackets) );
                end
                
                % catch bogus number of bytes in packet
                if( numberBytesInPacket == 0 )
                    % First check if actually could be no bytes
                    ncountTest = 1;[ tempHeader, ncount ] = fread( fidRead, [1,ncountTest], 'uint8' );
                    if( ncount ~= ncountTest );warning( 'NSxGetHeaders:readCountError', 'Unable to read correct number of elements' );if(~isempty(fidReadIn));fclose(fidRead);end;return;end
                    fseek( fidRead, -1, 'cof' );
                    if( tempHeader ~= defaulttempHeader )
                        % best guess at what's left
                        availableBytes = NSxbasicHeader.NumBytes - ...
                            ( basicHeaderPacketSize + ...
                            extentedHeaderPacketSize * ...
                            NSxbasicHeader.NumChannels  + ...
                            dataPacketHeaderSize );
                        NSxbasicHeader.NumTimeSamples(NSxbasicHeader.NumDataPackets) = ...
                            floor( ...
                            availableBytes / ...
                            ( sampleSize * ...
                            NSxbasicHeader.NumChannels ) );
                        numberBytesInPacket = ( sampleSize * ...
                            NSxbasicHeader.NumChannels * ...
                            NSxbasicHeader.NumTimeSamples(NSxbasicHeader.NumDataPackets) );
                        clear availableBytes
                    end;clear tempHeader
                end
                    
                NSxbasicHeader.TimeStampEnd(NSxbasicHeader.NumDataPackets) = ...
                    NSxbasicHeader.TimeStamp(NSxbasicHeader.NumDataPackets) + NSxbasicHeader.Period * ...
                    max( 0, (  NSxbasicHeader.NumTimeSamples(NSxbasicHeader.NumDataPackets) - 1 ) );
                
                if( fseek( fidRead, numberBytesInPacket, 'cof' ) == -1 )
                    warning( 'NSxGetHeaders:FilePositioningError', ...
                        'Invalid file positioning with message %s\n', ferror(fidRead,'clear') );
                    return;
                end
                
            else
                if( ftell( fidRead ) ~= NSxbasicHeader.NumBytes )
                    warning( 'NSxGetHeaders:badPacketHeader', ...
                        'Incomplete data header packet\n' );
                    return;
                else
                    break;
                end
            end
        end
    otherwise
        warning( 'NSxGetHeaders:badFileTypeID', ...
            'Unexpected file specification for %f\n', NSxbasicHeader.filespecDouble );
        return;
end%switch NSxbasicHeader.filespecDouble

if(~isempty(fidReadIn));
    fclose(fidRead);
end;

return;




