function summaryFig = NSxNEVSummaryPlot( NSxSourceFilename, SinglePlotChannels, AllPagesChannels, varargin )
% NSxNEVSummaryPlot - Plots a summary of the contents of a NS5 and NEV file
%
%   NSxNEVSummaryPlot( NSxSourceFilename, SinglePlotChannels, AllPagesChannels, varargin )
%
%   This function plots the contents of a NS5 and NEV file in an easily
%   interpretable manner. Given the large amount of NS5 data, these data
%   are plotted as a line between the maximum and minimum values over a
%   range of time.
%
%   Note: this subroutine is NS5-biased. That is, it will only plot data if
%   it is contained in the NS5 file.
%
%   Only NEV specification 2.1 and later are supported.
%
%   Inputs
%
%       NSxSourceFilename       -   Name of NS5 and NEV files to plot. The
%                                   file extension can be missing or
%                                   provided as NS5 ro NEV.
%
%       SinglePlotChannels      -   A list of the channels to be plotted in
%                                   a individual x-y axes. If not provided
%                                   or an empty matrix is provided, all
%                                   neural channels will be plotted.
%
%       AllPagesChannels        -   A list of the channels to be plotted in
%                                   a single x-y axes at the top of each
%                                   page. If not provided or an empty
%                                   matrix is provided, all external analog
%                                   input channels will be used.
%
%       singlePlotScale         -   Optional parameter pair of
%                                   'singlePlotScale' and a scalar value
%                                   that will be used to scale all channel
%                                   provided by SinglePlotChannels. The
%                                   value is the largest y-axis positive or
%                                   negative range that will be plotted. If
%                                   the data fix into a smaller range, then
%                                   the smaller range will be plotted. If
%                                   not provided, a default value of 100
%                                   will be used. Note: neural data will be
%                                   plotted scaled to the amplifier input,
%                                   and most likley will have units of uV.
%                                   External analog inputs typically have
%                                   units of mV.
%
%       singleSampleBin         -   The size of the bin (in units of time
%                                   samples that the maximum and minimum
%                                   values will be calculated. If not
%                                   provided, a default value of 1000
%                                   samples will be used.
%
%       timeRangeStart          -   Optional paramter pair of
%                                   'timeRangeStart' and the first point in
%                                   time to begin the filtering and saving
%                                   of data. If not provided, a default
%                                   value of 0 seconds will be used.
%
%       timeRangeEnd            -   Optional paramter pair of
%                                   'timeRangeEnd' and the last point in
%                                   time that filtering and saving of data
%                                   is performed. If not provided, a
%                                   default value of infinity will be used.
%
%       RippleBasedFile         -   Optional paramter pair of
%                                   'RippleBasedFile' and any value. If the
%                                   value provided is non-zero, it will be
%                                   assumed that this file uses Ripple's
%                                   version of the NSx specification, which
%                                   only affects the range of neuronal
%                                   channels. If not provided, a default
%                                   value of 0 will be used (and the file
%                                   will assume a Cerebus range of
%                                   channels).
%
%       DestinationDir          -   Optional paramter pair of
%                                   'DestinationDir' and any string. All
%                                   output files will be written to the
%                                   directory provided the string. If not
%                                   provided, all output will be writting
%                                   to the current directory. Do not
%                                   include the ending file separator.
%
%   Outputs - None
%
%   $author$, University of Utah
%   $date$
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Magic numbers
nRow = 9;nCol=6;
hAxis = NaN*ones( nRow, 2 );
hLine = NaN*ones( nRow, 4 );
hAxisSuplabel = NaN;numel(hAxisSuplabel);
SizeOfInt16 = 2;numel(SizeOfInt16);
SizeOfUint16 = 2;numel(SizeOfUint16);
SizeOfInt32 = 4;numel(SizeOfInt32);
SizeOfUint32 = 4;numel(SizeOfUint32);
singlePlotScaleDefault = 100;
singleSampleBinDefault = 1000;
timeRangeStartDefault = 0;
timeRangeEndDefault = Inf;
SinglePlotChannelsDefault = [];
AllPagesChannelsDefault = [];
nChannelsOfElectrodesCerebus = 128;
nChannelsOfElectrodesRipple = 5120;
ChannelsOfExternalInputsCerebus = 129:144;
ChannelsOfExternalInputsRipple = 10241:65531;
RippleBasedFileDefault = 0;
DestinationDirDefault = '.';

%% Create all outputs
summaryFig = [];

%% Parse inputs
pInput = inputParser();
pInput.addRequired( 'NSxSourceFilename', ...
    @(x) ( ischar(x) && isvector(x) && exist( x, 'file' ) ) );
pInput.addOptional( 'SinglePlotChannels', SinglePlotChannelsDefault, ...
    @(x) ( isempty(x) || ( isvector(x) && isnumeric(x) && all( x > 0 ) ) ) );
pInput.addOptional( 'AllPagesChannels', AllPagesChannelsDefault, ...
    @(x) ( isempty(x) || ( isvector(x) && isnumeric(x) && all( x > 0 ) ) ) );
pInput.addParamValue( 'singlePlotScale', singlePlotScaleDefault, ...
    @(x) ( isscalar(x) && isnumeric(x) && ( x > 0 ) ) && ( x <= Inf ) );
pInput.addParamValue( 'singleSampleBin', singleSampleBinDefault, ...
    @(x) ( isscalar(x) && isnumeric(x) && ( x > 0 ) ) && ( x <= Inf ) );
pInput.addParamValue( 'timeRangeStart', timeRangeStartDefault, ...
    @(x) ( isscalar(x) && isnumeric(x) && ( x >= 0 ) ) ); %#ok<*NVREPL>
pInput.addParamValue( 'timeRangeEnd', timeRangeEndDefault, ...
    @(x) ( isscalar(x) && isnumeric(x) && ( x <= Inf ) ) );
pInput.addParamValue( 'RippleBasedFile', RippleBasedFileDefault, ...
    @(x) ( isscalar(x) && isnumeric(x) ) );
pInput.addParamValue( 'DestinationDir', DestinationDirDefault, ...
    @(x) ( ischar(x) && isvector(x) && exist( x, 'dir' ) ) );
try
    pInput.parse( NSxSourceFilename, SinglePlotChannels, AllPagesChannels, varargin{:} );
catch mExp
    error( 'SummaryPlot:invalidInputParameter', ...
        'Error: %s', mExp.message );
end
SinglePlotChannels = pInput.Results.SinglePlotChannels;
AllPagesChannels = pInput.Results.AllPagesChannels;
singlePlotScale = pInput.Results.singlePlotScale;
singleSampleBin = pInput.Results.singleSampleBin;
timeRangeStart = pInput.Results.timeRangeStart;
timeRangeEnd = pInput.Results.timeRangeEnd;
RippleBasedFile = pInput.Results.RippleBasedFile;
DestinationDir = pInput.Results.DestinationDir;
clear pInput varargin
clear timeRangeStartDefault timeRangeEndDefault singlePlotScaleDefault
if( timeRangeEnd < timeRangeStart )
    warning( 'SummaryPlot:badTimeRange', ...
        'Time range requested (%f to %f) is badly defined\n\t%s is not filtered', ...
        timeRangeStart, timeRangeEnd, NSxSourceFilename );
    return;
end

%% Handle changes due to parameters
if( RippleBasedFile > 0 )
    ChannelsOfElectrodes = 1:nChannelsOfElectrodesRipple;
    ChannelsOfExternalInputs = ChannelsOfExternalInputsRipple;
else
    ChannelsOfElectrodes = 1:nChannelsOfElectrodesCerebus;
    ChannelsOfExternalInputs = ChannelsOfExternalInputsCerebus;
end
clear nChannelsOfElectrodesCerebus nChannelsOfElectrodesRipple
clear ChannelsOfExternalInputsRipple ChannelsOfExternalInputsCerebus

%% Create basic filenames and read NSx file headers
[ filenameNEV, filenameBase ] = forcefilenameextension( NSxSourceFilename, '.nev' );
[ pathStr, filenameStr, extensionStr ] = fileparts( NSxSourceFilename );numel( pathStr );clear pathStr
filenamePDF = forcefilenameextension( ...
    fullfile( DestinationDir, ...
    strcat( filenameStr, '_', extensionStr(2:end), '_SummaryPlot' ) ), ...
    '.pdf' );
clear pathStr filenameStr
fidRead = fopen( NSxSourceFilename, 'rb' );
[ NSxbasicHeader, NSxchannelHeader] = NSxGetHeaders( fidRead );
fclose( fidRead );clear fidRead ans

%% set up figure and subplots
summaryFig = figure;
set( summaryFig, 'PaperPositionMode', 'manual' ); 
orient landscape
for n = 1:nRow
    contFig = nCol*(n-1)+(1:(nCol-1));
    spikeFig = nCol*(n-1)+nCol;
    hAxis(n,1)=subtightplot( nRow, nCol, contFig, [ 0.025 0.05 ], ...
        [1/8.5 1/8.5], [1/11 1/11]  );
    hLine(n,1) = plot( hAxis(n,1), ...
        [0 1],[0 1], '.k-', ...
        'MarkerSize', 2 ...
        );
    hold on
    hLine(n,3) = plot( hAxis(n,1), ...
        [0 1], [0 1], 'r--', ...
        'LineWidth', 1 );
    hold off
    hAxis(n,2)=subtightplot( nRow, nCol, spikeFig, [ 0.025 0.05 ], ...
        [1/8.5 1/8.5], [1/11 1/11] );
    hLine(n,2) = plot( hAxis(n,2), ...
        [0 1], [0 1], 'k-' );
    hold on
    hLine(n,4) = plot( hAxis(n,2), ...
        [0 1], [0 1], 'r--', ...
        'LineWidth', 1);
    hold off
    if( n ~= nRow );
        set(hAxis(n,1:2),'XTickLabel','');
    else
        xlabel( hAxis(n,1), 'Time (s)' );
        xlabel( hAxis(n,2), 'Time (ms)' );
    end
end;clear n
hAxisSuplabel = suplabel( ...
    fixtexspecialcharacters( filenameBase ), 't' );
allAxes = [ hAxis(:); hAxisSuplabel ];

%% Convert desired time range into time stamps and indeces into NSx file
timeStartTimeStampNSx = floor( ...
    max( NSxbasicHeader.TimeStamp(1), ...
    floor( timeRangeStart * NSxbasicHeader.SampleRes ) ) );
timeStopTimeStampNSx = ceil( ...
    min( NSxbasicHeader.TimeStampEnd(end), ...
    ceil( timeRangeEnd * NSxbasicHeader.SampleRes ) ) );
for p = 1:NSxbasicHeader.NumDataPackets
    if( ...
            ( timeStartTimeStampNSx >= NSxbasicHeader.TimeStamp(p) ) && ...
            ( timeStartTimeStampNSx <= NSxbasicHeader.TimeStampEnd(p) ) ...
            )
        dataBlockNSx = p;
        if( timeStopTimeStampNSx >  NSxbasicHeader.TimeStampEnd(p) )
            fprintf( ...
                'Truncating ending time stamp to put all in one data block\n' );
            timeStopTimeStampNSx = NSxbasicHeader.TimeStampEnd(p);
        end
        break;
    end
    if( p == NSxbasicHeader.NumDataPackets )
        fprintf( 'Did not find valid time stamp range\n' )
        return;
    end
end;clear p
timeStartTimeStampNSxIndex = max( 0, ...
    floor( ...
    ( timeStartTimeStampNSx ...
    - NSxbasicHeader.TimeStamp(dataBlockNSx) ) ...
    / NSxbasicHeader.Period ) );
timeStopTimeStampNSxIndex = min( NSxbasicHeader.NumTimeSamples(dataBlockNSx), ...
    ceil( ...
    ( timeStopTimeStampNSx ...
    - NSxbasicHeader.TimeStamp(dataBlockNSx) ) ...
    / NSxbasicHeader.Period ) );
nTimeStampNSxIndeces = timeStopTimeStampNSxIndex - ...
    timeStartTimeStampNSxIndex + 1;

%% read NEV file headers and create memory map file
if( ~exist( filenameNEV, 'file' ) )
    fprintf( 'The NEV file %s does not exist, implying no waveform data is plotted\n', ...
        filenameNEV );
    NEVFileArray = [];
    NEVbasicHeader = [];
    NEVExists = 0;
else
    [ NEVFileArray, NEVbasicHeader, NEVwaveformHeader ] = NEVGetMemMapFile( filenameNEV );
    if( isempty( NEVFileArray ) )
        warning( 'SummaryPlot:unopenableNEVFile %s', ...
            filenameNEV );
        NEVExists = 0;
    else
        NEVExists = 1;
    end
end

%% Extract NEV times and thresholds
if( NEVExists == 1 )
    NEVWaveformSize = (NEVbasicHeader.datasize - 8)/2;
    NEVWaveformTime = ...
        single((((1:NEVWaveformSize)-1)/NEVbasicHeader.SampleRes)*1e3);
    clear NEVWaveformSize
    NEVThresholdHigh = [NEVwaveformHeader(:).highthreshold].*[NEVwaveformHeader(:).sf]/1000;
    NEVThresholdLow  = [NEVwaveformHeader(:).lowthreshold ].*[NEVwaveformHeader(:).sf]/1000;
    NEVThreshold = NaN*ones( size(  NEVThresholdHigh ) );
    NEVThreshold( NEVThresholdHigh ~= 0 ) = NEVThresholdHigh( NEVThresholdHigh ~= 0 );
    NEVThreshold( NEVThresholdLow ~= 0 ) = NEVThresholdLow( NEVThresholdLow ~= 0 );
    clear NEVThresholdLow NEVThresholdHigh
else
    NEVWaveformTime = ((1:48)-1)*1e3/30000;
    NEVThreshold = NaN;
end

%% Convert desired time range into time stamps in NEV file
if( NEVExists == 1 )
timeStartTimeStampNEV = floor( ...
    max( double( intmin( 'uint32' ) ), timeRangeStart ) ...
    * NEVbasicHeader.TimeRes );
timeStopTimeStampNEV = ceil( ...
    min( double( intmax( 'uint32' ) ), timeRangeEnd ) ...
    * NEVbasicHeader.TimeRes );
else
    timeStartTimeStampNEV = NaN;
    timeStopTimeStampNEV = NaN;
end

%% Get list of channels available & figure out what channels to plot
chanIdNSx = [NSxchannelHeader(:).id];
if( NEVExists == 1 )
    chanIdNEV = [NEVwaveformHeader(:).id]; 
    chanIdNSxNEV = intersect( chanIdNSx, chanIdNEV );
else
    chanIdNEV = [];
    NEVWaveformTime = ((1:48)-1)*1e3/30000;
    chanIdNSxNEV = chanIdNSx;
end

%% set up NSx memory map files
NSxmmf = NSxGetMemMapFile( NSxSourceFilename, ...
    'DataPacketIndex', dataBlockNSx, ...
    'SampleOffset', timeStartTimeStampNSxIndex, ...
    'WriteableFlag', false );

% Load pertenent data from NEV file
if( NEVExists == 1 )
    NEVTimeStamps = typecast( ...
        reshape( NEVFileArray.Data.Payload(1:2,:), [ (SizeOfUint32/SizeOfUint16)*NEVbasicHeader.NumPackets, 1] ), ...
        'uint32' );
    NEVTimeStampsIndeces = find( ...
        ( NEVTimeStamps >= timeStartTimeStampNEV ) & ...
        ( NEVTimeStamps <= timeStopTimeStampNEV ) ...
        );
    NEVTimeStamps = NEVTimeStamps( NEVTimeStampsIndeces );numel( NEVTimeStamps );
    NEVPacketID = typecast( ...
        reshape( NEVFileArray.Data.Payload(3,NEVTimeStampsIndeces), [ numel( NEVTimeStamps ), 1] ), ...
        'uint16' );numel( NEVPacketID );
    NEVUnitID = reshape( ...
        typecast( NEVFileArray.Data.Payload(4,NEVTimeStampsIndeces), 'uint8' ), ...
        [ 2, numel( NEVTimeStamps ) ] );
    NEVUnitID = NEVUnitID(1,:)';numel( NEVUnitID );
    NEVWaveforms = NEVFileArray.Data.Payload(5:end,NEVTimeStampsIndeces);
    clear NEVTimeStampsIndeces
else
    NEVTimeStamps = [];numel( NEVTimeStamps );
    NEVPacketID = [];
    NEVUnitID = [];
    NEVWaveforms = [];
end

%% Create vectors to hold NSx data
nTimeSamplesPlot = ceil(nTimeStampNSxIndeces/singleSampleBin);
allPlotData = NaN*ones(3,nTimeSamplesPlot,NSxbasicHeader.NumChannels,'single');
tNSx = NaN*ones(3,nTimeSamplesPlot,'single');
tNSx(1,:) = (((1:nTimeSamplesPlot)'-1)* ...
    singleSampleBin*NSxbasicHeader.Period + timeStartTimeStampNSx)/NSxbasicHeader.SampleRes;
tNSx(2,:)=tNSx(1,:);

%% Preload plot data
for n = 1:nRow;
    set( hLine(n,1), ...
        'XData', tNSx(:), ...
        'YData', reshape( allPlotData(:,:,1), [size(allPlotData,1)*size(allPlotData,2),1] ) ...
        );
    set( hLine(n,3), ...
        'XData', [tNSx(1,1) tNSx(1,end)], ...
        'YData', [ NaN NaN ] ...
        );
    axis( hAxis(n,1), ...
        [ tNSx(1,1) tNSx(1,end) -Inf Inf ] );
    set( hLine(n,2), ...
        'XData', NEVWaveformTime, ...
        'YData', NaN*NEVWaveformTime ...
        );
    set( hLine(n,4), ...
        'XData', [NEVWaveformTime(1) NEVWaveformTime(end)], ...
        'YData', [NaN NaN] ...
        );
    axis( hAxis(n,2), ...
        [ NEVWaveformTime(1) NEVWaveformTime(end) -Inf Inf ] );
end;clear n

%% Load pertenent data from NSx file
pc = 0;
for n = 1:singleSampleBin:nTimeStampNSxIndeces
    b = single( NSxmmf.Data.NSxData(:,n:min((n+(singleSampleBin-1)),nTimeStampNSxIndeces)) );
    allPlotData(1,1+(n-1)/singleSampleBin,:) = min( b, [], 2 );
    allPlotData(2,1+(n-1)/singleSampleBin,:) = max( b, [], 2 );
    if( floor( 100*(n-1)/nTimeStampNSxIndeces ) ~= pc )
        pc = floor( 100*(n-1)/nTimeStampNSxIndeces );
        fprintf( '.' );
        if( mod(pc-1,25) == 24 );fprintf( '\n' );end
    end
end;clear n b
if( pc < 100 )
    for n = (pc+1):100
        fprintf( '.' );
    end;clear n
    fprintf( '\n' )
end;clear pc

%% Rescale NSx data
allPlotData(1,:,:) = squeeze(allPlotData(1,:,:))*diag(single([NSxchannelHeader(:).resolution]));
allPlotData(2,:,:) = squeeze(allPlotData(2,:,:))*diag(single([NSxchannelHeader(:).resolution]));

%% Plot the data that is on all pages
if( isempty( AllPagesChannels ) )
    [ ~, AllPagesChannelsIndecesNSx ] = ...
        intersect( chanIdNSx(:), ChannelsOfExternalInputs(:) );
else
    [ ~, AllPagesChannelsIndecesNSx ] = ...
        intersect( chanIdNSx(:), AllPagesChannels(:) );
end
if( ~isempty( AllPagesChannelsIndecesNSx ) )
    plot( hAxis( 1, 1 ), ...
        tNSx(:), ...
        reshape( allPlotData(:,:,AllPagesChannelsIndecesNSx), ...
        [size(allPlotData,1)*size(allPlotData,2),length(AllPagesChannelsIndecesNSx)] ...
        ) );
    ylab = fixtexspecialcharacters( ...
        deblank( NSxchannelHeader(AllPagesChannelsIndecesNSx(1)).unitsOfData ) );
    set( ...
        get( hAxis(1,1), 'YLabel' ), ...
        'String', ylab, ...
        'FontSize', 10 );
    axis( hAxis( 1, 1 ), 'tight' );
    axis( hAxis( 1, 1 ), 'auto y' );
    
    set( hAxis(1,1), 'XTickLabel', '' );
end;clear AllPagesChannelsIndecesNSx ylab

%% Plot the data that is single use
if( isempty( SinglePlotChannels ) )
    SinglePlotChannelsWork = ...
        intersect( chanIdNSxNEV(:), ChannelsOfElectrodes(:) );
else
    SinglePlotChannelsWork = ...
        intersect( chanIdNSxNEV(:), SinglePlotChannels(:) );
end
for n = 1:length( SinglePlotChannelsWork )
    
    % Plotting indeces
    m = 1 + mod( n-1, (nRow-1) );
    pflag = 0;
    
    % Find NSx and NEV indeces
    zNSx = find( chanIdNSx == SinglePlotChannelsWork(n) );
    zNEV = find( chanIdNEV == SinglePlotChannelsWork(n) );
    
    % Plot NSx data
    if( length( zNSx ) == 1 )
        set( hLine(m+1,1), ...
            'YData', ...
            reshape( allPlotData(:,:,zNSx), ...
            [size(allPlotData,1)*size(allPlotData,2),1] ), ...
            'color', 'k' );
        ylab = fixtexspecialcharacters( [ ...
            deblank( NSxchannelHeader(zNSx).label ) ' ' ...
            deblank( NSxchannelHeader(zNSx).unitsOfData ) ] );
        set( ...
            get( hAxis(m+1,1), 'YLabel' ), ...
            'String', ylab, ...
            'FontSize', 8 );
    else
        set( hLine(m+1,1), ...
            'YData', ...
            reshape( NaN*allPlotData(:,:,1), ...
            [size(allPlotData,1)*size(allPlotData,2),1] ), ...
            'color', 'k' );
        ylab = ' ';
        set( ...
            get( hAxis(m+1,1), 'YLabel' ), ...
            'String', ylab, ...
            'FontSize', 8 );
    end;clear zNSx
    axis( hAxis(m+1,1), 'auto y' );
    if( max(abs( get( hAxis(m+1,1), 'YLim' ) ) ) > singlePlotScale )
        set( hAxis(m+1,1), 'YLim', [-singlePlotScale singlePlotScale] );
    end
    set( hAxis(m+1,1), 'YLim', get( hAxis(m+1,1), 'YLim' ) );
    
    % Plot NEV data
    if( ( length( zNEV ) == 1 ) && ( NEVExists == 1 ) )
        NEVWaveformsCh = ...
            (NEVwaveformHeader( zNEV ).sf/1000) * ...
            NEVWaveforms(:,((NEVUnitID<15)&(NEVPacketID==chanIdNEV(zNEV))) ...
            );
        if( ~isempty( NEVWaveformsCh ) )
            NEVWaveformsCh(end+1,:) = NaN; %#ok<AGROW>
            NEVWaveformsTm = ...
                repmat( [ NEVWaveformTime NaN ], 1, size( NEVWaveformsCh, 2 ) );
        else
            NEVWaveformsCh = NaN*ones( size( NEVWaveformTime ) );
            NEVWaveformsTm = NEVWaveformTime;
        end
    else
        NEVWaveformsCh = NaN*ones( size( NEVWaveformTime ) );
        NEVWaveformsTm = NEVWaveformTime;
    end
    set( hLine(m+1,2), ...
        'XData', NEVWaveformsTm(:), ...
        'YData', NEVWaveformsCh(:) ...
        );
    clear NEVWaveformsTm NEVWaveformsCh
    if( (m+1) ~= nRow )
        set( hAxis(m+1,2),'XTickLabel','');
    else
        xlabel( hAxis(m+1,2), 'Time (ms)' );
    end
    axis( hAxis(m+1,2), [ NEVWaveformTime(1) NEVWaveformTime(end) -Inf Inf ] );
    if( max(abs( get( hAxis(m+1,2), 'YLim' ) ) ) > singlePlotScale )
        set( hAxis(m+1,2), 'YLim', [-singlePlotScale singlePlotScale] );
    end
    set( hAxis(m+1,2), 'YLim', get( hAxis(m+1,1), 'YLim' ) );
    
    % Plot NEV threshold if exists in NSx data
     if( ( length( zNEV ) == 1 ) && ( NEVExists == 1 ) )
       set( hLine(m+1,3), ...
            'YData', [ NEVThreshold(zNEV) NEVThreshold(zNEV) ] ...
            );
        set( hLine(m+1,4), ...
            'YData', [ NEVThreshold(zNEV) NEVThreshold(zNEV) ] ...
            );
    end;clear zNEV
    
    % Print out if necessary
    if( n == length( SinglePlotChannelsWork ) )
        % clear unused axes
        for p = (m+1):8
            delete( hAxis(p+1,1) );
            delete( hAxis(p+1,2) );
        end;clear p
        % Set up correct printing process
        if( n <= 8 )
            pflag = 14; % Single page PDF
        else
            pflag = 17; % Last page of multi-page PDF
        end
    else
        if( m == 8 )
            if( n == 8 )
                pflag = 15; % First page of multi-page PDF
            else
                pflag = 16; % Middle page of multi-page PDF
            end
        end
    end
    if( pflag ~= 0 )
        PrintFigure( summaryFig, pflag, filenamePDF );
        % delate any extra axes, such as ones from stampit
        delete( setdiff( findall(summaryFig,'Type','Axes'), allAxes ) );
    end
end;clear n m pflag

%% Clean up
clear NSx tNSx NSxbasicHeader NSxchannelHeader NSxmmf NSxOffset NSxSourceFilename
clear allAxes hAxis hAxisSuplabel hLine spikeFig ylab
clear AllPagesChannels SinglePlotChannels singlePlotScale SinglePlotChannelsIndecesNSx
clear allPlotData
clear ans
clear chanIdNEV chanIdNSx NEVExists
clear contFig
clear dataBlockNSx
clear filenameBase filenameNEV filenameNSX filenamePDF
clear nCol nRow
clear NEVbasicHeader NEVwaveformHeader NEVWaveformSize NEVWaveformTime
clear nTimeSamplesPlot nTimeStampNSxIndeces
clear SizeOfInt16 SizeOfInt32 SizeOfUint16 SizeOfUint32
clear timeRangeEnd timeRangeEndDefault timeRangeStart timeRangeStartDefault
clear timeStartTimeStampNSx timeStartTimeStampNSxIndex timeStopTimeStampNSx timeStopTimeStampNSxIndex
clear singleSampleBin ans

return
