function [status, waveformdata]  = WFPlot( filename, NEVpathinput, unitsinput, pflaginput, Gridinput, NoiseFlag )
%
% function [status, waveformdata]  = WFPlot( filename, NEVpathinput, unitsinput, pflaginput, Gridinput, NoiseFlag );
%
% This plots a subsample of the waveforms associated with a unit
%
% Inputs
%
%    filename   - A character string containing the base filename to be plotted.
%                 Do not include the file extension. If filename is DATA, then it is assumed
%                 that the file DATA.NEV will be been opened and waveforms from this file plotted.
%
%    NEVpathinput - Path of the location the NEV file.
%
%    unitsinput - The units to plot waveforms. If empty or not specified, waveforms will be
%                 plotted for all units in the file. Noise units ignored.
%
%    pflaginput - Flag indicating what is plotted & printed
%                      < 0 - Don't plot anything
%                      = 0 - pause for brief period
%                      = 1 - pause for long period (DEFAULT)
%                      = 2 - pause until user inputs return
%                      = 3 - print on default printer
%                      = 4 - print as JPEG file
%                      = 5 - print as PDF file
%                      = 6 - print as multipage PDF file via Postscript file
%
%    Gridinput  - Grid used for plotting. If vector, then the first entry is the number of rows
%                 and second entry is number of columns. If scalar, then assumes a square grid.
%                 If the number of columns is greater than or equal to the number of row, then plotted landscape.
%                 Otherwise, plotted in tall mode (or portrait mode using all of the page).
%                 Default is [2,3]
%
%   NoiseFlag   - If not equal to zero, then noise units will be plotted.
%                 Otherwise, they are not plotted. Default is not to plot.
%
% Outputs
%
%    status     - Status of the results of the plotting, 1 if ok, 0 or empty if error occurred.
%
%    waveformdata - A structure of the unit number, time, and waveforms
%

% Make sure than the output exists
status = 0;
waveformdata = [];

% Check input arguments
if( nargin < 6 ), NoiseFlag = []; end;
if( nargin < 5 ), Gridinput = []; end;
if( nargin < 4 ), pflaginput = []; end;
if( nargin < 3 ), unitsinput = []; end;
if( nargin < 2 ), NEVpathinput = []; end;
if( nargin < 1 ), warning( 'WFPlot:argumentNumberError', ...
        '2 to 4 inputs required' ); return; end;

% Create local copies of inputs
if( isempty( NoiseFlag ) ), NoiseFlaginput  = 0; else NoiseFlaginput = NoiseFlag;end;
if( isempty( Gridinput ) ), Gridplot  = [2,3]; else Gridplot = Gridinput;end;
if( ( isempty( pflaginput ) ) || ( pflaginput > 6 ) ), pflag = 1; else pflag = pflaginput; end;
if( isempty( unitsinput ) ), unitstotest = []; else unitstotest = reshape( unitsinput, numel(unitsinput), 1 );end;
if( isempty( NEVpathinput ) ), NEVpath = '.'; else NEVpath = NEVpathinput;end;

if( numel(  Gridplot  ) == 1 )
    Gridplot = [ Gridplot(1), Gridplot(1) ];
else
    Gridplot = reshape( Gridplot, [ numel(Gridplot), 1 ] );
    Gridplot = [ Gridplot(1), Gridplot(2) ];
end;
Gridplotrows = Gridplot(1);
Gridplotcols = Gridplot(2);
Gridplottotal = Gridplotrows * Gridplotcols;

% Create filenames
[ foo, filename ] = fileparts( filename );clear foo %#ok<ASGLU>
filenameNEV = forcefilenameextension( ...
    fullfile( NEVpath, filename ), '.nev' );

% Open file
[fid, message] = fopen( filenameNEV, 'rb' );
if( fid == -1 ),
    warning( 'WFPlot:openError', ...
        ['Unable to open file: ' filename ', error message ' message ] );
    return;
end;
clear filenameNEV

% Read headers
[ NEVbasicHeader, NEVwaveformHeader ] = NEVGetHeaders( fid );
if( isempty( NEVbasicHeader ) )
    warning( 'WFPlot:readHeaderError', ...
        ['Unable to read headers in file: ' filename ', error message ' message ] );
    fclose(fid);
    return;
end;

% Calculate skip factors based on data size
wavesize = NEVbasicHeader.datasize - 8;
% tsskipsize = NEVbasicHeader.datasize - 4;
idskipsize = NEVbasicHeader.datasize - 2;
unitskipsize = NEVbasicHeader.datasize - 1;

% Read electrode but not unit
if( fseek( fid, NEVbasicHeader.dataptr+4, 'bof' ) == -1 ),
    warning( 'WFPlot:positionError', ...
        ['Unable to position file, error code ' ferror( fid ) ] );
    fclose(fid);
    return;
end;
idunit = fread( fid, inf, 'uint16', idskipsize );

% Read unit but not electrode
if( fseek( fid, (NEVbasicHeader.dataptr + 6), 'bof' ) == -1 ),
    warning( 'WFPlot:positionError', ...
    ['Unable to position file, error code ' ferror( fid ) ] );
    fclose(fid);
    return;
end;
unit = fread( fid, inf, 'uint8', unitskipsize );
unit = min( 15, unit );

% Combine electrode and unit to get unit number
if( size( unit, 1 ) ~= size( idunit, 1 ) ),
    warning( 'WFPlot:sizeError', ...
    'Sizes in error' );
    fclose(fid);
    return;
end;
idunit = round( 100 * idunit + unit ) / 100;
clear unit

% Select units to plot
unitsinfile = unique( idunit );
unitsinfile = unitsinfile( unitsinfile >= 1.00 );
if( NoiseFlaginput == 0 )
    unitsinfile = unitsinfile( mod( unitsinfile, 1 ) < 0.14 );
end
if( isempty( unitstotest  ) )
    unitstoplot = unitsinfile;
else
    unitstoplot = intersect( unitstotest, unitsinfile );
end;
unitstoplotsize = numel(  unitstoplot );

cnt = 1;

waveformdata = struct( ...
    'unit', cell( unitstoplotsize, 1), ...
    'time', cell( unitstoplotsize, 1), ...
    'waveforms', cell( unitstoplotsize, 1), ...
    'SNR', cell( unitstoplotsize, 1) ...
    );
clf reset;
if( Gridplotcols >= Gridplotrows );
    orient landscape;
    set( gcf, 'PaperUnits', 'inches' );
    set( gcf, 'PaperPosition', [ 0.25 0.5 10.5 7.0 ] );
else
    orient tall;
    set( gcf, 'PaperUnits', 'inches' );
    set( gcf, 'PaperPosition', [ 1.25, 1, 6.5, 9 ] );
end;
for n = 1:unitstoplotsize
    clear z iWFHeader
    z = find( idunit == unitstoplot(n) );
    iWFHeader = find( [NEVwaveformHeader(:).id] == fix( unitstoplot(n) ) );
    subplot( Gridplotrows, Gridplotcols, cnt );
    step = max( fix( size( z, 1 ) / 250 ), 1 );
    sf = 0.001 * NEVwaveformHeader( iWFHeader ).sf;
    if( isempty( sf ) );
        warning( 'WFPlot:scaleFactorError', ...
            [ 'Unable to get scale factor for unit ' num2str( unitstoplot(n), '%6.2f' ) ] );
        sf = 1;
    end;
    if( sf == 0 ); sf = 1; end;
    if( NEVwaveformHeader( iWFHeader ).numberbytes == 0 )
        wavesizech = wavesize;
        datasizecode = 'int8';
    else
        if( NEVwaveformHeader( iWFHeader ).numberbytes == 1 )
            wavesizech = wavesize;
            datasizecode = 'int8';
        else
            
            if( NEVwaveformHeader( iWFHeader ).numberbytes == 2 )
                wavesizech = wavesize / 2;
                datasizecode = 'int16';
            else
                if( NEVwaveformHeader( iWFHeader ).numberbytes == 4 )
                    wavesizech = wavesize / 4;
                    datasizecode = 'int32';
                else
                    warning( 'WFPlot:dataSizeError', ...
                        [ 'Unknown number of bytes for unit ' num2str( unitstoplot(n), '%6.2f' ) ] );
                    wavesizech = wavesize;
                    datasizecode = 'int8';
                end;
            end;
        end;
    end;
    t = 1.0e3 * ( (1:wavesizech) - 1 ) / NEVbasicHeader.SampleRes;
    waveformarray = NaN*ones( wavesizech, ( floor( size(z, 1) / step ) ) );
    wfn = 1;
    for m = 1:step:size(z, 1)
        if( fseek( fid, (NEVbasicHeader.dataptr + (z(m)-1) * NEVbasicHeader.datasize + 4), 'bof' ) == -1 ),
            warning( 'WFPlot:positionError', ...
                ['Unable to position file, error code ' ferror( fid ) ] );
            fclose(fid);
            return;
        end;
        id = fread( fid, 1, 'uint16' );
        unit = fread( fid, 2, 'uint8' );
        idunittest = id + unit( 1 ) / 100;
        if( idunittest ~= unitstoplot( n ) )
            if( ...
                    ( unit( 1 ) == 255 ) && ...
                    ( mod( idunittest, 1 ) >= 0.15 ) && ...
                    ( NoiseFlaginput ~= 0 ) ...
                    )
            else
                warning( 'WFPlot:logicError', ...
                    'Bad index' );
                fclose(fid);
                return;
            end
        end;
        waveform = sf * fread( fid, wavesizech, datasizecode );
        waveformarray(:,wfn)= waveform;
        wfn = wfn + 1;
    end;
    
    waveformdata(n).unit = unitstoplot(n);
    waveformdata(n).time = t;
    waveformdata(n).waveforms = waveformarray;
    waveformdata(n).range = mean( range( waveformdata(n).waveforms, 1 ) );
    waveformdata(n).noise = ...
        ( ...
        std( ...
        reshape( ...
        waveformdata(n).waveforms((end-7):end,:), ...
        (8*size(waveformdata(n).waveforms,2)), 1 ) ) ...
        );% Per RC Kelly (2007) J Neurosci 27:261
    waveformdata(n).SNR = waveformdata(n).range / ...
        ( 2 * waveformdata(n).noise ...
        );% Per RC Kelly (2007) J Neurosci 27:261
    
    plot( waveformdata(n).time, waveformdata(n).waveforms, 'k' );
    if( ...
            ( NEVwaveformHeader( iWFHeader ).lowthreshold ~= 0 ) ...
            && ...
            ( NEVwaveformHeader( iWFHeader ).lowthreshold > -6000 ) ...
            )
        hold on;
        plot( waveformdata(n).time, ...
            sf*NEVwaveformHeader( iWFHeader ).lowthreshold*ones( size( waveformdata(n).time ) ), ...
            'r--', 'linewidth', 1 );
        hold off;
    end
    if( ...
            ( NEVwaveformHeader( iWFHeader ).highthreshold ~= 0 ) ...
            && ...
            ( NEVwaveformHeader( iWFHeader ).highthreshold < 6000 ) ...
            )
        hold on;
        plot( waveformdata(n).time, ...
            sf*NEVwaveformHeader( iWFHeader ).highthreshold*ones( size( waveformdata(n).time ) ), ...
            'r-', 'linewidth', 1 );
        hold off;
    end
    titleleg = [ 'Ch ' int2str( fix( unitstoplot(n) ) ) ...
        ' Unit ' int2str(100*mod(unitstoplot(n),1)) ...
        ' Cnt ' int2str(size(z,1)) ...
        ' SNR ' num2str( waveformdata(n).SNR, '%5.1f' ) ];
    title( titleleg, 'FontSize', 8 );
    xlabel( 'Time (ms)', 'FontSize', 8 );
    ylabel( '\muV', 'FontSize', 8 );
    set( gca, 'FontSize', 8 );
    xlim( [ min(t) max(t) ] );
    
    cnt = cnt + 1;
    if( or( (cnt == (Gridplottotal+1)), (n==unitstoplotsize) ) ),
        
        cnt = 1;
        
        if( Gridplotcols >= Gridplotrows );
            orient landscape;
            set( gcf, 'PaperUnits', 'inches' );
            set( gcf, 'PaperPosition', [ 0.25 0.5 10.5 7.0 ] );
        else
            orient tall;
            set( gcf, 'PaperUnits', 'inches' );
            set( gcf, 'PaperPosition', [ 1.25, 1, 6.5, 9 ] );
        end;
        
        suptitle_withpatch( filename );
        stampit('d');
        
        switch ( pflag );
            
            case 0
                pause( 5 );
                
            case 1
                pause( 5 );
                
            case 2
                disp( 'Press return to continue... ' );
                pause;
                
            case 3
                print;
                
            case 4
                filenameOutput = forcefilenameextension( ...
                    fullfile( NEVpath, ...
                    [ filename '_Waveform_' num2str(n,'%3.3d') ] ), '.jpg' );
                print( '-djpeg', filenameOutput );
                
            case 5
                filenameOutput = forcefilenameextension( ...
                    fullfile( NEVpath, ...
                    [ filename '_Waveform_' num2str(n,'%3.3d') ] ), '.pdf' );
                print( '-dpdf', filenameOutput )
                
            case 6
                if( ~exist( 'filenameOutputPS', 'var' ) )
                    filenameOutputPDF = forcefilenameextension( ...
                        fullfile( NEVpath, ...
                        [ filename '_Waveforms' ] ), '.pdf' );
                    filenameOutputPS = forcefilenameextension( ...
                        fullfile( NEVpath, ...
                        [ filename '_Waveforms' ] ), '.ps' );
                    if( exist( filenameOutputPS, 'file' ) )
                        delete( filenameOutputPS )
                    end
                end
                print( '-dpsc2', '-append', '-r600', '-opengl', filenameOutputPS );
                
        end;% switch on pflag
        clear filenameOutput
        if( n ~= unitstoplotsize )
            clf reset;
        end
    end;% done with page
end;clear n
if( ( pflag == 6 ) && ( unitstoplotsize > 0 ) )
    ps2pdf( 'psfile', filenameOutputPS, 'pdffile', filenameOutputPDF, ...
        'gspapersize', 'letter', 'deletepsfile', 1, 'verbose', 0 );
end

fclose( fid );

status = 1;
return;
