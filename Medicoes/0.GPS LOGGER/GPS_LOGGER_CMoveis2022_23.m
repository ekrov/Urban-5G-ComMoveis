% ========================================================================
% C. Móveis 2022/23 - AARONIA GPS LOGGER - João R. Reis
% thanks to prof. Nuno Leonor for making this code available!
% ========================================================================
clc
clear
% ========================================================================
% Clear existing connections
% ========================================================================
if(~isempty(instrfind))
    fclose(instrfind);
end

% =======================================================================
% GPS Device Configurations
% =======================================================================
GPS_DEV.COM_PORT = 'COM8';             % <<< -- Define your COM Port (There are 2 COM Ports - One for GPS other for Accelerometer)
GPS_DEV.BAUDRATE = 625000;
GPS_DEV.TERMINATOR = 'CR/LF';
GPS_DEV.StopBits = 2;

% ========================================================================
% Open GPS Device
% ========================================================================
GPS_DEV.INSTRUMENT = serial(GPS_DEV.COM_PORT);
GPS_DEV.INSTRUMENT.Baudrate = GPS_DEV.BAUDRATE;
GPS_DEV.INSTRUMENT.Terminator = GPS_DEV.TERMINATOR;
GPS_DEV.INSTRUMENT.StopBits = GPS_DEV.StopBits;
fopen(GPS_DEV.INSTRUMENT);
fprintf(GPS_DEV.INSTRUMENT,'$PAAG,MODE,START');


% ========================================================================
% Config Spectral analyser
% ========================================================================
% ========================================================================
% OPEN VNA via USB
% ========================================================================
    
%%% FOR VNA VHZ8 R&S
%Open communication port
g = visa('RS','tcpip::172.16.10.10::inst0::INSTR'); %USB
%vna_port = visa('RS','tcpip0::192.168.10.2::inst0::INSTR');  %Ethernet
g.InputBufferSize = 22000;
fopen(g);

% ========================================================================
% CONFIGURE VNA
% ========================================================================

% CONFIGURE VNA

fprintf(g,[':FREQ:CENT ' '94MHz']);
fprintf(g,[':FREQ:SPAN ' '10MHz']);

fprintf(g,[':BAND:RES ' '100kHz']);
fprintf(g,[':BAND:VID ' '3kHz']);
% fprintf(g,[':DISP:WIND:TRAC:Y:RLEV' ref_level]);
    
% ========================================================================
% Start measurement
% ========================================================================

fileName = '5_fov_real.txt';
tic
while(1)
    tic
    % ========================================================================
    % GPS READING
    % ========================================================================
    GPS_String=fscanf(GPS_DEV.INSTRUMENT);
    
    if(contains(GPS_String,'GPGGA'))
        [GPS_READING] = Process_GPS_String(GPS_String);
        
        %%% NOTE %%%
        % SPECTRAL ANALYSER
        % stop continuous sweep 
        fprintf(g,['INIT:CONT ' 'OFF']);   

        %%% VNA WAIT FOR SWEEP
        fprintf(g,':INIT:IMM;*WAI');

        % There are two types of measurements we can perform:

        %%% READ MAGNITUDE TRACE - reads all trace (post processing is needed)
        %fprintf(g,'FORM ASC;:TRAC? TRACE1'); 
        %aux = fscanf(g);
        %MAGNITUDE(:) = str2double(split(aux,','));

        %%% SELECT AND READ MAX VALUE - recommended
        fprintf(g,':CALC:MARK:MAX');
        fprintf(g,':CALC:MARK:Y?');
        maxPtxValue = fscanf(g);
        ptxdBm =maxPtxValue
        
        % ptxdBm = write here a function to read from Spectrum Analyser
        %ptxdBm = 10*rand(); % dummy power

        %%% Display in Console line %%%
        disp([['Elapsed time ' num2str(toc)] 's,' ' Lat: ' num2str(GPS_READING.LAT) ' Long: ' num2str(GPS_READING.LON) ', Power: ' num2str(ptxdBm)]);
        
        %%% Save Data To File %%%
        fileID=fopen(fileName,'a');
        fprintf(fileID,'%s, %s, %s\n',GPS_READING.LAT, GPS_READING.LON, ptxdBm);
        fclose(fileID);
        
    end
end



