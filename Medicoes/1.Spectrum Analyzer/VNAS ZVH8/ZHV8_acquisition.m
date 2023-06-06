% ========================================================================
% C. Móveis 2022/23 - R&D ZVH8 - SPECTRUM ANALYZER
% João R. Reis
% ========================================================================

close all; clear all;  clc; instrreset;
  
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
% GET MEASUREMENT
% ========================================================================
    % stop continuous sweep 
    fprintf(g,['INIT:CONT ' 'OFF']);   

    %%% VNA WAIT FOR SWEEP
    fprintf(g,':INIT:IMM;*WAI');

    % There are two types of measurements we can perform:

        %%% READ MAGNITUDE TRACE - reads all trace (post processing is needed)
        fprintf(g,'FORM ASC;:TRAC? TRACE1'); 
        aux = fscanf(g);
        MAGNITUDE(:) = str2double(split(aux,','));

        %%% SELECT AND READ MAX VALUE - recommended
        fprintf(g,':CALC:MARK:MAX');
        fprintf(g,':CALC:MARK:Y?');
        maxPtxValue = fscanf(g);


% ========================================================================
% PLOT RESULTS
% ========================================================================
     

% Plot in real Time MAGNITUDE TRACE
h1 = figure (1);
p1 = plot(MAGNITUDE', 'linewidth',2);
hold on;
legend show;

% displays in console line the Max Value
disp (['Max. value (dBm): ' num2str(maxPtxValue)]);


fprintf(g,[':INIT:CONT ' 'ON']); % start continuous sweep
fclose(g);
fclose(instrfind);
