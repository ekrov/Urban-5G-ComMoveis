%% Localized coverage simulation (broadband simulation)
% This MATLAB Live Script performs antenna measurements and simulations, comparing the signal strengths obtained from real-world measurements with those predicted by simulation. The script utilizes MATLAB's built-in functions for data processing, visualization, and propagation modeling.

%% Antenna Measurements
% This section of the script performs measurements of antenna signal strengths. It begins by setting up the viewer object and loading the OpenStreetMap file represented by the variable `mapa_osm`. Then, it reads measurement data from the file "exemploFicheiroMediçãoCarlos.txt", which contains latitude, longitude, and signal strength values. The signal strength data is stored in the variable `signalStrength`. The script proceeds to set the transmitter position and height using the specified latitude, longitude, and antenna height values. The transmitter is displayed on the map using the `show` function. Next, the measurement data is processed and arranged in the format required by the site viewer. The processed data is plotted on the map using the `plot` and `contour` functions.

% ========================================================================
% Get file from measurement file
% ========================================================================

mapa_osm="map_estg.osm";    
viewer = siteviewer("Name","Medição ESTG","Basemap","openstreetmap","Buildings",mapa_osm);
viewer.Basemap = 'topographic';

fid = fopen('exemploFicheiroMediçãoCarlos.txt');
C = textscan(fid, '%f%f%f','Delimiter',',');
Latitude=C{1,1};
Longitude=C{1,2};
signalStrength=C{1,3};
% signalStrength= smoothdata(C{1,3},'gaussian',30);                   % filtered data if need be

% ========================================================================
% Set Transmitter Position and Height 
% ========================================================================
lat = 39.735143;
lon =  -8.820529;
tx = txsite("Latitude",lat,...
            "Longitude",lon,...
            "AntennaHeight",3, 'TransmitterFrequency', 2.6e9);
show(tx);
% ========================================================================
% Process data: arrange for siteViewer format
% ========================================================================

tbl = table(Latitude, Longitude, signalStrength);
pd = propagationData(tbl);

% ========================================================================
% PLOT data in site viewer
% ========================================================================
legendTitle = "RSSI" + newline + "(dBm)";                             % plot legend  

plot(pd, "LegendTitle", legendTitle, "Colormap", 'turbo');            % plot points
contour(pd, "LegendTitle", legendTitle, "Colormap", 'turbo');         % plot countour


%% ESTG Simulation
% This section simulates antenna signal strengths at the ESTG location. It starts by setting up the viewer object and loading the OpenStreetMap file represented by the variable `mapa_osm`. The measurement data is read from the file "exemploFicheiroMediçãoCarlos.txt" and stored in the variables `Latitude`, `Longitude`, and `signalStrength`. The transmitter parameters are defined based on the ITU-R M.[IMT-2020.EVAL] specification, including the carrier frequency, antenna height, and transmitter power. The transmitter site is created using the specified parameters, and the antenna pattern is defined using the provided azimuth and elevation values. The transmitter site and antenna pattern are displayed using the `show` and `pattern` functions, respectively. The script calculates the signal strength from the transmitter to multiple receiver locations specified by the latitude and longitude values. The received signal strengths are stored in the variable `sigStrDif`. The measurement data is processed and arranged for the site viewer format, and the resulting data is plotted on the map using the `plot` and `contour` functions.

% ========================================================================
% Get file from measurement file
% ========================================================================

mapa_osm="map_estg.osm";  

fid = fopen('exemploFicheiroMediçãoCarlos.txt');
C = textscan(fid, '%f%f%f','Delimiter',',');
Latitude=C{1,1};
Longitude=C{1,2};
signalStrength=C{1,3};
% signalStrength= smoothdata(C{1,3},'gaussian',30);                   % filtered data if need be

% ========================================================================
% Set Transmitter Position and Height 
% ========================================================================

% Define transmitter parameters using Table 8-2 (b) of Report ITU-R M.[IMT-2020.EVAL]
fq = 2.6e9; % Carrier frequency (4 GHz) for Dense Urban-eMBB
antHeight = 3; % m
G=5;
txPowerDBm = 0 + G; % Total transmit power in dBm
txPower = 10.^((txPowerDBm-30)/10); % Convert dBm to W

cellName='Sim antena';
cellAngles=-90;

%lat = 39.735143;
%lon =  -8.820529;

cellLat=39.735143;
cellLon=-8.820529;

% Create cell transmitter sites
tx = txsite('Name',cellName, ...
    'Latitude',cellLat, ...
    'Longitude',cellLon, ...
    'AntennaAngle',cellAngles, ...
    'AntennaHeight',antHeight, ...
    'TransmitterFrequency',fq, ...
    'TransmitterPower',txPower);


% Define pattern parameters
azvec = -180:180;
elvec = -90:90;

Am = 7; % Maximum attenuation (dB)

tilt = 0; % Tilt angle
az3dB = 60; % 3 dB bandwidth in azimuth
el3dB = 30; % 3 dB bandwidth in elevation

%Antenna factor é 32.69

% Define antenna pattern
[az,el] = meshgrid(azvec,elvec);
azMagPattern = -12*(az/az3dB).^2;
elMagPattern = -12*((el-tilt)/el3dB).^2;
combinedMagPattern = azMagPattern + elMagPattern;
combinedMagPattern(combinedMagPattern<-Am) = -Am; % Saturate at max attenuation
phasepattern = zeros(size(combinedMagPattern));

antennaElement = phased.CustomAntennaElement(...
    'AzimuthAngles',azvec, ...
    'ElevationAngles',elvec, ...
    'MagnitudePattern',combinedMagPattern, ...
    'PhasePattern',phasepattern);
   
tx.Antenna = antennaElement;


% ========================================================================
% Calculate signal strength
% ========================================================================
number_rx=length(Latitude);

% for i=1:number_rx
% 
%     %Casa 39.74120486553289, -8.814761942942786
%     rx = rxsite('Name','Carro 1', ...
%            'Latitude',Latitude(i,1), ...
%            'Longitude',Longitude(i,1));
%     %show(rx)
%     sigStrDif(i,1) = sigstrength(rx,tx,"raytracing");
%     fprintf('Recetor numero %d / %d feito\n', i, number_rx)
% 
% end

rx = rxsite('Latitude',Latitude, ...
       'Longitude',Longitude);
sigStrDif = sigstrength(rx,tx,"raytracing")';
%fid = fopen('exemploFicheiroMediçãoCarlos.txt');
%C = textscan(fid, '%f%f%f','Delimiter',',');
%Latitude=C{1,1};
%Longitude=C{1,2};
%signalStrength=C{1,3};
signalStrength=sigStrDif;
% signalStrength= smoothdata(C{1,3},'gaussian',30);                   % filtered data if need be

% ========================================================================
% Process data: arrange for siteViewer format
% ========================================================================

tbl = table(Latitude, Longitude, signalStrength);
pd = propagationData(tbl);

% ========================================================================
% PLOT data in site viewer
% ========================================================================
viewer = siteviewer("Name","Simulacao ESTG v2","Basemap","openstreetmap","Buildings",mapa_osm);
viewer.Basemap = 'topographic';
legendTitle = "RSSI" + newline + "(dBm)";                             % plot legend  

plot(pd, "LegendTitle", legendTitle, "Colormap", 'turbo');            % plot points
contour(pd, "LegendTitle", legendTitle, "Colormap", 'turbo');         % plot countour

show(tx);
pattern(tx);
%% Erro between Measurements and Simulations
% This section compares the measured and simulated antenna signal strengths at the ESTG location. It begins by setting up the viewer object and loading the OpenStreetMap file represented by the variable `mapa_osm`. The measurement data is read from the file "exemploFicheiroMediçãoCarlos.txt" and stored in the variables `Latitude`, `Longitude`, and `signalStrength`. The transmitter parameters and antenna pattern are defined as in the previous section. The transmitter site and antenna pattern are displayed using the `show` and `pattern` functions. The script calculates the difference between the simulated and measured signal strengths at multiple receiver locations specified by the latitude and longitude values. The differences are stored in the variable `signalStrength`. The measurement data is processed and arranged for the site viewer format, and the resulting data is plotted on the map using the `plot` and `contour` functions.

% ========================================================================
% Get file from measurement file
% ========================================================================

mapa_osm="map_estg.osm";    
viewer = siteviewer("Name","Medição vs Simulacao ESTG","Basemap","openstreetmap","Buildings",mapa_osm);
viewer.Basemap = 'topographic';


fid = fopen('exemploFicheiroMediçãoCarlos.txt');
C = textscan(fid, '%f%f%f','Delimiter',',');
Latitude=C{1,1};
Longitude=C{1,2};
signalStrength=C{1,3};
% signalStrength= smoothdata(C{1,3},'gaussian',30);                   % filtered data if need be

% ========================================================================
% Set Transmitter Position and Height 
% ========================================================================

% Define transmitter parameters using Table 8-2 (b) of Report ITU-R M.[IMT-2020.EVAL]
fq = 2.6e9; % Carrier frequency (4 GHz) for Dense Urban-eMBB
antHeight = 3; % m
%G=5.9;
G=5; 
txPowerDBm = 0 + G;% Total transmit power in dBm
txPower = 10.^((txPowerDBm-30)/10); % Convert dBm to W

cellName='Sim antena';
cellAngles=-90;

cellLat=39.735143;
cellLon=-8.820529;

% Create cell transmitter sites
tx = txsite('Name',cellName, ...
    'Latitude',cellLat, ...
    'Longitude',cellLon, ...
    'AntennaAngle',cellAngles, ...
    'AntennaHeight',antHeight, ...
    'TransmitterFrequency',fq, ...
    'TransmitterPower',txPower);


% Define pattern parameters
azvec = -180:180;
elvec = -90:90;

Am = 7; % Maximum attenuation (dB)

tilt = 0; % Tilt angle
az3dB = 60; % 3 dB bandwidth in azimuth
el3dB = 30; % 3 dB bandwidth in elevation

%Antenna factor é 32.69

% Define antenna pattern
[az,el] = meshgrid(azvec,elvec);
azMagPattern = -12*(az/az3dB).^2;
elMagPattern = -12*((el-tilt)/el3dB).^2;
combinedMagPattern = azMagPattern + elMagPattern;
combinedMagPattern(combinedMagPattern<-Am) = -Am; % Saturate at max attenuation
phasepattern = zeros(size(combinedMagPattern));

antennaElement = phased.CustomAntennaElement(...
    'AzimuthAngles',azvec, ...
    'ElevationAngles',elvec, ...
    'MagnitudePattern',combinedMagPattern, ...
    'PhasePattern',phasepattern);
   
tx.Antenna = antennaElement;

show(tx);
pattern(tx);
% ========================================================================
% Calculate difference in covid
% ========================================================================
number_rx=length(Latitude);
% 
% for i=1:number_rx
% 
%     Casa 39.74120486553289, -8.814761942942786
%     rx = rxsite('Name','Carro 1', ...
%            'Latitude',Latitude(i,1), ...
%            'Longitude',Longitude(i,1));
%     show(rx)
%     sigStrDif(i,1) = sigstrength(rx,tx);
%     fprintf('Recetor numero %d / %d feito\n', i, number_rx)
% 
% end

rx = rxsite('Latitude',Latitude, ...
       'Longitude',Longitude);
   
%    
% rtpm = propagationModel("raytracing", ...
%     "Method","sbr", ...
%     "MaxNumReflections",1, ...
%     "BuildingsMaterial","perfect-reflector", ...
%     "TerrainMaterial","perfect-reflector");   
   
sigStrDif = sigstrength(rx,tx,"raytracing")';
%sigStrDif = sigstrength(rx,tx,rtpm)';


%fid = fopen('exemploFicheiroMediçãoCarlos.txt');
%C = textscan(fid, '%f%f%f','Delimiter',',');
%Latitude=C{1,1};
%Longitude=C{1,2};
%signalStrength=C{1,3};
signalStrength=sigStrDif-signalStrength;
% signalStrength= smoothdata(C{1,3},'gaussian',30);                   % filtered data if need be

% ========================================================================
% Set Transmitter Position and Height 
% ========================================================================
lat = 39.735143;
lon =  -8.820529;
%tx = txsite("Latitude",lat,...
%            "Longitude",lon,...
%            "AntennaHeight",3, 'TransmitterFrequency', 2.6e9);
%show(tx);



% ========================================================================
% Process data: arrange for siteViewer format
% ========================================================================

tbl = table(Latitude, Longitude, signalStrength);
pd = propagationData(tbl);

% ========================================================================
% PLOT data in site viewer
% ========================================================================
legendTitle = "RSSI" + newline + "(dBm)";                             % plot legend  

plot(pd, "LegendTitle", legendTitle, "Colormap", 'turbo');            % plot points
contour(pd, "LegendTitle", legendTitle, "Colormap", 'turbo');         % plot countour


