%% RF Propagation and Visualization example

%Transmitter
%Biblioteca 39.73406426695488, -8.81966502069402
tx = txsite('Name','ESTG',...
       'Latitude',39.73406, ...
       'Longitude',-8.81966)
   
   
show(tx)

%Receiver
%Casa 39.74120486553289, -8.814761942942786
rx = rxsite('Name','Casa carlos', ...
       'Latitude',39.74120, ...
       'Longitude',-8.81476)
show(rx)

%Distancia entre antenas

dm = distance(tx,rx) % Unit: m
dkm = dm / 1000

%Angulo entre antenas
azFromEast = angle(tx,rx) % Unit: degrees counter-clockwise from East

azFromNorth = -azFromEast + 90 % Convert angle to clockwise from North

% Potencia recebida 
ss = sigstrength(rx,tx)

% Link margin ( uma métrica da robustes da comunicação, calculado
% subtraindo a sensividade do recetor e a potencia recebida
margin = abs(rx.ReceiverSensitivity - ss)

% Link de comunicações no mapa
link(rx,tx)

% Mapa de cobertura do transmissor
coverage(tx,'close-in', ...
       'SignalStrengths',-100:5:-60)
   
% Novo transmissor

[lat,lon] = location(tx,1000,90)

tx2 = txsite('Name','Transmitter2','Latitude',lat,'Longitude',lon,'AntennaHeight',30)

show(tx2)

% Calculate de SINR of the transmitter
% The second transmited signal is interfering with the first transmitted
% one
sinr([tx,tx2])

%% SINR Map for a 5G urban Macro-cells test environment

% Define a antena transmissora central
% Define center location site (cells 1-3)

% ------------- Fonte luminosa de leiria 39.74422011285419,
% -8.807632630982194 -------------
centerSite = txsite('Name','MathWorks Glasgow', ...
    'Latitude',39.74422,...
    'Longitude',-8.80763);

% ------------- Criar as antenas adjacentes -------------
numCellSites = 19;
siteDistances = zeros(1,numCellSites);
siteAngles = zeros(1,numCellSites);

% Definir distance e angulo entre antenas centrais
isd = 1000; % Inter-site distance
siteDistances(2:7) = isd;
siteAngles(2:7) = 30:60:360;

% Definir distancia e angulo entre antenas a uma distancia media
siteDistances(8:13) = 2*isd*cosd(30);
siteAngles(8:13) = 0:60:300;

% Definir distancia e angulo entre antenas a uma distancia alta

siteDistances(14:19) = 2*isd;
siteAngles(14:19) = 30:60:360;

% ------------- Definir parametros das células -------------

% Initialize arrays for cell transmitter parameters
numCells = numCellSites*3;
cellLats = zeros(1,numCells);
cellLons = zeros(1,numCells);
cellNames = strings(1,numCells);
cellAngles = zeros(1,numCells);

% Define cell sector angles
cellSectorAngles = [30 150 270];

% For each cell site location, populate data for each cell transmitter
cellInd = 1;
for siteInd = 1:numCellSites
    % Compute site location using distance and angle from center site
    [cellLat,cellLon] = location(centerSite, siteDistances(siteInd), siteAngles(siteInd));
    
    % Assign values for each cell
    for cellSectorAngle = cellSectorAngles
        cellNames(cellInd) = "Cell " + cellInd;
        cellLats(cellInd) = cellLat;
        cellLons(cellInd) = cellLon;
        cellAngles(cellInd) = cellSectorAngle;
        cellInd = cellInd + 1;
    end
end

% ------------- Criar os transmissores -------------

% Define transmitter parameters using Table 8-2 (b) of Report ITU-R M.[IMT-2020.EVAL]
fq = 4e9; % Carrier frequency (4 GHz) for Dense Urban-eMBB
antHeight = 25; % m
txPowerDBm = 44; % Total transmit power in dBm
txPower = 10.^((txPowerDBm-30)/10); % Convert dBm to W

% Create cell transmitter sites
txs = txsite('Name',cellNames, ...
    'Latitude',cellLats, ...
    'Longitude',cellLons, ...
    'AntennaAngle',cellAngles, ...
    'AntennaHeight',antHeight, ...
    'TransmitterFrequency',fq, ...
    'TransmitterPower',txPower);

% Launch Site Viewer
mapa_osm="map.osm";
viewer = siteviewer("Basemap","openstreetmap","Buildings",mapa_osm);
viewer.Basemap = 'topographic';
%viewer = siteviewer;

% Show sites on a map
show(txs);

%Receiver
%Casa 39.74120486553289, -8.814761942942786
rx = rxsite('Name','Casa carlos', ...
       'Latitude',39.74120, ...
       'Longitude',-8.81476)
show(rx)


%Distancia entre antenas

dm = distance(txs,rx) % Unit: m
dkm = dm / 1000

%Angulo entre antenas
azFromEast = angle(txs,rx); % Unit: degrees counter-clockwise from East

azFromNorth = -azFromEast + 90; % Convert angle to clockwise from North

% Potencia recebida 
ss = sigstrength(rx,txs);

% Link margin ( uma métrica da robustes da comunicação, calculado
% subtraindo a sensividade do recetor e a potencia recebida
margin = abs(rx.ReceiverSensitivity - ss);

% Link de comunicações no mapa
link(rx,txs)

% Mapa de cobertura do transmissor
coverage(txs,'longley-rice', ...
       'SignalStrengths',-100:5:-60)

viewer = siteviewer("Basemap","openstreetmap","Buildings",mapa_osm);
viewer.Basemap = 'topographic';
% Show sites on a map
show(txs);
% ------------- Criar elemento das antenas -------------

% Define pattern parameters
azvec = -180:180;
elvec = -90:90;
Am = 30; % Maximum attenuation (dB)
tilt = 0; % Tilt angle
az3dB = 65; % 3 dB bandwidth in azimuth
el3dB = 65; % 3 dB bandwidth in elevation

% Define antenna pattern
[az,el] = meshgrid(azvec,elvec);
azMagPattern = -12*(az/az3dB).^2;
elMagPattern = -12*((el-tilt)/el3dB).^2;
combinedMagPattern = azMagPattern + elMagPattern;
combinedMagPattern(combinedMagPattern<-Am) = -Am; % Saturate at max attenuation
phasepattern = zeros(size(combinedMagPattern));

% Create antenna element
antennaElement = phased.CustomAntennaElement(...
    'AzimuthAngles',azvec, ...
    'ElevationAngles',elvec, ...
    'MagnitudePattern',combinedMagPattern, ...
    'PhasePattern',phasepattern);
   
% Display radiation pattern
f = figure(4);
pattern(antennaElement,fq);

% ------------- Mapa de SINR para uma antena -------------
% Assign the antenna element for each cell transmitter
for tx = txs
    tx.Antenna = antennaElement;
end

% Define receiver parameters using Table 8-2 (b) of Report ITU-R M.[IMT-2020.EVAL] 
bw = 20e6; % 20 MHz bandwidth
rxNoiseFigure = 7; % dB
rxNoisePower = -174 + 10*log10(bw) + rxNoiseFigure;
rxGain = 0; % dBi
rxAntennaHeight = 1.5; % m

% Display SINR map
if isvalid(f)
    close(f)
end
%sinr(txs,'freespace', ...
sinr(txs,'longley-rice', ...
    'ReceiverGain',rxGain, ...
    'ReceiverAntennaHeight',rxAntennaHeight, ...
    'ReceiverNoisePower',rxNoisePower, ...    
    'MaxRange',isd, ...
    'Resolution',isd/20)

% Array retangular de 8 por 8 antenas


% Define array size
nrow = 8;
ncol = 8;

% Define element spacing
lambda = physconst('lightspeed')/fq;
drow = lambda/2;
dcol = lambda/2;

% Define taper to reduce sidelobes 
dBdown = 30;
taperz = chebwin(nrow,dBdown);
tapery = chebwin(ncol,dBdown);
tap = taperz*tapery.'; % Multiply vector tapers to get 8-by-8 taper values

% Create 8-by-8 antenna array
cellAntenna = phased.URA('Size',[nrow ncol], ...
    'Element',antennaElement, ...
    'ElementSpacing',[drow dcol], ...
    'Taper',tap, ...
    'ArrayNormal','x');
    
% Display radiation pattern
f = figure(5);
pattern(cellAntenna,fq);

% Mostrar mapa oara array antenas 8 por 8
% Launch Site Viewer
viewer = siteviewer("Basemap","openstreetmap","Buildings",mapa_osm);

% Show sites on a map
show(txs);
viewer.Basemap = 'topographic';

% Assign the antenna array for each cell transmitter, and apply downtilt.
% Without downtilt, pattern is too narrow for transmitter vicinity.
downtilt = 15;
for tx = txs
    tx.Antenna = cellAntenna;
    tx.AntennaAngle = [tx.AntennaAngle; -downtilt];
end

% Display SINR map
if isvalid(f)
    close(f)
end
sinr(txs,'longley-rice', ...
    'ReceiverGain',rxGain, ...
    'ReceiverAntennaHeight',rxAntennaHeight, ...
    'ReceiverNoisePower',rxNoisePower, ...    
    'MaxRange',isd, ...
    'Resolution',isd/20)

% Trocar a antena usada, de uma antena baseada em equações para um modelo
% de uma antena real

% Launch Site Viewer
viewer = siteviewer("Basemap","openstreetmap","Buildings",mapa_osm);

% Show sites on a map
show(txs);
viewer.Basemap = 'topographic';
% Design half-wavelength rectangular microstrip patch antenna
patchElement = design(patchMicrostrip,fq);
patchElement.Width = patchElement.Length;
patchElement.Tilt = 90;
patchElement.TiltAxis = [0 1 0];

% Display radiation pattern
f = figure(6);
pattern(patchElement,fq)

% Assign the patch antenna as the array element
cellAntenna.Element = patchElement;

% Display SINR map
if isvalid(f)
    close(f)
end
sinr(txs,'longley-rice',...
    'ReceiverGain',rxGain, ...
    'ReceiverAntennaHeight',rxAntennaHeight, ...
    'ReceiverNoisePower',rxNoisePower, ...    
    'MaxRange',isd, ...
    'Resolution',isd/20)

% Sumario
%This example shows how to construct a 5G urban macro-cell test environment consisting of a hexagonal network of 19 cell sites, 
%each containing 3 sectored cells. The signal-to-interference-plus-noise ratio (SINR) is visualized on a map for different antennas. 
%The following observations are made:

%A rectangular antenna array can provide greater directionality and therefore peak SINR values than use of a single antenna element.

%The outward-facing lobes on the perimeter of the SINR map represent areas where less interference occurs. A more realistic modelling 
%technique would be to replicate, or wrap around, cell sites to expand the geometry so that perimeter areas experience similar interference as interior areas.

%Using a rectangular antenna array, a propagation model that estimates increased path loss also results in higher SINR values due to less interference.

%Two antenna elements are tried in the antenna array: an equation-based element using Phased Array System Toolbox and a patch antenna element using Antenna Toolbox. 
%These produce similar SINR maps.

% Antenas direcionais

%% Display radiation pattern
%
f = figure(4);
pattern(antennaElement,fq);

f = figure(5);
pattern(cellAntenna,fq);

f = figure(6);
pattern(patchElement,fq)