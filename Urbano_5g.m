%%
% SINR Map for a 5G urban Macro-cells test environment

% Define a antena transmissora central
% Define center location site (cells 1-3)

% ------------- Fonte luminosa de leiria 39.74422011285419,
% -8.807632630982194 -------------
centerSite = txsite('Name','MathWorks Glasgow', ...
    'Latitude',39.74422,...
    'Longitude',-8.80763);

% ------------- Criar as antenas adjacentes -------------
% numCellSites = 19;
numCellSites = 7;
siteDistances = zeros(1,numCellSites);
siteAngles = zeros(1,numCellSites);

% Definir distance e angulo entre antenas centrais
isd = 1000; % Inter-site distance
siteDistances(2:7) = isd;
siteAngles(2:7) = 30:60:360;

% % Definir distancia e angulo entre antenas a uma distancia media
% siteDistances(8:13) = 2*isd*cosd(30);
% siteAngles(8:13) = 0:60:300;
% 
% % Definir distancia e angulo entre antenas a uma distancia alta
% 
% siteDistances(14:19) = 2*isd;
% siteAngles(14:19) = 30:60:360;

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
fq = 2.6e9; % Carrier frequency (4 GHz) for Dense Urban-eMBB
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
% viewer = siteviewer("Name","Coverage","Basemap","openstreetmap","Buildings",mapa_osm);
% viewer.Basemap = 'topographic';
% %viewer = siteviewer;
% % Show sites on a map
% show(txs);


%Receiver
%Casa 39.74120486553289, -8.814761942942786
rx = rxsite('Name','Casa carlos', ...
       'Latitude',39.74120, ...
       'Longitude',-8.81476)
% show(rx)


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
% link(rx,txs)
% 
%% Mapa de cobertura do transmissor
% coverage(txs,'longley-rice', ...
%        'SignalStrengths',-100:5:-60)

%viewer = siteviewer("Basemap","openstreetmap","Buildings",mapa_osm);
%viewer.Basemap = 'topographic';
% Show sites on a map
%show(txs);
%% RAIN PROPAGATION MODEL

% pm = propagationModel("longley-rice","TimeVariabilityTolerance",0.7);
% pm_rain=propagationModel('rain','RainRate',50);
% prop_model=pm+pm_rain;
% Mapa de cobertura do transmissor
% coverage(txs,"PropagationModel",pm)

%viewer = siteviewer("Name","RainPropModel","Basemap","openstreetmap","Buildings",mapa_osm);
%viewer.Basemap = 'topographic';

%pm = propagationModel("longley-rice","TimeVariabilityTolerance",0.7);
%pm_rain=propagationModel('rain','RainRate',50);

%prop_model=pm;
%Mapa de cobertura do transmissor
%coverage(txs,"PropagationModel",prop_model)


%% Display radiation pattern
%





%% Antennas

%All these microstrip antennas consisting of PCB with 6.6 mm thickness with dielectric
%constant EpsilionR 4.2, and loss-tangent of 0.02, and square ground plane of 100 mm x 100 mm,
%fed by 1.3 mm diameter coaxial probe, are designed to comply with the ISM band (2.4 - 2.5 GHz).

LGp = 100e-3;       % Ground plane length
WGp = 100e-3;       % Ground plane width
h = 6.6e-3;         % Height of the substrate

%Design a probe feed patchMicrostripElliptical antenna using a dimension of 33.5 mm major axis, 18.8 mm minor axis.
%The feed is offset by 11.6 mm from the origin along the X-axis.

a = 33.5e-3;       % Major Axis       
b = 18.8e-3;        % Minor Axis
f = 11.6e-3;        % Feed Offset     
d = dielectric('EpsilonR',4.2,'LossTangent',0.02);

%Create a patchMicrostripElliptical antenna using the defined parameters.

% p_Ellipse= patchMicrostripElliptical('MajorAxis',a,'MinorAxis',b,...
%     'Height',h,'Substrate',d,'GroundPlaneLength',LGp,'GroundPlaneWidth',... 
%     WGp,'FeedOffset',[-(a/2-f) 0]);
% figure; 
% show(p_Ellipse);

%Definir parametros de uma antenna circular
%Design a probe feed patchMicrostripCircular antenna using a dimension of 16 mm Radius.
%The feed is offset by 9.25 mm from the origin along the X-axis.

r = 16e-3;          % Radius  
f1 = 9.25e-3;       % Feed Offset

%Create a patchMicrostripCircular antenna using the defined parameters.

p_Circle= patchMicrostripCircular('Radius',r,'Height',h,'Substrate',d,...
    'GroundPlaneLength',LGp,'GroundPlaneWidth',WGp,'FeedOffset',[-(r-f1) 0]);
% figure; 
% show(p_Circle);

%Definir parametros de uma antena retangular
%Design a probe-fed rectangular patchMicrostrip antenna using a dimension of 28.20 mm length, 34.06 mm width.
%The feed is offset by 5.3 mm from the origin along the X-axis.

rect1 = 28.20e-3;  % length    
rect2 = 34.06e-3;  % width   
f2 = 5.3e-3;       % feedoffset

%Create a patchMicrostrip rectangular antenna using the defined parameters.

p_rect= patchMicrostrip('Length',rect1,'Width',rect2,'Height',h,'Substrate',d,...
    'GroundPlaneLength',LGp,'GroundPlaneWidth',WGp,... 
    'FeedOffset',[-(rect1/2-f2) 0]);
% figure; 
% show(p_rect);

%DEfinir parametros de antena triangular
%Design an equilateral patchMicrostripTriangular antenna using a dimension of 37.63 mm side.
%The feed is offset by 3.8 mm from the origin along the Y-axis.

side = 37.63e-3;    
f_off = 3.8e-3;     
p_triang= patchMicrostripTriangular('side',side,'Height',h,'Substrate',d,...
    'GroundPlaneLength',LGp,'GroundPlaneWidth',WGp,...
    'FeedOffset',[0 -side/2+f_off]);  
figure(1); 
show(p_triang);
title('Triangular antenna');

%Visualize reflection coefficient magnitute of patches

%Plot the reflection coefficient for these antennas over the band and a reference impedance of 50 ohms.
%Curves for the reflection coefficient magnitude are shown in the below figure.
%Manually mesh of patch antennas with different edge lengths. 

%Visualize Radiation Pattern
%The directivity of the antennas are around 6.37 dB for Elliptical patch, 7 dB for circular patch, 
%7.37 dB for rectangular patch and 6.16 dB for triangular patch.



%Mutually Coupling Between Rectangular and Triangular Patches
%This section is devoted to the study of cases, focusing on configurations with the lowest mutual couplings.
%The relative displacement 'd' is fixed as lambda/2, two patches placed side-by-side configuration and the patches positioned in the center of the rectangular 160 mm x 100 mm ground plane with a substrate of FR4.

LGp1 = 160e-3;       % Ground plane length
WGp1 = 100e-3;       % Ground plane width
Ground_plane1=antenna.Rectangle('Length',LGp1,'Width',WGp1);% Ground plane
dis = 0.0612;        % distance between two patches (d = lambda/2)

%Create a patchMicrostrip rectangular antenna using the defined parameters.

r_ant = pcbStack(p_rect);
rect_p = r_ant.Layers{1};
rect_p.Center = [-dis/2 0];

%Create a patchMicrostripTriangular antenna using the defined parameters. 
t_ant = pcbStack(p_triang);
triangle_p = t_ant.Layers{1};
triangle_p= rotateZ(triangle_p,180);
triangle_p= translate(triangle_p,[dis/2, 0, 0]);
patch = rect_p+triangle_p; % adding patches

%Define PCB Stack
%Use the pcbStack to define the metal and dielectric layers for mutually coupled patch antenna.
%the top-most layer is a patch layer, the second layer is dielectric layer, and the third layer is the ground plane. 

p_mc=pcbStack;
d4=dielectric('EpsilonR',4.2,'Thickness',h,'LossTangent',0.02);
p_mc.BoardThickness=d4.Thickness;
p_mc.BoardShape.Length=LGp1;
p_mc.BoardShape.Width=WGp1;
p_mc.Layers={patch,d4,Ground_plane1};
p_mc.FeedLocations=[-dis/2 -(rect1/2-f2) 1 3; dis/2 -5.3075e-3 1 3];
p_mc.FeedDiameter=1.3e-3;
% figure; 
% show(p_mc);

%Radiation Pattern of Mutually Coupled Patches
%The side-by -side configuration directivity is 7.4 dB

% figure; 
% pattern(p_mc,2.45e9);
% 
% %Pattern Magnitude
% 
% figure;
% [fm,~,t10] = pattern(p_mc,2.45e9,0,0:360); 
% polarpattern(t10,fm);

figure(2);
% [fmt,~,t10] = pattern(p_triang,2.60e9,0,0:360); 
% [fmr,~,t11] = pattern(p_rect,2.60e9,0,0:360); 
% [fmc,~,t12] = pattern(p_Circle,2.60e9,0,0:360); 
% [fm,~,t10] = pattern(p_mc,2.60e9,0,0:360);
% 
% p=polarpattern(t10,fmt,t10,fmr,t10,fmc,t10,fm);
% legend('Antena triangular', 'Antena rectangular', 'Antena circular','Antena rect + triang');
% title('Polar pattern of antennas');

% figure; 
% pattern(p_triang,2.60e9);
% title('Radiation diagram antena triangular');


%%
% ------------- Criar elemento das antenas -------------
% show(txs);
% Define pattern parameters


%azvec = -180:180;
%elvec = -90:90;
%Am = 30; % Maximum attenuation (dB)
%tilt = 0; % Tilt angle
%az3dB = 65; % 3 dB bandwidth in azimuth
%el3dB = 65; % 3 dB bandwidth in elevation

% Define antenna pattern
%[az,el] = meshgrid(azvec,elvec);
%azMagPattern = -12*(az/az3dB).^2;
%elMagPattern = -12*((el-tilt)/el3dB).^2;
%combinedMagPattern = azMagPattern + elMagPattern;
%combinedMagPattern(combinedMagPattern<-Am) = -Am; % Saturate at max attenuation
%phasepattern = zeros(size(combinedMagPattern));

% % Create antenna element
%antennaElement = phased.CustomAntennaElement(...
%    'AzimuthAngles',azvec, ...
%    'ElevationAngles',elvec, ...
%    'MagnitudePattern',combinedMagPattern, ...
%    'PhasePattern',phasepattern);

%antennaElement=p_triang;
% Display radiation pattern

%f = figure(4);
%pattern(antennaElement,fq);

% Define array size
%nrow = 8;
%ncol = 8;

% Define element spacing
%lambda = physconst('lightspeed')/fq;
%drow = lambda/2;
%dcol = lambda/2;

% Define taper to reduce sidelobes 
%dBdown = 30;
%taperz = chebwin(nrow,dBdown);
%tapery = chebwin(ncol,dBdown);
%tap = taperz*tapery.'; % Multiply vector tapers to get 8-by-8 taper values

% Create 8-by-8 antenna array
%cellAntenna = phased.URA('Size',[nrow ncol], ...
%    'Element',p_triang, ...
%    'ElementSpacing',[drow dcol], ...
%    'Taper',tap, ...
%    'ArrayNormal','x');

%f = figure(5);
%pattern(cellAntenna,fq);
%title('Diagrama radiacao antena triangular no array');
% Display radiation pattern
% f = figure(6);
% pattern(patchElement,fq)

% % Design half-wavelength rectangular microstrip patch antenna

%patchElement = design(patchMicrostrip,fq);
%patchElement.Width = patchElement.Length;
%patchElement.Tilt = 90;
%patchElement.TiltAxis = [0 1 0];


%Assign the patch antenna as the array element
%cellAntenna.Element = patchElement;


%f = figure(6);
%pattern(patchElement,fq)


%%
% Trocar a antena usada, de uma antena baseada em equações para um modelo
% de uma antena real

% Launch Site Viewer
% viewer = siteviewer("Basemap","openstreetmap","Buildings",mapa_osm);
% 
% % Show sites on a map
% show(txs);
% viewer.Basemap = 'topographic';
% % Design half-wavelength rectangular microstrip patch antenna
% patchElement = design(patchMicrostrip,fq);
% patchElement.Width = patchElement.Length;
% patchElement.Tilt = 90;
% patchElement.TiltAxis = [0 1 0];
% 
% % Display radiation pattern
% f = figure(6);
% pattern(patchElement,fq)

% Assign the patch antenna as the array element
%cellAntenna.Element = patchElement;


%%
% ------------- Mapa de SINR para uma antena -------------

% viewer = siteviewer("Name","SINR uma antena por cel","Basemap","openstreetmap","Buildings",mapa_osm);
% viewer.Basemap = 'topographic';


% Assign the antenna element for each cell transmitter
for tx = txs
    tx.Antenna = p_triang;
end

% Define receiver parameters using Table 8-2 (b) of Report ITU-R M.[IMT-2020.EVAL] 
bw = 20e6; % 20 MHz bandwidth
rxNoiseFigure = 7; % dB
rxNoisePower = -174 + 10*log10(bw) + rxNoiseFigure;
rxGain = 0; % dBi
rxAntennaHeight = 1.5; % m

% % Display SINR map
% if isvalid(f)
%     close(f)
% end
% %sinr(txs,'freespace', ...
% sinr(txs,'longley-rice', ...
%     'ReceiverGain',rxGain, ...
%     'ReceiverAntennaHeight',rxAntennaHeight, ...
%     'ReceiverNoisePower',rxNoisePower, ...    
%     'MaxRange',isd, ...
%     'Resolution',isd/20)




%%
% Array retangular de 8 por 8 antenas
% E SINR

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
    'Element',p_triang, ...
    'ElementSpacing',[drow dcol], ...
    'Taper',tap, ...
    'ArrayNormal','x');
    
% Display radiation pattern
f = figure(5);
pattern(cellAntenna,fq);
title('Diagrama radiacao antena no array');

% Mostrar mapa oara array antenas 8 por 8
% Launch Site Viewer
viewer = siteviewer("Name","SINR para array 8x8 antenas","Basemap","openstreetmap","Buildings",mapa_osm);
viewer.Basemap = 'topographic';
% Show sites on a map
show(txs);


%cellAntenna.Element =p_triang;

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




% % Display SINR map
% if isvalid(f)
%     close(f)
% end
% sinr(txs,'longley-rice',...
%     'ReceiverGain',rxGain, ...
%     'ReceiverAntennaHeight',rxAntennaHeight, ...
%     'ReceiverNoisePower',rxNoisePower, ...    
%     'MaxRange',isd, ...
%     'Resolution',isd/20)

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

%% 
mapa_osm="map.osm";
viewer = siteviewer("Name","BestTx","Basemap","openstreetmap","Buildings",mapa_osm);
viewer.Basemap = 'topographic';
%Show sites on a map

%Receiver
%Casa 39.74120486553289, -8.814761942942786
rx = rxsite('Name','Casa carlos', ...
       'Latitude',39.74120, ...
       'Longitude',-8.81476)
show(rx)

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
fq = 2.6e9; % Carrier frequency (4 GHz) for Dense Urban-eMBB
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

show(txs);
%Distancia entre antenas

dm = distance(txs,rx) % Unit: m
dkm = dm / 1000

%Angulo entre antenas
azFromEast = angle(txs,rx); % Unit: degrees counter-clockwise from East

azFromNorth = -azFromEast + 90; % Convert angle to clockwise from North

% Potencia recebida 
max_value=-200;

for tx = txs 
    ss = sigstrength(rx,tx);
    if(ss>max_value)
        best_tx=tx;
        max_value=ss;
    end
end

% Link margin ( uma métrica da robustes da comunicação, calculado
% subtraindo a sensividade do recetor e a potencia recebida
margin = abs(rx.ReceiverSensitivity - ss);

% Link de comunicações no mapa
%link(rx,txs)
link(rx,best_tx)

% Mapa de cobertura do transmissor
%coverage(txs,'longley-rice','SignalStrengths',-100:5:-60)

%viewer = siteviewer("Basemap","openstreetmap","Buildings",mapa_osm);
%viewer.Basemap = 'topographic';
% Show sites on a map
%show(txs);

%% HANDOVER
number_steps=17;
latitudes_mobile=[39.739006,39.738719,39.738315,39.737893,39.737529,39.737018,39.736676,39.736125,39.73564,39.735,39.7344,39.73404,39.73364,39.733560,39.733470,39.733526,39.733707,39.7341]
longitudes_mobile=[-8.817274,-8.815051,-8.814523,-8.814,-8.813384,-8.812783,-8.812316,-8.811591,-8.810848,-8.810,-8.809421,-8.808913,-8.8077,-8.807259,-8.80663,-8.806145,-8.805375,-8.8043]


viewer = siteviewer("Name","HandOver","Basemap","openstreetmap","Buildings",mapa_osm);
viewer.Basemap = 'topographic';
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
fq = 2.6e9; % Carrier frequency (4 GHz) for Dense Urban-eMBB
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

%Show sites on a map
show(txs);
%Receiver

handover_pot=zeros(1,14);

for i=1:number_steps

    %Casa 39.74120486553289, -8.814761942942786
    rx = rxsite('Name','Carro 1', ...
           'Latitude',latitudes_mobile(i), ...
           'Longitude',longitudes_mobile(i))
    show(rx)

    %Distancia entre antenas

    dm = distance(txs,rx) % Unit: m
    dkm = dm / 1000

    %Angulo entre antenas
    azFromEast = angle(txs,rx); % Unit: degrees counter-clockwise from East

    azFromNorth = -azFromEast + 90; % Convert angle to clockwise from North

    % Potencia recebida 
    max_value=-200;

    tx_to_test=[txs(13), txs(16)];
    j=0;
    for tx = tx_to_test 
        j=j+1;
        ss = sigstrength(rx,tx);
        handover_pot(i,j)=ss;
        if(ss>max_value)
            best_tx=tx;
            max_value=ss;
            
            %handover_pot(i)=max_value;
            
        end
    end
    
    % Link margin ( uma métrica da robustes da comunicação, calculado
    % subtraindo a sensividade do recetor e a potencia recebida
    margin = abs(rx.ReceiverSensitivity - ss);

    % Link de comunicações no mapa
    %link(rx,txs)
    link(rx,best_tx)

    % Mapa de cobertura do transmissor
    %coverage(txs,'longley-rice','SignalStrengths',-100:5:-60)

    %viewer = siteviewer("Basemap","openstreetmap","Buildings",mapa_osm);
    %viewer.Basemap = 'topographic';
    % Show sites on a map
    %show(txs);
end
f = figure(10);
title("Handover - mobile receptor");
x=1:1:number_steps;
plot(x,handover_pot(:,1))

%% Simulacao ESTG

mapa_osm="map_estg.osm";    
viewer = siteviewer("Name","Simulacao ESTG","Basemap","openstreetmap","Buildings",mapa_osm);
viewer.Basemap = 'topographic'

% Define transmitter parameters using Table 8-2 (b) of Report ITU-R M.[IMT-2020.EVAL]
fq = 2.6e9; % Carrier frequency (4 GHz) for Dense Urban-eMBB
antHeight = 3; % m
G=5.9;
txPowerDBm = 0 + G; % Total transmit power in dBm
txPower = 10.^((txPowerDBm-30)/10); % Convert dBm to W

cellName='Sim antena';
cellAngles=0;

cellLat=39.735220;
cellLon=-8.820404;

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
el3dB = 70; % 3 dB bandwidth in elevation

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

%f = figure(4);
%pattern(antennaElement,fq);

%f = figure(5);
%[fmt,~,t10]=pattern(antennaElement,fq,0,0:360);

%p=polarpattern(t10,fmt);
%patternAzimuth(antennaElement,fq)
%figure;
%[fmt,~,t10] = pattern(p_triang,2.60e9,0,0:360); 

rtpm = propagationModel("raytracing", ...
    "Method","sbr", ...
    "MaxNumReflections",0, ...
    "BuildingsMaterial","perfect-reflector", ...
    "TerrainMaterial","perfect-reflector");

coverage(tx,rtpm, ...
    "SignalStrengths",-120:1:txPower, ...
    "MaxRange",250, ...
    "Resolution",3, ...
    "Transparency",0.6,...
    "SignalStrengths",-120:0.1:txPower,...
    'ColorLimits',[-120,txPower])



%Receiver
%Casa 39.74120486553289, -8.814761942942786
rx = rxsite('Name','Carlos', ...
       'Latitude',39.735213, ...
       'Longitude',-8.819667)
show(rx)

% Link margin ( uma métrica da robustes da comunicação, calculado
% subtraindo a sensividade do recetor e a potencia recebida
ss = sigstrength(rx,tx);
margin = abs(rx.ReceiverSensitivity - ss);

% Link de comunicações no mapa
%link(rx,txs)
link(rx,tx)


%show(tx);
%% Medicoes das antenas
% ========================================================================
% Get file from measurement file
% ========================================================================

mapa_osm="map_estg.osm";    
viewer = siteviewer("Name","Medição ESTG","Basemap","openstreetmap","Buildings",mapa_osm);
viewer.Basemap = 'topographic'


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


%% Erro das medições e simulados
% ========================================================================
% Get file from measurement file
% ========================================================================

mapa_osm="map_estg.osm";    
viewer = siteviewer("Name","Medição ESTG","Basemap","openstreetmap","Buildings",mapa_osm);
viewer.Basemap = 'topographic'


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
