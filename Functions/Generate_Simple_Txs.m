function [txs] = Generate_Simple_Txs(numCellSites, isd)
% GENERATE_SIMPLE_TXS generates cell transmitter sites for a cellular network.
%
% Syntax:
%   [txs] = Generate_Simple_Txs(numCellSites, isd)
%
% Description:
%   The function Generate_Simple_Txs creates cell transmitter sites for a
%   cellular network based on the given parameters. It computes the
%   locations of adjacent antennas based on the distance and angle from the
%   center site. It then defines the parameters of each cell transmitter
%   and creates the transmitter sites.
%
% Input Arguments:
%   - numCellSites: Number of cell sites in the network.
%   - isd: Inter-site distance between adjacent antennas.
%
% Output Argument:
%   - txs: An array of txsite objects representing the cell transmitter
%     sites.
%
% Usage Example:
%   numCellSites = 7;
%   isd = 500; % meters
%   txs = Generate_Simple_Txs(numCellSites, isd);
%
% References:
%   The implementation of this function is based on the MATLAB 5G Toolbox
%   documentation for txsite and location.
%
% See also txsite, location

centerSite = txsite('Name','MathWorks Glasgow', ...
    'Latitude',39.74422,...
    'Longitude',-8.80763);
% ------------- Criar as antenas adjacentes -------------
siteDistances = zeros(1,numCellSites);
siteAngles = zeros(1,numCellSites);

% Definir distance e angulo entre antenas centrais
siteDistances(2:numCellSites) = isd;
siteAngles(2:numCellSites) = linspace(30, 330, numCellSites - 1);

% Definir distancia e angulo entre antenas a uma distancia media
%siteDistances(8:13) = 2*isd*cosd(30);
%siteAngles(8:13) = 0:60:300;

% Definir distancia e angulo entre antenas a uma distancia alta

%siteDistances(14:19) = 2*isd;
%siteAngles(14:19) = 30:60:360;

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
end

