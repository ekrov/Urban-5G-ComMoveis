
% ========================================================================
% C. Móveis 2022/23 - AARONIA GPS LOGGER - João R. Reis
% thanks to Tiago Oliveira for making this code available!
% ========================================================================

% ========================================================================
% Get file from measurement file
% ========================================================================
fid = fopen('5_fov_real.txt');
C = textscan(fid, '%f%f%f','Delimiter',',');
Latitude=C{1,1};
Longitude=C{1,2};
signalStrength=C{1,3};
% signalStrength= smoothdata(C{1,3},'gaussian',30);                   % filtered data if need be

% ========================================================================
% Set Transmitter Position and Height 
% ========================================================================
lat = 39.82852;
lon =  -8.85835;
tx = txsite("Latitude",lat,...
            "Longitude",lon,...
            "AntennaHeight",3, 'TransmitterFrequency', 2.4e9);
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





