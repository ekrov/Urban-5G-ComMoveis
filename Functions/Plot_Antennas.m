function [] = Plot_Antennas(p_triang, p_rect, p_Circle, p_mc)
% PLOT_ANTENNAS plots the polar patterns and radiation diagrams of different antennas.
%
% Syntax:
%   [] = Plot_Antennas(p_triang, p_rect, p_Circle, p_mc)
%
% Description:
%   The function Plot_Antennas plots the polar patterns and radiation
%   diagrams of different antennas. It takes the antenna patterns as input
%   and generates visualizations for comparison. The function uses the
%   pattern function to compute the patterns and polarpattern function to
%   create the polar plots.
%
% Input Arguments:
%   - p_triang: Antenna pattern for the triangular antenna.
%   - p_rect: Antenna pattern for the rectangular antenna.
%   - p_Circle: Antenna pattern for the circular antenna.
%   - p_mc: Antenna pattern for the combined rectangular and triangular
%            antenna.
%
% Usage Example:
%   p_triang = ...; % Define the antenna pattern for the triangular antenna
%   p_rect = ...; % Define the antenna pattern for the rectangular antenna
%   p_Circle = ...; % Define the antenna pattern for the circular antenna
%   p_mc = ...; % Define the antenna pattern for the combined rectangular
%                % and triangular antenna
%
%   Plot_Antennas(p_triang, p_rect, p_Circle, p_mc);
%
% See also pattern, polarpattern

figure();
[fmt,~,t10] = pattern(p_triang,2.60e9,0,0:360); 
[fmr,~,t11] = pattern(p_rect,2.60e9,0,0:360); 
[fmc,~,t12] = pattern(p_Circle,2.60e9,0,0:360); 
[fm,~,t10] = pattern(p_mc,2.60e9,0,0:360);

p=polarpattern(t10,fmt,t10,fmr,t10,fmc,t10,fm);
legend('Antena triangular', 'Antena rectangular', 'Antena circular','Antena rect + triang');
title('Polar pattern of antennas');

figure; 
pattern(p_triang,2.60e9);
title('Radiation diagram antena triangular');
end

