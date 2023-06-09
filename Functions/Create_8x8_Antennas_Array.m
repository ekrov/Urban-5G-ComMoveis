function [cellAntenna] = Create_8x8_Antennas_Array(single_antenna,fq)
% CREATE_8X8_ANTENNAS_ARRAY creates an 8-by-8 rectangular array of antennas.
%
% Syntax:
%   cellAntenna = Create_8x8_Antennas_Array(single_antenna)
%
% Description:
%   The function Create_8x8_Antennas_Array creates an 8-by-8 antenna array
%   with the specified single_antenna element.
%
% Input Arguments:
%   - single_antenna: An instance of a phased antenna element representing
%     a single antenna in the array.
%   - fq: The central frequency used on this simulation
%    
% Output Argument:
%   - cellAntenna: A phased URA object representing the 8-by-8 antenna
%     array.
%
% Usage Example:
%   % Create a single antenna element
%   singleAntenna = phased.CosineAntennaElement;
%
%   % Create the 8-by-8 antenna array
%   arrayAntenna = Create_8x8_Antennas_Array(singleAntenna);
%
% References:
%   The implementation of this function is based on the MATLAB Phased Array
%   System Toolbox documentation for phased.URA.
%
% See also phased.URA

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
    'Element',single_antenna, ...
    'ElementSpacing',[drow dcol], ...
    'Taper',tap, ...
    'ArrayNormal','x');
end

