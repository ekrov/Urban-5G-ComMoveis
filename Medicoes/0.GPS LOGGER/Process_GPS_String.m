function [GPS_READING] = Process_GPS_String(GPS_String)

GPS_READING.STRING = GPS_String;
GPS_READING.MATLAB_TIME_STAMP = now;

Aux_Indexes = find(GPS_String==',');

%% TIME_STAMP
if( (Aux_Indexes(2)-Aux_Indexes(1)) > 1)
    GPS_READING.TIME_STR    =            GPS_String(Aux_Indexes(1)+1:Aux_Indexes(2)-1);
else
    GPS_READING.TIME_STR    = NaN;
end


%% LATITUDE VALUE
if( (Aux_Indexes(2)-Aux_Indexes(1)) > 1)
    GPS_READING.LAT         = str2double(GPS_String(Aux_Indexes(2)+1:Aux_Indexes(2)+2)) + str2double(GPS_String(Aux_Indexes(2)+3:Aux_Indexes(3)-1))/60;
else
    GPS_READING.LAT         = NaN;
end

%% LATITUDE LABEL
if( (Aux_Indexes(2)-Aux_Indexes(1)) > 1)
    GPS_READING.LAT_LABEL   =            GPS_String(Aux_Indexes(3)+1:Aux_Indexes(4)-1);
else
    GPS_READING.LAT_LABEL   = NaN;
end

%% LONGITUDE VALUE
if( (Aux_Indexes(2)-Aux_Indexes(1)) > 1)
%     GPS_READING.LON         = str2double(GPS_String(Aux_Indexes(4)+1:Aux_Indexes(5)-1))/100;
    GPS_READING.LON         = str2double(GPS_String(Aux_Indexes(4)+1:Aux_Indexes(4)+3)) + str2double(GPS_String(Aux_Indexes(4)+4:Aux_Indexes(5)-1))/60;
else
    GPS_READING.LON         = NaN;
end

%% LONGITUDE LABEL
if( (Aux_Indexes(2)-Aux_Indexes(1)) > 1)
    GPS_READING.LON_LABEL   =            GPS_String(Aux_Indexes(5)+1:Aux_Indexes(6)-1);
else
    GPS_READING.LON_LABEL   = NaN;
end

%% Process Coords Labels
if(GPS_READING.LAT_LABEL == 'S')
    GPS_READING.GPS_LAT = GPS_READING.LAT * (-1);
end
if(GPS_READING.LON_LABEL == 'W')
    GPS_READING.LON = GPS_READING.LON * (-1);
end

end