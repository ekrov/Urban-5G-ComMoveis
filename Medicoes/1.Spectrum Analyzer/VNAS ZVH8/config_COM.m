% COM configuration
function [serialport] = config_COM(port_COM)

     s = serial(horzcat('COM',port_COM));        % set COM para windows
   % s = serial('/dev/cu.usbserial-FTEPL004');        % set COM par MAC

    
    s.terminator = 'CR';                          % set Terminator ENTER
    fopen(s);                                   % open COM
    s.RecordDetail = 'verbose';                 % start REc
    s.RecordName = 'reports/serialFile.txt';
    record(s,'off');
    
    serialport = s;
end

