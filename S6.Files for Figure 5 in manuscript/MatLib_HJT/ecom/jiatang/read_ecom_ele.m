function [ELE,StopInfo] = read_ecom_ele(FileToRead,GridInfo,iStep)
% function read_ecom_ele a function to read water elevations obtained from
% ECOM simulations
% Created by Jiatang, SYSU

ELE = [];
% for iStep = iniStep:endStep
try
    m = memmapfile(...
        FileToRead,   ...
        'Format', {   ...
        'uint32'  [1]  'tag1';
        'single'  [1]  'time';
        'uint32'  [1]  'tag2';
        'uint32'  [1]  'tag3';
        'single' [GridInfo.nx GridInfo.ny] 'ELE';
        'uint32'  [1]  'tag4';
        'uint32'  [1]  'tag5';
        'single' [GridInfo.nx GridInfo.ny] 'UA';
        'uint32'  [1]  'tag6';
        'uint32'  [1]  'tag7';
        'single' [GridInfo.nx GridInfo.ny] 'VA';
        'uint32'  [1]  'tag8';
        }, ...
        'Repeat',1,  ...
        'offset',(iStep-1)*4*(3+3*(2+GridInfo.nx*GridInfo.ny)) ... %数据块的大小
        );
    mdata = m.data;
catch ME
    StopInfo.StopStep = iStep-1;
    StopInfo.EndOfFile = 1;
    StopInfo.ME = ME;
    return;
end
StopInfo.StopStep = iStep;
StopInfo.EndOfFile = 0;

ELE = mdata(1).ELE(:,:);
%
ELE = double(ELE);
% end

return
end