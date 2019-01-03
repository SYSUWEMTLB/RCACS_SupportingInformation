function [KMKHRHO,StopInfo] = read_ecom_KMKHRHO2(FileToRead,GridInfo,iStep)
% function read_ecom_KMKHRHO a function to read 3D density fields obtained from
% ECOM simulations
% Created by Jiatang, SYSU

KMKHRHO = [];
% for iStep = iniStep:endStep
try
    m = memmapfile(...
        FileToRead,   ...
        'Format', {   ...
        'uint32'  [1]  'tag1';
        'single'  [1]  'time';
        'uint32'  [1]  'tag2';
        'uint32'  [1]  'tag3';
        'single' [GridInfo.nx GridInfo.ny GridInfo.nz] 'AAM';
        'uint32'  [1]  'tag4';
        'uint32'  [1]  'tag5';
        'single' [GridInfo.nx GridInfo.ny GridInfo.nz] 'KH';
        'uint32'  [1]  'tag6';
        'uint32'  [1]  'tag7';
        'single' [GridInfo.nx GridInfo.ny GridInfo.nz] 'WTSURF';
        'uint32'  [1]  'tag8';
        }, ...
        'Repeat',1,  ...
        'offset',(iStep-1)*4*(3+3*(2+GridInfo.nx*GridInfo.ny*GridInfo.nz)) ... %数据块的大小
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

KMKHRHO.time = double(mdata(1).time);
KMKHRHO.AAM = double(mdata(1).AAM(:,:,:));
KMKHRHO.KH = double(mdata(1).KH(:,:,:));
KMKHRHO.RHO = double(mdata(1).WTSURF(:,:,:));
% end

return
end