function [Vel,StopInfo] = read_ecom_uv(FileToRead,GridInfo,iStep)
% function read_ecom_uv a function to read 3D velocity fields obtained from
% ECOM simulations
% Created by Jiatang, SYSU

Vel = [];
% for iStep = iniStep:endStep
try
    m = memmapfile(...
        FileToRead,   ...
        'Format', {   ...
        'uint32'  [1]  'tag1';
        'single'  [1]  'time';
        'uint32'  [1]  'tag2';
        'uint32'  [1]  'tag3';
        'single' [GridInfo.nx GridInfo.ny GridInfo.nz] 'u';
        'uint32'  [1]  'tag4';
        'uint32'  [1]  'tag5';
        'single' [GridInfo.nx GridInfo.ny GridInfo.nz] 'v';
        'uint32'  [1]  'tag6';
        'uint32'  [1]  'tag7';
        'single' [GridInfo.nx GridInfo.ny GridInfo.nz] 'w';
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
Vel.time = double(mdata(1).time);
Vel.u = double(mdata(1).u(:,:,:));
Vel.v = double(mdata(1).v(:,:,:));
Vel.w = double(mdata(1).w(:,:,:));
% end

return
end