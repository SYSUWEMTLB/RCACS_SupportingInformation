function [Fields2D,StopInfo] = read_ecom_2dFields(FileToRead,GridInfo,iStep)
% function read_ecom_2dFields a function to read 2D fields obtained from
% ECOM simulations
% Created by Jiatang, SYSU
% Updated by Jiatang on 2014/07/02

Fields2D = [];
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

Fields2D.time = double(mdata(1).time);
Fields2D.EL = double(mdata(1).ELE(:,:));
Fields2D.UA = double(mdata(1).UA(:,:));
Fields2D.VA = double(mdata(1).VA(:,:));
% end

return
end