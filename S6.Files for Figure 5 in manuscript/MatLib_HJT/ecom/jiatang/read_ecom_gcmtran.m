function [TS,StopInfo] = read_ecom_gcmtran(FileToRead,GridInfo,iStep)

TS = [];
% for iStep = iniStep:endStep
try
    m = memmapfile(...
        FileToRead,   ...
        'Format', {   ...
        'uint32'  [1]  'tag1';
        'single'  [1]  'time';
        'uint32'  [1]  'tag2';
        'uint32'  [1]  'tag3';
        'single' [GridInfo.nx GridInfo.ny GridInfo.nz-1] 'T';
        'uint32'  [1]  'tag4';
        'uint32'  [1]  'tag5';
        'single' [GridInfo.nx GridInfo.ny GridInfo.nz-1] 'S';
        'uint32'  [1]  'tag6';
        'uint32'  [1]  'tag7';
        'single' [GridInfo.nx GridInfo.ny GridInfo.nz-1] 'U';
        'uint32'  [1]  'tag8';
        'uint32'  [1]  'tag9';
        'single' [GridInfo.nx GridInfo.ny GridInfo.nz-1] 'V';
        'uint32'  [1]  'tag10';
        'uint32'  [1]  'tag11';
        'single' [GridInfo.nx GridInfo.ny] 'ES';
        'uint32'  [1]  'tag12';
        'uint32'  [1]  'tag13';
        'single' [GridInfo.nx GridInfo.ny] 'ED';
        'uint32'  [1]  'tag14';
        }, ...
        'Repeat',1,  ...
        'offset',(iStep-1)*4*(3+4*(2+GridInfo.nx*GridInfo.ny*(GridInfo.nz-1))+2*(2+GridInfo.nx*GridInfo.ny)) ... %数据块的大小
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

TS.time = double(mdata(1).time);
TS.T = double(mdata(1).T(:,:,:));
TS.S = double(mdata(1).S(:,:,:));
TS.U = double(mdata(1).U(:,:,:));
TS.V = double(mdata(1).V(:,:,:));
TS.ES = double(mdata(1).ES(:,:));
TS.ED = double(mdata(1).ED(:,:));
% end

return
end