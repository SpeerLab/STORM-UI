%%%%
% bin2mat(binfile, matfile)
% Wenqin Wang 2009 & 2010
% MList will be one array of (structures of arrays).
% MList will be one array of structures (== struct array) if it's masterlist-only from either
% the mlist itself or the input argument
function MList=bin2mat(binfile, matfile, masterlistonly)



sizeofminfo=72;      % 72 bytes per minfo

fid = fopen(binfile);
fseek(fid,0,'eof'); 
file_length = ftell(fid);% size of the file in bytes
fseek(fid,0,'bof');

% read header
version = char(fread(fid,4,'*char')');% M425
frames = fread(fid,1,'*int32');% number of frames. real frames.
status = fread(fid,1,'*int32');% identified = 2, stormed = 6, traced =3, tracked = 4
header_length = ftell(fid);






nmol=zeros(frames+1,1);  % nmol matrix stores the number of molecules in each frame;
nmol(1)= fread (fid, 1,'int32');% number of molecules in the 0th frame (master list)

fnames = {'x','y','xc','yc','h','a','w','phi','ax','bg','i','c','density',...
    'frame','length','link','z','zc','selfframe'}; %length=index. density = valid = fit iterations
%x, y, xc, yc in pixels.
%z and zc in nanometers

ftypes = {'*single','*single','*single','*single','*single','*single','*single',...
    '*single','*single','*single','*single','*int32','*int32','*int32','*int32',...
    '*int32','*single','*single','*single'}; % an asterisk keeps the precision the same. e.g., '*single' is equiv. to 'single=>single'

ftypes2 = {'single','single','single','single','single','single','single',...
    'single','single','single','single','int32','int32','int32','int32',...
    'int32','single','single', 'single'};
    
lengthfnames=max(size(fnames));
for k=1:lengthfnames
    ftypes2func{k}=str2func(ftypes2{k}); % function handles
end



for f=1:frames
    fseek(fid, sizeofminfo*nmol(f),'cof');
    nmol(f+1)=fread(fid,1,'*int32');
end

nmolcum=cumsum(nmol);

if nmolcum(end)==nmolcum(1) % this means molecule lists in 1, 2, ... frames do not exist
    keepframeinfo=0;
else
    keepframeinfo=1;
end

% the byte offset of the last molecule 
%testoffset= header_length  + (nmolcum(frames)+nmol(frames+1)-1)*sizeofminfo + (frames+1)*4; 


if nargin ==3 && masterlistonly~=0
    keepframeinfo=0;
end


if keepframeinfo==0 % in this case, allocate memory only once!
    
    for k=1:lengthfnames-1 %this will make Matlab allocate memory needed for this big struct array! major time-saving step
        MList(nmol(1)).(fnames{k})=[];
    end
    for index=1:nmol(1)
        fseek (fid, header_length+4 + (index-1)*sizeofminfo, 'bof');
    
        for k = 1:lengthfnames - 1 % disregard selfframe because it's always 0 in this case
            MList(index).(fnames{k}) = fread(fid,1,ftypes{k});
        end
    end
    
else
    
    for k=1:lengthfnames 
        MList(nmol(1)).(fnames{k})=[];
    end
    
    MList = struct;
    for index=1:nmol(1)
        fseek (fid, header_length+4 + (index-1)*sizeofminfo + 14*4, 'bof');
        length = fread(fid,1,'*int32');
%       if ~keepframeinfo
%           length=0;
%       end
        fseek (fid, header_length+4 + (index-1)*sizeofminfo, 'bof');

 %       for k = 1:lengthfnames
 %           MList(index).(fnames{k}) = zeros(length+1,1,ftypes2{k});
 %       end
        
        
        for k = 1:lengthfnames - 1 % disregard selfframe for now
            MList(index).(fnames{k})(length+1) = ftypes2func{k}(0); % to save time, fill the last point with 0
            MList(index).(fnames{k})(1) = fread(fid,1,ftypes{k});
        end
        MList(index).selfframe(length+1) = int32(0); % take care of selfframe
        MList(index).selfframe(1) = int32(0); % take care of selfframe

        fr = MList(index).frame(1);
        lk = MList(index).link(1);

        f=2;
        while lk ~= -1 % link = -1 means there is no "next appearance"
            offset = header_length  + (nmolcum(fr)+lk)*sizeofminfo + (fr+1)*4; % from Insight3: fr is for real. link = 3 means its next appearance is the 4-th molecule in the fr-th frame.
            fseek (fid, double(offset), 'bof');
            for k = 1:lengthfnames - 1 % disregard selfframe for now
                MList(index).(fnames{k})(f) = fread(fid,1,ftypes{k});
            end
            MList(index).selfframe(f) = fr;
            fr = MList(index).frame(f);
            lk = MList(index).link(f);
            f=f+1;
        end
    end
    
    
end





fclose(fid);
save(matfile);
