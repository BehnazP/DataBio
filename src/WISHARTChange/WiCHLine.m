function [when,prob] = WiCHLine(files,pol_name,time_name,noL,pval,pol,ROI)
% WICHLINE identifies the changes in a sequence of SAR images line by line
%
% [when,prob] = WiCHLine(files,pol_name,time_name,noL,pval,pol,ROI)
% Input:
% files     -    a cell structure containing all the images path and name
%                Example:
%                {'C:\Download\201712.VV_0'},{'C:\Download\201712.VV_1'}
%                {'C:\Download\201712.VH_0'},{'C:\Download\201712.VV_1'} and
%                so on, the program support geotiff, envi and emi images
% pol_name  -    a char vector containing the name of polarization appear
%                after the name of image following a dot
%                Example: .VV .VH
% time_name -    a vector containing either a sequence of number 0:10 or a
%                personalize numbers such as 062, 064, 067 at the end of the
%                name of the files
% nol       -    number of looks
% pval      -    a significant threshold for p-test
%	pol	      -    9 = 'full',
%                5 = 'azim',
%                4 = 'dual',
%                3 = 'fdiag',
%                2 = 'ddiag',
%                1 = 'single'
% ROI 	    -	   is a region of interest
%
% Output:
% when 	    -	   a 3D matrix (nrow x ncol x ntimes) showing if lnR shows change after each time point for each pixel
% prob   	  -	   the average P-value of the ROI
%
% Behnaz Pirzamanbein
% bepi@dtu.dk
% Image Analysis and Computer Graphics section
% Department of Applied Mathematics and Computer Science
% Technical University of Denmark
% First version 17 July 2018
% last version 2018-11-03

if nargin < 7
    ROI = [];
end

sz = size(files,2);
if isempty(sz)
    error('Exit without selecting any file')
end

sz_time = size(time_name,2);
sz_pol = size(pol_name,2);

wc = struct('P',nan);
ntime_1 = sz_time-1;
pP = [];
count = 0;
h = waitbar(0,'Processing...','Name','WiCHLine');

[pathstr, name, ext] = fileparts(files{1});
switch ext
    case '.tif'
        info = geotiffinfo(files{1});
        nrow = info.Height;
        ncol = info.Width;
        when = nan(nrow,ncol,ntime_1);
        for row  = 1:nrow
            waitbar(row/(nrow),h)
            lines = nan(1,ncol,sz_pol,sz_time);
            for i = 1:sz
                for j = 1:sz_pol
                    if contains(name,pol_name{j})
                        for k = 1:sz_time
                            if endsWith(name,['_',time_name{k}])
                                lines(1,:,j,k) = imread(files{i},'PixelRegion',{[row row],[1,ncol]});
                                break
                            end
                        end
                    end
                end
            end
            for i = 1:ntime_1
                wc(i).P = WishartChange(lines(1,:,:,i:sz_time),noL,pol);
            end
            when(row,:,:) = WiCHWhen(wc,pval);
            if isempty(ROI)
                pP(:,:,row) = WiCHProb(wc);
            else
                if any(ROI(row,:))
                    count = count+1;
                    pP(:,:,count) = WiCHProb(wc,ROI(row,:));
                end
            end
        end
    otherwise
        fname = [name,'.hdr'];
        hdrfile = fullfile(pathstr,fname);
        info = envihdrread(hdrfile);
        nrow = info.lines;
        ncol = info.samples;
        when = nan(nrow,ncol,ntime_1);
        fileINDX = zeros(sz,1);
        position = zeros(sz,1);
        for row  = 1:nrow
            waitbar(row/(nrow),h)
            lines = nan(1,ncol,sz_pol,sz_time);
            for i = 1:sz
                fileINDX(i) = fopen(files{i},'r');
                [pathstr, name, ~] = fileparts(files{i});
                fname = [name,'.hdr'];
                hdrfile = fullfile(pathstr,fname);
                info = envihdrread(hdrfile);
                if (info.data_type == 6 || info.data_type == 9)
                    iscx =1;
                else
                    iscx = 0;
                end
                [precision,machineformat] = envInfo(info);
                hdr_sz = info.header_off_set;
                if row == 1
                    fseek(fileINDX(i),hdr_sz,-1);
                else
                    fseek(fileINDX(i),position(i),-1);
                end
                for j = 1:sz_pol
                    if contains(name,pol_name{j})
                        for k = 1:sz_time
                            if endsWith(name,['_',time_name{k}])
                                 D = fread(fileINDX(i),info.samples,precision,0, machineformat);
                                 position(i) = ftell(fileINDX(i));
                                 if iscx
                                     D = complex(D(1:2:end),D(2:2:end));
                                 end
                                 lines(1,:,j,k) = D;
                                break
                            end
                        end
                    end
                end
                fclose(fileINDX(i));
            end
            for i = 1:ntime_1
                wc(i).P = WishartChange(lines(1,:,:,i:sz_time),noL,pol);
            end
            when(row,:,:) = WiCHWhen(wc,pval);
            if isempty(ROI)
                pP(:,:,row) = WiCHProb(wc);
            else
                if any(ROI(row,:))
                    count = count+1;
                    pP(:,:,count) = WiCHProb(wc,ROI(row,:));
                end
            end
        end
end
prob = mean(pP,3);
fclose('all');
close(h)
end                                   
