function [D,info]=emiread(datafile,hdrfile)
%emi2Envi Reads emi image file.
%[D,INFO]=emi2EnviSAR(DATAFILE,HDRFILE) provides images data (D) and a
%structure (INFO) whose fields correspond to header items.
%[D,INFO]=ENVIREAD(DATAFILE) assumes the header file is named
%"DATAFILE.hdr" and exists in the same directory.
%D will be in whatever number format (double, int, etc.) as in the emi
%file.

%Original version by Ian Howat, Ohio State Universtiy, ihowat@gmail.com
%Heavily modified by Behnaz Pirzamanbein, bepi@dtu.dk.

if nargin < 1
    error('You must specify at least one input');
elseif nargin <2
    hdrfile=[deblank(datafile),'.hdr']; %implicit name
end
info  = emihdrread(hdrfile);
D = emidataread(datafile,info);


