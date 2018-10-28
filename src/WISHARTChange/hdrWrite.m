function hdrWrite(fname,nrow,ncol,nband,class,interleave,endian,iscx,isemi)

filehdr = fopen([fname,'.hdr'],'w');

switch class
    case 'uint8';   dt = 1;      % byte
           
    case 'int16';   dt = 2;      % integer
           
    case 'int32';   dt = 3;      % long int

    case 'single'  
        if  ~iscx,  dt = 4;      % float
        else
                    dt = 6;
        end
               
    case 'double'
        if ~iscx,   dt = 5;      % double
        else
                    dt = 9;
        end
    
    case 'uint16';  dt = 12;     % unsigned int
            
    case 'uint32';  dt = 13;     % unsigned long
        
    case 'int64';   dt = 14;     % unsigned long
        
    case 'uint64';  dt = 15;     % unsigned long
             
    otherwise
        error('Data type not recognized');
end

switch isemi
    case 1
        header_offset = 4*ncol;
        NAME = 'EMI';
        ft = 'emi file';
    case 0 
        header_offset = 0;
        NAME = 'ENVI';
        ft = 'ENVI Standard';
end

fprintf(filehdr,'%s \n',NAME);
fprintf(filehdr,'%s \n','description = {');
fprintf(filehdr,'%s \n','Exported from MATLAB}');
fprintf(filehdr,'%s %i \n','samples =',ncol);   
fprintf(filehdr,'%s %i \n','lines   =',nrow); 
fprintf(filehdr,'%s %i \n','bands   =',nband);   
fprintf(filehdr,'%s %i \n','header off set =',header_offset);
fprintf(filehdr,'%s %s \n','file type =',ft);
fprintf(filehdr,'%s %i \n','data type =',dt);
fprintf(filehdr,'%s %s \n','interleave =',interleave);
fprintf(filehdr,'sensor type = Unknown\n');
fprintf(filehdr,'%s %i \n','byte order =',endian);
fprintf(filehdr,'wavelength units = Unknown\n');
fclose(filehdr);

end

