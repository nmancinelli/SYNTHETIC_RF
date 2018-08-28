function [depth, vpv, vsv, rho] = load_from_card_file ( path2cardfile )
%load_from_card_file Load velocity model from a card file.
%   Just give it the path.
    radius = zeros(1,750);
    vpv = zeros(1,750);
    vsv = zeros(1,750);
    rho = zeros(1,750);
    

    fprintf('  loading %s ... \n', path2cardfile)
    fid = fopen(path2cardfile);

    %read header lines
    for i = 1:3 
        hdr = fgets(fid);
    end

    format = '%f';

    for iline = 1:750;
        nfo = textscan(fid,format,9);
        radius(iline) = nfo{1}(1)/1000.;
        rho(iline)    = nfo{1}(2)/1000.;
        vpv(iline)    = nfo{1}(3)/1000.;
        vsv(iline)    = nfo{1}(4)/1000.;
    end
    
    depth = max(radius) - radius;
    
    zcutoff = 350.;
    
    vpv   = vpv(depth <= zcutoff);
    vsv   = vsv(depth <= zcutoff);
    rho   = rho(depth <= zcutoff);
    depth = depth(depth <= zcutoff);

    fclose(fid);

end

