% Written and developed by                                                %
% Patrizio Graziosi, patrizio.graziosi@warwick.ac.uk, during the          %  
% Marie Curie - Individual Fellowships  GENESIS - project ID 788465       %
% Generic transport simulator for new generation thermoelectric materials %
% ----------------------------------------------------------------------- %
% This file is distributed under the terms of the GNU                     %
% General Public License. See the file `LICENSE' in  the root directory   %
% of the present distribution.                                            %
% ----------------------------------------------------------------------- %
%                                                                         %
% Please cite the code source when publishing results obtained            %
% using the present code                                                  %
%                                                                         %
% ----------------------------------------------------------------------- %


function Data_from_bxsf(fileName,num_of_bands,num_of_points,num_lines_to_skip,num_elements_each_row)

inputFile=strcat(fileName);
fid=fopen(inputFile);


% num_lines_to_skip = 19; % number of line to be skipped in order to get the first ilne of energy values, 19 is the typical values in bsxf files
% skip the lines
for i=1:(num_lines_to_skip-4)
    temp=fgetl(fid);
end
% the a, b, c vectors that ate in from 3 to 1 lines above the number of
% lines to skip
temp = fgetl(fid); 
a = str2num(temp);
temp = fgetl(fid); 
b = str2num(temp);
temp = fgetl(fid); 
c = str2num(temp);
%
temp = fgetl(fid); % thus we arrive to the num_line_to_skip

if rem(num_of_points,num_elements_each_row) ~= 0
    num_lines_to_read = floor(num_of_points/num_elements_each_row)+1;
elseif rem(num_of_points,num_elements_each_row) == 0
    num_lines_to_read = floor(num_of_points/num_elements_each_row);
end

for id_band = 1:num_of_bands
    if id_band == 1
        starting_line = num_lines_to_skip+1;
    else
        starting_line = num_lines_to_skip+(num_lines_to_read+1)*(id_band-1)+1;
    end
    end_line = starting_line+num_lines_to_read-1;

    
    Band_temp=[];
    for id_line = starting_line:end_line
        
        temp = fgetl(fid);
        
        Band_temp = [Band_temp temp];

    end
    BandMatrix(:,id_band) = str2num(Band_temp);
    num_lines_inbetween = 1; % number of line in between two bands, 1 is the typical values in bsxf files
    % skip the lines
    for i=1:num_lines_inbetween
        temp=fgetl(fid);
    end
    
end

%these are only the axes, regardless of the angles between them
% kx_array = ( -cell_vector_width:(2*cell_vector_width/round(num_of_points_kx-1)):cell_vector_width ) * 2*pi/(alat*1e-9) ; 
% ky_array = ( -cell_vector_width:(2*cell_vector_width/round(num_of_points_ky-1)):cell_vector_width ) * 2*pi/(alat*1e-9) ; 
% kz_array = ( -cell_vector_width:(2*cell_vector_width/round(num_of_points_kz-1)):cell_vector_width ) * 2*pi/(alat*1e-9) ; 
%

save E_temp % it contains the BandMatrix, one column each band