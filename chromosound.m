function [output] = chromosound(infile)

% assign parameters
time = 'test'; % this is used to find the last time step
colors = 'test'; % used to assign mass colors
condensin_mass_color = 2; % color of condensin beads
DNA_mass_color = [1 3 4]; % color of DNA beads
cohesin_mass_color = 5; % color of cohesin beads
mass_condensin = []; % this is used to look for condensin beads
mass_cohesin = []; % this is used to look for cohesin beads
springs_condensin = []; % this is used to look for springs bound to condensin
springs_cohesin = []; % this is used to look for springs bound to cohesin

% open up the file
fid_in = fopen(infile);

% assign tline so that the lines can be looped through
tline = fgetl(fid_in);

m = 1; % used to assign mass coords
n = 1; % used to assign springs

% loop through all the lines to find the masses and appropriate springs
while ischar(tline)
    % find the mass coordinates
    if size(strfind(tline,'mass '),1) ~= 0
        % split the string into pieces to parse coordinates
        b = strsplit(tline);
        % record the color of each of the beads
        if size(b,2) == 7
            % rows without an eighth entry correpond to a color of 1 (red)
            mass_info(m,1) = 1;
        elseif size(b,2) == 8
            % rows with an eighth entry have their color logged appropriately
            mass_info(m,1) = str2double(b{8});
        else
        end
        mass_info(m,2) = str2double(b{3}); % puts the mass number in the table
        % increase the counter by 1
        m = m+1;
    else
    end
    if size(strfind(tline,'spring '),1) ~= 0
        % split the string into pieces to parse coordinates
        b = strsplit(tline);
        
        % save the desired springs into a cell array
        springs(n,1) = str2double(b{3});
        springs(n,2) = str2double(b{4});
        
        % increase the counter by 1
        n = n+1;
    else
    end
    if size(strfind(tline,'Time '),1) ~= 0
        % assign the time line to the variable time, this will always pick the last one
        if strcmp(time,'test')==1
            time_init = tline;
        else
        end
        
        time = tline;
    else
    end
    % loop to the next line
    tline = fgetl(fid_in);
end

% close the file
fclose('all');

m = 1; % used to track coords over time
mass_coords = []; % set the empty matrix for checks later

% if this is not the original file, it looks for the last time step
if size(strfind(time,'Time '),1) ~= 0
    % open up the file (again)
    fid_in = fopen(infile);
    % assign tline so that the lines can be looped through (starts it over)
    tline = fgetl(fid_in);
    % loop through all the lines
    while ischar(tline)
        % looks for the beginning of the last time step
        if size(strfind(tline,'Time '),1)>0
            % loop to the next line
            tline = fgetl(fid_in);
            % this assigns the coordinates (which are listed backwards for some reason)
            for d = 0:size(mass_info,1)-1
                % split the line into components
                e = strsplit(tline);
                % log all of the coordinates into a cell matrix
                mass_coords(size(mass_info,1)-d,1,m) = str2double(e{1});
                mass_coords(size(mass_info,1)-d,2,m) = str2double(e{2});
                mass_coords(size(mass_info,1)-d,3,m) = str2double(e{3});
                % loop to the next line
                tline = fgetl(fid_in);
            end
            m = m + 1; % increase the counter by 1
        else
            % loop to the next line if it isn't the right one
            tline = fgetl(fid_in);
        end
    end
else
end

% close the file
fclose('all');

% assign counters for the following section
m = 1; % used to assign DNA
n = 1; % used to assign condensin
p = 1; % used to assign cohesin

% separate the condensin into a separate file
for z = 1:size(mass_info)
    if max(mass_info(z,1) == DNA_mass_color) == 1
        % assign everything with DNA colors to the DNA matrix
        mass_DNA(m,1:2) =  mass_info(z,1:2);
        m = m+1;
    elseif max(mass_info(z,1) == condensin_mass_color) == 1
        % assign everything with the condensin color to the condensin matrix
        mass_condensin(n,1:2) =  mass_info(z,1:2);
        n = n+1;
    elseif max(mass_info(z,1) == cohesin_mass_color) == 1
        % assign everything with the cohesin color to the condensin matrix
        mass_cohesin(p,1:2) =  mass_info(z,1:2);
        p = p+1;
    else
    end
end

m = 1; % used to count springs

% create a list of DNA-DNA springs
for z = 1:size(springs,1)
    if max(springs(z,1) == mass_DNA(:,2)) == 1 && max(springs(z,2) == mass_DNA(:,2)) == 1
        % if the spring contains a condensin bead, log it in a table
        springs_DNA(m,:) = springs(z,:);
        % increase the counter by 1
        m = m + 1;
    else
    end
end

m = 1; % used to count springs

% create a list of condensin-condensin springs
if size(mass_condensin,1)>0 % make sure we have condensin
    for z = 1:size(springs,1)
        if max(springs(z,1) == mass_condensin(:,2)) == 1 && max(springs(z,2) == mass_condensin(:,2)) == 1
            % if the spring contains a condensin bead, log it in a table
            springs_condensin(m,:) = springs(z,:);
            % increase the counter by 1
            m = m + 1;
        else
        end
    end
else
end

m = 1; % used to count springs

% create a list of cohesin-cohesin springs
if size(mass_cohesin,1)>0 % make sure we have condensin
    for z = 1:size(springs,1)
        if max(springs(z,1) == mass_cohesin(:,2)) == 1 && max(springs(z,2) == mass_cohesin(:,2)) == 1
            % if the spring contains a condensin bead, log it in a table
            springs_cohesin(m,:) = springs(z,:);
            % increase the counter by 1
            m = m + 1;
        else
        end
    end
else
end

m = 1; % used to count springs

% create a list of DNA_condensin springs
if size(mass_cohesin,1)>0 % make sure we have condensin
    for z = 1:size(springs,1)
        if max(springs(z,1) == mass_DNA(:,2)) == 1 && max(springs(z,2) == mass_condensin(:,2)) == 1
            % if the spring contains a condensin bead, log it in a table
            springs_DNA_condensin(m,:) = springs(z,:);
            % increase the counter by 1
            m = m + 1;
        else
        end
    end
else
end

% loop through the time steps
for z = 1:size(mass_coords,3)
    for n = 1:size(springs_DNA,1) % loop through all DNA springs
        springs_DNA_dist(n,z) = distance_between_3D_chromoshake(...
            [mass_coords(springs_DNA(n,1)+1,1,z) mass_coords(springs_DNA(n,2)+1,1,z)],...
            [mass_coords(springs_DNA(n,1)+1,2,z) mass_coords(springs_DNA(n,2)+1,2,z)],...
            [mass_coords(springs_DNA(n,1)+1,3,z) mass_coords(springs_DNA(n,2)+1,3,z)]);
    end
    for n = 1:size(springs_cohesin,1) % loop through all cohesin springs
        springs_cohesin_dist(n,z) = distance_between_3D_chromoshake(...
            [mass_coords(springs_cohesin(n,1)+1,1,z) mass_coords(springs_cohesin(n,2)+1,1,z)],...
            [mass_coords(springs_cohesin(n,1)+1,2,z) mass_coords(springs_cohesin(n,2)+1,2,z)],...
            [mass_coords(springs_cohesin(n,1)+1,3,z) mass_coords(springs_cohesin(n,2)+1,3,z)]);
    end
    for n = 1:size(springs_condensin,1) % loop through all condensin springs
        springs_condensin_dist(n,z) = distance_between_3D_chromoshake(...
            [mass_coords(springs_condensin(n,1)+1,1,z) mass_coords(springs_condensin(n,2)+1,1,z)],...
            [mass_coords(springs_condensin(n,1)+1,2,z) mass_coords(springs_condensin(n,2)+1,2,z)],...
            [mass_coords(springs_condensin(n,1)+1,3,z) mass_coords(springs_condensin(n,2)+1,3,z)]);
    end
    for n = 1:size(springs_DNA_condensin,1) % loop through all DNA-condensin springs
        springs_DNA_condensin_dist(n,z) = distance_between_3D_chromoshake(...
            [mass_coords(springs_DNA_condensin(n,1)+1,1,z) mass_coords(springs_DNA_condensin(n,2)+1,1,z)],...
            [mass_coords(springs_DNA_condensin(n,1)+1,2,z) mass_coords(springs_DNA_condensin(n,2)+1,2,z)],...
            [mass_coords(springs_DNA_condensin(n,1)+1,3,z) mass_coords(springs_DNA_condensin(n,2)+1,3,z)]);
    end
end

% each row in these tables represents a single spring
% each column in these tables represents a single time step

output.springs_DNA_dist = springs_DNA_dist; % distanced for DNA springs
output.springs_cohesin_dist = springs_cohesin_dist; % distanced for cohesin springs
output.springs_condensin_dist = springs_condensin_dist; % distanced for condensin springs
output.springs_DNA_condensin_dist = springs_DNA_condensin_dist; % distanced for DNA-condensin springs