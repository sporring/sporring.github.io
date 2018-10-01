function result = mapRegions(map,regions);

% result = mapRegions(map,regions);
% mapRegions finds new values for regions using the map.
% the map is assumed to be a many-to-one
% that is one value per column
% the new value is founded by using the old value as 
% col index in the map.

% example mapRegions(links4to1, regions1)
% gives the same number to all regions at level 1
% which maps to the same region at level 4

%                         November 1997, Ole Fogh Olsen

[Y I] = max(map);

result = I(regions);



