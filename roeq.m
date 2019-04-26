function [ro_eq] = roeq(ro1,ro2,d1,d2)
%Calculation of equivalent resistivity of mixed materials coaxial line
%   Function resistivities of materials and correspondent diameters
%   to calculate resistivity of line with single equivalent material
K = ro2/ro1;
ro_eq = ro1*((d1*sqrt(K)+d2)^2/(d1+d2)^2);
end