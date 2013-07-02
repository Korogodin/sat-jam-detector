function [ECI, V_ECI] = ECEF2ECI( ECEF, GST, V_ECEF )
% Note: Richard Rieber's code in base of this function

%Checking number of inputs for errors
if nargin < 2 || nargin > 3
    error('Incorrect number of inputs.  See help eci2ecef.')
end

[b,n] = size(ECEF);
%Checking to see if length of ECEF matrix is the same as the length of the GST vector
if n ~= length(GST)
    error('Size of ECEF vector not equal to size of GST vector.  Check inputs.')
end

%Checking to see if the ECEF vector has 3 elements
if b ~= 3
    error('ECEF vector must have 3 elements to the vector (X,Y,Z) coordinates')
end

ECI = zeros(3,n);

if nargin == 3
    V_ECI = zeros(3,n);
end

const_OrbitConverter;
for j = 1:n  %Iterating thru the number of positions provided by user
    % Rotating the ECEF vector into the ECEF frame via the GST angle about the Z-axis
    ECI(:,j) = U3(-GST(j))*ECEF(:,j);
    dT = 1;
    ECI2(:,j) = U3(-(GST(j) + dT*w_earth))*(ECEF(:,j) + V_ECEF*dT);
        
    if nargin == 3
%         V_ECI(:,j) = V_ECEF(:,j) + cross(ECEF(:,j), [0; 0; w_earth]);
        V_ECI(:,j) = ECI2(:,j) - ECI(:, j);
    end
end

end

