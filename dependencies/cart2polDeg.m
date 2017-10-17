function [ theta, r ] = cart2polDeg( x,y )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
[theta,r] = cart2pol(x,y);
theta = rad2deg(theta);
end

