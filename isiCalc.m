function [ isiArray ] = isiCalc( spikeTimes )
%UNTITLED2 Summary of this function goes here
%   Inputs :
%           spikeTimes

% Calculate ISI
isiArray = diff( spikeTimes );
end

