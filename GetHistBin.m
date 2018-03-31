function [ hist_bin_edge, hist_bin_center ] = GetHistBin( hist_bin_center_lo, hist_bin_center_hi, bin_size )
%   GetHistBin takes in the first bin center, the last bin center and bin
%   size, and return hist bin location in format of edges and centers

hist_bin_center = hist_bin_center_lo:bin_size:hist_bin_center_hi;
hist_bin_edge = [hist_bin_center - bin_size/2, hist_bin_center(end) + bin_size/2];
%   the left side of the comma is an array of the lower bound of each bin,
%   while the right side of the comma is only one value -- the higher bound
%   of the last bin.


end

