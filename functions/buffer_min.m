function val_min=buffer_min(A,rows,cols,r_pix)
% Minimum Value Within Circular Buffer
% Description: Finds the minimum value of raster A within a circular buffer
% (radius r_pix pixels) around each point (row, col). Returns NaN for
% points outside of the raster or buffers containing only NaNs. If r_pix
% is 0, the value at the single pixel is returned.
%
% Inputs:
%   A     - raster (e.g., factor of safety)
%   rows  - row of each point
%   cols  - column of each point
%   r_pix - buffer radius (pixels)

% Create Circular Window Offsets
[dr,dc]=ndgrid(-r_pix:r_pix,-r_pix:r_pix);
mask=(dr.^2+dc.^2)<=r_pix^2;
dr=dr(mask)';
dc=dc(mask)';

% Rows and Columns of All Pixels in Each Buffer (Points x Pixels in Window)
rows_buff=rows(:)+dr;
cols_buff=cols(:)+dc;

% Remove Pixels Outside of Raster
in_raster=rows_buff>=1 & rows_buff<=size(A,1) & cols_buff>=1 & cols_buff<=size(A,2);

% Extract Values Within Buffer
vals=nan(size(rows_buff));
vals(in_raster)=A(sub2ind(size(A),rows_buff(in_raster),cols_buff(in_raster)));

% Minimum Value Within Buffer
val_min=min(vals,[],2,'omitnan');

end
