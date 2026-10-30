function fs=inf_slope_fs(c,gamma,phi,slope,depth,mw,gamma_w,z,fs_cap)
% Infinite Slope Factor of Safety
% Description: Computes factor of safety using the infinite slope method
% (same form as the initiation model).

% Compute Stability Using Infinite Slope
fs=(c./(gamma.*depth.*sind(slope)))+((cosd(slope).*tand(phi))./sind(slope)).*(1-mw.*(gamma_w./gamma));

% Make Sure that NaNs are Reflected in FS Raster
fs(isnan(z))=NaN;

% Replace Infinite FS
fs(isinf(fs))=fs_cap;

end
