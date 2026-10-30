function [slope,aspect,dx,dy] = gradient_prince(Z,x_cellsize,y_cellsize)
% Descendant of the mighty gradient king!
[m,n]=size(Z);
slope=zeros(m,n);
aspect=zeros(m,n);
dx=zeros(m,n);
dy=zeros(m,n);

% Calculate gradient using (dz/dx, dz/dy)
for j=2:m-1
    for k=2:n-1
        dz_dx=(((Z(j-1,k+1))+2*(Z(j,k+1))+(Z(j+1,k+1)))-((Z(j-1,k-1))...
            +2*(Z(j,k-1))+(Z(j+1,k-1))))/(8*x_cellsize);
        dz_dy=(((Z(j+1,k-1))+2*(Z(j+1,k))+(Z(j+1,k+1)))-((Z(j-1,k-1))...
            +2*(Z(j-1,k))+(Z(j-1,k+1))))/(8*y_cellsize);
        
        % Compute slope
        rise_run=sqrt(dz_dx^2+dz_dy^2);
        slope_degrees=atan(rise_run)*57.29578;
%         slope(j,k)=round(slope_degrees);
        slope(j,k)=slope_degrees;
        
        % Compute aspect (http://desktop.arcgis.com/en/arcmap/10.3/tools/spatial-analyst-toolbox/how-aspect-works.htm)
        asp=57.29578*atan2(dz_dy,-dz_dx); % CLOCKWISE from NORTH
%         asp=(atan2d(dz_dy,dz_dx)-90)*-1;
%         asp=asp+(asp<0)*360;

        cell=asp;
        if asp<0
            cell=90-asp;
        elseif asp>90
            cell=360-asp+90;
        else
            cell=90-asp;
        end
        aspect(j,k)=round(cell);
        
        dy(j,k)=dz_dy;
        dx(j,k)=dz_dx;
    end
end

% Convert arrays to single
slope=single(slope);
aspect=single(aspect);
dx=single(dx);
dy=single(dy);

end

