tissue_z=z_interface1;
tissue_x_back=x(x_index)-4e-2;
tissue_x_after=x(x_index)+4e-2;
z_new=tissue_z:dz:zmax5;
x_new=tissue_x_back:dx:tissue_x_after;
y_new=tissue_x_back:dy:tissue_x_after;
xnewindex_back=find(x==tissue_x_back);
xnewindex_after=find(x==tissue_x_after);
znewindex_back=find(z==tissue_z);
znewindex_after=find(z==zmax5);
P3D=p(xnewindex_back:xnewindex_after,xnewindex_back:xnewindex_after,znewindex_back:znewindex_after);
nx_new=length(x_new);
ny_new=length(y_new);

%