function g=element_geometry(m,k)
g=struct('center',m.center(k,:),'scale',m.scale(k),'volume',m.volume(k),'C',m.C(:,:,k));
end
