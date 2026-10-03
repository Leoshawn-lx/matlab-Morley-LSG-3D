function [reduced,full]=element_dofs(m,k)
s=double(m.scalarMap(k,:)); a=double(m.scalarFreeMap(s));a=a(:)';
reduced=[a,a+m.nFreeScalar,a+2*m.nFreeScalar];
reduced(repmat(a==0,1,3))=0;
if nargout>1,full=[s,s+m.scalarNdof,s+2*m.scalarNdof];end
end
