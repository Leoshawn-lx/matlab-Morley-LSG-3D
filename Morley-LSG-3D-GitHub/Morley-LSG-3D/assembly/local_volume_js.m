function K=local_volume_js(v,g,lambda,iota,mu)
% Exact degree-two integration for the P2 volume bilinear form.
[b,w]=tet_quadrature(2);[~,G,H]=morley_basis3d(b*v,g);
GG=zeros(10,10);for j=1:3,GG=GG+G(:,:,j)'*bsxfun(@times,G(:,:,j),w);end
Hr=reshape(H,10,9);HH=Hr*Hr';K=zeros(30,30);
for a=1:3
    ia=(a-1)*10+(1:10);Ha=reshape(H(:,a,:),10,3);
    for c=1:3
        ic=(c-1)*10+(1:10);Hc=reshape(H(:,c,:),10,3);
        V=mu*(G(:,:,c)'*bsxfun(@times,G(:,:,a),w))+lambda*(G(:,:,a)'*bsxfun(@times,G(:,:,c),w));
        W=mu*(Hc*Ha')+lambda*(Ha*Hc');
        if a==c,V=V+mu*GG;W=W+mu*HH;end
        K(ia,ic)=g.volume*(V+iota^2*W);
    end
end
end
