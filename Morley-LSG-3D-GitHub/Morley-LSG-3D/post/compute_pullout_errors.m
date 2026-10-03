function err = compute_pullout_errors(m,uFull,ref,iota,cfg)
%COMPUTE_PULLOUT_ERRORS Approximation errors for the 3-D Morley solution.
%
% Reports the three direct analogues of the reference-paper quantities:
%   displacement      ||u-u_h||_L2,
%   strain            ||eps(u)-eps_h(u_h)||_L2,
%   strain gradient   ||grad eps(u)-grad_h eps_h(u_h)||_L2,
% plus the manuscript three-bar norm
%   |||u-u_h|||_{iota,h} = sqrt(|e|_1,h^2+iota^2 |e|_2,h^2).

[b,wq]=tet_quadrature(cfg.errorOrder);
nt=size(m.elem,1);

accL2=0; accStrain=0; accGradStrain=0;
accH1=0; accH2=0; accDiv=0; accTrans=0;
refL2=0; refStrain=0; refGradStrain=0; refH1=0; refH2=0;
outsideVol=0; totalVol=0;

for t=1:nt
    verts=m.node(double(m.elem(t,:)),:);
    g=element_geometry(m,t);
    x=b*verts;
    [phi,G,H]=morley_basis3d(x,g);
    [~,fd]=element_dofs(m,t);
    coeff=reshape(uFull(fd(:)),10,3);

    uh=phi*coeff;
    nq=size(x,1);
    guh=zeros(nq,3,3);
    for j=1:3
        guh(:,:,j)=G(:,:,j)*coeff;
    end

    Huh=zeros(3,3,3);
    for c=1:3
        for j=1:3
            for k=1:3
                Huh(c,j,k)=sum(coeff(:,c).*H(:,j,k));
            end
        end
    end

    ex=pullout_reference_fields(x,ref);
    dU=uh-ex.u;
    l2p=sum(dU.^2,2);
    refL2p=sum(ex.u.^2,2);

    h1p=zeros(nq,1); refH1p=zeros(nq,1);
    strainp=zeros(nq,1); refStrainp=zeros(nq,1);
    divp=zeros(nq,1);
    for i=1:3
        divp=divp+guh(:,i,i);
        for j=1:3
            dg=guh(:,i,j)-ex.grad(:,i,j);
            h1p=h1p+dg.^2;
            refH1p=refH1p+ex.grad(:,i,j).^2;

            epsh=0.5*(guh(:,i,j)+guh(:,j,i));
            epse=0.5*(ex.grad(:,i,j)+ex.grad(:,j,i));
            strainp=strainp+(epsh-epse).^2;
            refStrainp=refStrainp+epse.^2;
        end
    end

    h2p=zeros(nq,1); refH2p=zeros(nq,1);
    gradStrainp=zeros(nq,1); refGradStrainp=zeros(nq,1);
    for c=1:3
        for j=1:3
            for k=1:3
                hn=Huh(c,j,k);
                he=reshape(ex.hess(:,c,j,k),nq,1);
                h2p=h2p+(hn-he).^2;
                refH2p=refH2p+he.^2;
            end
        end
    end
    for i=1:3
        for j=1:3
            for k=1:3
                gsn=0.5*(Huh(i,j,k)+Huh(j,i,k));
                gse=0.5*(reshape(ex.hess(:,i,j,k),nq,1)+reshape(ex.hess(:,j,i,k),nq,1));
                gradStrainp=gradStrainp+(gsn-gse).^2;
                refGradStrainp=refGradStrainp+gse.^2;
            end
        end
    end

    transp=uh(:,1).^2+uh(:,2).^2;
    wt=g.volume*wq;
    accL2=accL2+sum(wt.*l2p);
    accStrain=accStrain+sum(wt.*strainp);
    accGradStrain=accGradStrain+sum(wt.*gradStrainp);
    accH1=accH1+sum(wt.*h1p);
    accH2=accH2+sum(wt.*h2p);
    accDiv=accDiv+sum(wt.*divp.^2);
    accTrans=accTrans+sum(wt.*transp);

    refL2=refL2+sum(wt.*refL2p);
    refStrain=refStrain+sum(wt.*refStrainp);
    refGradStrain=refGradStrain+sum(wt.*refGradStrainp);
    refH1=refH1+sum(wt.*refH1p);
    refH2=refH2+sum(wt.*refH2p);

    outsideVol=outsideVol+sum(wt.*double(ex.rActual<ref.r(1)));
    totalVol=totalVol+sum(wt);

    if t==1 || t==nt || mod(t,1000)==0
        fprintf('  error integration %d/%d\n',t,nt);
    end
end

err.L2=sqrt(accL2);
err.strain=sqrt(accStrain);
err.gradStrain=sqrt(accGradStrain);
err.H1=sqrt(accH1);
err.H2=sqrt(accH2);
err.threebar=sqrt(accH1+iota^2*accH2);
err.divUh=sqrt(accDiv);
err.transverseL2=sqrt(accTrans);

err.refL2=sqrt(refL2);
err.refStrain=sqrt(refStrain);
err.refGradStrain=sqrt(refGradStrain);
err.refH1=sqrt(refH1);
err.refH2=sqrt(refH2);
err.refThreebar=sqrt(refH1+iota^2*refH2);
err.relL2=err.L2/max(err.refL2,realmin);
err.relStrain=err.strain/max(err.refStrain,realmin);
err.relGradStrain=err.gradStrain/max(err.refGradStrain,realmin);
err.relThreebar=err.threebar/max(err.refThreebar,realmin);
err.geometryOutsideReferenceVolumeFraction=outsideVol/max(totalVol,realmin);
end
