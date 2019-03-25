function [Kcc_a,Kce_a,Kee_a,Fe,Fc] = area(nelx,nely,E,A,L,flip)

Kt=E/L;
Kf = 0;
Keo = [Kt 0       0        -Kt 0       0
        0   12*Kf   6*L*Kf   0   -12*Kf  6*L*Kf
        0   6*L*Kf  4*L^2*Kf 0   -6*L*Kf 2*L^2*Kf
        -Kt 0       0        Kt  0       0
        0   -12*Kf  -6*L*Kf  0   12*Kf   -6*L*Kf
        0   6*L*Kf  2*L^2*Kf 0   -6*L*Kf 4*L^2*Kf];
Pov = [0 1 0 0 0 0
       1 0 0 0 0 0
       0 0 1 0 0 0
       0 0 0 0 1 0
       0 0 0 1 0 0
       0 0 0 0 0 1];
Kev=Pov'*Keo*Pov;
node_list=1:(2*(nelx+1)+3*(nely-1));
node_dofs_map=reshape((1:3*(2*(nelx+1)+3*(nely-1)))',3,[])';
ELEMENT(1:nelx,:)=[(1:nelx)',(1:nelx)'+1];
ELEMENT(nelx+(1:nelx),:)=[nelx+1+(1:nelx)',(1:nelx)'+nelx+2];
ELEMENT(2*nelx+(1:(nely)),:)=[1,2*nelx+2+1
                              ((2*nelx+2+1):(2*nelx+2+nely-2))',((2*nelx+2+1):(2*nelx+2+nely-2))'+1
                                                    (2*nelx+2+nely-1), nelx+2];
ELEMENT(2*nelx+nely+(1:(nely)),:)=[fix(nelx/2),2*nelx+2+nely
                                                    ((2*nelx+2+nely):(2*nelx+2+2*nely-3))',((2*nelx+2+nely):(2*nelx+2+2*nely-3))'+1
                                                    (2*nelx+2+2*nely-2), nelx+1+fix(nelx/2)];
ELEMENT(2*nelx+2*nely+(1:(nely)),:)=[nelx+1,2*nelx+2+2*nely-2+1
                                                    ((2*nelx+2+2*nely-2+1):(2*nelx+2+3*(nely-1)-1))',((2*nelx+2+2*nely-2+1):(2*nelx+2+3*(nely-1)-1))'+1
                                                    (2*nelx+2+3*(nely-1)), 2*nelx+2];
ELEMENT_DOF=[node_dofs_map(ELEMENT(:,1),:),node_dofs_map(ELEMENT(:,2),:)];
[io,jo,ko]=find(Keo);
[iv,jv,kv]=find(Kev);
Io=reshape(ELEMENT_DOF(1:2*nelx,io)',[],1);
Jo=reshape(ELEMENT_DOF(1:2*nelx,jo)',[],1);
% Ko=repmat(ko,2*nelx,1);
Ko=[repmat(ko,nelx,1)/4;4*repmat(ko,nelx,1)];
Iv=reshape(ELEMENT_DOF(2*nelx+(1:3*nely),iv)',[],1);
Jv=reshape(ELEMENT_DOF(2*nelx+(1:3*nely),jv)',[],1);
Kv=4*repmat(kv,3*nely,1);
K=sparse([Io;Iv],[Jo;Jv],[Ko;Kv],max(node_dofs_map(:)),max(node_dofs_map(:)));
K=(K+K')/2;
F=sparse(3*((nelx+1))+1,1,-1,max(node_dofs_map(:)),1);
coupling_dofs=[(1:3:(3*(nelx+1)-2));2:3:(3*(nelx+1)-1)];
coupling_dofs=coupling_dofs(:);
reduced_dofs=setdiff(node_dofs_map(:),coupling_dofs);

% Total Compliance
u0=zeros(size(K,1),1);
u0(reduced_dofs)=K(reduced_dofs,reduced_dofs)\F(reduced_dofs);
Kcc_a=K(coupling_dofs,coupling_dofs);
Kce_a=K(coupling_dofs,reduced_dofs);
Kee_a=K(reduced_dofs,reduced_dofs);
Fc=F(coupling_dofs);
Fe=F(reduced_dofs);

