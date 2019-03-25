function x=top88_plustruss_mma_cond_internship_by_set(nelx,nely,L,volfrac,tarmass_A1,tarmass_A2,tarmass_A3,tarmass_A4,tarmass_A5,penal,rmin,ft)
filename = 'Topology_Optimization.gif';
% NEED TO DO: 
% 1. OPTIMIZATION BY GROUPING AREAS OF BEAMS
% 2. GUASS-JORDAN ELIMINATION (NEW CONDENSATION TECHNIQUE)
% 3. COMPARISON OF TECHNIQUES
close all
%% MATERIAL PROPERTIES
E0 = 1;
Emin = 1e-9;
nu = 0.3;
P=4;
Sl=1;
rho = 1;
%% PREPARE FINITE ELEMENT ANALccYSIS
A11 = [12  3 -6 -3;  3 12  3  0; -6  3 12 -3; -3  0 -3 12];
A12 = [-6 -3  0  3; -3 -6 -3 -6;  0 -3 -6  3;  3 -6  3 -6];
B11 = [-4  3 -2  9;  3 -4 -9  4; -2 -9 -4 -3;  9  4 -3 -4];
B12 = [ 2 -3  4 -9; -3  2  9 -2;  4  9  2  3; -9 -2  3  2];
% D=E0/(1-nu^2)*[1 nu 0;nu 1 0; 0 0 (1-nu)/2];
% B=1/2*[-1 0 1 0 1 0 -1 0;0 -1 0 -1 0 1 0 1;-1 -1 -1 1 1 1 1 -1];
% DB=D*B;
% B=1/2*[-1 0 1 0 1 0 -1 0;0 -1 0 -1 0 1 0 1;-1 -1 -1 1 1 1 1 -1];
% DB=D*B;
% Cvm=[1 -0.5 0;-0.5 1 0;0 0 3];
% Sel=DB'*Cvm*DB;
KE = 1/(1-nu^2)/24*([A11 A12;A12' A11]+nu*[B11 B12;B12' B11]);
nodenrs = reshape(1:(1+nelx)*(1+nely),1+nely,1+nelx);
edofVec = reshape(2*nodenrs(1:end-1,1:end-1)+1,nelx*nely,1);
edofMat = repmat(edofVec,1,8)+repmat([0 1 2*nely+[2 3 0 1] -2 -1],nelx*nely,1);
iK = reshape(kron(edofMat,ones(8,1))',64*nelx*nely,1);
jK = reshape(kron(edofMat,ones(1,8))',64*nelx*nely,1);
%fixeddofs = union([1:2:2*(nely+1)],[2*(nelx+1)*(nely+1)]);
fixeddofs = [2*(fix(nelx/2)*(nely+1)+1:(nely+1):((nelx)*(nely+1)+1))-1; 2*(fix(nelx/2)*(nely+1)+1:(nely+1):((nelx)*(nely+1)+1))];
fixeddofs = fixeddofs(:);
%% PREPARE FILTER
iH = ones(nelx*nely*(2*(ceil(rmin)-1)+1)^2,1);
jH = ones(size(iH));
sH = zeros(size(iH));
k = 0;
for i1 = 1:nelx
    for j1 = 1:nely
        e1 = (i1-1)*nely+j1;
        for i2 = max(i1-(ceil(rmin)-1),1):min(i1+(ceil(rmin)-1),nelx)
            for j2 = max(j1-(ceil(rmin)-1),1):min(j1+(ceil(rmin)-1),nely)
                e2 = (i2-1)*nely+j2;
                k = k+1;
                iH(k) = e1;
                jH(k) = e2;
                sH(k) = max(0,rmin-sqrt((i1-i2)^2+(j1-j2)^2));
            end
        end
    end
end
H = sparse(iH,jH,sH);
Hs = sum(H,2);
%% INITIALIZE ITERATION
% x = repmat(1,nely,nelx);
% xPhys = x;
% loop = 0;
% change = 1;

%% INITIALIZE ITERATION
x = ones(nely,nelx);
A_max = 2;
xPhys = x;
loop = 0;
change = 1;
m = 3;
n = length(xPhys(:))+2;
epsimin = 0.0000001;
eeen    = ones(n-2,1);
eeem    = ones(m,1);
zeron   = zeros(n,1);
zerom   = zeros(m,1);
xval    = [xPhys(:); 5; 5];
xold1   = xval;
xold2   = xval;
xmin    = zeron;
xmin(end-4:end) = 0.3*ones(2,1);
xmax    = cat(1,eeen, A_max*ones(2,1));
low     = xmin;
upp     = xmax;
C       = 1000*eeem;
d       = 0*eeem;
a0      = 1;
a       = zerom;
outeriter = 0;
maxoutit  = 1500;
kkttol  =0.001;

%%%% The iterations start:
kktnorm = kkttol+10;
% kktnorm = kkttol;
outit = 0;
change=1;
%% START ITERATION
while kktnorm > kkttol && outit < maxoutit && change > 0.01
%    if outit == 0
%        A_old = A_max;
%    else
%        A_old = A;
%    end
   I1 = 100*(sqrt(xval(end-1))^4)/12;
   I2 = 100*(sqrt(xval(end))^4)/12;   
%    I3 = 100*(sqrt(xval(end-2))^4)/12;
%    I4 = 100*(sqrt(xval(end-1))^4)/12;
%    I5 = 100*(sqrt(xval(end))^4)/12;
   
   %Condensation of Structure
   %if A_old - A > 0
        [Kcc,Kce, Kce_pr,Kee, Kee_pr, Fe, Fe_pr,Fc,K_cnd,F_cnd,compliance_truss, u_c, coupling_dofs_truss] = truss_stiffness_condensation_by_set_t(nelx,nely,E0,xval(end-1),xval(end), I1, I2, L);
        
       
     %% Generate the Projection matrix, load vector and stiffness matrix ready for assembly
        coupling_nodes = [nely+1 (nely+1)*fix((nelx+1)/2) (nely+1)*(nelx+1)];
        coupling_dofs = [coupling_nodes*2-1; coupling_nodes*2];
        coupling_dofs = coupling_dofs(:);
        
        coupling_nodes_pr = [nely+1:nely+1:(nely+1)*(nelx+1)];
        coupling_dofs_pr = [coupling_nodes_pr*2-1;coupling_nodes_pr*2];
        coupling_dofs_pr = coupling_dofs_pr(:);
        
        coupling_dofs_truss = sort(coupling_dofs_truss(:));
        
        coupling_dofs_sp1 = repmat(coupling_dofs_truss(:),6);
        coupling_dofs_sp1 = coupling_dofs_sp1(:,1);
        
        coupling_dofs_sp2 = sort(coupling_dofs_sp1);
                
        K_cnd_vec = K_cnd(:);
        K_cnd = sparse(coupling_dofs_sp1',coupling_dofs_sp2',K_cnd_vec,2*(nelx+1),2*(nelx+1));
        
        Pr = sparse(coupling_dofs,[1 2 11 12 21 22],ones(1,6),2*(nelx+1)*(nely+1),2*(nelx+1));
        Pr_pr = sparse(coupling_dofs_pr,1:2*(nelx+1),ones(1,2*(nelx+1)),2*(nelx+1)*(nely+1),2*(nelx+1));
        
        Kt = [Pr*K_cnd*Pr' Pr_pr*zeros(size(Kce_pr)); zeros(size(Kce_pr))'*Pr_pr' zeros(size(Kee_pr))];
        
        % Truss Displacement Recovery
%         U_c = (K_cnd^-1)*F_cnd;
%         R_g = -(Kee^-1)*Kce';
%         U_e = R_g*U_c;
%         U_t = [U_c; U_e];
        
        F_cnd_s = sparse([1 2 11 12 21 22]',ones(1,6)',F_cnd,2*(nelx+1),1);
        
        Kt = (Kt+Kt')/2;
        Ft = [Pr*F_cnd_s; zeros(size(Fe_pr))];
        U = zeros(2*(nely+1)*(nelx+1)+length(Fe_pr),1);
        F = U; %
        F = F + Ft;
        alldofs = [1:length(F)];
        freedofs = setdiff(alldofs,fixeddofs); %For Design Zone
  % end
   % New Condensed Stiffness Matrix (For new Area)
    outit   = outit+1;
    outeriter = outeriter+1;
    if ft == 1
    xPhys(:)=xval(end-5);
    elseif ft == 2
    xPhys(:) = (H*xval(1:end-5))./Hs;
    end
    
    %% FE-ANALYSIS
    sK = reshape(KE(:)*(Emin+xPhys(:)'.^penal*(E0-Emin)),64*nelx*nely,1);
    K0 = sparse(iK,jK,sK,size(Kt,1),size(Kt,2)); 
    K = K0 + Kt; % Question: For this step, we need to keep Kt the same size as K0, 
                 % so I have replaced the Kt slaves with zeros. Will this
                 % cancel the speed-up gain we get from condensation.               
    K = (K+K')/2;
    U(freedofs) = K(freedofs,freedofs)\F(freedofs); % We get a near singular matrix because of the zeros in the matrix
%      cndK=log(condest(K(freedofs,freedofs)))/log(10);
    %% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
    %Compliance of DZ only
    dv = ones(nely,nelx);
    %Total Compliance (DZ + Truss)
    c = abs(F'*U);% + compliance_truss); %Recovery of Truss Compliance
    ce=reshape(sum((U(edofMat)*KE).*U(edofMat),2),nely,nelx);
    dc = -penal*(E0-Emin)*xPhys.^(penal-1).*ce;
    
    %% FILTERING/MODIFICATION OF SENSITIVITIES
    if ft == 1
        dc(:) = H*(x(:).*dc(:))./Hs./max(1e-3,x(:));
    elseif ft == 2
        dc(:) = H*(dc(:)./Hs);
        dv(:) = H*(dv(:)./Hs);
    end
    %% Gather info for MMA
    f0val = c; %Full Truss + DZ Compliance

    % Constraint Function Definition
    vol_frac = (mean(xPhys(:))-volfrac)/volfrac; %Volume Fraction
   % mass_constraints = (xval(end-4)*L*nelx*rho +  xval(end-3)*nelx*L*rho +... 
   % xval(end-2)*nely*L*rho + xval(end-1)*nely*L*rho + xval(end)*nely*L*rho + sum(rho*xval(1:nelx*nely))- tarmass)/tarmass; %Check Mass Constraint
    
    mass_constraints_A1 = (xval(end-4)*L*nelx*rho - tarmass_A1)/tarmass_A1;
    mass_constraints_A2 = (xval(end-4)*L*nelx*rho - tarmass_A2)/tarmass_A2;
    mass_constraints_A3 = (xval(end-4)*L*nelx*rho - tarmass_A3)/tarmass_A3;
    mass_constraints_A4 = (xval(end-4)*L*nelx*rho - tarmass_A4)/tarmass_A4;
    mass_constraints_A5 = (xval(end-4)*L*nelx*rho - tarmass_A5)/tarmass_A5;


    mass =xval(end-4)*L*nelx*rho +  xval(end-3)*nelx*L*rho +... 
    xval(end-2)*nely*L*rho + xval(end-1)*nely*L*rho + xval(end)*nely*L*rho + sum(rho*xval(1:nelx*nely));
    fval = [vol_frac; mass_constraints_A1; mass_constraints_A2; mass_constraints_A3; mass_constraints_A4; mass_constraints_A5];
    
    % Construction of Sensitivity of Area of the Truss for Optimization (Truss)
    [dK_cnd_dA, Kee_a, B_a] = area_cond_by_set(nelx,nely,E0,L,Kce,Kee);
    dK_cnd_vec_dA = dK_cnd_dA(:);
    
    dK_cnd_dA = sparse(coupling_dofs_sp1',coupling_dofs_sp2',dK_cnd_vec_dA,2*(nelx+1),2*(nelx+1));
    K_area = [Pr*dK_cnd_dA*Pr' Pr_pr*zeros(size(Kce_pr)); zeros(size(Kce_pr))'*Pr_pr' zeros(size(Kee_pr))];
    K_area = (K_area+K_area')/2;
    
    [Kee_a1, B_a1] = area_cond_by_set_dA1(nelx,nely,E0,L,Kce,Kee);
    [Kee_a2, B_a2] = area_cond_by_set_dA2(nelx,nely,E0,L,Kce,Kee);
    [Kee_a3, B_a3] = area_cond_by_set_dA3(nelx,nely,E0,L,Kce,Kee);
    [Kee_a4, B_a4] = area_cond_by_set_dA4(nelx,nely,E0,L,Kce,Kee);
    [Kee_a5, B_a5] = area_cond_by_set_dA5(nelx,nely,E0,L,Kce,Kee);
    [dtruss_dA1, dtruss_dA2, dtruss_dA3, dtruss_dA4, dtruss_dA5] = truss_compliance_sensitivity_by_area_by_set(Kee_a1,Kee_a2,Kee_a3,Kee_a4,Kee_a5, Fe, B_a1,B_a2,B_a3,B_a4,B_a5, u_c);    
    
    dc_dA1 = U'*K_area*U + dtruss_dA1;
    dc_dA2 = U'*K_area*U + dtruss_dA2;
    dc_dA3 = U'*K_area*U + dtruss_dA3;
    dc_dA4 = U'*K_area*U + dtruss_dA4;
    dc_dA5 = U'*K_area*U + dtruss_dA5;

    
    %Construction of Sensitivity of Objective Functions
    df0dc_dx = dc(:);
    df0dx = [df0dc_dx; dc_dA1; dc_dA2; dc_dA3; dc_dA4; dc_dA5]; 
    
    dfdx_dv = dv(:)'/length(dv(:))/volfrac;%;-1*dfandisp(:)'/0.38
    
    dfdx_A1 = L*nelx*rho;
    dfdx_A2 = L*nelx*rho;
    dfdx_A3 = nely*L*rho;
    dfdx_A4 = nely*L*rho;
    dfdx_A5 = nely*L*rho;
  
    dfdx = [dfdx_dv 0 0 0 0 0; zeros(1,length(dfdx_dv)) dfdx_A1 0 0 0 0; zeros(1,length(dfdx_dv)) 0 dfdx_A2 0 0 0;...
        zeros(1,length(dfdx_dv)) 0 0 dfdx_A3 0 0;...
        zeros(1,length(dfdx_dv)) 0 0 0 dfdx_A4 0; zeros(1,length(dfdx_dv)) 0 0 0 0 dfdx_A5];
  
    innerit=0;
    outvector1 = [outeriter innerit f0val fval'];
    outvector2 = xval;
    %% MMA code optimization
    [xmma,ymma,zmma,lam,xsi,eta,mu,zet,S,low,upp] = ...
        mmasub(m,n,outeriter,xval,xmin,xmax,xold1,xold2, ...
            f0val,df0dx,fval,dfdx,low,upp,a0,a,C,d);
    xold2 = xold1;
    xold1 = xval;
    xval  = xmma;
    change=norm(xval-xold1);
    
    %% PRINT RESULTS
    fprintf(' It.:%5i Obj.:%11.4f Vol.:%7.3f kktnorm.:%7.3f Area.:%7.3f\n Mass: %7.3f\n',outit,c, ...
        mean(xPhys(:)),kktnorm, xval(end-4), mass);
    figure(2)
    hold on
    plot(outit,c,'bo','MarkerFaceColor','b')
    plot(outit,mean(xPhys(:))*100,'ro','MarkerFaceColor','r')    
    %plot(outit,cndK,'go','MarkerFaceColor','g')
    title(['Convergence V = ',num2str(mean(xPhys(:))*100),',compliance =',num2str(c),', iter = ', num2str(outit)])
    grid on
    legend('compliance','V %')
    xlabel('iter')
    
     %% PLOT Areas
    figure(3)
    hold on
    plot(outit,xval(end),'ro','MarkerFaceColor','r')
   % plot(outit,xval(end-3),'go','MarkerFaceColor','g')
    plot(outit,xval(end-1),'yo','MarkerFaceColor','y')
    %plot(outit,xval(end-1),'co','MarkerFaceColor','c')
    %plot(outit, xval(end),'mo','MarkerFaceColor','m')
    title(['iter = ', num2str(outit), ' A1 =' , num2str(xval(end-1)),' A3 =' , num2str(xval(end))])
    grid on
    legend('A1','A2')
    xlabel('iter')
    
    %% %% The residual vector of the KKT conditions is calculated:
    [residu,kktnorm,residumax] = ...
        kktcheck(m,n,xval,ymma,zmma,lam,xsi,eta,mu,zet,S, ...
        xmin,xmax,df0dx,fval,dfdx,a0,a,C,d);
    outvector1 = [outeriter innerit f0val fval(:)'];
    outvector2 = xval;
    
    %% PLOT DENSITIES
    h = figure(1); 
    set(h,'Color',[1 1 1]);
    [Yy,Xx]=find(nodenrs);
    Yy=nely+1-Yy;
    Xx=Xx-1;
    hold on; patchplot2 = patch('Vertices',[Xx,Yy],'Faces',edofMat(:,[2,4,6,8])/2,'FaceVertexCData',(1-xval(1:end-5))*[1 1 1],'FaceColor','flat','EdgeColor','none'); axis equal; axis off; drawnow;hold off
    title(['Design zone at iteration ',num2str(outit)])
    colormap(gray);
end