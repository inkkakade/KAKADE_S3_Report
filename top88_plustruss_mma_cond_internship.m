function x=top88_plustruss_mma_cond_internship(nelx,nely,L,volfrac,tarmass,penal,rmin,ft)
filename = 'Topology_Optimization.gif1';
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
B11 = [-4  3 -2  9;  3 -4 -9  4; -2 -9 -4 -3;  9  6 -3 -4];
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
A_max = 3;
xPhys = x;
loop = 0;
change = 1;
m = 2;
n = length(xPhys(:))+1;
epsimin = 0.0000001;
eeen    = ones(n-1,1);
eeem    = ones(m,1);
zeron   = zeros(n,1);
zerom   = zeros(m,1);
xval    = [xPhys(:); 2];
xold1   = xval;
xold2   = xval;
xmin    = zeron;
xmin(end) = 0.3;
xmax    = cat(1,eeen,A_max);
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
   A = xval(end);
   I = 10*(sqrt(A)^4)/12;
   %Condensation of Structure
        [Kcc,Kce,Kee,Fe,Fc,K_cnd,F_cnd,compliance_truss, u_c,B] = truss_stiffness_condensation(nelx,nely,E0,A,I,L);
   
     %% Generate the Projection matrix, load vector and stiffness matrix ready for assembly
        coupling_nodes = [nely+1:nely+1:(nely+1)*(nelx+1)];
        coupling_dofs = [coupling_nodes*2-1;coupling_nodes*2];
        coupling_dofs = coupling_dofs(:);
        Pr = sparse(coupling_dofs,1:2*(nelx+1),ones(1,2*(nelx+1)),2*(nelx+1)*(nely+1),2*(nelx+1));
        Kt = [Pr*K_cnd*Pr' Pr*zeros(size(Kce)); zeros(size(Kce))'*Pr' zeros(size(Kee))];
        
        % Truss Displacement Recovery
%         U_c = (K_cnd^-1)*F_cnd;
%         R_g = -(Kee^-1)*Kce';
%         U_e = R_g*U_c;
%         U_t = [U_c; U_e];
        
        Kt = (Kt+Kt')/2;
        Ft = [Pr*F_cnd; zeros(size(Fe))];
        U = zeros(2*(nely+1)*(nelx+1)+length(Fe),1);
        F = U; %
        F = F + Ft;
        alldofs = [1:length(F)];
        freedofs = setdiff(alldofs,fixeddofs); %For Design Zone
%    end
   % New Condensed Stiffness Matrix (For new Area)
    outit   = outit+1;
    outeriter = outeriter+1;
    if ft == 1
    xPhys(:)=xval(end-1);
    elseif ft == 2
    xPhys(:) = (H*xval(1:end-1))./Hs;
    end
    
    %% FE-ANALYSIS
    sK = reshape(KE(:)*(Emin+xPhys(:)'.^penal*(E0-Emin)),64*nelx*nely,1);
    K0 = sparse(iK,jK,sK,size(Kt,1),size(Kt,2)); 
    K = K0 + Kt; % Question: For this step, we need to keep Kt the same size as K0, 
                 % so I have replaced the Kt slaves with zeros.              
    K = (K+K')/2;
    U(freedofs) = K(freedofs,freedofs)\F(freedofs); % We get a near singular matrix because of the zeros in the matrix
%      cndK=log(condest(K(freedofs,freedofs)))/log(10);
    %% OBJECTIVE FUNCTION AND SENSITIVITY ANALYSIS
    %Compliance of DZ only
    dv = ones(nely,nelx);
    %Total Compliance (DZ + Truss)
    DZ_c = F'*U;
    c = DZ_c;% + compliance_truss; 
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
    f0val = c;   
    % Constraint Function Definition
    vol_frac = (mean(xPhys(:))-volfrac)/volfrac; %Volume Fraction
    mass_constraints = (xval(end)*L*rho + mean(xPhys(:))*rho - tarmass)/tarmass; %Check Mass Constraint
    mass = xval(end)*L*rho + mean(xPhys(:))*rho;
    mass_dz = mean(xPhys(:))*rho/mass;
    mass_t =  xval(end)*L*rho/mass; 
    fval = [vol_frac; mass_constraints];
    
    % Construction of Sensitivity of Area of the Truss for Optimization (Truss)
    [dK_cnd_dA, Kee_a, B_a] = area_cond(nelx,nely,E0,L,Kce,Kee,B);
    K_area = [Pr*dK_cnd_dA*Pr' Pr*zeros(size(Kce)); zeros(size(Kce))'*Pr' zeros(size(Kee))];
    K_area = (K_area+K_area')/2;
    %[dtruss_compliance_sensitivity_dA] = truss_compliance_sensitivity(Kee_a, Fe, B_a, u_c);
    dc_dA = (U'*K_area*U)/10e4;% + dtruss_compliance_sensitivity_dA;
    
    %Construction of Sensitivity of Objective Functions
    df0dc_dx = dc(:);
    df0dx = [df0dc_dx; dc_dA]; 
    
    dfdx_dv = dv(:)'/length(dv(:))/volfrac;%;-1*dfandisp(:)'/0.38
    dfdx_A = L*rho;
    dfdx = [dfdx_dv 0; zeros(1,length(dfdx_dv)) dfdx_A];
  
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
        mean(xPhys(:)),kktnorm, A, mass);
    fprintf(' It.:%5i DZ C.:%11.4f T C.:%7.3f Mass Pylon.:%7.3f Mass Engine.:%7.3f\n',outit,DZ_c, ...
       compliance_truss,mass_dz, mass_t);
    figure(2)
    hold on
    plot(outit,c,'bo','MarkerFaceColor','b')
    plot(outit,mean(xPhys(:))*100,'ro','MarkerFaceColor','r')    
    %plot(outit,cndK,'go','MarkerFaceColor','g')
    title(['Convergence V = ',num2str(mean(xPhys(:))*100),',compliance =',num2str(c),', iter = ', num2str(outit)])
    grid on
    legend('compliance','V %')
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
    hold on; patchplot2 = patch('Vertices',[Xx,Yy],'Faces',edofMat(:,[2,4,6,8])/2,'FaceVertexCData',(1-xval(1:end-1))*[1 1 1],'FaceColor','flat','EdgeColor','none'); axis equal; axis off; drawnow;hold off
    title(['Design zone at iteration ',num2str(outit)])
    colormap(gray);
    drawnow
    % Capture the plot as an image 
      frame = getframe(h); 
      im = frame2im(frame);     
      [imind,cm] = rgb2ind(im,256); 
    % Write to the GIF File 
      if outit == 1 
          imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
      else 
          imwrite(imind,cm,filename,'gif','WriteMode','append'); 
      end 
end