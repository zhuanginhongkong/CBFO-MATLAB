classdef CBFO % A CLASS OF COLORED BODY-FITTED OPTIMIZATION (CBFO)
    properties
        Param
        Input
        Output
    end
    methods
        function obj = CBFO(nelx,nely,volfrac)
            obj.Param.penal          =   3 ;
            obj.Param.rmin           =   2.5 ;
            obj.Param.maxedge        =  10 ;
            obj.Param.minedge        =   2.5 ;
            obj.Param.d1             =   0.6 ;
            obj.Param.d2             =   1.2 ;
            obj.Param.E0             = 0.1 ;
            obj.Param.Emin           = 1e-4 ;
            obj.Param.nu             = 0.3 ;
            obj.Input.nelx           = nelx ;
            obj.Input.nely           = nely ;
            obj.Input.volfrac        = volfrac ;
            obj.Input.BoundaryPoints = [ 0, 0; nelx, nely] ;
            obj = obj.Optimize ;
        end
        function obj = Optimize(obj)
            %% READ INPUTS AND PARAMETERS
            nelx     = obj.Input.nelx ;
            nely     = obj.Input.nely ;
            BDY      = obj.Input.BoundaryPoints ;
            volfrac  = obj.Input.volfrac;
            maxedge  = obj.Param.maxedge ;
            minedge  = obj.Param.minedge ;
            d1       = obj.Param.d1 ;
            d2       = obj.Param.d2 ;
            penal    = obj.Param.penal ;
            rmin     = obj.Param.rmin ;
            E0       = obj.Param.E0 ;
            Emin     = obj.Param.Emin ;
            nu       = obj.Param.nu ;
            %% Function Citation
            OcUpdate = @ obj.OptimalityCriteriaUpdate ; 
            GenMesh = @ obj.GenerateMesh ;              
            Filter = @ obj.FilterIndex ;                
            ConPt = @ obj.ContourPoints ;               
            FEA  = @ obj.FiniteElementAnalysis ;        
            %% INITIAL REMESHING
            [xn,yn] = meshgrid(0: nelx, 0: nely) ;
            [p, t, ~, ~, Ve, pmid,tp,tnp,Passive] = GenMesh( 1, xn, yn, BDY) ;
            xphy = repmat( volfrac./0.935, length(t), 1) ; figure(1) ;
            xphy(1:length(tp))=0;
            clf ; colormap summer ;
            patch('Faces', tnp,'Vertices', p, 'FaceVertexCData',...
                xphy(length(tp)+1:end),'FaceColor','flat') ; axis off equal tight
            saveas(gcf,'./initial.png') ;
            %% PREPARE FILTER
            [H, Hs] = Filter( t, pmid, rmin) ;
            %% INITIALIZE ITERATION
            loop = 0; change = 1;
            %% START ITERATION
            while change > 0.11
                loop = loop + 1 ;
                %% FEA AND SENSITIVITY ANALYSIS
                [ ~, dCe, J] = FEA( t, p, BDY, xphy, E0, Emin, nu,...
                    penal) ;
                dCe(:) = H*(xphy(:).* dCe(:))./Hs./max(1e-3, xphy(:)) ;
                vol = sum(Ve.* xphy)/( 0.935* nelx * nely);
                Obj(loop) = J;   volt(loop) = vol;
                %% OPTIMALITY CRITERIA UPDATE OF DESIGN VARIABLES
                [xphy , change] = OcUpdate( xphy, dCe, Ve, volfrac,...
                    tp,tnp, p, loop, nelx, nely,Passive);
                fprintf(' It.:%5i Obj.:%8.4f Vol.:%7.3f ch.:%7.3f\n',...
                    loop,J,vol,change);
            end
            %% SECOND REMESHING
            xBF = imgaussfilt(griddata(pmid(:,1), pmid(:,2), xphy , xn,...
                yn,'nearest'),2) ;
            clf; contour(xn,yn,xBF,[volfrac volfrac],'linewidth',2,'color',[62 43 109]/255);  hold on
            [c] = ConPt(contour( xn, yn, xBF, [volfrac volfrac]), d1, d2) ;
            [p,t,~,~,Ve,pmid,tp,tnp,Passive,tv] = GenMesh(2,xn,yn,BDY,c,xBF, maxedge,...
                minedge, 80) ;
            xphy = interp2(xn,yn,xBF,pmid(:,1),pmid(:,2),'cubic') ;
            colormap summer;
            patch('Faces',tnp,'Vertices',p,'FaceVertexCData',xphy(length(tp)+1:end),'FaceColor','flat');
            colorbar; hold on; axis off equal tight ;
            %% PREPARE FILTER
            [H,Hs] = Filter(t, pmid, rmin);
            %% START ITERATION
            while change > 0.02
                loop = loop + 1;
                %% FEA AND SENSITIVITY ANALYSIS
                [~,dCe,J] = FEA( t, p, BDY, xphy, E0, Emin, nu,...
                    penal);
                dCe(:) = H*(xphy(:).*dCe(:))./Hs./max(1e-3, xphy(:));
                vol = sum( Ve.* xphy)/( 0.935* nelx * nely);
                Obj(loop) = J;   volt(loop) = vol;
                %% OPTIMALITY CRITERIA UPDATE OF DESIGN VARIABLES AND PHYSICAL DENSITIES
                [ xphy, change] = OcUpdate( xphy, dCe, Ve, volfrac,...
                    tp,tnp, p, loop, nelx, nely,Passive);
                fprintf(' It.:%5i Obj.:%8.4f Vol.:%7.3f ch.:%7.3f\n',...
                    loop, J, vol, change);
            end
            save('History.mat','Obj','volt');          % SAVE OPTIMIZATION DATA
            save('Structure.mat','xn','yn','c','xphy','maxedge','minedge','BDY','p','t','xBF','pmid','volfrac','nelx','nely','E0','Emin','nu','penal'); 
            %% FINAL RESULT PLOTTTING
            xBF = imgaussfilt(griddata( pmid(:,1), pmid(:,2), xphy, xn,...
                yn,'nearest'), 2);
            clf; contour(xn,yn,xBF,[0.5 0.5],'linewidth',2,'color',...
                [62 43 109]/255);  hold on
            [c] = ConPt(contour(xn,yn,xBF,[0.5 0.5]), d1, d2);
            [p,t,t1,t2,Ve,pmid,tp,tnp,Passive,tv] = GenMesh(2,xn,yn,BDY,c,xBF, maxedge,...
                minedge, 600) ;
            clf; patch('Faces', tv, 'Vertices', p, 'EdgeColor', 'k',...
                'FaceColor',[0 127 102]/255); hold on; axis off equal tight
            patch( 'Faces', t2, 'Vertices', p, 'EdgeColor', 'k', ...
                'FaceColor',[255 255 102]/255) ; 
            saveas(gcf,'./final.png');
            t = [t1; t2];
            xphy = ones( length(t), 1); xphy(1:length(t1)) = 0 ;
            [~,~,J] = FEA( t, p, BDY, xphy, E0, Emin, nu, penal);
            vol = sum( Ve.* xphy)/( nelx* nely) ;
            fprintf(' Optimization %5s Obj.:%8.4f Vol.:%7.3f,',...
                'Result', J, vol); disp('Script Run Successfuly!!!!!!!!!!')
            p = p'; t = t';
            q = mean(pdetriq(p,t));
            obj.Output = xphy;
        end
    end
    %% SUBFUNCTIONS:
    methods (Static)
        %% FIND POINTS ON THE CONTOUR
        function [c] = ContourPoints(c,d1,d2)
            tol = 1e-12; num = 1; col = 1;
            while col < size(c,2)
                idx = col+1:col+c(2,col);
                s(num).x = c(1,idx);     s(num).y = c(2,idx);
                s(num).isopen = abs(diff(c(1,idx([1 end])))) >tol || abs(diff(c(2,idx([1 end])))) >tol;
                num = num+1; col = col+c(2,col)+1;
            end
            c = [];
            for k = 1:num-1
                ct = [s(k).x; s(k).y];
                if length(ct)>4
                    ct1 = ct; ndist = sqrt(sum(diff(ct,1,2).^2,1));
                    for i = 1:size(ndist,2)
                        if  ndist(i) < d1
                            ct1(:,i) = [0;0];  ct1(:,(i+1)) = 0.5*(ct(:,i)+ct(:,(i+1)));
                        end
                    end
                    ct1(:,all(ct1==0,1)) = [];
                    if  sqrt(sum((ct1(:,1)-ct1(:,end)).^2,1)) < d1
                        ct1(:,end) = [];
                    end
                    if size(ct1,2)>2
                        if s(k).isopen == 0
                            ct1 = [ct1 ct1(:,1)];
                        end
                        ndist = sqrt(sum(diff(ct1,1,2).^2,1)); ct=ct1(:,1);
                        for i = 1:size(ndist,2)
                            if  ndist(i) > d2
                                ct = [ct 0.5*(ct1(:,i)+ct1(:,(i+1))) ct1(:,i+1)];
                            else
                                ct = [ct ct1(:,i+1)];
                            end
                        end
                        if s(k).isopen == 1
                            c = [c; ct'];
                        else
                            c = [c; ct(:,1:end-1)'];
                        end
                    end
                end
            end
        end
        %% BODY FITTED MESH GENERATOR
        function [p,t,t1,t2,Ve,pmid,tp,tnp,Passive,tv] = GenerateMesh( meshtime, xn, yn, BDY, c, dN,maxedge,minedge,maxiter)
            Uniqnode = @ CBFO.UniqueNode ;
            x = xn; x(2:2:end,:) = x(2:2:end,:)+0.5; pi = [x(:),yn(:)];
            px = 0.75 * BDY(2,1); py = 0.74 * BDY(2,2);
            [xp,yp] = meshgrid(px,py:BDY(2,2)); [xp2,yp2] = meshgrid(px:BDY(2,1),py);
            P = [xp(:),yp(:); xp2(:),yp2(:)];
            Forcepts = [BDY(2,1) py];
            if meshtime ==1
                pi(:,1) = min(BDY(2,1),max(BDY(1,1),pi(:,1))); pi(:,2) = min(BDY(2,2),max(BDY(1,2),pi(:,2)));
                p=[P; Forcepts;  BDY(2,1) BDY(1,2); BDY(1,1) BDY(2,2); BDY(1,1) BDY(1,2); BDY(2,1) BDY(2,2); BDY(2,1),(BDY(1,2)+BDY(2,2))/2; pi];
                p = unique(p,'rows','stable');
                t = delaunayn(p);
                [p, t] = Uniqnode(p,t) ;
                pmid = (p(t(:,1),:)+p(t(:,2),:)+p(t(:,3),:))/3 ;
                t1 = [] ; t2 = [] ;
                tnp = t(~(pmid(:,1)>px & pmid(:,2)>py),:);
                tp = t((pmid(:,1)>px & pmid(:,2)>py),:);
                Passive = length(tp);
                t = [tp ; tnp];
                pmid = (p(t(:,1),:)+p(t(:,2),:)+p(t(:,3),:))/3;
            elseif  meshtime == 2
                d = zeros(size(pi,1),1);
                for i = 1:size(pi,1)
                    d(i) = sqrt(min((pi(i,1)-c(:,1)).^2+(pi(i,2)-c(:,2)).^2));
                end
                r0 = 1./min(max(minedge,d),maxedge).^2;
                pfix=[ c; P; Forcepts; BDY(2,1) BDY(1,2); BDY(1,1) BDY(2,2); BDY(1,1) BDY(1,2); BDY(2,1) BDY(2,2);];
                p = [pfix; pi(r0./max(r0)>0.5,:)];
                p = unique(p,'rows','stable');
                p1 = 0; warning off; t = delaunayn(p);
                for i = 1:maxiter
                    if max(sum((p-p1).^2,2))>1e-6
                        t = delaunayn(p);
                        edges = unique(sort([t(:,[1,2]);t(:,[1,3]);t(:,[2,3])],2),'rows') ;
                        p1 = p;
                    end
                    midpoint = (p(edges(:,1),:)+p(edges(:,2),:))/2;
                    d = zeros(size(midpoint,1),1);
                    for j = 1:size(midpoint,1)
                        d(j) = sqrt(min((midpoint(j,1)-c(:,1)).^2+(midpoint(j,2)-c(:,2)).^2));
                    end
                    L = sqrt(sum((p(edges(:,1),:)-p(edges(:,2),:)).^2,2));
                    L = min(max(0.1,L),8);
                    L1 = min(max(minedge,d),maxedge);
                    L0 = 1.2*L1*sqrt(sum(L.^2)/sum(L1.^2));
                    Fb = max(L0-L,0)./L *[1,1].*(p(edges(:,1),:)-p(edges(:,2),:));
                    Fp = full(sparse(edges(:,[1,1,2,2]),ones(size(d))*[1,2,1,2],[Fb,-Fb],size(p,1),2));
                    Fp(1:size(pfix,1),:) = 0 ;
                    p = p + 0.2*Fp ;
                    p(:,1) = min(BDY(2,1),max(BDY(1,1),p(:,1))) ;
                    p(:,2) = min(BDY(2,2),max(BDY(1,2),p(:,2))) ;
                end
                [p, t] = Uniqnode( p, t);
                pmid = (p(t(:,1),:)+p(t(:,2),:)+p(t(:,3),:))/3;
                dE = interp2(xn,yn,dN,pmid(:,1),pmid(:,2),'cubic');
                tnp = t(~(pmid(:,1)>px & pmid(:,2)>py),:);
                tp = t((pmid(:,1)>px & pmid(:,2)>py),:);
                Passive = length(tp);
                dEnp = dE(~(pmid(:,1)>px & pmid(:,2)>py),:);
                tv=tnp(dEnp<0.5,:);
                t1=[tp; tnp(dEnp<0.5,:)]; t2=tnp(dEnp>=0.5,:); t=[tp;tnp];
                pmid = (p(t(:,1),:)+p(t(:,2),:)+p(t(:,3),:))/3;
            end
            Ve = zeros(length(t),1) ;
            for kk = 1:length(t)
                Ve(kk) = 0.5.*det([ones(3,1) p(t(kk,:),:)]);
            end
        end
        %% REMOVE REPEATED NODES
        function [P,T] = UniqueNode(p,t)
            Ve=zeros(length(t),1) ;
            for kk = 1:length(t)
                Ve(kk) = 0.5.*det([ones(3,1) p(t(kk,:),:)]);
            end
            t((Ve==0),:) = [];
            P = unique(p(unique(sort(t(:))),:),'rows','stable');
            for i = 1:length(t)
                for j = 1:3  T(i,j) = find(P(:,1)==p(t(i,j),1)& P(:,2)==p(t(i,j),2)); end
            end
        end
        %% PREPARE ADAPTIVE FILTER
        function [H,Hs] = FilterIndex(t,pmid,rmin)
            iH = ones(length(t)^2,1) ;
            jH = ones(size(iH)) ;
            sH = zeros(size(iH)) ;
            k = 0 ;
            % Fullscale mesh filter
                for i1 = 1:length(t)
                    for i2 = 1:length(t)
                        k = k+1;
                        iH(k) = i1;
                        jH(k) = i2;
                        sH(k) = max(0,rmin-sqrt((pmid(i1,1)-pmid(i2,1))^2+(pmid(i1,2)-pmid(i2,2))^2));
                    end
                end
            H = sparse(iH,jH,sH);
            Hs = sum(H,2);
        end
        %% ELEMENT STIFFNESS MATRIX
        function [Ke] = ElementMatrixKe(X,Y,E0,nu)
            D = E0/(1-nu^2)*[1  nu  0;  nu  1  0;  0  0  (1-nu)/2];
            J = [X(1)-X(3) Y(1)-Y(3);X(2)-X(3) Y(2)-Y(3)];
            Be = 1/det(J)*[J(2,2) 0 -J(1,2) 0 -J(2,2)+J(1,2) 0; 0 -J(2,1) 0 J(1,1) 0 J(2,1)-J(1,1); -J(2,1) J(2,2) J(1,1) -J(1,2) J(2,1)-J(1,1) -J(2,2)+J(1,2)];
            Ke = 1/2*det(J)*Be'*D*Be;
        end
        %% FINITE ELEMENT ANALYSIS
        function [Ce,dCe,J] = FiniteElementAnalysis(t,p,BDY,x,E0,Emin,nu,penal)
            EM = @CBFO.ElementMatrixKe ;
            NT = length(t); KK = zeros(6,6*NT);
            for i = 1:NT
                KK(:,6*i-5:6*i) = (Emin+x(i).^penal*(E0-Emin)) * EM(p(t(i,:),1),p(t(i,:),2),1,nu);
            end
            elemDof = zeros(NT,6); elemDof(:,[1 3 5]) = 2*t-1; elemDof(:,[2 4 6]) = 2*t;
            iK = reshape(kron(elemDof,ones(6,1))',36*NT,1);
            jK = reshape(kron(elemDof,ones(1,6))',36*NT,1);
            sK = reshape(KK,36*NT,1);
            NK = sparse(iK,jK,sK,2*length(p),2*length(p));
            kin = 0.1;            % Definition of the stiffness of spring at input port
            kout = 0.1;           % Definition of the stiffness of spring at output port
            %% LOAD CASES
            fixedNodes1=find(p(:,1)==BDY(1,1) & p(:,2)==BDY(1,2)); % The node on the bottom left corner is fixed
            fixedDof1=[2*fixedNodes1-1; 2*fixedNodes1]; % Definiton of fixed DOF
            fixedNodes2=find(p(:,2)==BDY(2,2)); % Symmetric constraint
            fixedDof2= 2*fixedNodes2; % Definiton of fixed DOF
            fixedDof=[fixedDof1;fixedDof2];
            forceNodes=find(p(:,1)==BDY(1,1) & p(:,2)==BDY(2,2)); % The node on the middle left
            OutputNodes=find(p(:,1)==BDY(2,1) & p(:,2)==0.74*BDY(2,2)); 
            din = 2*forceNodes-1;
            dout = 2*OutputNodes;
            AllDof = 1:2*length(p);
            freeDofs = setdiff(AllDof,fixedDof);
            U1 = zeros(2*length(p),1);
            U2 = zeros(2*length(p),1);
            F = zeros(2*length(p),2);
            F(din,1) = 100;
            F(dout,2) = 1;
            NK(din,din) = NK(din,din) + kin;
            NK(dout,dout) = NK(dout,dout) + kout;
            NK = (NK+NK')/2;
            U1(freeDofs,:) = NK(freeDofs,freeDofs) \ F(freeDofs,1);   %Calculate displacement matrix
            U2(freeDofs,:) = NK(freeDofs,freeDofs) \ F(freeDofs,2);   %Calculate displacement matrix
            J =U1(dout,1);
            for i = 1:NT
                Ce(i) = sum((U2(elemDof(i,:))'*KK(:,6*i-5:6*i)).*U1(elemDof(i,:))',2);
            end
            dCe= Ce';
            dCe(dCe>0)=-dCe(dCe>0)*0.01;
        end
        %% OPTIMALITY CRITERIA UPDATE OF DESIGN VARIABLES AND PHYSICAL DENSITIES
        function [xnew,change] = OptimalityCriteriaUpdate(x,dCe,Ve,volfrac,tp,tnp,p,loop,nelx,nely,Passive)
            l1 = 0; l2 = 1e9;  move = 0.5;
            while (l2-l1)/(l1+l2) > 1e-3
                lmid = 0.5*(l2+l1);
                xnew = max(0,max(x-move,min(1,min(x+move,x.*sqrt(-dCe/lmid)))));
                xnew(1:Passive) = 0;
                xPhys = xnew;
                if sum(xPhys(:).*Ve(:)) > 0.935* volfrac*nelx*nely, l1 = lmid; else l2 = lmid; end
            end
            change = max(abs(xnew(:)-x(:)));
            colormap summer;
            clf; patch('Faces',tnp,'Vertices',p,'FaceVertexCData',xnew(length(tp)+1:end),'FaceColor','flat');
            colorbar; hold on; axis off equal tight ;
            saveas(gcf,['./fig', int2str(loop) '.png']);
        end
    end
end
