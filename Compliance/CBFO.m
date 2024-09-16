classdef CBFO_OC  % A CLASS OF COLORED BODY-FITTED OPTIMIZATION (CBFO)
    properties
        Param
        Input
        Output
    end
    methods
        function obj = CBFO_OC( nelx, nely, volfrac, loadcase)
            obj.Param.penal          =   3 ;
            obj.Param.rmin           =   6 ;
            obj.Param.maxedge        =  10 ;
            obj.Param.minedge        =   2 ;
            obj.Param.d1             = 0.4 ;
            obj.Param.d2             =   1 ;
            obj.Param.E0             = 1e5 ;
            obj.Param.Emin           = 1e-4 ;
            obj.Param.nu             = 0.3  ;
            obj.Input.nelx           = nelx ;
            obj.Input.nely           = nely ;
            obj.Input.volfrac        = volfrac ;
            obj.Input.loadcase       = loadcase ;
            obj.Input.Boundary = [ 0, 0; nelx, nely] ;
            obj = obj.Optimize;
            %         obj.GUIDemo(obj) ; % GUI DISPLAY
        end
        function obj = Optimize(obj)
            %% READ INPUTS AND PARAMETERS
            nelx     = obj.Input.nelx ;
            nely     = obj.Input.nely ;
            BDY      = obj.Input.Boundary ;
            volfrac  = obj.Input.volfrac;
            loadcase = obj.Input.loadcase ;
            maxedge  = obj.Param.maxedge ;
            minedge  = obj.Param.minedge ;
            d1       = obj.Param.d1 ;
            d2       = obj.Param.d2 ;
            penal    = obj.Param.penal ;
            rmin     = obj.Param.rmin ;
            E0       = obj.Param.E0 ;
            Emin     = obj.Param.Emin ;
            nu       = obj.Param.nu ;

            %% FUNCTION CITATION
            OcUpdate = @ obj.OptimalityCriteriaUpdate ;
            GenMesh = @ obj.GenerateMesh ;
            Filter = @ obj.FilterIndex ;
            ConPt = @ obj.ContourPoints ;
            FEA  = @ obj.FiniteElementAnalysis ;

            %% INITIAL MESHING
            [ xn, yn] = meshgrid(0: nelx, 0: nely) ;
            [ p, t, ~, ~, Ve, pmid] = GenMesh( 1, xn, yn, BDY) ; %p 节点坐标 dim*NodeNumber t 单元连接性 (nodeNumber4
            xphy = repmat( volfrac, length(t), 1) ; figure(1) ;
            clf ; colormap summer ;
            patch( 'Faces', t,'Vertices', p, 'FaceVertexCData', xphy, 'FaceColor', 'flat') ; colorbar ;
            %         saveas( gcf, './initial.png') ;

            % PREPARE FILTER
            [ H, Hs] = Filter( pmid, rmin) ;

            %% START ITERATION
            % INITIALIZE ITERATION
            loop = 0 ; change = 1;
            while change > 0.1
                loop = loop + 1 ;

                % FEA AND SENSITIVITY ANALYSIS
                [ dCe, J] = FEA( t, p, BDY, xphy, E0, Emin, nu, penal, loadcase) ;
                dCe(:) = H*( xphy(:).* dCe(:))./Hs./max( 1e-3, xphy(:)) ;
                vol = sum( Ve.* xphy)/( nelx * nely) ;

                % OPTIMALITY CRITERIA UPDATE OF DESIGN VARIABLES
                [ xphy, change] = OcUpdate( xphy, dCe, Ve, volfrac, t, p, loop, nelx, nely) ;
                fprintf( 'It.:%5i Obj.:%8.4f Vol.:%7.3f ch.:%7.3f\n', loop, J, vol, change) ;
                obj.Output.FirstMeshing.Xphy{loop} = xphy ;
            end
            obj.Output.FirstMeshing.p = p ;
            obj.Output.FirstMeshing.t = t ;
            obj.Output.Designcheck = loop;

            %% FIRST REMESHING
            xBF = imgaussfilt( griddata(pmid(:,1), pmid(:,2), xphy , xn, yn, 'nearest'), 2) ;
            [c] = ConPt(contour( xn, yn, xBF, [0.5 0.5]), d1, d2) ;
            [ p, t, ~, ~, Ve, pmid] = GenMesh( 2, xn, yn, BDY, c, xBF, maxedge, minedge, 80) ;
            xphy = interp2( xn, yn, xBF, pmid(:,1), pmid(:,2), 'cubic') ;

            % PREPARE FILTER
            [ H, Hs] = Filter(pmid, rmin) ;

            %% START ITERATION
            while change > 0.01
                loop = loop + 1 ;

                % FEA AND SENSITIVITY ANALYSIS
                [ dCe, J] = FEA( t, p, BDY, xphy, E0, Emin, nu, penal, loadcase) ;
                dCe(:) = H*(xphy(:).*dCe(:))./Hs./max(1e-3, xphy(:));
                vol = sum(Ve.*xphy)/(nelx*nely) ;

                % OPTIMALITY CRITERIA UPDATE OF DESIGN VARIABLES
                [ xphy, change] = OcUpdate( xphy, dCe, Ve, volfrac, t, p, loop, nelx, nely) ;
                fprintf('It.:%5i Obj.:%8.4f Vol.:%7.3f ch.:%7.3f\n', loop, J, vol, change) ;
                obj.Output.SecondMeshing.Xphy{loop- obj.Output.Designcheck } = xphy ;
            end
            obj.Output.SecondMeshing.p = p ;
            obj.Output.SecondMeshing.t = t ;

            %% SECOND REMESHING AND FINAL RESULT PLOTTTING
            xBF = imgaussfilt( griddata(pmid(:,1), pmid(:,2), xphy, xn, yn,'nearest'), 1) ;
            clf; contour(xn, yn, xBF, [0.5 0.5],'linewidth',2,'color', [62 43 109]/255); set(gca,'YDir','reverse') ; hold on
            [c] = ConPt(contour(xn, yn, xBF, [0.5 0.5]), d1, d2) ;
            [p,t,t1,t2,Ve] = GenMesh(2, xn, yn, BDY, c, xBF, maxedge, minedge, 80) ;
            patch('Faces', t1, 'Vertices', p, 'EdgeColor', 'k', 'FaceColor',[0 127 102]/255) ;  axis off equal tight
            patch('Faces', t2, 'Vertices', p, 'EdgeColor', 'k', 'FaceColor',[255 255 102]/255) ; set(gca,'YDir','reverse') ;
            %         saveas(gcf,'./final.png') ;
            xphy = ones(length(t), 1); xphy(1:length(t1)) = 1e-9 ;
            [~,J] = FEA(t, p, BDY, xphy, E0, Emin, nu, penal, loadcase) ;
            vol = sum(Ve.* xphy)/(nelx* nely) ;
            fprintf(' Optimization %5s Obj.:%8.4f Vol.:%7.3f,', 'Result', J, vol) ;
            disp('Script Run Successfuly!!!!!!!!!!');
        end
    end

    % STATIC METHODS
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
        function [ p, t, t1, t2, Ve, pmid] = GenerateMesh( meshtime, xn, yn, BDY, c, dN, maxedge, minedge, maxiter)
            Uniqnode = @ CBFO_OC.UniqueNode ;
            x = xn ; x( 2:2:end, :) = x( 2:2:end, :)+0.5 ; pi = [ x(:), yn(:)] ;
            if meshtime == 1
                pi(:,1) = min(BDY(2,1),max(BDY(1,1),pi(:,1))); pi(:,2) = min(BDY(2,2),max(BDY(1,2),pi(:,2)));
                p = [ BDY(2,1) BDY(1,2); BDY(1,1) BDY(2,2); BDY(1,1) BDY(1,2); BDY(2,1) BDY(2,2); BDY(2,1),(BDY(1,2)+BDY(2,2))/2; pi];
                p = unique( p, 'rows', 'stable');
                t = delaunayn(p);
                [ p, t] = Uniqnode( p, t) ;
                pmid = (p(t( :, 1), :) + p(t( :, 2), :) + p(t( :, 3), :))/3 ;
                t1 = [] ; t2 = [] ;
            elseif  meshtime == 2
                d = zeros( size(pi,1), 1) ;
                for i = 1:size( pi, 1)
                    d(i) = sqrt(min((pi(i,1)-c(:,1)).^2+(pi(i,2)-c(:,2)).^2)) ;
                end
                r0 = 1./min( max(minedge,d), maxedge).^2 ;
                pfix=[c; BDY(2,1) BDY(1,2); BDY(1,1) BDY(2,2); BDY(1,1) BDY(1,2); BDY(2,1) BDY(2,2); BDY(2,1),(BDY(1,2)+BDY(2,2))/2] ;
                p = [pfix; pi(r0./max(r0)>0.5,:)] ;
                p = unique(p,'rows','stable') ;
                p1 = 0 ; warning off ; t = delaunayn(p) ;
                for i = 1 : maxiter
                    if max( sum((p-p1).^2, 2)) > 1e-6
                        t = delaunayn(p) ;
                        edges = unique( sort([t(:,[1,2]);t(:,[1,3]);t(:,[2,3])],2), 'rows') ;
                        p1 = p ;
                    end
                    midpoint = ( p(edges(:,1),:) + p( edges(:,2), :))/2 ;
                    d = zeros( size(midpoint,1), 1) ;
                    for j = 1 : size( midpoint, 1)
                        d(j) = sqrt(min((midpoint(j,1)-c(:,1)).^2+(midpoint(j,2)-c(:,2)).^2)) ;
                    end
                    L = sqrt(sum((p(edges(:,1),:)-p(edges(:,2),:)).^2,2)) ;
                    L1 = min(max(minedge,d),maxedge) ;
                    L0 = 1.2*L1*sqrt(sum(L.^2)/sum(L1.^2)) ;
                    Fb = max(L0-L,0)./L *[1,1].*(p(edges(:,1),:)-p(edges(:,2),:)) ;
                    Fp = full(sparse( edges(:,[1,1,2,2]), ones(size(d))*[1,2,1,2], [Fb,-Fb], size(p,1), 2)) ;
                    Fp( 1 : size(pfix,1), :) = 0 ;
                    p = p + 0.2*Fp ;
                    p(:,1) = min(BDY(2,1), max(BDY(1,1), p(:,1))) ;
                    p(:,2) = min(BDY(2,2), max(BDY(1,2), p(:,2))) ;
                end
                [ p, t] = Uniqnode( p, t);
                pmid = (p(t(:,1),:)+p(t(:,2),:)+p(t(:,3),:))/3 ;
                dE = interp2( xn, yn, dN, pmid(:,1), pmid(:,2), 'cubic') ;
                t1 = t( dE<0.5, :); t2 = t( dE>=0.5, :); t = [ t1; t2]; hold on
                pmid = (p(t(:,1),:)+p(t(:,2),:)+p(t(:,3),:))/3 ;
            end ;   Ve = zeros(length(t),1) ;
            for kk = 1 : length(t)
                Ve(kk) = 0.5.*det([ones(3,1) p(t(kk,:), :)]) ;
            end
        end

        %% REMOVE REPEATED NODES
        function [ P, T] = UniqueNode( p, t)
            Ve = zeros( length(t), 1) ;
            for kk = 1 : length(t)
                Ve(kk) = 0.5.*det([ones(3,1) p(t(kk,:),:)]) ;
            end
            t((Ve==0),:) = [] ;
            P = unique(p(unique(sort(t(:))),:), 'rows', 'stable') ;
            for i = 1 : length(t)
                for j = 1:3  T(i,j) = find(P(:,1)==p(t(i,j),1)& P(:,2)==p(t(i,j),2)); end
            end
        end

        %% PREPARE ADAPTIVE FILTER
        function [ H, Hs] = FilterIndex(pmid, rmin)
            maxNumNeighbor = 200;
            iH = ones(length(pmid)*maxNumNeighbor, 1) ;
            jH = ones(length(pmid)*maxNumNeighbor, 1) ;
            sH = zeros(length(pmid)*maxNumNeighbor, 1) ;
            indexCount = 1;
            for pii = 1 : length(pmid)
                centreEle = pmid(pii,:);
                neighborEleIndex = find(abs(pmid(:,1)-centreEle(1))<=rmin & abs(pmid(:,2)-centreEle(2))<=rmin);
                neighborEleCoord = pmid(neighborEleIndex,:);
                neighborWeight = max((rmin - pdist2(centreEle,neighborEleCoord)),0);
                indexCountNext = indexCount + length(neighborWeight);

                iH(indexCount:indexCountNext-1) = pii;
                jH(indexCount:indexCountNext-1) = neighborEleIndex;
                sH(indexCount:indexCountNext-1) = neighborWeight;
                indexCount = indexCountNext;
            
            end
            H = sparse( iH, jH, sH) ;
            Hs = sum( H, 2) ;
        end

            %% ELEMENT STIFFNESS MATRIX
            function [Ke] = ElementMatrixKe(X,Y,E0,nu)
                D = E0/(1-nu^2)*[1  nu  0;  nu  1  0;  0  0  (1-nu)/2] ;
                J = [X(1)-X(3) Y(1)-Y(3) ; X(2)-X(3) Y(2)-Y(3)] ;
                Be = 1/det(J)*[J(2,2) 0 -J(1,2) 0 -J(2,2)+J(1,2) 0; 0 -J(2,1) 0 J(1,1) 0 J(2,1)-J(1,1); -J(2,1) J(2,2) J(1,1) -J(1,2) J(2,1)-J(1,1) -J(2,2)+J(1,2)];
                Ke = 1/2*det(J)*Be'*D*Be ;
            end

            %% FINITE ELEMENT ANALYSIS
            function [ dCe, J] = FiniteElementAnalysis( t, p, BDY, x, E0, Emin, nu, penal, loadcase)
                EM = @CBFO_OC.ElementMatrixKe ;
                NT = length(t); KK = zeros( 6, 6*NT) ;
                for pi = 1:NT
                    KK(:,6*pi-5:6*pi) = (Emin+x(pi).^penal*(E0-Emin)) * EM(p(t(pi,:),1),p(t(pi,:),2),1,nu) ;
                end
                elemDof = zeros(NT,6); elemDof(:,[1 3 5]) = 2*t-1; elemDof(:,[2 4 6]) = 2*t;
                iK = reshape( kron( elemDof, ones(6,1))', 36*NT,1) ;
                jK = reshape( kron( elemDof, ones(1,6))', 36*NT,1) ;
                sK = reshape( KK, 36*NT, 1) ;
                NK = sparse( iK, jK, sK, 2*length(p), 2*length(p)) ;
                NK = (NK+NK')/2 ;

                %% LOAD CASES
                if loadcase == 1
                    fixedNodes = find(p(:,1)==BDY(1,1)) ;
                    forceNodes = find(p(:,1)==BDY(2,1) & p(:,2)==(BDY(1,2)+BDY(2,2))/2) ;
                    fixedDof = [2*fixedNodes-1; 2*fixedNodes] ;
                    AllDof = 1:2*length(p) ;
                    freeDofs = setdiff( AllDof, fixedDof) ;
                    F = sparse( 2*forceNodes, 1, -100, 2*length(p), 1) ;
                elseif loadcase == 2
                    fixedNodes1 = find(p( :, 1) == BDY( 1, 1)) ;
                    fixedDof1 = 2*fixedNodes1 - 1 ;
                    fixedNodes2 = find(p(:,1) == BDY(2,1) & p(:,2) == BDY(2,2)) ;
                    fixedDof2 = 2*fixedNodes2 ;
                    fixedDof = [fixedDof1;fixedDof2] ;
                    forceNodes = find(p(:,1)==BDY(1,1) & p(:,2)==BDY(2,2)) ;
                    AllDof = 1:2*length(p) ;
                    freeDofs = setdiff( AllDof, fixedDof) ;
                    F = sparse( 2*forceNodes, 1, -50, 2*length(p), 1) ;
                end
                U = zeros( 2*length(p), 1) ;
                U( freeDofs, :) = NK(freeDofs,freeDofs) \ F(freeDofs,1) ;
                for pi = 1 : NT
                    Ce(pi) = sum((U(elemDof(pi,:))'*KK(:,6*pi-5:6*pi)).*U( elemDof(pi,:))', 2) ;
                    dCe(pi) = -penal*(E0-Emin)*x(pi).^(penal-1).*(Ce(pi)./(Emin+x(pi).^penal*(E0-Emin)));
                end
                J = 0.5.*sum(Ce); dCe = dCe';
            end

            %% OPTIMALITY CRITERIA UPDATE OF DESIGN VARIABLES AND PHYSICAL DENSITIES
            function [ xnew, change] = OptimalityCriteriaUpdate( x, dCe, Ve, volfrac, t, p, loop, nelx, nely)
                l1 = 0 ; l2 = 1e9 ;  move = 0.5 ;
                while (l2-l1)/(l1+l2) > 1e-3
                    lmid = 0.5*(l2+l1) ;
                    xnew = max( 0, max( x-move, min( 1, min( x+move, x.*sqrt(-dCe/lmid))))) ;
                    xPhys = xnew ;
                    if sum(xPhys(:).*Ve(:)) > volfrac*nelx*nely, l1 = lmid; else l2 = lmid ; end
                end
                change = max(abs(xnew(:)-x(:))) ;
                colormap summer ;
                clf; patch( 'Faces', t, 'Vertices', p, 'FaceVertexCData', xnew, 'FaceColor', 'flat') ;
                colorbar ; hold on ; axis off equal tight ;set( gca, 'YDir', 'reverse');
                drawnow ;
                %         saveas(gcf,['./fig', int2str(loop) '.png']);
            end
        end
    end
