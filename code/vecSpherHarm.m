function vecSpherHarmOut = vecSpherHarm(type,N,phi, t)

t = t(:);
phi = phi(:);
phi = phi';

idx = find(abs(abs(t)-1) > 1e-10);
idx2 = find(abs(abs(t)-1) <= 1e-10);
[legDiff, leg] = assLegendreDerivative(N,t);

t = t';
vecSpherHarmOut = zeros((N+1)^2, 3, length(t));
   
    if type == 1
        x = sqrt(1-t.^2).*cos(phi);
        y = sqrt(1-t.^2).*sin(phi);
        z = t;
        
        spH = zeros((N+1)^2,length(t));
        counter = 0;
        for n = 0:N
            clear j;
            if n > 0
                j = -n:n;
                spH(counter+(1:n),:) = shiftdim(leg(n+1,n+1:-1:2,:),1).*sqrt(2).*cos(j(1:n)'*phi);    
                spH(counter+(n+2:(2*n+1)),:) = shiftdim(leg(n+1,2:n+1,:),1).*sqrt(2).*sin(j(n+2:end)'*phi);    
            end
            spH(counter+n+1,:) = leg(n+1,1,:);    
            counter = counter + 2*n+1;
        end
        vecSpherHarmOut(:,1,:) = repmat(x,(N+1)^2,1).*spH;
        vecSpherHarmOut(:,2,:) = repmat(y,(N+1)^2,1).*spH;
        vecSpherHarmOut(:,3,:) = repmat(z,(N+1)^2,1).*spH;
    else
        d2 = zeros((N+1)^2,length(idx2));
        c2 = zeros((N+1)^2,length(idx));
        c3 = zeros((N+1)^2,length(idx));
        d3 = zeros((N+1)^2,length(idx2));
        counter = 1;
        for n = 1:N 
            clear j;
            j = -n:n;

            % The case of t^2 = 1
            if ~isempty(idx2)
                %Calculate for t^2=1 the values in e_t direction. First for
                %the two j=-n, j=n, and j = 0.
                %and then for -n < j < n, with j ~=0.
            
                % j = -n
                % d3(counter+1,:) = zeros(size(idx2));
                % j = n
                % d3(counter+2*n+1,:) = zeros(size(idx2));
                %j =-1
                d3(counter+n,:) = -1*(sign(t(idx2))).^(n+1).*sqrt((2*n+1)*n*(n+1)/(2*pi)).*t(idx2).*cos(j(n)*phi(idx2))/2; 
                % j = 0
                d3(counter+n+1,:) = zeros(size(idx2)); 
                % j = 1
                d3(counter+n+2,:) = -1*(sign(t(idx2))).^(n+1).*sqrt((2*n+1)*n*(n+1)/(2*pi)).*t(idx2).*sin(j(n+2)*phi(idx2))/2; 
                if n > 2
                    % j = -n+1, ... , -2
                    d3(counter+(2:n-1),:) = (repmat(sqrt((n-j(2:n-1)'+1).*(n+j(2:n-1)')),1,length(idx2)).*shiftdim(leg(n+1,n-1:-1:2,idx2),1)).*sqrt(2).*cos(j(2:n-1)'*phi(idx2));
                    % j = 2,...,n-1
                    d3(counter+(n+3:2*n),:) = (repmat(sqrt((n-j(n+3:2*n)'+1).*(n+j(n+3:2*n)')),1,length(idx2)).*shiftdim(leg(n+1,2:n-1,idx2),1)).*sqrt(2).*cos(j(n+3:2*n)'*phi(idx2));
                end
                %For j =0 c2 and d2 are always equals to zero.
                %Calculate for t^2=1 the values in e_z direction for j = 0
                
                % For j = 0 d2=0
                % For |j| > 1 is d2 =0
                %For t^2=1 and j = -1 the e_phi direction
                d2(counter+n,:) = (sign(t(idx2))).^(n+1).*sqrt((2*n+1)*n*(n+1)/(2*pi))*j(n).*sin(-phi(idx2))/2;
                %And for j = 1
                d2(counter+(n+2),:) = (sign(t(idx2))).^(n+1).*sqrt((2*n+1)*n*(n+1)/(2*pi))*j(n+2).*cos(phi(idx2))/2;
            end
            
            %Calculate for t^2~=1 and j < 0 the values
            c2(counter+(1:n),:) = -repmat(j(1:n)',1,length(idx))*sqrt(2).*shiftdim(leg(n+1,(n+1):-1:2,idx),1).*sin(j(1:n)'*phi(idx))./repmat(sqrt(1-t(idx).^2),length(j(1:n)),1);
            c3(counter+(1:n),:) = repmat(sqrt(1-t(idx).^2),length(j(1:n)),1).*shiftdim(legDiff(n+1,n+1:-1:2,idx),1)*sqrt(2).*cos(j(1:n)'*phi(idx));
            
            %Calculate for t^2~=1 and j > 0 the values
            c2(counter+(n+2:2*n+1),:) = repmat(j(n+2:end)',1,length(idx))*sqrt(2).*shiftdim(leg(n+1,2:n+1,idx),1).*cos(j(n+2:end)'*phi(idx))./repmat(sqrt(1-t(idx).^2),length(j(n+2:end)),1);
            c3(counter+(n+2:2*n+1),:) = repmat(sqrt(1-t(idx).^2),length(j(n+2:end)),1).*shiftdim(legDiff(n+1,2:(n+1),idx),1)*sqrt(2).*sin(j(n+2:end)'*phi(idx));
            
            %For j =0 c2 and d2 are always equals to zero.
            %Calculate for j=0 and t^2~=1
            c3(counter+n+1,:) = sqrt(1-t(idx).^2).*shiftdim(legDiff(n+1,1,idx),1);
            
            counter = counter + 2*n+1;
        end
        
        counter = 1;
        c1 = zeros((N+1)^2, length(t));
        for n=1:N
           c1(counter+(1:2*n+1),:) = sqrt(1/(n*(n+1)))*ones(2*n+1,length(t));
           counter = counter + 2*n+1;
        end
        if type == 2
            vecSpherHarmOut(:,1,idx) = c1(:,idx).*((c2).*repmat(-sin(phi(idx)),(N+1)^2,1) - (c3).*repmat((t(idx).*cos(phi(idx))),(N+1)^2,1));
            vecSpherHarmOut(:,2,idx) = c1(:,idx).*((c2).*repmat(cos(phi(idx)),(N+1)^2,1)   - (c3).*repmat(t(idx).*sin(phi(idx)),(N+1)^2,1));
            vecSpherHarmOut(:,3,idx) = c1(:,idx).*((c3).*repmat(sqrt(1-t(idx).^2), (N+1)^2,1));

            if ~isempty(idx2)
                vecSpherHarmOut(:,1,idx2) = c1(:,idx2).*((d2).*repmat(-sin(phi(idx2)),(N+1)^2,1) - (d3).*repmat((t(idx2).*cos(phi(idx2))),(N+1)^2,1));
                vecSpherHarmOut(:,2,idx2) = c1(:,idx2).*((d2).*repmat(cos(phi(idx2)),(N+1)^2,1)   - (d3).*repmat(t(idx2).*sin(phi(idx2)),(N+1)^2,1));
                vecSpherHarmOut(:,3,idx2) = c1(:,idx2).*((d3).*repmat(sqrt(1-t(idx2).^2), (N+1)^2,1));     
            end
            vecSpherHarmOut(1,:,:) = NaN;
        elseif  type == 3
            vecSpherHarmOut(:,1,idx) = c1(:,idx).*(-(c3).*repmat(-sin(phi(idx)),(N+1)^2,1) - (c2).*repmat((t(idx).*cos(phi(idx))),(N+1)^2,1));
            vecSpherHarmOut(:,2,idx) = c1(:,idx).*(-(c3).*repmat(cos(phi(idx)),(N+1)^2,1)   - (c2).*repmat(t(idx).*sin(phi(idx)),(N+1)^2,1));
            vecSpherHarmOut(:,3,idx) = c1(:,idx).*((c2).*repmat(sqrt(1-t(idx).^2), (N+1)^2,1));

            if ~isempty(idx2)
                vecSpherHarmOut(:,1,idx2) = c1(:,idx2).*(-(d3).*repmat(-sin(phi(idx2)),(N+1)^2,1) - (d2).*repmat((t(idx2).*cos(phi(idx2))),(N+1)^2,1));
                vecSpherHarmOut(:,2,idx2) = c1(:,idx2).*(-(d3).*repmat(cos(phi(idx2)),(N+1)^2,1)   - (d2).*repmat(t(idx2).*sin(phi(idx2)),(N+1)^2,1));
                vecSpherHarmOut(:,3,idx2) = c1(:,idx2).*((d2).*repmat(sqrt(1-t(idx2).^2), (N+1)^2,1));     
            end
            vecSpherHarmOut(1,:,:) = NaN;
        else
            error('Type must be an integer between 1 and 3.');
        end
    end
end