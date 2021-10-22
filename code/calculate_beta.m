function beta = calculate_beta(device, N, type)
    arguments
        device (1,1) struct
        N (1,1) {mustBeInteger, mustBePositive} = 200
        type {mustBeMember(type,{'homo', 'regul', 'LGS', 'triag', 'matica', 'resive'})} = 'resive'
    end
    
    n=1:N;
    
    switch type
        case 'homo'
            beta = 1./(device.sigmaCerebrum.*(2*n+1));
        case 'regul'
            beta = zeros(size(n));
            for k=1:N  
                % Matrix des LGS
                A(1,:)=[device.radiusCerebrum^(2*k+1) -device.radiusCerebrum^(2*k+1) -1 0 0 0 0];
                A(2,:)=[k*device.radiusCerebrum^(2*k+1)*device.sigmaCerebrum -k*device.radiusCerebrum^(2*k+1)*device.sigmaFluid (k+1)*device.sigmaFluid 0 0 0 0];
                A(3,:)=[0 device.radiusFluid^(2*k+1) 1 -device.radiusFluid^(2*k+1) -1 0 0];
                A(4,:)=[0 k*device.radiusFluid^(2*k+1)*device.sigmaFluid -(k+1)*device.sigmaFluid -k*device.radiusFluid^(2*k+1)*device.sigmaBone (k+1)*device.sigmaBone 0 0];
                A(5,:)=[0 0 0 device.radiusBone^(2*k+1) 1 -device.radiusBone^(2*k+1) -1];
                A(6,:)=[0 0 0 k*device.radiusBone^(2*k+1)*device.sigmaBone -(k+1)*device.sigmaBone -k*device.radiusBone^(2*k+1)*device.sigmaScalp (k+1)*device.sigmaScalp];
                A(7,:)=[0 0 0 0 0 k*device.radiusScalp^(2*k+1) -(k+1)];
                b=[-1/(device.sigmaCerebrum*(2*k+1)) (k+1)/(2*k+1) 0 0 0 0 0]';
                A=A+diag(max(max(abs(A)))*0.1*ones(1,7));
                x=A\b;
                beta(k)=x(7);
            end   
        case 'LGS'
            beta = zeros(size(n));
            for k=1:N  
                % Matrix des LGS
                A(1,:)=[device.radiusCerebrum^(2*k+1) -device.radiusCerebrum^(2*k+1) -1 0 0 0 0];
                A(2,:)=[k*device.radiusCerebrum^(2*k+1)*device.sigmaCerebrum -k*device.radiusCerebrum^(2*k+1)*device.sigmaFluid (k+1)*device.sigmaFluid 0 0 0 0];
                A(3,:)=[0 device.radiusFluid^(2*k+1) 1 -device.radiusFluid^(2*k+1) -1 0 0];
                A(4,:)=[0 k*device.radiusFluid^(2*k+1)*device.sigmaFluid -(k+1)*device.sigmaFluid -k*device.radiusFluid^(2*k+1)*device.sigmaBone (k+1)*device.sigmaBone 0 0];
                A(5,:)=[0 0 0 device.radiusBone^(2*k+1) 1 -device.radiusBone^(2*k+1) -1];
                A(6,:)=[0 0 0 k*device.radiusBone^(2*k+1)*device.sigmaBone -(k+1)*device.sigmaBone -k*device.radiusBone^(2*k+1)*device.sigmaScalp (k+1)*device.sigmaScalp];
                A(7,:)=[0 0 0 0 0 k*device.radiusScalp^(2*k+1) -(k+1)];
                b=[-1/(device.sigmaCerebrum*(2*k+1)) (k+1)/(2*k+1) 0 0 0 0 0]';
                x=A\b;
                beta(k)=x(7);
            end 
        case 'triag'
            beta = zeros(size(n));
            device.sigma = [device.sigmaCerebrum, device.sigmaFluid, device.sigmaBone, device.sigmaScalp];
            device.radius = [device.radiusCerebrum, device.radiusFluid, device.radiusBone, device.radiusScalp];
            i = 1;
            for k=1:N 
                A(1,:) = [(k+1)*device.sigma(i+1) + k*device.sigma(i+0), -(2*k+1)*device.sigma(i+1), 0,0,0,0,0]; 
                A(2,:) = [k*device.sigma(i+0), -k*device.sigma(i+1), (k+1)*device.sigma(i+1)/device.radius(i+0)^(2*k+1), 0, 0, 0, 0];
                A(3,:) = [0, (k+1)*device.sigma(i+2) + k*device.sigma(i+1), (k+1)*(device.sigma(i+2)-device.sigma(i+1))/device.radius(i+1)^(2*k+1), -(2*k+1)*device.sigma(i+2), 0,0,0];
                A(4,:) = [0, 0, -device.sigma(i+1)*(2*k+1)/device.radius(i+1)^(2*k+1), k*(device.sigma(i+1)-device.sigma(i+2)), ((k+1)*device.sigma(i+2)+k*device.sigma(i+1))/device.radius(i+1)^(2*k+1), 0,0];
                A(5,:) = [0, 0, 0, (k+1)*device.sigma(i+3) + k*device.sigma(i+2), ((k+1)*(device.sigma(i+3)-device.sigma(i+2)))/device.radius(i+2)^(2*k+1), -(2*k+1)*device.sigma(i+3), 0];
                A(6,:) = [0, 0, 0, 0, -device.sigma(i+2)*(2*k+1)/device.radius(i+2)^(2*k+1), k*(device.sigma(i+2)-device.sigma(i+3)), ((k+1)*device.sigma(i+3) + k*device.sigma(i+2))/device.radius(i+2)^(2*k+1)];
                A(7,:) = [0, 0, 0, 0, 0, device.radius(i+3)^(2*k+1)*k, -(k+1)];
                b = [(k+1)*(device.sigma(i+0) - device.sigma(i+1))/(device.radius(i+0)^(2*k+1)*(2*k+1)*device.sigma(i+0)), (k+1)/(device.radius(i+0)^(2*k+1)*(2*k+1)), 0,0,0,0,0]';
                x=A\b;
                beta(k)=x(7);
            end
        case 'matica'
            beta = device.radiusBone.*device.sigmaBone.*device.radiusFluid.*device.sigmaFluid.*(1+2.*n).^2.*device.radiusScalp.*(device.radiusBone.^2.*device.radiusCerebrum.*device.radiusFluid.^2.*device.radiusScalp).^(2.*n).*(device.radiusBone.^(2+4.* ...
            n).*(1+n).*(device.radiusCerebrum.*(device.radiusBone.*device.radiusCerebrum.^2.*device.radiusFluid).^(2.*n).*(device.sigmaBone+(-1).*device.sigmaFluid).*((-1).*device.sigmaCerebrum+device.sigmaFluid).*n.* ...
            (1+n)+device.radiusFluid.*(device.radiusBone.*device.radiusCerebrum.*device.radiusFluid.^2).^(2.*n).*(device.sigmaBone+(device.sigmaBone+device.sigmaFluid).*n).*(device.sigmaFluid+(device.sigmaCerebrum+device.sigmaFluid).*n)).*(device.sigmaBone+( ...
            -1).*device.sigmaScalp)+(-1).*device.radiusFluid.*n.*(1+n).*((-1).*device.radiusCerebrum.*(device.radiusCerebrum.*device.radiusFluid).^(4.*n).*(device.sigmaCerebrum+(-1).*device.sigmaFluid).* ...
            (device.sigmaFluid+(device.sigmaBone+device.sigmaFluid).*n)+device.radiusCerebrum.^(2.*n).*device.radiusFluid.^(1+6.*n).*(device.sigmaBone+(-1).*device.sigmaFluid).*(device.sigmaFluid+(device.sigmaCerebrum+device.sigmaFluid).*n)).* ...
            device.radiusScalp.*(device.radiusBone.*device.radiusScalp).^(2.*n).*(device.sigmaBone+(-1).*device.sigmaScalp)+device.radiusBone.*device.radiusCerebrum.^(2.*n).*device.radiusFluid.*(device.radiusBone.*device.radiusFluid).^(4.*n).*( ...
            1+n).*(device.radiusCerebrum.^(1+2.*n).*(device.sigmaCerebrum+(-1).*device.sigmaFluid).*(device.sigmaFluid+(device.sigmaBone+device.sigmaFluid).*n)+device.radiusFluid.^(1+2.*n).*((-1).* ...
            device.sigmaBone+device.sigmaFluid).*(device.sigmaFluid+(device.sigmaCerebrum+device.sigmaFluid).*n)).*(device.sigmaBone+device.sigmaBone.*n+n.*device.sigmaScalp)+device.radiusBone.*(device.radiusCerebrum.*(device.radiusBone.*device.radiusCerebrum.^2.*device.radiusFluid).^(2.*n).*( ...
            device.sigmaBone+(-1).*device.sigmaFluid).*((-1).*device.sigmaCerebrum+device.sigmaFluid).*n.*(1+n)+device.radiusFluid.*(device.radiusBone.*device.radiusCerebrum.*device.radiusFluid.^2).^(2.*n).*(device.sigmaBone+(device.sigmaBone+ ...
            device.sigmaFluid).*n).*(device.sigmaFluid+(device.sigmaCerebrum+device.sigmaFluid).*n)).*device.radiusScalp.*(device.radiusBone.*device.radiusScalp).^(2.*n).*(device.sigmaScalp+n.*(device.sigmaBone+device.sigmaScalp))).^(-1);
        case 'resive'
            beta = zeros(size(n));
            sigma = [device.sigmaCerebrum, device.sigmaFluid, device.sigmaBone, device.sigmaScalp];
            radius = [device.radiusCerebrum, device.radiusFluid, device.radiusBone, device.radiusScalp];
            for k=1:N
               M = eye(2);
               for l=1:3
                   Mtmp = [k+1+k*sigma(l)/sigma(l+1), (k+1)*(1-sigma(l)/sigma(l+1))*radius(l)^(-(2*k+1));...
                      k*(1-sigma(l)/sigma(l+1))*radius(l)^(2*k+1), k+(k+1)*sigma(l)/sigma(l+1)];
                   M = Mtmp*M;
               end
               M = M./(2*k+1)^3;
               beta(k) = k/(sigma(end)*(2*k+1)*(k*M(1,1) - (k+1)*M(2,1)*radius(end)^(-(2*k+1)))); 
            end
    end
    
    n=1:sum(not(isnan(beta)));
    beta = beta(1:length(n));
end