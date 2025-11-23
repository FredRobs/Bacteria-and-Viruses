%parametes
rho=28;
sigma=10;
stepsize=0.001;
init=[0.0,0.1,0.0];

%range of bifurcation
range= 0.4:0.0005:0.8;

%initalise poincare matrix
intersection=NaN*zeros(1000,length(range));
position=0;

%ODE options relative and absolute tolerance
options=odeset('RelTol',1e-5,'AbsTol',1e-5);

function model=lorenz(t,x,sigma,rho,beta)
    model=[sigma*(x(2)-x(1));
        x(1)*(rho-x(3))-x(2);
        x(1)*x(2)-beta*x(3)];
end

%Bifurcation loop
for beta=range
    %print value of beta + iterate position of matrix
    beta
    position=position+1;
    temp=1;

    %obtain model
    modelFunction=@(t,x)lorenz(t,x,sigma,rho,beta);
    
    %simulate system for given beta
    [t,x]=ode45(modelFunction,0:stepsize:1000,init,options);
    
    %consider time above transient time
    cutoff=t>200;
    Xout=x(cutoff,:);

    for i=2:length(Xout)
        %has it crossed section when x1=0?
        if Xout(i,1)<0 && Xout(i-1,1)>0
             %save intersection value
            intersection(temp,position)=Xout(i,3);
            temp=temp+1;
        end
    end
end


hold on
plot(range,intersection,'.k','MarkerSize',2)
xlabel('\beta')
ylabel('x_1')
set(gca,'fontsize',12)
set(gca,'fontweight','bold')
box on


% Code adapted from Lazaros Moysis (2025). 
% Bifurcation diagram for the Lorenz Chaotic system 
% (https://uk.mathworks.com/matlabcentral/fileexchange/156752-bifurcation-diagram-for-the-lorenz-chaotic-system),
% MATLAB Central File Exchange. Retrieved November 23, 2025.


%[appendix]
%---
%[metadata:view]
%   data: {"layout":"onright"}
%---
