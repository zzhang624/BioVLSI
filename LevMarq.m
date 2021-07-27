function [x_out] = LevMarq(Phi,dt,NumExp,x_guess,plotflag)

%NumExp=1;
y=Phi;
t=[0:length(y)-1].';
%t=dt*[0:length(y)-1].';
%t=t/t(2);
%x_guess=[1;1E-4]; %[A tau]
x_guess=repmat(x_guess, NumExp,1);
tau=1;
k_max=300;
eps_1=1E-35;
eps_2=1E-35;

k=0;
nu=2;
x=x_guess;
M=zeros(length(t),1);
for ii=1:NumExp
    M=M+x(2*ii-1)*exp(-x(2*ii)*t);
end
%M=x(1)*exp(-x(2)*t);
f=y-M;
J=zeros(length(t),length(x));
for ii=1:NumExp
    J(:,2*ii-1)=-exp(-x(2*ii)*t);
    J(:,2*ii)=x(2*ii-1)*t.*exp(-x(2*ii)*t);
end
%J=[-exp(-x(2)*t) x(1)*t.*exp(-x(2)*t)];
A=J'*J;
g=J'*f;
flag=(norm(g,Inf)<=eps_1);
mu=tau*max(max(A.*eye(size(A))));

while ((flag==0)&&(k<k_max))
    k=k+1;
    h=-inv(A+mu*eye(size(A)))*g;
    if (norm(h)<= eps_2*(norm(x)+eps_2))
        flag=1;
    else
        x_new=x+h;
        F_x=0.5*f'*f;
        %M_new=x_new(1)*exp(-x_new(2)*t);
        M_new=zeros(length(t),1);
        for ii=1:NumExp
            M_new=M_new+x_new(2*ii-1)*exp(-x_new(2*ii)*t);
        end
        f_x_new=y-M_new;
        %f_x_new=f+J*h;
        F_x_new=0.5*(f_x_new')*f_x_new;
        L_diff=0.5*h'*(mu*h-g);
        rho=(F_x-F_x_new)/L_diff;
        if rho>0
            x=x_new;
            M=zeros(length(t),1);
            for ii=1:NumExp
                M=M+x(2*ii-1)*exp(-x(2*ii)*t);
            end
            %             M=x(1)*exp(-x(2)*t);
            f=y-M;
            J=zeros(length(t),length(x));
            for ii=1:NumExp
                J(:,2*ii-1)=-exp(-x(2*ii)*t);
                J(:,2*ii)=x(2*ii-1)*t.*exp(-x(2*ii)*t);
            end
            %J=[-exp(-x(2)*t) x(1)*t.*exp(-x(2)*t)];
            A=J'*J;
            g=J'*f;
            flag=(norm(g,Inf)<=eps_1);
            mu=tau*max(1/3,abs(1-(2*rho-1)^3));
            nu=2;
        else
            mu=mu*nu;
            nu=2*nu;
        end
        
    end
end

x_out=zeros(size(x));
for ii=1:NumExp
    x_out(2*ii-1)=x(2*ii-1);
    x_out(2*ii)=x(2*ii)/dt;
end
%x_out=[x(1),x(2)/Data.dt];

if plotflag==true
figure
plot(t/1E-12,y,'LineWidth',2)
hold on
plot(t/1E-12,y-M,'--','LineWidth',2)
plot(t/1E-12,M,'-','LineWidth',2)
xlabel('Time (ps)')
ylabel("Normalized Autocorrelation Function for 1ubq")
title('Levenberg-Marquardt Fit of Autocorrelation')
legend({'Original Signal','Residual','Best Exponential Fit'})
ACLfig('fullslide')
grid on
saveas(gcf, 'auto_1ubq_corrected.png')
end
end