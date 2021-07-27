Data=load('PhiTest.mat');
y=Data.Phi0;
t=Data.dt*[0:length(y)-1].';
t=t/t(2);
x_guess=[1;1E-3]; %[A tau]

tau=1;
k_max=200;
eps_1=1E-35;
eps_2=1E-35;

k=0;
nu=2;
x=x_guess;
M=x(1)*exp(-x(2)*t);
f=y-M;
J=[-exp(-x(2)*t) x(1)*t.*exp(-x(2)*t)];
A=J'*J;
g=J'*f;
flag=(norm(g,Inf)<=eps_1);
mu=tau*max(max(A.*eye(size(A))));

while ((flag==0)&&(k<k_max))
    k=k+1;
    h=-inv(A+mu*eye(size(A)))*g;
    if norm(h)<= eps_2*(norm(x)+eps_2)
        flag=1;
    else
        x_new=x+h;
        F_x=0.5*f'*f;
        M_new=x_new(1)*exp(-x_new(2)*t);
        %f_x_new=y-M_new;
        f_x_new=f+J*h;
        F_x_new=0.5*(f_x_new')*f_x_new;
        L_diff=0.5*h'*(mu*h-g);
        rho=(F_x-F_x_new)/L_diff;
        if rho>0
            x=x_new;
            M=x(1)*exp(-x(2)*t);
            f=y-M;
            J=[-exp(-x(2)*t) x(1)*t.*exp(-x(2)*t)];
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

figure
plot(t*Data.dt/1E-12,y,'LineWidth',2)
hold on
plot(t*Data.dt/1E-12,y-M,'--','LineWidth',2)
plot(t*Data.dt/1E-12,M,'-','LineWidth',2)
xlabel('Time (ps)')
ylabel("Normalized Autocorrelation Function")
title('Levenberg-Marquardt Fit of Autocorrelation')
legend({'Original Signal','Residual','Best Exponential Fit'})
ACLfig('fullslide')
grid on
x_out=[x(1),x(2)/Data.dt];