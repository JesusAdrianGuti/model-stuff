
U_Tilda = zeros(n+1,1);

%Adams-Moulton Method

for i=5:n+1

    U_Tilda(i)=u(i-1)+(h/24)*(55*f(t(i-1),u(i-1))-59*f(t(i-2),u(i-2))+37*f(t(i-3),u(i-3))-9*f(t(i-4),u(i-4)));

    u(i)=u(i-1)+(h/24)*(9*f(t(i),U_Tilda(i))+19*f(t(i-1),u(i-1))-5*f(t(i-2),u(i-2))+f(t(i-3),u(i-3)));

end

 

plot(t,exact,'g',t,u,'b--o');

title('when h=0.5');

 

%when h=0.1

h=0.1;

t=lower:h:upper;

n=(upper-lower)/h;

exact = zeros(n+1,1);

u = zeros(n+1,1);

%exact equation

for i=1:n+1

    exact(i)=(1/4)+(1/4)*exp(-4*t(i));

end

 

%Runge-Kutta Method

u(1)=0.5;

for i=1:3

    k1=h*f(t(i),u(i));

    k2=h*f(t(i)+(h/2),u(i)+(k1/2));

    k3=h*f(t(i)+(h/2),u(i)+(k2/2));

    k4=h*f(t(i)+(h),u(i)+(k3));

    

    u(i+1)=u(i)+(1/6)*(k1+2*k2+2*k3+k4);

end

 

U_Tilda = zeros(n+1,1);

%Adams-Moulton Method

for i=5:n+1

    U_Tilda(i)=u(i-1)+(h/24)*(55*f(t(i-1),u(i-1))-59*f(t(i-2),u(i-2))+37*f(t(i-3),u(i-3))-9*f(t(i-4),u(i-4)));

    u(i)=u(i-1)+(h/24)*(9*f(t(i),U_Tilda(i))+19*f(t(i-1),u(i-1))-5*f(t(i-2),u(i-2))+f(t(i-3),u(i-3)));

end

 

figure;

plot(t,exact,'p',t,u,'b--o');

title('when h=0.1');

 

%when h=0.01;

h=0.01;

t=lower:h:upper;

n=(upper-lower)/h;

exact = zeros(n+1,1);

u = zeros(n+1,1);

%exact equation

for i=1:n+1

    exact(i)=(1/4)+(1/4)*exp(-4*t(i));

end

 

%Runge-Kutta Method

u(1)=0.5;

for i=1:3

    k1=h*f(t(i),u(i));

    k2=h*f(t(i)+(h/2),u(i)+(k1/2));

    k3=h*f(t(i)+(h/2),u(i)+(k2/2));

    k4=h*f(t(i)+(h),u(i)+(k3));

    

    u(i+1)=u(i)+(1/6)*(k1+2*k2+2*k3+k4);

end

 

U_Tilda = zeros(n+1,1);

%Adams-Moulton Method

for i=5:n+1

    U_Tilda(i)=u(i-1)+(h/24)*(55*f(t(i-1),u(i-1))-59*f(t(i-2),u(i-2))+37*f(t(i-3),u(i-3))-9*f(t(i-4),u(i-4)));

    u(i)=u(i-1)+(h/24)*(9*f(t(i),U_Tilda(i))+19*f(t(i-1),u(i-1))-5*f(t(i-2),u(i-2))+f(t(i-3),u(i-3)));

end

 

figure;

plot(t,exact,'r',t,u,'b--o');

title('when h=0.01');

 

%when h=0.001

h=0.001;

t=lower:h:upper;

n=(upper-lower)/h;

exact = zeros(n+1,1);

u = zeros(n+1,1);

%exact equation

for i=1:n+1

    exact(i)=(1/4)+(1/4)*exp(-4*t(i));

end

 

%Runge-Kutta Method

u(1)=0.5;

for i=1:3

    k1=h*f(t(i),u(i));

    k2=h*f(t(i)+(h/2),u(i)+(k1/2));

    k3=h*f(t(i)+(h/2),u(i)+(k2/2));

    k4=h*f(t(i)+(h),u(i)+(k3));

    

    u(i+1)=u(i)+(1/6)*(k1+2*k2+2*k3+k4);

end

 

U_Tilda = zeros(n+1,1);

%Adams-Moulton Method

for i=5:n+1

    U_Tilda(i)=u(i-1)+(h/24)*(55*f(t(i-1),u(i-1))-59*f(t(i-2),u(i-2))+37*f(t(i-3),u(i-3))-9*f(t(i-4),u(i-4)));

    u(i)=u(i-1)+(h/24)*(9*f(t(i),U_Tilda(i))+19*f(t(i-1),u(i-1))-5*f(t(i-2),u(i-2))+f(t(i-3),u(i-3)));

end

 

figure;

plot(t,exact,'o',t,u,'b--o');

title('when h=0.001');