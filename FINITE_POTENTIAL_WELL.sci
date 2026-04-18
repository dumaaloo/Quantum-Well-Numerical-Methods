clc
clear
me=9.1*((10)^(-31))
e=1.6*(10^(-19))
h=1.054*(10^(-34))
n=input('enter the no. of intervals:')
a1=-2*(10^-10)
a2=-a1//input('enter the ending point of finite well:')
x1=2*a1//input('enter initial value of x:')
x2=2*a2//input('enter final value of x:')
dx=(x2-x1)/(n-1)
c=-(h^2)/(2*me*(dx^2)*(e))
a=ones(n-3,1)
ke=diag(a,1)+diag(a,-1)+(-2)*(eye(n-2,n-2))
x=linspace(x1,x2,n)
v=zeros(n,1)
for i=1:n
    if abs(x(i))<=a2
       v(i)=-100
   end 
end
pe=diag(v(2:n-1))
H=(ke.*c)+pe
[A,B]=spec(H)
count=0
for i=1:n-2
    if B(i,i)<0
        count=count+1
    else
        break
    end
end
//normalisation-----------------------
figure
for i=1:count
    Y=A(:,i).^2
    N(i)=inttrap(x(2:n-1),Y)
    wfn(:,i)=A(:,i)/sqrt(N(i))//------------------normalised wavefunction
    subplot(2,2,i)
    plot(x(2:n-1),wfn(:,i))
    xtitle("normalised wavefunctions psi ("+string(i)+string(")"))
end

figure
for i=1:count
    subplot(2,2,i)
    plot(x(2:n-1),(wfn(:,i).*wfn(:,i)))
    xtitle("wavefunctions squared |psi^2| ("+string(i)+string(")"))
end

//expectation values-------------------
for i=1:count
    Y=wfn(:,i).*(x(2:n-1)').*wfn(:,i)
    Y1=wfn(:,i).*(x(2:n-1)').*(x(2:n-1)').*wfn(:,i)
    Y4=wfn(:,i).*(v(2:n-1)).*wfn(:,i)
    exp_x(i)=inttrap(x(2:n-1),Y)
    exp_x2(i)=inttrap(x(2:n-1),Y1)
    exp_v(i)=inttrap(x(2:n-1),Y4)
end
xnew=x(2:n-1)
for j=1:count
for i=1:n-3
	xmid(i)=(xnew(i)+xnew(i+1))/2
	dA(i,j)=(wfn(i+1,j)-wfn(i,j))/dx
	Amid(i,j)=(wfn(i+1,j)+wfn(i,j))/2
end
end
xmid2=x(3:n-2)
for j=1:count
for i=1:n-4
	d2A(i,j)=(wfn(i+2,j)-2*wfn(i+1,j)+wfn(i,j))/(dx^2)
    Amid2(i,j)=(Amid(i+1,j)+Amid(i,j))/2
end
end
for j=1:count
	Y3=Amid(:,j).*dA(:,j)
    Y2=Amid2(:,j).*d2A(:,j)
	exp_p(j)=(-%i)*h*(inttrap(xmid,Y3))
    exp_p2(j)=-1*(h^2)*(inttrap(xmid2,Y2))
end
//verification of energy conservation--------------------
exp_ke=(exp_p2./(2*me))/e
E=exp_ke+exp_v
//Uncertainity principle---------------------------------
sdx=sqrt(exp_x2-(exp_x).^2)
sdp=sqrt(exp_p2-(exp_p).^2)
ucp=(sdx.*sdp*2)/h
//potential well plot------------------------------------
eval=diag(B)
figure
t=gca()
t.data_bounds=[-0.2*(10^-9),-101;0.2*(10^-9),1]
plot(x,v)
for i=1:count
    xstring(x2-0.3*10^(-10),eval(i),string(eval(i)))
    plot(x,eval(i),'r-')
end
title("Finite Potential Well")
//errror in energy eigen values
err=(abs(E(1:count)-eval(1:count))./eval(1:count))*100
