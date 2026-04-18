clear
me=9.1*((10)^(-31))
e=1.6*(10^(-19))
h=1.054*(10^(-34))
a1=1*(10^-10)
v0=500*e
z0=a1*sqrt(2*me*v0/(h^2))
E0=input('enter the initial guess for energy eigen value:')
E0_j=E0*e
function a=f(x)
    a=tan(x)-sqrt(((z0/x)^2)-1)
endfunction
function a=g(x)
    a=-(1/tan(x))-sqrt(((z0/x)^2)-1)
endfunction
function a=f1(x)
    a=((1/cos(x))^2)+z0*z0/(sqrt(((z0/z)^2)-1)*(z^3))
endfunction
function a=g1(x)
    a=((1/sin(x))^2)+z0*z0/sqrt(((z0/z)^2)-1)*(z^3)
endfunction
z=a1*sqrt(2*me*abs(v0-E0_j)/(h*h))
while(1)
    i=1
    while(1)
        zmin=(i-1)*%pi/2
        zmax=i*%pi/2
        if z>zmin && z<zmax
            break
        else
            i=i+1
        end
    end
    if modulo(i,2)==0 then
        znew=z-(g(z)/g1(z)) 
        //disp(znew)
    else
        znew=z-(f(z)/f1(z))
    end
    E=v0-((h*h*z*z)/(a1*a1*2*me))
    Enew=v0-((h*h*znew*znew)/(a1*a1*2*me))
    diff_ev=abs(E-Enew)/e
    if diff_ev > 10^-5 then
        z=znew
    else
        disp(Enew/e)
        break
    end
end
k=1
while(1)
    zmin=(k-1)*%pi/2
    zmax=k*%pi/2
    if z0>zmin && z0<zmax
        break
    else
        k=k+1
    end
end
figure
scf()
x=linspace(0,k*%pi/2,500)
dx=(k*%pi/2)/500
j=1

while j<=k
    if modulo(j,2)==0
        i=(((j-1)*%pi/2)):dx:((j*%pi/2))
        plot(i,-cotg(i),'b-*')
    else
        i=(((j-1)*%pi/2)):dx:((j*%pi/2))
        plot(i,tan(i),'r-o')
    end
    j=j+1 
end
plot(x(1:20:500),sqrt(((z0./x(1:20:500)).^2)-1),'g-sq')
t=gca()
t.data_bounds=[0,0;k*%pi/2,10]
curve=t.children()
h1=curve(1)
h2=curve(k)
h3=curve(k+1)
legend([h3,h2,h1],['tan(x)','cot(x)','sqrt((z0/z)^2-1)'])

/*for i=1:k
    if modulo(i,2)==0 && x<i*%pi && x>(i-1)*%pi
        plot(x,-cotg(x))
    elseif modulo(i,2)==1 && x<i*%pi && x>(i-1)*%pi
        plot(x,tan(x))
    end
end


