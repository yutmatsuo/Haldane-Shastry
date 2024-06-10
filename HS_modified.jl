using LinearAlgebra
using StatsBase

function spin(a,p,n)
   s=zeros(Int8,n)
    for i=1:n
        s[i]=(rem(a,p^i)-rem(a,p^(i-1)))/p^(i-1)
    end
    return s
end

function k(i,j,a,p,n)
    return a+(spin(a,p,n)[j]-spin(a,p,n)[i])*(p^(i-1)-p^(j-1))
end

function ssum(p,n)
    s=0.
    for j=2:n
        for i=1:j-1
            s=s+1/(4*(sin(pi*(i-j)/n))^2)
        end
    end
    return s
end

function hs(p,n)
    h=zeros(p^n,p^n)
    for a=0:p^n-1
        for b=0:p^n-1
            if a==b
                h[a+1,b+1]=h[a+1,b+1]+ssum(p,n)
            end
            for i=2:n
                for j=1:i-1
                    if a==k(i,j,b,p,n)
                        h[a+1,b+1]=h[a+1,b+1]-1/(4*(sin(pi*(i-j)/n))^2)
                    end
                end
            end
        end
    end
    return h
end

function hse(p,n)
    return map(x->round(Int64,2*x),eigvals(hs(p,n))) 
end

function frq(x)
    m=maximum(x)
    v=zeros(m+1)
    for i in x
        v[i+1]+=1
    end
    for j=0:m
        if v[j+1]!=0
            print(j,": ",round(Int64,v[j+1]),"\n")
        end
    end
end
