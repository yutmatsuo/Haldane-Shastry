using LinearAlgebra

"""
This program computes a generalization of Haldane-Shastry for the vector case
namely, su(2) k=2.
"""
function spin(a,n)
   s=zeros(Int8,n)
    for i=1:n
        s[i]=(rem(a,3^i)-rem(a,3^(i-1)))/3^(i-1)
    end
    return s
end
"""
spin configuration of nxn site is described by 3^n dim vector.
spin(a,n) maps a number in a=1 to 3^n to spin configuration (s_1,..., s_n) with s_i=0,1,2.
s_i is obtained as spin(a,n)[i] namely i-th component of spin(a,n).
"""

jj=[4 0 0 0 0 0 0 0 0;0 2 0 2 0 0 0 0 0;0 0 0 0 4 0 0 0 0; 0 2 0 2 0 0 0 0 0;0 0 1 0 2 0 1 0 0;0 0 0 0 0 2 0 2 0;0 0 0 0 4 0 0 0 0;0 0 0 0 0 2 0 2 0;0 0 0 0 0 0 0 0 4]

"""
jj is the matrix representing J_{ab}J_{ba} for a-th and b-th spin state.
It is obtained in a separate Mathematica file "spin 1 case.nb"
"""

function dd(a,n,i,j)
    return a-spin(a,n)[i]*3^(i-1)-spin(a,n)[j]*3^(j-1)
end

"""
dd(a,n,i,j) gives the spin state obtained from a by putting 0 into the i-th and j-th component.
In Haldane-Shastry Hamiltonian, this is necessary to give a nonvanishing matrix component between
i-th and j-th component.
"""

function hs(n)
    h=zeros(3^n,3^n)
    for a=0:3^n-1
        for b=0:3^n-1
            for i=2:n
                for j=1:i-1
                    if dd(a,n,i,j)==dd(b,n,i,j)
                        h[a+1,b+1]=h[a+1,b+1]+jj[spin(a,n)[i]+3*spin(a,n)[j]+1,spin(b,n)[i]+3*spin(b,n)[j]+1,]/(4*(sin(pi*(i-j)/n))^2)
                    end
                end
            end
        end
    end
    return h
end

"""
hs(n) gives the (generalized) Haldane-Shastry Hamiltonian where 1-K_{ab} is replaced by J_{ab}J_{ba}
"""

function hse(n)
    return eigvals(hs(n))
end
"""
hse(n) gives the integer part of the eigenvalue multiplied by 2.
"""

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

"""
frq(x) count the number of each integer in x.
frq(hse(n)) gives the eigenvalue and its degeneracy.
"""

function printval(v)
    j=1
    for i in v
        print(j,": ",i,"\n")
        j=j+1
    end
end
