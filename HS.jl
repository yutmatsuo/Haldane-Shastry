using LinearAlgebra

function spin(a,n)
   s=zeros(Int8,n)
    for i=1:n
        s[i]=(rem(a,2^i)-rem(a,2^(i-1)))/2^(i-1)
    end
    return s
end


function hs(n)
    h=zeros(2^n,2^n)
    for a=0:2^n-1
        for b=0:2^n-1
            spa=spin(a,n)
            spb=spin(b,n)

            for i=2:n
                for j=1:i-1
                    if spa[i]==spb[i]&& spa[j]==spb[j]
                        h[a+1,b+1]=h[a+1,b+1]+1/(4*sin(pi*(i-j)/n)^2);
                    elseif spa[i]==spb[j]&& spa[j]==spb[i]
                        h[a+1,b+1]=h[a+1,b+1]-1/(4*sin(pi*(i-j)/n)^2);
                    end
                end
            end
        end
    end
    return h
end

function hse(n)
    return eigvals(hs(n))
end
