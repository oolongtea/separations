function sr3=muthree(a, N, M, G, threemu, Kgprime)
yp=(a(1).*N.*Kgprime.*a(2))./((a(3)+N).*(Kgprime+G).*(a(2)+M));
sr3=sum((threemu-yp).^2);
