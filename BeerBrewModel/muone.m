function sr=muone(a, G, onemu)
yp=(a(1).*G)./(a(2)+G);
sr=sum((onemu-yp).^2);

