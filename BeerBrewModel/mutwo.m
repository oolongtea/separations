function sr2=mutwo(a, M, G, twomu)
yp2=(a(1).*M.*a(3))./((a(2)+M).*(a(3)+G));
sr2=sum((twomu-yp2).^2);