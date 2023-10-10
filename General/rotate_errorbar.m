function rotate_errorbar(M,S,T,C,LW,D)

R = [cos(T) sin(T); -sin(T) cos(T)];
RT = [cos(pi/2) sin(pi/2);-sin(pi/2) cos(pi/2)];
RT

D = [1 0]*D;
DR = D*R*RT;

SR = [S 0]*R;




plot(M(1)+SR(1)*[-1 1],M(2)+SR(2)*[-1 1],'color',C,'LineWidth',LW)
plot(M(1)+SR(1)+DR(1)*[-1 1],M(2)+SR(2)+DR(2)*[-1 1],'color',C,'LineWidth',LW)
plot(M(1)-SR(1)+DR(1)*[-1 1],M(2)-SR(2)+DR(2)*[-1 1],'color',C,'LineWidth',LW)
