function dydt=Lindbladcav(t,y)
global sdownc supc annic creac nn gamma kappa gammap H sz
rho = reshape(y,[2*nn,2*nn]);
drhodt = -1i*(H*rho-rho*H)+0.5*gammap*(sz*rho*sz-0.5*(sz*sz*rho+rho*sz*sz))+gamma*(sdownc*rho*supc-0.5*(supc*sdownc*rho+rho*supc*sdownc))+kappa*(annic*rho*creac-0.5*(creac*annic*rho+rho*creac*annic));
dydt = reshape(drhodt,[4*nn*nn,1]);
end
