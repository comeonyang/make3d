function test_rodr
% TEST_RODR Regression test of RODR() and IRODR()
%   TEST_RODR() performs a regression tests of the RODR() and IRODR()
%   functions. It compares RODR(w) with EXMP(HAT(w)) and IRODR(R) with
%   IHAT(LOGM(R)). It also checks the derivatives numerically.
%
%   The function might isse a buch of warning message. This is
%   normal as some singular values are tested for completeness.
%
%   See also RODR(), IRODR().

% Special cases
fprintf('Testing [0 0 0]\n') ;
[er, edr, eir, edir] = test([0 0 0]') ;
report(er, edr, eir, edir) ;

% General cases
T=100;
er=[] ;
edr=[] ;
eir=[] ;
edir=[] ;
for t=1:T
  fprintf('Testing random vector %d of %d\r',t,T) ;
  v=rand(3,1) ;
  [er_, edr_, eir_, edir_] = test(v) ;  
  er   = [er   er_   ] ;
  edr  = [edr  edr_  ] ;
  eir  = [eir  eir_  ] ;
  edir = [edir edir_ ] ;
end
fprintf('\n') ;
report(er,edr,eir,edir) ;

% Around the clock
% Note that -pi and +pi are singular. Since they give
% the same rotation R, the inverse Rodrigues formula will map both
% cases back to +pi (by convention), while
% the logm function (used for the regression test here)
% is not defined for such cases, for which the error is big.
% Also, the derivatives of the inverse formula are not computed 
% there (even if it might be possible to do so on the local branch)
T=100;
er=[] ;
edr=[] ;
eir=[] ;
edir=[] ;
th_range=linspace(-pi,pi,T) ;
for t=1:T
  fprintf('Testing angle %f [rad]\r',th_range(t));
  v = th_range(t)*[1;0;0] ;
  [er_, edr_, eir_, edir_] = test(v) ;  
  er   = [er   er_   ] ;
  edr  = [edr  edr_  ] ;
  eir  = [eir  eir_  ] ;
  edir = [edir edir_ ] ;
end
report(er,edr,eir,edir) ;
fprintf('n') ;


% --------------------------------------------------------------------
%                                                     Helper functions
% --------------------------------------------------------------------
% Computes the rodrigues function and its inverse and compare the
% result with `numerical ground truth'. Compares the derivatives as well.
function [er, edr, eir, edir] = test(v)

 rodrigues  = @rodr ;
irodrigues = @irodr ;

% Numerical derivations
e1 = [1;0;0] ;
e2 = [0;1;0] ;
e3 = [0;0;1] ;
E{1} = [0  0 0 ; 0 0 -1 ;  0 1 0] ;
E{2} = [0  0 1 ; 0 0  0 ; -1 0 0] ;
E{3} = [0 -1 0 ; 1 0  0 ;  0 0 0] ;

step=1e-5;
RR = expm(hat(v)) ;
dR1 = (expm(hat( v + step * e1)) - RR)/step ;
dR2 = (expm(hat( v + step * e2)) - RR)/step ;
dR3 = (expm(hat( v + step * e3)) - RR)/step ;
dRR = [ dR1(:), dR2(:), dR3(:) ] ;

E{1} = RR*E{1} ;
E{2} = RR*E{2} ;
E{3} = RR*E{3} ;

duu=zeros(3,3) ;
step=1e-5;
for k=1:3
  duu(:,k) = ...
      (ihat(logm( RR + step * E{k})) - v)/step ;
end

% Rodrigues and inverse formulas
[R,dR] =  rodrigues( v ) ;
[u,du] = irodrigues( R ) ;

du = du * reshape([ E{:} ],9,3) ;

% Errors
er   = max(abs( R(:)- RR(:)))/max(abs(  R(:))+eps) ;
edr  = max(abs(dR(:)-dRR(:)))/max(abs(dRR(:))+eps) ;
eir  = max(abs( v   - u    ))/max(abs(v(:))+eps) ;
edir = max(abs(du(:)-duu(:)))/max(abs(duu(:))+eps) ; 

% Print statistics
function report(er, edr, eir, edir)
%keyboard
rng=1:length(er) ;
fprintf('rel(x,y) = max|x-y|/(max|y|+eps)\n') ;
fprintf('er   = rel(  rodrigues(u),  expm(hat(u)))\n') ;
fprintf('edr  = rel(  rodrigues(u), dexpm(hat(u)))\n') ;
fprintf('eir  = rel( irodrigues(R),  ihat(Log(R))\n') ;
fprintf('edir = rel( irodrigues(R), dihat(Log(R))\n') ;

fprintf('t     | er         |edr         |eir         |edir       \n') ;
fprintf('%5d | %10.3e | %10.3e | %10.3e | %10.3e \n', [rng;er;edr;eir;edir]) ;
