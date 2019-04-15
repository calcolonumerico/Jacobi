function [x,niter,resrel] = jacobi(A,b,TOL,MAXITER)
%jacobi Summary of this function goes here
%   Detailed explanation goes here

 %Controlli numero input
   if(nargin<2)
        error("Inserire come input almeno A e b.")
   end

  
 %Controlli sulla matrice
 if(nargin==2)
    if(~ismatrix(A))
       error("Il primo input deve essere una matrice.")
    elseif(~issparse(A)||isempty(A))
       error("La matrice deve essere sparsa.")
    elseif(size(A,1)~=size(A,2))
        error("La matrice deve essere quadrata.")
    elseif(~isnumeric(A) ||~isreal(A)||any(any(~isfinite(A)))|| any(any(~isa(A,'double')))||any(any(isnan(A))))
      error("Gli elementi della matrice devono essere numeri reali double finiti.")
    end
    
    %Controlli sul vettore 
   if(~isvector(b)||isscalar(b)||isrow(b))
       error('Errore, b deve essere un vettore colonna.')
   elseif(length(b)~=length(A))
       error("Il vettore b e la matrice A hanno dimensioni diverse.")
   elseif(any(ischar(b))||any(~isfinite(b)) ||  any(~isreal(b))||any(~isa(b,'double'))||any(isnan(b)))
       error('Errore, b deve contenere reali finiti double.')
   end
    
     TOL=10^-6;
     MAXITER=500;
 end
  

  if(nargin==3) 
  %Verifica se TOL è uno scalare, e se lo è verifica se è ammissibile
   if(~isscalar(TOL))
       error('Errore, TOL deve essere uno scalare')
   elseif (~isfinite(TOL) ||  ~isreal(TOL) || ischar(TOL)||any(isnan(TOL)))
        error('Errore, TOL deve essere settato come un numero reale.')
   end
   MAXITER=500;
  end
   
if(nargin==4)
  %Verifica se MAXITER è uno scalare, e se lo è verifica se è ammissibile;
   if(~isscalar(MAXITER))
       error('Errore, MAXITER deve essere uno scalare')
   elseif (~isfinite(MAXITER) ||  ~isreal(MAXITER) || ischar(MAXITER)||isnan(MAXITER))
        error('Errore, MAXITER deve essere settato come un numero reale finito.')
   end
end
  
   %Controlli ipotesi Jacobi
   
   
      
   %Verifica se TOL appartiene a dei valori corretti
   if(sign(TOL)<=0||TOL<eps)
      TOL=10^-6;
      warning('Valore TOL errato. Utilizzo valore di default.') 
   end
    %Verifica se NMAX appartiene a dei valori corretti
   if(MAXITER<2||sign(MAXITER)<=0|| MAXITER>500)
      MAXITER=500;
      warning('Valore MAXITER errato. Utilizzo valore di default.') 
   end

   
   %Metodo Jacobi
   
   
   
   %Parametro output numero iterazioni
   if(nargout==2)
      niter=n;
   end
   
    %Parametro output residuo relativo
   if(nargout==3)
      resrel=norm(b-A*x)/norm(b);
   end


end

