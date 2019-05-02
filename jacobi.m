function [x,niter,resrel] = jacobi(A,b,TOL,MAXITER)

 %Controlli numero input
 if(nargin<2)
      error("Inserire come input almeno A e b.")
 end
  %Controlli sulla matrice
    if(~ismatrix(A))
       error("Il primo input deve essere una matrice.")
    elseif(~issparse(A)||isempty(A))
       error("La matrice deve essere sparsa e non vuota.")
    elseif(size(A,1)~=size(A,2))
        error("La matrice deve essere quadrata.")
    elseif(length(A)<100)
       error("La matrice deve essere di grandi dimensioni(almeno 100x100)")
    elseif(~isnumeric(A) ||~isreal(A)||any(~isfinite(A(:)))|| any(~isa(A(:),'double'))||any(isnan(A(:))))
      error("Gli elementi della matrice devono essere numeri reali double finiti.")
    end
    
   %Controllo che la diagonale non abbia elementi nulli
   if any(find(abs(spdiags(A,0))<=eps(norm(A,Inf))))==1
       error("Sono presenti elementi nulli all'interno della diagonale principale.")
   end
    
    
    %Controlli sul vettore 
   if(~isvector(b)||isscalar(b)||isrow(b))
       error('Errore, b deve essere un vettore colonna.')
   elseif(length(b)~=length(A))
       error("Il vettore b e la matrice A hanno dimensioni diverse.")
   elseif(any(ischar(b))||any(~isfinite(b)) ||  any(~isreal(b))||any(~isa(b,'double'))||any(isnan(b)))
       error('Errore, b deve contenere reali finiti double.')
   end
   

 if(nargin<3 ||isempty(TOL))
     TOL=10^-6;
  %Verifica se TOL è uno scalare, e se lo è verifica se è ammissibile
 elseif(~isscalar(TOL))
       error('Errore, TOL deve essere uno scalare')
   elseif (~isfinite(TOL) ||  ~isreal(TOL) || ischar(TOL)||isnan(TOL))
        error('Errore, TOL deve essere settato come un numero reale.')
 end

   
if(nargin<4 ||isempty(MAXITER))
    MAXITER=500;
  %Verifica se MAXITER è uno scalare, e se lo è verifica se è ammissibile;
elseif(~isscalar(MAXITER))
       error('Errore, MAXITER deve essere uno scalare')
   elseif (~isfinite(MAXITER) ||  ~isreal(MAXITER) || ischar(MAXITER)||isnan(MAXITER))
        error('Errore, MAXITER deve essere settato come un numero reale finito.')
end

     
   %Verifica se TOL appartiene a dei valori corretti
   if(sign(TOL)<=0||TOL<eps)
      TOL=10^-6;
      warning('Valore TOL errato. Utilizzo valore di default.') 
   end
    %Verifica se MAXITER appartiene a dei valori corretti
   if(MAXITER<2||sign(MAXITER)<=0|| MAXITER>500)
      MAXITER=500;
      warning('Valore MAXITER errato. Utilizzo valore di default.') 
   end

  
   
   %Metodo Jacobi
   TOLX = TOL;     %Tolleranza gestita dinamicamente         
   itr=1; 
   n=length(A); 
   d=spdiags(A,0);    
   %Prealloco vettori
   x0=sparse(zeros(n,1));
   x=(b-((sum(A.*x0,2)-d)))./d;
   val=TOLX*norm(x,Inf);
   
   while(itr<MAXITER && (norm(x-x0,Inf)>norm(x0,Inf)*TOLX))
      
      if(val>realmin)
                TOLX=val;
            else
                TOLX=realmin;
      end
      x0=x;
      x=(b-((sum((A.*x0),2)-d)))./d;
      itr=itr+1;
   end
   residuo=norm(b-A*x)/norm(b);
    
   %Parametro output numero iterazioni
   if(nargout>=2)
      niter=itr;
   end
   
    %Parametro output residuo relativo
   if(nargout==3)
      resrel=residuo;
   end
  
   %Se non convergo stampo a video warning
   if(itr >= MAXITER)
    warning('L ''algoritmo si è arrestato all''iterazione %d, senza convergere. Il residuo relativo è %f. ',itr,residuo)
   end

end

