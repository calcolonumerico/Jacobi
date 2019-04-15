function [x,niter,resrel] = jacobi(A,b,TOL,MAXITER)
%jacobi Summary of this function goes here
%   Detailed explanation goes here

 %Controlli numero input
   if(nargin<2)
        error("Inserire come input almeno A e b.")
   elseif(nargin==2)
       TOL=10^-6;
       MAXITER=500;
   elseif(nargin==3)
       MAXITER=500;
   end
  
 %Controlli sulla matrice

    if(~ismatrix(A))
       error("Il primo input deve essere una matrice.")
    elseif(~issparse(A)||isempty(A))
       error("La matrice deve essere sparsa.")
    elseif(size(A,1)~=size(A,2))
        error("La matrice deve essere quadrata.")
  %  elseif(length(A)<100)
     %   error("La matrice deve essere di grandi dimensioni(almeno 100x100)")
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

  %Controllo che la diagonale non abbia elementi nulli
   if any(find(abs(diag(A))<=eps(norm(A,Inf))))==1
       error("Sono presenti elementi nulli all'interno della diagonale principale.")
   end
   %Controllo sulla convergenza
  if(abs(diag(A))<=abs(sum(A,2)-diag(A)))
       error("La matrice non ha diagonale strettamante dominante. Il metodo non converge.")
  end
 
  %Verifica se TOL è uno scalare, e se lo è verifica se è ammissibile
   if(~isscalar(TOL))
       error('Errore, TOL deve essere uno scalare')
   elseif (~isfinite(TOL) ||  ~isreal(TOL) || ischar(TOL)||isnan(TOL))
        error('Errore, TOL deve essere settato come un numero reale.')
   end

   

  %Verifica se MAXITER è uno scalare, e se lo è verifica se è ammissibile;
   if(~isscalar(MAXITER))
       error('Errore, MAXITER deve essere uno scalare')
   elseif (~isfinite(MAXITER) ||  ~isreal(MAXITER) || ischar(MAXITER)||isnan(MAXITER))
        error('Errore, MAXITER deve essere settato come un numero reale finito.')
   end

     
   %Verifica se TOL appartiene a dei valori corretti
   if(sign(TOL)<=0||TOL<10^-6)
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
                   %Prealloco vettori
   x0=sparse(zeros(n,1));
   x= (speye(n)-spdiags(1./spdiags(A,0),0,n,n)*A)*x0 + spdiags(1./spdiags(A,0),0,n,n)*b;%? preso da corrado

   while(itr<MAXITER && (norm(x-x0)>norm(x0)*TOLX))
    x0=x;
    for i=1:n
        sigma=0;
        for j=1:n
            if j~=i
                sigma=sigma+A(i,j)*x(j);
            end
        end
         x(i)=(1/A(i,i))*(b(i)-sigma);
    end
    
    if(TOL*norm(x,Inf)>realmin)
                TOLX=TOL*norm(x,Inf);
            else
                TOLX=realmin;
            end
    
    
    itr=itr+1;
    end

    x=x(n);
   %Parametro output numero iterazioni
   if(nargout==2)
      niter=itr;
   end
   
    %Parametro output residuo relativo
   if(nargout==3)
      resrel=norm(b-A*x)/norm(b);
   end


end

