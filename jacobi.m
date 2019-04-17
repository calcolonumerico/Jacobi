function [x,niter,resrel] = jacobi(A,b,TOL,MAXITER)
%jacobi La funzione risolve un sistema sparso di grandi dimensioni attraverso il
% metodo iterativo di Jacobi
%
% Sintassi:
% [x,niter,resrel]=jacobi(A,b)
% [x,niter,resrel]=jacobi(A,b,TOL)
% [x,niter,resrel]=jacobi(A,b,TOL,MAXITER)
%
% Descrizione:
% x=jacobi(A,b), trova e restituisce la soluzione x di un sistema sparso.Tolleranza e numero massimo di iterazioni sono impostate con valori di default.
% x=jacobi(A,b,TOL), permette all'utente di specificare la Tolleranza desiderata.
% x=jacobi(A,b,TOL,MAXITER), permette all'utente di specificare la Tolleranza desiderata 
% e il numero massimo di iterazioni.
% [x, niter]=jacobi(___), resituisce anche il numero di iterazioni
% effettuate per calcolare la soluzione.
% [x,niter,resrel]=jacobi(___) restituisce anche il valore del residuo
% relativo
%
% Parametri di ingresso:
%   A               = matrice sparsa di grandi dimensioni.
%   b               = vettore dei termini noti.
%   TOL (facolativo) = tolleranza gestita dall'utente.In caso di omissione,
%                        � posta pari a 10^-6.
%   MAXITER(facoltativo) = numero massimo di iterazioni. Se
%                       omesso � posto pari a 500.
%
% Parametri di uscita:
%   x               = valore della soluzione del sistema.
%   niter (facoltativo)  = numero di iterazioni effettuate. 
%   resrel (facoltativo) = residuo relativo.
%
% Diagnostica:
% Il programma si arresta mostrando un messaggio di errore nelle seguenti situazioni:
%   -Se i parametri di input sono meno di due.
%  - Se il primo parametro d'ingresso non � una matrice congrua al problema.
%  - Se il secondo parametro d'ingresso non � un vettore congruo al problema.
%  - Se sulla diagonale ci sono elementi nulli.
%  - Se TOL o MAXITER sono valori non ammissibili per il problema.
% 
% Accuratezza:
%L'accuratezza dipende dal numero massimo di iterazioni(MAXITER) e dalla tolleranza 
%(TOL) specificata.
%
% Algoritmo
% La funzione implementa l'algoritmo di jacobi.
%
% Esempi di utilizzo:
%  A=gallery('poisson',10);
%  x=ones(100,1);
%  b=A*x;
%  x=jacobi(A,b,10^-5,400);
% 
%  
% x =
% 
%    0.999999999985448
%
%----------------------------------------------------------------------
%
%  A=gallery('poisson',10);
%  x=ones(100,1);
%  b=A*x;
% [x niter]=jacobi(A,b,10^-5,400);
% x =
% 
%    0.999999999985448
% niter =
% 
%     36
% 
%
%----------------------------------------------------------------------
%  A=gallery('poisson',15);
%  x=ones(225,1);
%  b=A*x;
% [x niter resrel]=jacobi(A,b,10^-5,400)
% x =
% 
%    0.999999999992724
% 
% 
% niter =
% 
%     37
% 
% 
% resrel =
% 
%    1.283330914010793
%
%   Autori:
%       Iodice Ivano
%       Vincenzo De Francesco

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
    elseif(~isnumeric(A) ||~isreal(A)||any(any(~isfinite(A)))|| any(any(~isa(A,'double')))||any(any(isnan(A))))
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
  %Verifica se TOL � uno scalare, e se lo � verifica se � ammissibile
 elseif(~isscalar(TOL))
       error('Errore, TOL deve essere uno scalare')
   elseif (~isfinite(TOL) ||  ~isreal(TOL) || ischar(TOL)||isnan(TOL))
        error('Errore, TOL deve essere settato come un numero reale.')
 end

   
if(nargin<4 ||isempty(MAXITER))
    MAXITER=500;
  %Verifica se MAXITER � uno scalare, e se lo � verifica se � ammissibile;
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

  %Controllo sulla convergenza
  if(abs(spdiags(A,0))<=abs(sum(A,2)-spdiags(A,0)))
       warning("La matrice non ha diagonale strettamante dominante. Il metodo potrebbe non convergere.")
  end
   
   %Metodo Jacobi
   TOLX = TOL;     %Tolleranza gestita dinamicamente         
   itr=1; 
   n=length(A); 
                 %Prealloco vettori
   x0=sparse(zeros(n,1));
   x= (speye(n)-spdiags(1./spdiags(A,0),0,n,n)*A)*x0 + spdiags(1./spdiags(A,0),0,n,n)*b;

   while(itr<MAXITER && (norm(x-x0)>norm(x0)*TOLX))
      if(TOL*norm(x,Inf)>realmin)
                TOLX=TOL*norm(x,Inf);
            else
                TOLX=realmin;
      end
      x0=x;
      x=(b-((sum(A,2)-spdiags(A,0)).*x0))./spdiags(A,0);
      itr=itr+1;
   end
   residuo=norm(b-A*x)/norm(b);
   x=x(n);   
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
    warning('L ''algoritmo si � arrestato all''iterazione %d, senza convergere. Il residuo relativo � %f. ',itr,residuo)
   end

end
