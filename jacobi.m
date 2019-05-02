function [x,niter,resrel] = jacobi(A,b,TOL,MAXITER)
% jacobi La funzione risolve un sistema sparso di grandi dimensioni attraverso il metodo iterativo di Jacobi
% 
%Sintassi:
%[x,niter,resrel] = jacobi(A,b,TOL,MAXITER)
% 
%Descrizione:
% x=jacobi(A,b), trova e restituisce la soluzione x di un sistema sparso.Tolleranza e numero massimo di iterazioni sono impostate con valori di default. 
% x=jacobi(A,b,TOL), permette all'utente di specificare la Tolleranza desiderata. 
% x=jacobi(A,b,TOL,MAXITER), permette all'utente di specificare la Tolleranza desiderata e il numero massimo di iterazioni. 
% [x, niter]=jacobi(_), resituisce anche il numero di iterazioni effettuate per calcolare la soluzione.
% [x,niter,resrel]=jacobi(_) restituisce anche il valore del residuo relativo
% 
% Parametri di ingresso:
% A = matrice sparsa di grandi dimensioni. 
% b = vettore colonna dei termini noti.
% TOL (facoltativo) = tolleranza.Numero di cifre significative richieste dall'utente. Se non specificato è posta pari a 10^-6.
% MAXITER (facoltativo)= numero massimo di iterazioni.Limita il numero di iterazioni che la funzione può compiere. Se omesso è posto pari a 500.
% 
% Parametri di uscita:
% x = valore della soluzione del sistema. 
% 
% Diagnostica:
% Il programma si arresta mostrando un messaggio di errore nelle seguenti situazioni:
% Se i parametri di input sono meno di due. 
% Se il primo parametro d'ingresso non è una matrice congrua al problema. 
% Se il secondo parametro d'ingresso non è un vettore congruo al problema.
% Se sulla diagonale ci sono elementi nulli. 
% Se TOL o MAXITER sono valori non ammissibili per il problema.
% 
% Accuratezza:
% L'accuratezza dipende dal numero massimo di iterazioni(MAXITER) e dalla tolleranza (TOL) specificata e dal condizionamento della matrice.
% 
% Algoritmo:
% La funzione implementa l'algoritmo di jacobi.
% 
% Esempi di utilizzo:
% A=gallery('poisson',10); 
% x=ones(100,1); 
% b=A*x; 
% x=jacobi(A,b,10^-5,400)
% 
% x =
% 
%    1.000000476837158
%    .
%--------------------------------
% A=gallery('poisson',15); 
% x=ones(225,1);
%  b=A*x; 
% [x niter resrel]=jacobi(A,b,10^-5,400)
% 
% x =
% 
%   1.000000000931323
%   .
% 
% niter =
% 
%     18
% 
% resrel =
% 
%      3.925279456634614e-06
% Autori:
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

