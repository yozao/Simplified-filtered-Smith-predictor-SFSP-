function [Fr, kr, V, V_, S_]=SFSPd(Gd,Ts,d,p,bV,pV)
%
%
%   Essa função retorna os parametros do SFSP a tempo discreto
%
%   Gd é o modelo do processo discretizado e SISO
%
%   bV indica a velocidade de rejeição de perturbações, polos de V
%   
%   pV são os zeros do controlador primário equivalente, que permite
%   incluir integradores e também cancelar a dinamica de malha aberta
% 
%   O Fr é o filtro de referencia, que por padrão será unitário
%   kr é o controlador primário, trata-se de um ganho ou um vetor que
%   garante erro zero e regime estácionario
%
%   O filtro de robustez V tem o formato: numV/denV
%   numV=[v1, v2, v3, v4, ... ,vn]; 
%   denV=[b1, b2, b3, b5, ... ,bs, 1];
%   Onde n = s+1
%   numV são os indices do Polinômio v0*z^n+v1*z^(n-1)+...+vn
%   denV são os indices do Polinômio (z-b1)*(z-b2)*...*(z-bs)
%
%   V_ é o resultado da decomposição em frações parciais da expressão
%   ((zI-A)^-1)*V, referente ao termo que mantem o denominador denV.
%   S_ é a implementação dos termos faltantes da implementação estável S,
%   que é o proprio atraso de transporte distribuido.
%
%   Como boa pratica deve-se incluir em 'pV' os polos do processo em malha
%   aberta e pelo menos um integrador para rejeitar perturbações constantes
%
%   Os argumentos de entrada são: (Gd,d,p,bV,pV)
%   Os argumentos de saida são: [Fr, kr, V, V_, S_]
%
%   Se estamos sintonizando um controlador baseado em modelo, então, todos
%   os parametros utilizados são do modelo identificado.
%
%% modelo em espaço de estados
[nG,dG] =  tfdata(Gd,'v'); %obtém numerador e denominador do modelo
[A,C,B,~] = tf2ss(nG,dG);
A = A'; B = B'; C = C'; D=0;

%   
pV=[1 pV]; %Acrescenta o integrador padrão do controlador equivalente
nb=length(bV); %ordem do Filtro Beta
nr=length(pV); %Número de variaveis para o numerador
pV=sort(pV, 'descend'); %%organiza de forma decrescente
z = tf('z',Ts);

K = acker(A,B,p); %Ganho de realimentação
kr = eye(size(B,2))/(C/(eye(size(A)) + (-A + B*K))*B); %ganho de referencia

Gk = tf(ss(A,B,K,0,Ts)); 
[nGk,~] = tfdata(Gk,'v');
dGk=dG;

%% Filtro de referencia
Fr = 1;


%% Filtro de robustez
d=round(d/Ts);
if (nb+1 ~= nr) %verifica se a dimensões estão corretas
     return
end
% Constroi as  strings e os simbolicos
for k=0:(nb-1)
    Simbes(k+1)=string(strcat('v',num2str(k)));
end
syms(Simbes)
syms w
numGk = poly2sym(nGk,w);
numG = poly2sym(nG,w);
denG = poly2sym(dG,w);

denV=1;
numV=sym(strcat('v',num2str(nb)));
for k=1:(nb)
numV = numV+sym(strcat('v',num2str(k-1)))*w^((nb+1)-k);
denV = denV*(w-bV(k));
end

S_1 = simplify((denG*denV + numGk*denV - numG*numV*w^-d));
%% criando as matrizes do sistema linear aff*v=bff

%vetor de variaveis
for k=0:(nb)
simbol(1,nb+1-k)=[sym(strcat('v',num2str(k)))];
end

%% criando as equações para forma o sistema linear
pp=1e-12;
nn=1;
for j=1:nr
    if (pV(j)==pp)
        af(j)=subs(simplify(diff(S_1,w,nn)),w,pV(j));
        aux=double(coeffs(af(j),simbol));
        aff(j,:)=aux(2:end);
        bff(j,:)=-aux(1);
        nn=nn+1;
    else
        af(j)=subs(simplify(diff(S_1,w,0)),w,pV(j));
        aux=double(coeffs(af(j),simbol));
        aff(j,:)=aux(2:end);
        bff(j,:)=-aux(1);
        nn=1;
        pp=pV(j);
    end
end

%polinomio do denominador de V
aux1=1;
for j=1:nb
aux1 = aux1*(z-bV(j));
end
[dV,~,~]=tfdata(aux1,'v');

x = real(inv(aff)*bff);% resolve o sistema linear e obtém as variaveis
nV = [x']; %numerador do V

V=zpk(tf(nV, dV, Ts)); %constroi o V
%% Implementação estável de S(z)
aux1 =  zpk(ss(V*Gd));
aux2 = ss(A,B,K*A^d,0,Ts);
V_ = minreal(zpk(aux1-aux2),1e-2);

S_ = 0; %filtro FIR
for i = 1:d
    S_ = S_ + K*A^(i-1)*B*z^-i;
end

end

