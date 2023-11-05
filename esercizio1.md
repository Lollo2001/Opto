# TITOLO
## LABORATORIO 1

### Punto 1
Per una guida d'onda in silicio di tipo slab avente l'indice di rifrazione del core pari a $3.48$ e l'indice di rifrazione del cludding pari a $1.44$, è stato riportato sul grafico sottostante l'indice di rifrazione efficace $n_{eff}$ in funzione dello spessore $d$ del core per il modo fondamentale $TE$ a $1.55\mu m$. I risultati, ottenuti applicando la condizione di continuita del campo magnetico lungo l'asse z per la soluzione dell'equazione d'onda, mostrano i soli primi tre modi $TE$ della guida. I modi con indice superiore al secondo hanno un andamento molto simile a quelli visualizzati nel grafico, ovvero tenderebbero anch'essi a $n_{core}$ per lo spessore del core che tende ad infinito, ma avrebbero l'intersezione con l'asse delle ascisse a valori sempre maggiori di $d$. Osservando il grafico è possibile valutare quale sia lo spessore più opportuno per il core in base alle nostre specifiche di progetto. In particolare prenendo uno spessore minore di circa $0.25 \mu m$ la guida è monomodale, ovvero presenta solo il modo $TE_0$. Prendendo invece uno spessore compreso tra circa $0.25 \mu m$ e $0.5 \mu m$ la guida è multimodale, ovvero presenta anche il modo $TE_1$ oltre al modo $TE_0$. Tale ragionamento può essere esteso per tutti i modi $TE$ ma è importante notare come al crescere delle dimensioni dello spessore del core scelto, cresce anche il numero di modi che si propagano nella guida. Come suggerito dal nome stesso una guida monomodale supporta un solo modo di propagazione e ciò significa che all'interno della guida è possibile propagare una sola frequenza. Per una guida multimodale, invece, è possibile propagare più modi ovvero sono presenti più onde aventi frequenze diverse. 

![n_eff(d)](figure/es1/es1_1.jpg)

### Punto 2
Per fare in modo che la guida d'onda in silicio specificata nel punto 1 abbiamo i soli modi $TE_0$ e $TE_1$, è necessario scegliere uno spessore del core compreso tra circa $0.25 \mu m$ e $0.5 \mu m$. Tutte le analisi seguenti sono condotte considerando uno slab avente lo spessore del core pari a $0.45 \mu m$. La scelta non è stata casuale in quanto l'intenzione era quella di avere una guida d'onda in cui il campo fosse il più confinato possibile dentro al core ($\Gamma$ più vicina possibile a $1$ come verrà descritto nel Punto 3) senza però avere anche il modo $TE_2$ nella guida.
Supponendo che la guida sia infinitivamente lunga in tutte le direzioni e supponendo che la discontinuità del materiale sia lungo l'asse x, dire che una guida sia $TE$ significa che essa ha il campo elettrico perpendicolare alla direzione di propagazione  dell'onda ($z$) e che ha componente del campo elettrico solo lungo l'asse $y$. Il campo magnetico, dovendo essere perpendicolare al campo elettrico, ha componente lungo l'asse $x$ e $z$. Nei grafici sottostanti sono state calcolate in maniera analitica le intensità dei campi $E_y$, $H_z$ e $H_x$ sia per il modo  $TE_0$ che per il modo $TE_1$ per mezzo delle equazioni di Maxwell (le linee verticali tratteggiate delimitano la porzione di core lungo l'asse $x$).
E' utile notare dai grafici che i campi hanno maggior intensità dentro al core della guida e che man mano che ci inoltriamo nel cludding allontanandoci dal core, l'intensità dei campi tende progressivamente a zero.

<div style="display: flex;">
    <img src="figure/es1/es1_2_Ey_modo0.jpg" alt="Ey_0" style="width: 50%;">
    <img src="figure/es1/es1_2_Ey_modo1.jpg" alt="Ey_1" style="width: 50%;">
</div>

<div style="display: flex;">
    <img src="figure/es1/es1_2_Hz_modo0.jpg" alt="Hz_0" style="width: 50%;">
    <img src="figure/es1/es1_2_Hz_modo1.jpg" alt="Hz_1" style="width: 50%;">
</div>
<div style="display: flex;">
    <img src="figure/es1/es1_2_Hx_modo0.jpg" alt="Hx_0" style="width: 50%;">
    <img src="figure/es1/es1_2_Hx_modo1.jpg" alt="Hx_1" style="width: 50%;">
</div>


### Punto 3

![n_eff(d)](figure/es1/es1_3.jpg)